import os
import polars as pl
import polars.selectors as cs
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn3_circles, venn3
from itertools import product
import pandas as pd
import scipy.stats as stat
import numpy as np

from src.uniprot_utils import create_entry_cache, get_function, get_go_terms

# Embed text as editable TrueType fonts (type 42) rather than the default
# Type 3 fonts, which Illustrator opens as outlines instead of text boxes.
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

ip_protein = "Fy2"
analysis_directory = Path(
    '/Users/henrysanford/Dropbox @RU Dropbox/Vinogradova Laboratory/Macrophage project/Manuscript 1/01_Data Analysis/01_Proteomics/04_IP-MS/03_analysis'
)
data_directory = Path(
    analysis_directory /"census-out-processing" 
)
protein_id_columns = ["uniprot", "protein", "description"]

# CRAPome key rebuilt for the macrophage protein list (see build_crapome_key.py)
CRAPOME_KEY_PATH = Path(__file__).parent / "crapome" / "macrophage_crapome_key.csv"

FERRY_MEMBERS = ["TBCK", "CRYZL1", "C12orf4", "PPP1R21", "GATD1"]
ENRICHMENT_FC_CUTOFF = 1.5
P_VALUE_CUTOFF = 0.05
APPLY_2_PEP_REQ = True

def process_ip_protein(ip_protein):
    ip_ms_files = Path(
        data_directory / ip_protein.lower() / "output" / ip_protein.lower()
    )
    results_dir = ip_ms_files / "04_results"
    results_dir.mkdir(exist_ok=True)


    donors = {
        "processed_census-out_20260507_TZ3-154": "d1",
        "processed_census-out_20260517_TZ3-156": "d2",
        "processed_census-out_20260525_TZ3-158": "d3",
        'processed_census-out_20260517_TZ3-155':"d1", 
        'processed_census-out_20260428_TZ3-152_IP':"d2", 
        'processed_census-out_20260525_TZ3-157':"d3"
    }

    def convert_to_long_format(df, id_vars, values_name):
        condition_renames = {
            "M0_IgG": "M0-IgG",
            "LPS_IgG": "LPS-IgG",
            f"M0_{ip_protein}": f"M0-{ip_protein}",
            f"LPS_{ip_protein}": f"LPS-{ip_protein}",
        }
        df = (
            # convert to long format
            df.unpivot(
                index=id_vars,
                value_name=values_name,
                variable_name="sample",
            ).drop_nans(subset=values_name)
            # metadata cleaning
            .with_columns(pl.col("sample").str.replace_many(condition_renames))
        )

        print(df)

        df = df.with_columns(  # extract donor from sample name
            pl.col("sample")
            .str.splitn("_", 3)
            .struct.rename_fields(["condition", "technical_replicate", "donor"])
            .alias("parts")
        ).unnest("parts").with_columns(
            pl.col("donor").replace(donors)
        )
        print(set(df["donor"]))

        return df


    # read replicate level data
    df = pl.read_csv(
        ip_ms_files
        / f"03_combined_files/{ip_protein.lower()}/03_combfiles_forpca_channelratio_or_rawsignal_wp.csv"
    )

    # peptide number logic
    pep_num = (
        df
        .select(cs.by_name(protein_id_columns), cs.contains("pepNum"))
        .unpivot(on=cs.contains("pepNum"), index=protein_id_columns, value_name="pep_num")
        .with_columns(pl.col("variable").str.replace("pepNum_", "").replace(donors).alias("donor")).drop("variable")
    )

    df = df.filter(
        pl.col("uniprot").str.contains("contaminant").not_(),
        pl.col("description").str.contains("Keratin").not_(),
    ).drop(cs.contains("pepNum", "empty"), "")

    df = convert_to_long_format(
        df, protein_id_columns, "signal_intensity"
    )


    def filter_to_min_2_reps(
        df: pl.DataFrame,
        id_vars: list,
    ) -> pl.DataFrame:
        """Filters dataframe to targets quantified in
        2 biological replicates"""
        # filter to min two biological replicates
        return (
            df.group_by(id_vars)
            .agg(pl.col("donor").n_unique().alias("n_replicates"))
            .filter(pl.col("n_replicates") > 1)
            .drop("n_replicates")
            .join(other=df, on=id_vars)
        )


    n_proteins = len(set(df["protein"]))
    print(f"{n_proteins} before filtering")

    df = (
        df.with_columns(
            pl.concat_str(
                ["condition", "donor", "technical_replicate"], separator="_"
            ).alias("clean_sample_name")
        )
        .with_columns(
            pl.col("condition").str.split("-").list.get(0).alias("super_condition"),
            pl.col("condition").str.split("-").list.get(1).alias("enrichment"),
        )
        .join(other=pep_num, on=protein_id_columns + ["donor"], how="left")
    )

    df.write_csv(data_directory / f"unfiltered_{ip_protein}.csv")

    if APPLY_2_PEP_REQ:
        df = df.filter(pl.col("pep_num").gt(1))

        n_proteins = len(set(df["protein"]))
        print(f"{n_proteins} after filtering to 2+ peptides")

    df = filter_to_min_2_reps(
        df,
        protein_id_columns,
    )

    n_proteins = len(set(df["protein"]))
    print(f"{n_proteins} after filtering to 2+ reps")


    if ip_protein == "Fy2":
        protein_name_to_normalize_to = "PPP1R21"
    elif ip_protein == "Tbck":
        protein_name_to_normalize_to = "TBCK"
    else:
        AssertionError("not valid protein")

    fy2_flag_si = df.filter(
        # filter to protein with FLAG
        pl.col("protein").eq(protein_name_to_normalize_to),
        # filter to conditions with FLAG
        pl.col("condition").str.contains_any([ip_protein]),
    )

    #display(fy2_flag_si.select("protein", "sample","signal_intensity",))

    normalization_grouping = ["donor"]

    median_si = (fy2_flag_si.group_by(normalization_grouping)
        .agg(pl.col("signal_intensity").median().alias("median_si")))

    #display(median_si)

    fy2_normalization_factors = (
        median_si
        # add fy2 median signal intensity values
        .join(other=fy2_flag_si, on=normalization_grouping, how="left")
        # calculate normalization factors
        .with_columns(
            pl.col("median_si").truediv("signal_intensity").alias("normalization_factor")
        )
    )

    #display(fy2_normalization_factors.select("condition","donor","technical_replicate", "normalization_factor"))

    normalized_df = (
        df.join(other=fy2_normalization_factors.select("sample", "normalization_factor"), on="sample", how="left")
        .with_columns(pl.col("normalization_factor").fill_null(1))
        .with_columns(pl.col("signal_intensity").mul("normalization_factor"))
    )

    # view data in a boxplot
    boxplot_data = pl.concat(
        items=[
            normalized_df.with_columns(
                pl.lit(f"normalized to {ip_protein} median").alias("normalization")
            ),
            df.with_columns(pl.lit("unnormalized").alias("normalization")),
        ],
        how="diagonal",
    )

    g = sns.FacetGrid(
        boxplot_data.sort(by = "enrichment"),
        row="normalization",
        hue="condition",
        sharex=False,
        height=4,
        row_order=["unnormalized", f"normalized to {ip_protein} median"],
        aspect=1.2,
        margin_titles=True,
    )
    g.map_dataframe(
        sns.boxplot,
        x="clean_sample_name",
        y="signal_intensity",
        fliersize=1,
        log_scale=True,
    )
    g.fig.set_size_inches(7, 4) 
    g.set_titles(col_template="{col_name}", row_template="{row_name}")

    for ax in g.axes[0, :]:  # First row (index 0)
        ax.set_xticklabels([])
    g.add_legend()
    g.set_axis_labels("Sample", "Signal Intensity")
    # plt.tight_layout()
    g.set_xlabels("")
    g.set_xticklabels(rotation=90, ha="right")
    #plt.show()

    channel_ratio_df = normalized_df.rename({"signal_intensity":"channel_ratio"})

    sns.kdeplot(
        data = channel_ratio_df,
        x = "channel_ratio",
        hue = "condition",
    )



    # taking median across technical replicates within each donor
    donor_condition_medians = (
        channel_ratio_df.filter(pl.col("condition").str.contains(f"M0-{ip_protein}"))
        .group_by(["donor", "condition"] + protein_id_columns)
        .agg(pl.col("channel_ratio").median().alias("control_median"))
        .select([ "control_median", "donor"] + protein_id_columns)
    )

    #display(donor_condition_medians)

    normalized_to_controls = channel_ratio_df.join(
        other=donor_condition_medians, on=["donor"] + protein_id_columns
    ).with_columns(
        pl.col("channel_ratio").truediv("control_median").alias("ratio_to_control")
    )


    sns.set_style("whitegrid")
    ax = sns.boxplot(
        data=normalized_to_controls.filter(
            pl.col("condition").str.contains("ptwist").not_()
        ).sort(by="clean_sample_name"),
        y="clean_sample_name",
        x="ratio_to_control",
        fliersize=1,
        log_scale=True
    )
    plt.xlabel("log2 FC from IgG")
    plt.ylabel("")






    def get_p_value(row, cond_1, cond_2):
        ttest_result = stat.ttest_ind(row[cond_1], row[cond_2], nan_policy="omit")
        return ttest_result[1]


    def get_expr_unadj(row):
        if (row["log2_FC"] > np.log2(ENRICHMENT_FC_CUTOFF)) & (
            row["-log10_pval"] > -np.log10(P_VALUE_CUTOFF)
        ):
            return "Significant Up"
        elif (row["log2_FC"] < -np.log2(ENRICHMENT_FC_CUTOFF)) & (
            row["-log10_pval"] > -np.log10(P_VALUE_CUTOFF)
        ):
            return "Significant Down"
        else:
            return "Stable"


    def volcano_plot(df, conditions_list):

        # without p_adj

        list_cond_1 = df.filter(like=conditions_list[0]).columns.tolist()
        list_cond_2 = df.filter(like=conditions_list[1]).columns.tolist()

        important_col_df = df[list_cond_1 + list_cond_2]

        idx_cond_1 = important_col_df.columns.get_indexer(list_cond_1)
        idx_cond_2 = important_col_df.columns.get_indexer(list_cond_2)

        important_col_df["p_value"] = important_col_df.apply(
            get_p_value,
            axis=1,
            args=(idx_cond_1, idx_cond_2),
        )

        for condition in conditions_list:
            list_columns_cond = important_col_df.filter(like=condition).columns.tolist()
            cond_df = important_col_df[list_columns_cond]

            important_col_df["mean_" + condition] = cond_df.mean(axis=1)

        important_col_df["FC"] = (
            important_col_df["mean_" + conditions_list[0]]
            / important_col_df["mean_" + conditions_list[1]]
        )
        important_col_df["log2_FC"] = np.log2(important_col_df["FC"])  #

        volcano_df = important_col_df[["p_value", "log2_FC"]]
        volcano_df = volcano_df.dropna()

        volcano_df["p_value"] = volcano_df["p_value"].astype("float")

        volcano_df["-log10_pval_adj"] = -1 * np.log10(
            len(volcano_df) * volcano_df["p_value"]
        )
        volcano_df["-log10_pval"] = -1 * np.log10(volcano_df["p_value"])

        volcano_df["Regulation"] = volcano_df.apply(get_expr_unadj, axis=1)
        volcano_df["Regulation"] = volcano_df["Regulation"].astype("category")

        volcano_df = volcano_df.reset_index()
        volcano_df = volcano_df.set_index(["uniprot", "description"])
        return volcano_df


    ttest_input = (
        channel_ratio_df.select(
            ["clean_sample_name", "channel_ratio"] + protein_id_columns
        )
        .with_columns(pl.col("channel_ratio"))
        .pivot(index=protein_id_columns, on="clean_sample_name")
        .to_pandas()
        .set_index(protein_id_columns)
    )

    comparisons = [
        [f"LPS-{ip_protein}", "LPS-IgG"],
        [f"M0-{ip_protein}", "M0-IgG"],
        [f"LPS-{ip_protein}", f"M0-{ip_protein}"],
    ]
    significance_test_results = []
    for comparison in comparisons:
        significance_test = volcano_plot(df=ttest_input, conditions_list=comparison)

        significance_test["comparison"] = comparison[0] + " vs. " + comparison[1]

        significance_test_results.append(pl.from_pandas(significance_test))


    # represent inf p-value as zero
    # this is a result of normalizing to Fy2
    # in scipy zero variance = inf t-test result
    volcano_df = pl.concat(significance_test_results).with_columns(pl.col("-log10_pval").replace(np.inf, 0))


    #volcano_df.write_csv(results_dir / "volcano_plot/volcano_data_long.csv")

    volcano_df

    plot_data = volcano_df.filter(pl.col("comparison") == f"LPS-{ip_protein} vs. M0-{ip_protein}")

    sns.scatterplot(
        data = plot_data,
        x = "log2_FC",
        y = "-log10_pval"
    )

    labeled_proteins = [ "PPP1R21", "TBCK",'CRYZL1', "PURA", "PURB"]
    labeled_data = plot_data.filter(
        (pl.col("protein").is_in(labeled_proteins))
    )

    for row in labeled_data.iter_rows(named=True):
        plt.text(
            row["log2_FC"], 
            row["-log10_pval"], 
            f"  {row['protein']}", 
            fontsize=8, 
            ha='left', 
            va='center'
        )

    plot_data.sort(by = "log2_FC")

    wide_volcano_df = (
        volcano_df.drop("-log10_pval_adj")
        .rename({"-log10_pval": "neg_log10_pval"})
        .pivot(on=["comparison"], values=["p_value", "log2_FC", "neg_log10_pval", "Regulation"])
    )

    # add comparison to original ferry paper
    schumacher_protein = set(
        pl.read_excel(
            '/Users/henrysanford/Dropbox @RU Dropbox/Vinogradova Laboratory/Macrophage project/Manuscript 1/01_Data Analysis/01_Proteomics/00_reference lists/Ferry_paper/mmc3.xlsx', sheet_name="GST-FERRY Proteomics"
        ).filter("FERRY-complex_enriched")["gene_name"]
    )
    formatted_df = pl.from_pandas(ttest_input.reset_index()).with_columns(
        pl.col("protein").is_in(schumacher_protein).alias("GST-FERRY IP protein")
    )

    # add uniprot function/go-terms
    cache = create_entry_cache(formatted_df)

    formatted_df = formatted_df.with_columns([
        pl.col("uniprot").map_elements(
            lambda x: get_function(x, cache=cache),
        ).alias("uniprot_function"),
        pl.col("uniprot").map_elements(
            lambda x: get_go_terms(x, cache=cache),
        ).alias("uniprot_goterms")
    ])

    formatted_df = formatted_df.join(other=wide_volcano_df, on="protein")

    crapome = pl.read_csv(
        CRAPOME_KEY_PATH
    ).with_columns(
        pl.col("Num of Expt. (found/total)")
        .str.split(" / ")
        .list.get(0)
        .alias("Crapome Num of Expt. (out of 716 total)")
    ).rename({"User Input": "uniprot", "contaminant": "Crapome contaminant"}).select(["uniprot", "Crapome contaminant", "Crapome Num of Expt. (out of 716 total)"])


    set(
        formatted_df.with_columns(pl.col("uniprot_goterms").str.split("|"))
        .explode("uniprot_goterms")
        .filter(pl.col("uniprot_goterms").str.contains("ribosom"))["uniprot_goterms"]
    )


    ferry_member_names = {
        "TBCK": "Fy1",
        "PPP1R21": "Fy2",
        "C12orf4": "Fy3",
        "CRYZL1": "Fy4",
        "GATD1": "Fy5",
    }

    rbps = set(
        pl.read_excel(
            "/Users/henrysanford/Dropbox @RU Dropbox/Vinogradova Laboratory/Macrophage project/Manuscript 1/01_Data Analysis/01_Proteomics/00_reference lists/RBP/Census_of_RBP_2014.xlsx",
            sheet_name="RBP table",
        )["gene name"]
    )


    formatted_df = formatted_df.with_columns(
        # pl.col("protein").is_in(ip_enriched_proteins).alias("IP enriched"),
        pl.col("uniprot_goterms").str.contains("ribosom").alias("ribosomal"),
        pl.col("protein").replace(ferry_member_names, default=None).alias("FERRY subunits"),
        pl.col("protein").is_in(rbps).alias("RBP"),
    ).with_columns(
        pl.when(pl.col("FERRY subunits").is_not_null())
        .then(pl.lit("FERRY"))
        .when(pl.col("ribosomal"))
        .then(pl.lit("Ribosomal"))
        .when(pl.col("RBP"))
        .then(pl.lit("RBP"))
        .otherwise(pl.lit("Other"))
        .alias("category")
    ).select(
        ~cs.contains("-"),
        cs.contains("-"),
    ).join(
        other=crapome, on="uniprot", how="left"
    )

    formatted_df.write_csv(
        data_directory / f"ipms_table_{ip_protein}.csv"
    )

    volcano_df.join(
        formatted_df.select("protein", "uniprot", "category").unique(),
        on="protein",
        how="left",
    ).drop("-log10_pval_adj").with_columns(pl.lit(ip_protein).alias("dataset")).write_csv(data_directory / f"long_volcano_data_{ip_protein}.csv")

for ip_protein in ["Fy2", "Tbck"]:
    process_ip_protein(ip_protein)
pl.concat([pl.read_csv(data_directory / f"long_volcano_data_{p}.csv") for p in ["Tbck", "Fy2"]], how = "vertical").write_csv(data_directory / "long_volcano_data.csv")
pl.concat([pl.read_csv(data_directory / f"unfiltered_{p}.csv") for p in ["Tbck", "Fy2"]], how = "vertical").write_csv(data_directory / "unfiltered_data.csv")


ip_ms_files = Path(data_directory / ip_protein.lower() / "output" / ip_protein.lower())
results_dir = ip_ms_files / "04_results"
results_dir.mkdir(exist_ok=True)

venn_diagram_dir = (
    analysis_directory / "endogenous_pulldown_analysis/2_plus_reps_2_plus_peps/venn_diagram/"
)

venn_diagram_dir.mkdir(exist_ok=True)
wang_supplement_map = {
    "Tbck": "Table. S11",
    "Fy2": "Table. S12",
    "C12orf4": "Table. S13",
}

# membership column labels (also used as venn set labels)
SCHUHMACHER = "Schuhmacher et. al.\nGST-FERRY IP, 2023"
FY2 = "Macrophage\n Fy2 IP"
TBCK = "Macrophage\n Tbck IP"
wang_label = lambda ap: f"Wang et. al.\n{ap} AP-MS"


def read_wang_table(ap_protein):
    return (
        pl.read_excel(
            analysis_directory
            / "reference_data"
            / "1-s2.0-S2095927325010576-mmc1.xlsx",
            sheet_name=wang_supplement_map[ap_protein],
            read_options={"header_row": 1},
        )
        .filter(
            pl.col("Protein FDR Confidence: Combined").eq("High"),
            pl.col("^Abundance: F\d+: Sample$").gt(5000),
            pl.col("# Peptides").gt(1),
            pl.col("Description").str.contains("GN="),
        )
        .with_columns(
            pl.col("Description")
            .str.split("GN=")
            .list.get(1)
            .str.split(" ")
            .list.get(0)
            .alias("protein")
        )
        .drop(
            [
                "Checked",
                "Protein FDR Confidence: Combined",
                "Master",
                "# Protein Groups",
                r"^Found in Sample: \[S\d+\] F\d+: Sample$",
                "Accession",
                "Description",
            ]
        )
    )


# --- enriched protein sets ---
df = pl.read_csv(data_directory / "long_volcano_data.csv")


def enriched(name):
    return set(
        df.filter(
            pl.col("Regulation").eq("Significant Up"),
            pl.col("comparison").str.contains("IgG"),
            pl.col("comparison").str.contains(name),
        )["protein"]
    )


fy2_enriched_proteins = enriched("Fy2")
tbck_enriched_proteins = enriched("Tbck")

# --- GO enrichment (BP/MF/CC) on proteins enriched in Tbck and Fy2 vs IgG ---
# run each IP protein separately rather than pooling the two gene sets
# restrict to non-crapome proteins (Crapome.org contaminant flag = found in > 200
# of 716 experiments). The key is rebuilt for the macrophage protein list by
# ipms/build_crapome_key.py so every quantified protein is covered.

crapome = (
    pl.read_csv(CRAPOME_KEY_PATH)
    .rename({"User Input": "uniprot"})
    .select(["uniprot", "contaminant"])
)

joined = df.join(crapome, on="uniprot", how="left")

# surface proteins absent from the key instead of silently treating them as
# clean (the original bug: null contaminant passed the filter)
unmatched = sorted(set(joined.filter(pl.col("contaminant").is_null())["protein"]))
if unmatched:
    print(
        f"WARNING: {len(unmatched)} protein(s) not found in the Crapome key "
        f"(treated as non-crapome): {unmatched}"
    )

# gene names flagged as crapome contaminants by Crapome.org
crapome_proteins = set(
    joined.filter(pl.col("contaminant") == True)["protein"]
)


def non_crapome(proteins):
    return sorted(p for p in proteins if p not in crapome_proteins)


go_gene_sets = {
    #"Fy2": non_crapome(fy2_enriched_proteins),
    "Tbck": non_crapome(tbck_enriched_proteins),
}

# Export the crapome-filtered enriched gene sets. Fisher-exact GO enrichment runs in R via go_enrich()
# with the full human background 
pl.DataFrame(
    {
        "dataset": [dataset for dataset, genes in go_gene_sets.items() for _ in genes],
        "gene": [gene for dataset, genes in go_gene_sets.items() for gene in genes],
    }
).write_csv(data_directory / "go_gene_set.csv")

# original ferry paper comparison
schuhmacher_protein = set(
    pl.read_excel(
        "/Users/henrysanford/Dropbox @RU Dropbox/Vinogradova Laboratory/Macrophage project/Manuscript 1/01_Data Analysis/01_Proteomics/00_reference lists/Ferry_paper/mmc3.xlsx",
        sheet_name="GST-FERRY Proteomics",
    ).filter("FERRY-complex_enriched")["gene_name"]
)

# --- build one master binary table (with Wang metadata) ---
wang_tables = {ap: read_wang_table(ap) for ap in wang_supplement_map}

membership = {
    SCHUHMACHER: schuhmacher_protein,
    FY2: fy2_enriched_proteins,
    TBCK: tbck_enriched_proteins,
    **{wang_label(ap): set(tbl["protein"]) for ap, tbl in wang_tables.items()},
}

all_elements = sorted(set.union(*membership.values()))
master = pl.DataFrame(
    {
        "protein": all_elements,
        **{
            label: [int(el in s) for el in all_elements]
            for label, s in membership.items()
        },
    }
)

for ap, tbl in wang_tables.items():
    master = master.join(
        tbl.rename(
            {c: f"Wang et. al. {ap}: {c}" for c in tbl.columns if c != "protein"}
        ),
        on="protein",
        how="left",
    )

master.write_excel(
    venn_diagram_dir / "intersection_table.xlsx",
    column_widths=150,
    freeze_panes=(1, 1),
)


def venn_3_wrapper(columns, name):
    sets = [set(master.filter(pl.col(c).eq(1))["protein"]) for c in columns]
    fig, ax = plt.subplots(figsize=(2, 2))
    v = venn3(
        sets,
        set_labels=columns,
        ax=ax,
        set_colors=["#FFCC31", "#4869B2", "green"],
    )
    venn3_circles(sets, ax=ax, linewidth=0.25)
    for text in [*v.set_labels, *v.subset_labels]:
        if text:
            text.set_fontsize(6)
            text.set_color("black")
            text.set_fontfamily("Arial")

    plt.tight_layout()
    plt.savefig(venn_diagram_dir / f"{name}_venn_diagram.pdf")
    #plt.show()


venn_3_wrapper([SCHUHMACHER, FY2, TBCK], "Schumacher")
for ap in wang_supplement_map:
    venn_3_wrapper([wang_label(ap), FY2, TBCK], f"Wang_{ap}")
venn_3_wrapper([wang_label(ap) for ap in wang_supplement_map], "Wang_three_pulldowns")