import polars as pl
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib_venn import venn3_circles, venn3

ip_protein = "Tbck"
analysis_directory = Path(
    '/Users/henrysanford/Dropbox @RU Dropbox/Vinogradova Laboratory/Macrophage project/Manuscript 1/01_Data Analysis/01_Proteomics/04_IP-MS/03_analysis'
)
data_directory = Path(
    analysis_directory /"census-out-processing"
)

ip_ms_files = Path(
    data_directory / ip_protein.lower() / "output" / ip_protein.lower()
)
results_dir = ip_ms_files / "04_results"
results_dir.mkdir(exist_ok=True)

venn_diagram_dir = analysis_directory / "endogenous_pulldown_analysis/venn_diagram_test/"

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
    return((
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
            ]
        )
    ))


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

# original ferry paper comparison
schuhmacher_protein = set(
    pl.read_excel(
        '/Users/henrysanford/Dropbox @RU Dropbox/Vinogradova Laboratory/Macrophage project/Manuscript 1/01_Data Analysis/01_Proteomics/00_reference lists/Ferry_paper/mmc3.xlsx', sheet_name="GST-FERRY Proteomics"
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
master = pl.DataFrame({
    "protein": all_elements,
    **{label: [int(el in s) for el in all_elements] for label, s in membership.items()},
})

for ap, tbl in wang_tables.items():
    master = master.join(
        tbl.rename({c: f"{ap}: {c}" for c in tbl.columns if c != "protein"}),
        on="protein",
        how="left",
    )

master.write_csv(venn_diagram_dir / "master_intersection_table.csv")


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
    plt.show()


venn_3_wrapper([SCHUHMACHER, FY2, TBCK], "Schumacher")
for ap in wang_supplement_map:
    venn_3_wrapper([wang_label(ap), FY2, TBCK], f"Wang_{ap}")
venn_3_wrapper([wang_label(ap) for ap in wang_supplement_map], "Wang_three_pulldowns")
