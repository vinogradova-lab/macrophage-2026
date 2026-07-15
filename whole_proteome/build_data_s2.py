"""Assemble Data S2 (whole proteome): a contents index + seven flat sheets.

  S2-1 Unenriched proteomics  per protein: per-replicate abundances, per-stimulus-vs-M0 stats,
                              functional-category flags.
  S2-2 Bulk RNA vs WP         per gene/protein: per-comparison RNA (DESeq2) and WP fold changes.
  S2-3 Reactivity changes     per protein: per-stimulus IA-DTB reactivity blocks.
  S2-4 Immunoprecipitation MS the two FERRY-subunit IP-MS pulldowns merged on uniprot.
  S2-5 GO-term enrichment     per significant GO term: Fisher-exact enrichment across the WP,
                              reactivity, and phospho analyses, tagged with dataset/background.
  S2-6 GSEA of PC loadings    per significant GO:BP term: pre-ranked GSEA of the whole-proteome
                              PCA loadings, one block per principal component.
  S2-7 Phosphoproteomics      per phosphosite: expression-normalized per-replicate ratios,
                              TLR4-vs-M0 stats, and UniProt function.

Run after wp_downstream_analysis.ipynb has generated the volcano/wp_vs_rnaseq/GSEA tables and
the three visualization .Rmd files have written their go_enrich() CSVs (see the S2-5 loaders).
S2-7 additionally needs phosphoproteomics.ipynb to have written the expression-normalized
phosphorylation table to the manuscript tree (see PHOSPHO_TABLE_CSV).

    conda run -n polars python whole_proteome/build_data_s2.py
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl

# make `import src.*` work when run from anywhere in the repo
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
from src.cross_ref import (  # noqa: E402
    AM_BENIGN,
    AM_PATHOGENIC,
    CLINVAR_RESIDUE_TOLERANCE,
    am_class,
    conservation_bin,
    load_complex_membership,
    load_cys_cross_ref,
)
from src.macrophage import (  # noqa: E402
    CONDITIONS,
    FIGURE1_PANELS,
    MANUSCRIPT,
    PROTEOMICS,
    RC_CONDITIONS,
    require,
)
from src.reactivity import (  # noqa: E402
    CONDITION_BLOCK,
    has_reactivity_change,
    load_rc_long,
    pivot_rc_wide,
    rc_change_positions,
)
from src.supp_data import SuppSheet, write_supplementary_workbook  # noqa: E402

# Same data root the notebook uses (input_folder_path in wp_downstream_analysis.ipynb).
INPUT_FOLDER = Path(
    "/Users/henrysanford/dev/test_data/macrophage/01_whole_proteome/03_results/3_reps"
)
ABUNDANCE_CSV = (
    INPUT_FOLDER
    / "03_combined_files/1/03_combfiles_forpca_channelratio_or_rawsignal_wp.csv"
)
VOLCANO_CSV = INPUT_FOLDER / "04_results/volcano_plots/volcano_data.csv"
RNA_VS_WP_DIR = INPUT_FOLDER / "04_results/wp_vs_rnaseq"

# Cysteine reactivity changes come from src.reactivity, which pivots the current run's canonical
# long-format table (src.macrophage.RC_LONG_FORMAT). The old test_data rc_df.csv is NOT used: it
# predates the passes_ratio_in_two_replicates filter and overcounts reactivity changes.

# Endogenous IP-MS pulldowns of the two FERRY subunits (already supplement-ready tables).
IPMS_DIR = PROTEOMICS / "04_IP-MS/03_analysis/census-out-processing"
# Functional-category flags: the notebook's canonical volcano_data_long_format.csv (2026.1.Hs
# MSigDB query_gene_ontology_msigdb gene_dicts + HALLMARK_GLYCOLYSIS; same file
# multi_modal_significance.py reads). NOT the stale volcano_data_with_group_annotation.csv.
ANNOTATION_CSV = FIGURE1_PANELS / "volcano_plot/2026_ontology/volcano_data_long_format.csv"

OUTPUT_XLSX = REPO_ROOT / "Data S2.xlsx"
ID_COLS = ["uniprot", "protein", "description"]
# Comparisons follow the shared CONDITIONS order, matching the WP volcano column order.

WP_SHEET_NAME = "S2-1 Unenriched proteomics"
WP_TITLE = (
    "Whole proteome data from M0 and TLR/STING stimulated macrophages "
    "(related to Figure 1 and Extended Data Fig. 1)"
)
WP_DESCRIPTION = (
    "TMT-exp data showing protein expression changes in macrophages stimulated with TLR/STING "
    "agonists compared to M0 macrophages. Data are from n = 3 donors."
)

RNA_SHEET_NAME = "S2-2 Bulk RNA vs WP"
RNA_TITLE = (
    "Whole proteome and transcriptome data from M0 and TLR/STING stimulated macrophages "
    "(related to Figure 1 and Extended Data Fig. 1)"
)
RNA_DESCRIPTION = (
    "TMT-exp and bulk RNA-seq data showing gene and protein expression changes in macrophages "
    "stimulated with TLR/STING agonists compared to M0 macrophages. TMT-exp and bulk RNA-Seq are "
    "from n = 3 donors each (different donors for each assay)."
)

RC_SHEET_NAME = "S2-3 Reactivity changes"
RC_TITLE = (
    "Cysteine reactivity changes identified from IA-DTB reactivity TMT-ABPP data "
    "(related to Figures 2, 3, and 4, and Extended Data Fig. 3, 4, 5, and 6)."
)
RC_DESCRIPTION = (
    "Proteins with identified cysteine reactivity changes in at least one condition (TLR1-2, "
    "TLR3, TLR4, TLR7, TLR8, STING) compared to M0 macrophages. Data are median from n = 3 donors. "
    "Per-condition columns list one value per residue in the 'residue' column, '|'-separated; "
    "'(C1;C2)' marks a change that could not be attributed to a single residue and counts once. "
    "The cross-reference columns on the right describe only the cysteines that change: complex "
    "membership (CORUM, ComplexPortal), conservation depth (out of 102 organisms), AlphaMissense "
    f"pathogenicity (pathogenic > {AM_PATHOGENIC}, benign < {AM_BENIGN}), and pathogenic ClinVar "
    f"variants within {CLINVAR_RESIDUE_TOLERANCE} residues."
)

RC_LEAD_COLS = ["uniprot", "protein", "description", "rc_n_peptides", "sequence", "residue"]
# Cross-reference block, to the right of the per-condition blocks, in funnel order:
# complex membership -> conservation -> AlphaMissense -> ClinVar. Each annotation gets a
# protein-level summary (filterable) plus a residue-keyed detail string (lossless).
RC_CROSS_REF_COLS = [
    "complex_membership",
    "corum_complexes",
    "complexportal_complexes",
    "max_conservation_depth",
    "conservation_bin",
    "conservation_by_residue",
    "max_am_pathogenicity",
    "am_class",
    "am_pathogenicity_by_residue",
    "clinvar_pathogenic_variants",
    "clinvar_phenotypes",
    "clinvar_by_residue",
]

# IP-MS: the two endogenous FERRY-subunit pulldowns, merged into one sheet on uniprot.
IPMS_SHEET_NAME = "S2-4 Immunoprecipitation MS"
IPMS_FILES = [("ipms_table_Fy2.csv", "Fy2"), ("ipms_table_Tbck.csv", "Tbck")]
# Shared identifier/annotation columns (front) and CRAPome columns (back); the remaining,
# bait-specific columns already carry the bait name (e.g. M0-Fy2_d1_1, p_value_LPS-Tbck ...).
IPMS_ID_COLS = [
    "uniprot",
    "protein",
    "description",
    "uniprot_function",
    "uniprot_goterms",
    "ribosomal",
    "FERRY subunits",
    "RBP",
    "category",
    "GST-FERRY IP protein",
]
IPMS_CRAPOME_COLS = ["Crapome contaminant", "Crapome Num of Expt. (out of 716 total)"]
IPMS_TITLE = (
    "Immunoprecipitation mass spectrometry data "
    "(related to Figure 5 and Extended Data Fig. 7)"
)
IPMS_DESCRIPTION = (
    "Immunoprecipitation mass spectrometry showing proteins co-enriched with TBCK or PPP1R21 in "
    "LPS stimulated macrophages (TLR4) or unstimulated (M0) macrophages. Data are from n = 3 "
    "donors (different donors for each immunoprecipitation assay)."
)

# --- S2-5: GO-term enrichment across modalities -------------------------------------------
# Each per-analysis CSV is a go_enrich() output from the three visualization .Rmd files
# (macrophage_figure_setting.R: clusterProfiler::enricher over MSigDB 2026.1.Hs C5 GO),
# normalized to one flat table of four tag columns + the standard go_enrich() schema. Only
# significant terms are kept, at each figure's own cutoff (all also require n_proteins > 5).
GO_STD_COLS = [
    "ontology", "GO_id", "Term", "Overlap", "n_proteins",
    "p_value", "adjusted_p_value", "neg_log10_fdr", "Genes",
]
GO_TAG_COLS = ["Dataset", "Analysis", "Subset", "Background"]
GO_OUT_COLS = GO_TAG_COLS + GO_STD_COLS

# Per-test background (the clusterProfiler `universe`), spelled out in the table per the caption.
GO_BG_ALL = "All annotated genes (MSigDB 2026.1.Hs C5 GO)"
GO_BG_WP_DETECTED = "Per-comparison detected proteome"

# Source CSVs. The WP and reactivity results live under the absolute Manuscript figure tree;
# the two phospho results are written into the repo's phosphoproteomics/ folder.
WP_GO_CSV = FIGURE1_PANELS / "go_wp_vs_RNA" / "rna_vs_wp_correlation_go_term_analysis_BP.csv"
RC_TLR4_GO_CSV = (
    MANUSCRIPT / "02_Figures/6-figure format_EV"
    / "Figure 2_Reactivity_part1_GTPase regulation/Panels/go_enrichment/go_results.csv"
)
CROSS_OMICS_GO_CSV = (
    MANUSCRIPT / "02_Figures/6-figure format_EV"
    / "Figure 3_Reactivity_part2_endocytosis/Panels/go_term_venn_diagram"
    / "cross_omics_go_enrichment.csv"
)
PHOSPHO_DIR = REPO_ROOT / "phosphoproteomics"
PHOSPHO_BP_GO_CSV = PHOSPHO_DIR / "go_bp_phospho.csv"
PHOSPHO_CC_GO_CSV = PHOSPHO_DIR / "go_cc_reactivity.csv"

GO_SHEET_NAME = "S2-5 GO-term enrichment"
GO_TITLE = (
    "GO-term enrichment analysis of whole proteome, reactivity, and phosphoproteomic data"
)
GO_DESCRIPTION = (
    "Fisher's exact test and the 2026 MSigDB gene sets were used for all analyses, and the "
    "background for each test is specified in the table. Only significant terms are shown."
)

# --- S2-6: GSEA of the WP PCA loadings ----------------------------------------------------
# One gseapy.prerank report per PC, written by wp_downstream_analysis.ipynb
# (perform_gsea_proteomics): that PC's protein loadings ranked and tested against MSigDB
# 2026.1.Hs c5.go.bp, 1000 permutations, BH correction, seed 42. Kept at the same FDR cutoff
# the panel uses (the "GSEA of PC loadings" chunk in wp_visualization.Rmd).
PCA_GSEA_DIR = INPUT_FOLDER / "04_results/pca_plots"
PCA_GSEA_PCS = ["PC1", "PC2"]
PCA_GSEA_FDR_CUTOFF = 0.01
# gseapy's "Name" column is dropped: it is the constant run label "prerank".
PCA_GSEA_COLS = [
    "Term", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val",
    "Tag %", "Gene %", "Lead_genes",
]

GSEA_SHEET_NAME = "S2-6 GSEA of PC loadings"
GSEA_TITLE = (
    "Gene set enrichment analysis of whole proteome PCA loadings "
    "(related to Figure 1 and Extended Data Fig. 1)"
)
GSEA_DESCRIPTION = (
    "Pre-ranked GSEA of the protein loadings on each principal component of the whole proteome "
    "PCA, against the GO biological process gene sets of the 2026 MSigDB. A positive NES marks a "
    "term enriched at the positive end of that principal component. Only significant terms "
    "(FDR q-value < 0.01) are shown."
)

# --- S2-7: expression-normalized phosphoproteomics -----------------------------------------
# Written by phosphoproteomics.ipynb: median-normalized channel ratios divided by that protein's
# whole-proteome fold change (an inner join on uniprot/condition, which is what trims 11,690
# sites to 9,591), then an unpaired t-test of all TLR4 channels vs all M0 channels. This is the
# file rc_visualization.Rmd and the Figure 3 panels read. NOT the sibling `_annotated.csv` (same
# data plus three figure-specific GO flags, written by reactivity_polars.ipynb) and NOT the
# repo's gitignored normalized_phosphorylation_table.csv.
PHOSPHO_TABLE_CSV = (
    PROTEOMICS
    / "05_phosphoproteomics/processed_results"
    / "expression-normalized_phosphorylation_table.csv"
)
# Phospho has no adjusted p-value, so STAT_METRICS (which lists -log10_pval_adj) is not reused.
PHOSPHO_STAT_METRICS = ["log2_FC", "FC", "p_value", "-log10_pval", "Regulation"]
PHOSPHO_LEAD_COLS = [
    "uniprot",
    "protein",
    "residue",
    "sequence",
    "description",
    "uniprot_function",
]
# "<metric>_phospho - TLR4 vs. M0 (N Residues)" -> ("<metric>", "TLR4 vs. M0"). The trailing
# residue count shifts whenever the table is regenerated, so it is parsed off, mirroring
# _STAT_SUFFIX. `index_phospho` (a reset_index artifact) drops out: "index" is not a metric.
_PHOSPHO_STAT_SUFFIX = re.compile(r"^(?P<metric>.+?)_phospho - (?P<comp>.+? vs\. M0)")
# 4 donors (d1, d2, d3, d5) x 3 technical replicates x 2 conditions.
_PHOSPHO_CHANNEL = re.compile(r"^(M0|TLR4)_d\d+_\d+$")

PHOSPHO_SHEET_NAME = "S2-7 Phosphoproteomics"
PHOSPHO_TITLE = (
    "Expression-normalized phosphoproteomics data from M0 and LPS stimulated macrophages "
    "(related to Figure 3 and Extended Data Fig. 4 and 5)"
)
PHOSPHO_DESCRIPTION = (
    "Phosphosite quantification in macrophages stimulated with LPS (TLR4) compared to M0 "
    "macrophages. Per-replicate channel ratios were normalized to protein expression by dividing "
    "each phosphosite's channel ratio by that protein's whole proteome fold change, so changes "
    "reflect phosphorylation rather than protein abundance; phosphosites on proteins without "
    "whole proteome coverage are not reported. Data are from n = 4 donors x 3 technical "
    "replicates, restricted to phosphosites quantified in at least 2 donors. Significance is an "
    "unpaired two-sample t-test of all TLR4 channels against all M0 channels; p-values are raw "
    "and are not corrected for multiple testing. 'Regulation' marks a phosphosite Significant Up "
    "or Significant Down at p < 0.05 and a fold change greater than 2."
)


# Per-comparison statistics kept from volcano_data.csv (helper cols index/name dropped).
STAT_METRICS = ["log2_FC", "p_value", "-log10_pval", "-log10_pval_adj", "Regulation"]

# wp_vs_rnaseq per-condition columns -> supplementary metric names (order = column order).
RNA_METRIC_RENAME = {
    "log2FoldChange": "RNA_log2FoldChange",
    "pvalue": "RNA_pvalue",
    "padj": "RNA_padj",
    "logFC_wp": "WP_log2FoldChange",
    "expression": "regulation_protein_rna",
}

# Functional-category flags are every column after "condition"; detect them rather than
# hardcode, so Data S2 tracks the notebook's category set.
_LONG_FORMAT_LEADING_COLS = "condition"

# "<metric>_Macrophage WP - <COND> vs. M0 (N Proteins)" -> ("<metric>", "<COND> vs M0")
_STAT_SUFFIX = re.compile(r"^(?P<metric>.+?)_Macrophage WP - (?P<comp>.+? vs\. M0)")


def _require(path):
    return require(path, "Run whole_proteome/wp_downstream_analysis.ipynb first to generate it.")


def load_abundance():
    """Per-replicate channel ratios keyed by uniprot (pepNum cols dropped, names cleaned).

    protein/description are dropped here and taken from the stats spine instead: the
    ``description`` strings differ between source tables by trailing whitespace, so uniprot
    is the only reliable join key.
    """
    df = pd.read_csv(_require(ABUNDANCE_CSV), index_col=0)
    channel_cols = [
        c for c in df.columns if c not in ID_COLS and not c.startswith("pepNum")
    ]
    df = df[["uniprot"] + channel_cols].copy()
    # "M0_D5_1_processed_census-out_20250118_TZ3-104" -> "M0_D5_1"
    df = df.rename(columns={c: re.sub(r"_processed.*", "", c) for c in channel_cols})
    return df


def load_stats():
    """Per-comparison differential-expression stats, columns renamed to '<COND> vs M0 <metric>'.

    Comparison order follows the volcano_data.csv column order.
    """
    df = pd.read_csv(_require(VOLCANO_CSV))
    keep, rename, seen = list(ID_COLS), {}, []
    for col in df.columns:
        m = _STAT_SUFFIX.match(col)
        if not m:
            continue
        metric, comp = m.group("metric"), m.group("comp").replace("vs.", "vs")
        if metric not in STAT_METRICS:
            continue
        keep.append(col)
        rename[col] = f"{comp} {metric}"
        if comp not in seen:
            seen.append(comp)
    ordered = list(ID_COLS) + [
        c for comp in seen for c in keep if rename.get(c, "").startswith(f"{comp} ")
    ]
    return df[ordered].rename(columns=rename)


def load_flags():
    """Boolean functional-category flags (columns after 'condition') per uniprot."""
    df = pd.read_csv(_require(ANNOTATION_CSV))
    flag_cols = list(df.columns[df.columns.get_loc(_LONG_FORMAT_LEADING_COLS) + 1 :])
    flags = df.groupby("uniprot", as_index=False)[flag_cols].any()
    return flags


def build_whole_proteome():
    """S2-1: identifiers + per-replicate abundances + per-comparison stats + category flags."""
    stats = load_stats()  # spine: identifiers + the analyzed protein set
    abundance = load_abundance()
    flags = load_flags()

    # Join on uniprot (the only reliable key, see load_abundance); identifiers from the spine.
    return (
        stats[ID_COLS]
        .merge(abundance, on="uniprot", how="left")
        .merge(stats.drop(columns=["protein", "description"]), on="uniprot", how="left")
        .merge(flags, on="uniprot", how="left")
    )


def build_rna_vs_wp():
    """S2-2: per-comparison RNA (DESeq2) and WP fold changes with the protein/RNA quadrant.

    Each ``<COND> vs. M0.csv`` merges DESeq2 RNA results with WP means; the notebook's
    concat of two merge strategies leaves near-identical duplicate uniprot rows, so keep the
    first per uniprot. Identifiers come from the union of all comparisons.
    """
    id_frames, blocks = [], []
    for cond in CONDITIONS:
        df = pd.read_csv(_require(RNA_VS_WP_DIR / f"{cond} vs. M0.csv"))
        df = df.drop_duplicates(subset="uniprot", keep="first")
        id_frames.append(df[ID_COLS])
        block = df[["uniprot", *RNA_METRIC_RENAME]].rename(
            columns={k: f"{cond} vs M0 {v}" for k, v in RNA_METRIC_RENAME.items()}
        )
        blocks.append(block)

    final = pd.concat(id_frames).drop_duplicates(subset="uniprot", keep="first")
    for block in blocks:
        final = final.merge(block, on="uniprot", how="left")
    return final


def _residue_detail(positions, values, fmt="{}"):
    """'C145:102|C273:34' over the residues that have a value, in position order."""
    parts = [
        f"C{pos}:{fmt.format(values[pos])}" for pos in sorted(positions) if pos in values
    ]
    return "|".join(parts)


def build_cross_ref_block(rc_wide):
    """The funnel's evidence per protein, over the cysteines that change in any condition.

    Each annotation gets a protein-level summary and a residue-keyed detail string. max() before
    binning is monotonic, so each summary answers the funnel's question for that filter on its own
    ("is any changing cysteine conserved in 51+ organisms?"). It does not force one *single*
    cysteine to clear every filter the way the funnel's residue-level cascade does -- on the
    current data the two agree, but the residue columns are what settle it.

    Counts here run higher than the funnel's: it annotates behind a .drop_nulls() that discards a
    residue missing any of conservation/pPSE/AlphaMissense, and its stored Complex call is
    CORUM-only, whereas ``complex_membership`` follows cell 38's code and ORs in ComplexPortal.
    """
    changed = rc_change_positions(rc_wide)
    cache = load_cys_cross_ref()
    conservation, pathogenicity, clinvar_n, clinvar_pheno = {}, {}, {}, {}
    for row in cache.iter_rows(named=True):
        key = (row["uniprot"], row["pos"])
        if row["conservation_depth"] is not None:
            conservation[key] = row["conservation_depth"]
        if row["am_pathogenicity"] is not None:
            pathogenicity[key] = row["am_pathogenicity"]
        if row["clinvar_pathogenic_variants"]:
            clinvar_n[key] = row["clinvar_pathogenic_variants"]
            clinvar_pheno[key] = row["clinvar_phenotypes"]

    rows = []
    for uniprot in rc_wide["uniprot"]:
        positions = changed.get(uniprot, set())
        cons = {p: conservation[(uniprot, p)] for p in positions if (uniprot, p) in conservation}
        am = {p: pathogenicity[(uniprot, p)] for p in positions if (uniprot, p) in pathogenicity}
        cv = {p: clinvar_n[(uniprot, p)] for p in positions if (uniprot, p) in clinvar_n}
        phenotypes = sorted(
            {clinvar_pheno[(uniprot, p)] for p in positions if (uniprot, p) in clinvar_pheno}
        )
        max_cons = max(cons.values(), default=None)
        max_am = max(am.values(), default=None)
        rows.append(
            {
                "uniprot": uniprot,
                "max_conservation_depth": max_cons,
                "conservation_bin": conservation_bin(max_cons),
                "conservation_by_residue": _residue_detail(positions, cons, "{:.0f}"),
                "max_am_pathogenicity": max_am,
                "am_class": am_class(max_am),
                "am_pathogenicity_by_residue": _residue_detail(positions, am, "{:.4f}"),
                "clinvar_pathogenic_variants": sum(cv.values()),
                "clinvar_phenotypes": "|".join(phenotypes),
                "clinvar_by_residue": _residue_detail(positions, cv),
            }
        )
    block = pl.DataFrame(rows)
    complexes = load_complex_membership(rc_wide["uniprot"]).select(
        "uniprot", "complex_membership", "corum_complexes", "complexportal_complexes"
    )
    return block.join(complexes, on="uniprot", how="left")


def build_reactivity():
    """S2-3: per-stimulus IA-DTB reactivity blocks + the funnel's cross-reference evidence.

    One row per protein with a reactivity change in at least one condition (as the sheet's
    description says), from the current run's long-format table via src.reactivity.pivot_rc_wide.
    """
    rc_wide = pivot_rc_wide(load_rc_long()).filter(has_reactivity_change())
    out = rc_wide.join(build_cross_ref_block(rc_wide), on="uniprot", how="left")
    cols = list(RC_LEAD_COLS)
    for cond in RC_CONDITIONS:
        cols += [f"{cond}_{metric}" for metric in CONDITION_BLOCK]
    cols += RC_CROSS_REF_COLS
    missing = [c for c in cols if c not in out.columns]
    if missing:
        sys.exit(f"ERROR: reactivity table missing expected columns: {missing}")
    return out.select(cols).to_pandas()


def build_ipms():
    """S2-4: the two FERRY-subunit IP-MS tables merged on uniprot into one sheet.

    Shared identifier/annotation columns lead and the CRAPome columns trail; each bait's own
    signal and stat columns (which already carry the bait name) sit between them, Fy2 then Tbck.
    Rows are the union of proteins across the two pulldowns.
    """
    shared = IPMS_ID_COLS + IPMS_CRAPOME_COLS
    spine_parts, bait_blocks = [], []
    for csv_name, bait in IPMS_FILES:
        df = pd.read_csv(_require(IPMS_DIR / csv_name))
        spine_parts.append(df[shared])
        bait_cols = [c for c in df.columns if c not in shared]
        # The bait-signal and stat columns already carry the bait name; the IgG-control signal
        # columns (e.g. M0-IgG_d1_1) are identically named across both pulldowns but are distinct
        # experiments, so tag them with the bait to keep them separate after the merge.
        block = df[["uniprot"] + bait_cols].rename(
            columns={c: c if bait in c else f"{c}_{bait}" for c in bait_cols}
        )
        bait_blocks.append(block)

    spine = pd.concat(spine_parts).drop_duplicates(subset="uniprot", keep="first")
    out = spine[IPMS_ID_COLS]
    for block in bait_blocks:
        out = out.merge(block, on="uniprot", how="outer")
    return out.merge(spine[["uniprot"] + IPMS_CRAPOME_COLS], on="uniprot", how="left")


def _tag_go(df, dataset, analysis, background, subset=""):
    """Add the four tag columns and return the unified S2-5 column order.

    ``subset`` is a per-row facet (Series aligned on ``df``'s index) or a scalar label; any
    standard column the source lacks (e.g. Overlap for the cross-omics frame) is filled with NA.
    """
    out = df.copy()
    out["Dataset"] = dataset
    out["Analysis"] = analysis
    out["Subset"] = subset
    out["Background"] = background
    for col in GO_STD_COLS:
        if col not in out.columns:
            out[col] = pd.NA
    return out[GO_OUT_COLS]


def load_wp_go():
    """S2-5: WP RNA-vs-WP correlation-quadrant GO:BP (per-comparison detected-proteome background).

    Source is already BP-only and quadrant-labeled (condition/direction/color); the figure cutoff
    is re-applied defensively.
    """
    df = pd.read_csv(require(WP_GO_CSV, "Re-run the go_enrich chunk in wp_visualization.Rmd."))
    df = df[(df["adjusted_p_value"] < 0.05) & (df["n_proteins"] > 5) & (df["ontology"] == "BP")]
    subset = (
        df["condition"].astype(str) + " " + df["direction"].astype(str)
        + " (" + df["color"].astype(str) + ")"
    )
    return _tag_go(df, "Whole proteome", "RNA-vs-WP correlation quadrants (GO:BP)",
                   GO_BG_WP_DETECTED, subset)


def load_reactivity_tlr4_go():
    """S2-5: TLR4 cysteine-reactive proteins, GO:BP over-representation (all-annotated background)."""
    path = require(RC_TLR4_GO_CSV, "Re-run reactivity/rc_visualization.Rmd (go_enrich chunk).")
    df = pd.read_csv(path)
    if "adjusted_p_value" not in df.columns:
        sys.exit(
            f"ERROR: {path} is the old enrichGO schema (pre go_enrich() refactor).\n"
            "Re-run the go_enrich() chunk in reactivity/rc_visualization.Rmd to refresh it."
        )
    df = df[(df["adjusted_p_value"] < 0.05) & (df["n_proteins"] > 5)]
    return _tag_go(df, "Reactivity", "TLR4-reactive proteins (GO:BP)", GO_BG_ALL)


def load_cross_omics_go():
    """S2-5: cross-omics enrichment of the TLR4 shared/omic-specific changes.

    The source is one row per (GO term, omic) and is NOT pre-filtered (the figure shows the same
    terms across omics for contrast), so the standard cutoff is applied here to keep only the
    significant cells. Schema differs: ``ont`` -> ontology, ``omic`` -> Subset; Overlap is absent
    (filled NA) and neg_log10_fdr is derived from the BH p-value.
    """
    path = require(CROSS_OMICS_GO_CSV, "Re-run reactivity/rc_visualization.Rmd (cross-omics chunk).")
    df = pd.read_csv(path).rename(columns={"ont": "ontology"})
    df = df[(df["adjusted_p_value"] < 0.05) & (df["n_proteins"] > 5)]
    df["neg_log10_fdr"] = -np.log10(df["adjusted_p_value"])
    return _tag_go(df, "Cross-omics", "Cross-omics significant changes", GO_BG_ALL, df["omic"])


def load_phospho_bp_go():
    """S2-5: up/down-phosphorylated proteins, GO:BP over-representation (all-annotated background)."""
    path = require(PHOSPHO_BP_GO_CSV, "Re-run phosphoproteomics/visualization.Rmd (go_enrich chunk).")
    df = pd.read_csv(path)
    df = df[(df["adjusted_p_value"] < 0.01) & (df["n_proteins"] > 5)]
    return _tag_go(df, "Phosphoproteomics", "Phosphosite up/down (GO:BP)", GO_BG_ALL,
                   df["direction"])


def load_phospho_cc_go():
    """S2-5: phospho-substrate localization, GO:CC over-representation (all-annotated background)."""
    path = require(PHOSPHO_CC_GO_CSV, "Re-run phosphoproteomics/visualization.Rmd (go_cc chunk).")
    df = pd.read_csv(path)
    df = df[(df["adjusted_p_value"] < 0.2) & (df["n_proteins"] > 5)]
    return _tag_go(df, "Phosphoproteomics", "Phosphosite localization (GO:CC)", GO_BG_ALL)


def build_go_enrichment():
    """S2-5: significant GO-term enrichment across WP, reactivity, and phospho, one flat table."""
    blocks = [
        load_wp_go(),
        load_reactivity_tlr4_go(),
        load_cross_omics_go(),
        load_phospho_bp_go(),
        load_phospho_cc_go(),
    ]
    return pd.concat(blocks, ignore_index=True)


def build_pc_gsea():
    """S2-6: significant pre-ranked GSEA hits for each PCA loading vector, one flat table.

    One gseapy report per PC, concatenated behind a leading ``PC`` tag column. Within a PC rows
    run from the most positive to the most negative NES, so they read from the +PC end of the
    axis to the -PC end, matching the panel's orientation.
    """
    blocks = []
    for pc in PCA_GSEA_PCS:
        df = pd.read_csv(_require(PCA_GSEA_DIR / pc / "gseapy.gene_set.prerank.report.csv"))
        df = df[df["FDR q-val"] < PCA_GSEA_FDR_CUTOFF].sort_values("NES", ascending=False)
        df["PC"] = pc
        blocks.append(df[["PC"] + PCA_GSEA_COLS])
    return pd.concat(blocks, ignore_index=True)


def build_phospho():
    """S2-7: per phosphosite identifiers + UniProt function + TLR4-vs-M0 stats + channels."""
    path = require(
        PHOSPHO_TABLE_CSV,
        "Run phosphoproteomics/phosphoproteomics.ipynb first to generate it "
        "(the expression-normalization and volcano cells).",
    )
    df = pd.read_csv(path)

    rename = {}
    for col in df.columns:
        m = _PHOSPHO_STAT_SUFFIX.match(col)
        if m and m.group("metric") in PHOSPHO_STAT_METRICS:
            rename[col] = f'{m.group("comp").replace("vs.", "vs")} {m.group("metric")}'
    # Metric order follows PHOSPHO_STAT_METRICS, not the source column order.
    stats = [
        c
        for metric in PHOSPHO_STAT_METRICS
        for c in rename
        if rename[c].endswith(f" {metric}")
    ]
    channels = sorted(c for c in df.columns if _PHOSPHO_CHANNEL.match(c))

    # ID (redundant with protein+residue, and stale with respect to the sorted residue column)
    # and identifier (redundant with uniprot+protein+description+residue+sequence) are dropped.
    cols = PHOSPHO_LEAD_COLS + stats + channels
    missing = [c for c in cols if c not in df.columns]
    if missing:
        sys.exit(f"ERROR: phospho table missing expected columns: {missing}")
    return df[cols].rename(columns=rename)


def main():
    sheets = [
        SuppSheet(1, WP_SHEET_NAME, WP_TITLE, WP_DESCRIPTION, build_whole_proteome()),
        SuppSheet(2, RNA_SHEET_NAME, RNA_TITLE, RNA_DESCRIPTION, build_rna_vs_wp()),
        SuppSheet(3, RC_SHEET_NAME, RC_TITLE, RC_DESCRIPTION, build_reactivity()),
        SuppSheet(4, IPMS_SHEET_NAME, IPMS_TITLE, IPMS_DESCRIPTION, build_ipms()),
        SuppSheet(5, GO_SHEET_NAME, GO_TITLE, GO_DESCRIPTION, build_go_enrichment()),
        SuppSheet(6, GSEA_SHEET_NAME, GSEA_TITLE, GSEA_DESCRIPTION, build_pc_gsea()),
        SuppSheet(
            7, PHOSPHO_SHEET_NAME, PHOSPHO_TITLE, PHOSPHO_DESCRIPTION, build_phospho()
        ),
    ]
    write_supplementary_workbook(OUTPUT_XLSX, sheets)
    print(f"Wrote {OUTPUT_XLSX}")
    for s in sheets:
        print(f"  '{s.sheet_name}': {s.df.shape[0]} rows x {s.df.shape[1]} columns")


if __name__ == "__main__":
    main()
