"""Assemble Data S3 (whole-cell metabolomics): a contents index + one flat sheet, one row per
metabolite (name/HMDB/KEGG/class annotations, per-sample raw and cell-volume-normalized
QRILC-imputed intensities, and M0-vs-LPS stats).

Targeted polar dMRM assay: one comparison (M0 vs LPS), 179 metabolites, 4 donors. Run after
metabolomics_analysis.ipynb has generated the combined-files and volcano tables.

    conda run -n polars python metabolomics/build_data_s3.py
"""

import sys
from pathlib import Path

import pandas as pd

# make `import src.*` work when run from anywhere in the repo
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
from src.macrophage import require  # noqa: E402
from src.supp_data import SuppSheet, write_supplementary_workbook  # noqa: E402

# Same data root the notebook uses (data_dir in metabolomics_analysis.ipynb).
METABO_DIR = Path("/Users/henrysanford/dev/test_data/macrophage/metabolomics")
CROSS = METABO_DIR / "cross_replicate_analysis"
CHANNEL_RATIO_CSV = CROSS / "03_combined_files/combreplicates_channel_ratio.csv"
RAW_SIGNAL_CSV = CROSS / "03_combined_files/combreplicates_raw_signal.csv"
VOLCANO_CSV = CROSS / "04_results/volcano_plots/volcano_data.csv"
ANNOTATION_CSV = (
    METABO_DIR
    / "20240924_list of polar metabolites for the targeted assay_annotations_LSP_v3.csv"
)

OUTPUT_XLSX = REPO_ROOT / "Data S3.xlsx"
ID_COL = "Compound"

# Annotation columns (order) from the annotation CSV, joined Compound == NAME (NAME dropped).
ANNOTATION_COLS = [
    "HMDB",
    "KEGG",
    "Alias",
    "Pathway",
    "Metabolite Class (Vardhana lab)",
    "Super Pathway (MSK metabolomics Core dMRM method",
    "Chemical Taxonomy, Super Class (HMDB)",
    "Chemical Taxonomy, Sub Class (HMDB)",
]
# Lead identifiers: Compound + 4 database IDs up front (rest of annotation follows).
LEAD_COLS = [ID_COL, "HMDB", "KEGG", "Alias"]

# Per-sample value blocks (channel-ratio + raw signal), tagged with these suffixes.
NORMALIZED_SUFFIX = "_cell-volume-normalized-QRILC-imputation"
RAW_SUFFIX = "_raw-signal-intensity"

# The single comparison and its volcano metrics (the "index_*" helper column is dropped).
COMPARISON = "LPS vs M0"
STAT_METRICS = ["log2_FC", "p_value", "-log10_pval", "-log10_pval_adj", "Regulation"]

SHEET_NAME = "S3-1 Polar metabolites"
TITLE = (
    "Abundance of polar metabolites in M0 and TLR4 stimulated macrophages "
    "(related to Figure 4 and Extended Data Fig. 3)"
)
DESCRIPTION = (
    "Polar metabolite abundance values in macrophages stimulated with TLR/STING agonists compared "
    "to M0 macrophages as determined by LC-MS/MS analysis. The channel ratio values after quantile "
    "regression imputation of left censored data (QRILC) are shown alongside raw signal intensity "
    "values and differential expression calculations. Data are from n = 4 donors."
)


def load_values(csv, suffix):
    """Per-sample intensities keyed by Compound, sample columns tagged with ``suffix``."""
    df = pd.read_csv(require(csv))
    sample_cols = [c for c in df.columns if c != ID_COL]
    return df[[ID_COL, *sample_cols]].rename(columns={c: f"{c}{suffix}" for c in sample_cols})


def load_stats():
    """M0-vs-LPS volcano metrics keyed by Compound, renamed to '<COMPARISON> <metric>'.

    Volcano columns are '<metric>_metabolomics - M0 vs. LPS (N Metabolites)'; match each by its
    '<metric>_metabolomics' prefix.
    """
    df = pd.read_csv(require(VOLCANO_CSV))
    keep, rename = [ID_COL], {}
    for metric in STAT_METRICS:
        prefix = f"{metric}_metabolomics"
        matches = [c for c in df.columns if c.startswith(prefix)]
        if not matches:
            sys.exit(f"ERROR: volcano table missing a '{prefix}' column")
        col = matches[0]
        keep.append(col)
        rename[col] = f"{COMPARISON} {metric}"
    return df[keep].rename(columns=rename)


def load_annotation():
    """Metabolite annotations keyed by Compound (annotation ``NAME`` -> ``Compound``)."""
    df = pd.read_csv(require(ANNOTATION_CSV))
    df = df.rename(columns={"NAME": ID_COL})
    return df[[ID_COL, *ANNOTATION_COLS]].drop_duplicates(subset=ID_COL, keep="first")


def build_metabolomics():
    """S3-1: annotations + per-sample raw/normalized intensities + M0-vs-LPS stats.

    Spine is the channel-ratio Compound set; annotation is a left join (unmatched metabolites
    keep null annotation fields, as in the notebook).
    """
    normalized = load_values(CHANNEL_RATIO_CSV, NORMALIZED_SUFFIX)
    raw = load_values(RAW_SIGNAL_CSV, RAW_SUFFIX)
    stats = load_stats()
    annotation = load_annotation()

    merged = (
        normalized[[ID_COL]]
        .merge(annotation, on=ID_COL, how="left")
        .merge(normalized, on=ID_COL, how="left")
        .merge(raw, on=ID_COL, how="left")
        .merge(stats, on=ID_COL, how="left")
    )

    # Order: lead identifiers, remaining annotation, normalized block, raw block, stats.
    normalized_cols = [c for c in normalized.columns if c != ID_COL]
    raw_cols = [c for c in raw.columns if c != ID_COL]
    stat_cols = [f"{COMPARISON} {m}" for m in STAT_METRICS]
    trailing_annotation = [c for c in ANNOTATION_COLS if c not in LEAD_COLS]
    column_order = LEAD_COLS + trailing_annotation + normalized_cols + raw_cols + stat_cols
    return merged[column_order]


def main():
    sheets = [SuppSheet(1, SHEET_NAME, TITLE, DESCRIPTION, build_metabolomics())]
    write_supplementary_workbook(OUTPUT_XLSX, sheets, freeze_cols=len(LEAD_COLS))
    print(f"Wrote {OUTPUT_XLSX}")
    for s in sheets:
        print(f"  '{s.sheet_name}': {s.df.shape[0]} rows x {s.df.shape[1]} columns")


if __name__ == "__main__":
    main()
