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

# The annotation CSV's super-pathway header is missing its closing paren; the sheet ships the
# closed form.
SUPER_PATHWAY_SRC = "Super Pathway (MSK metabolomics Core dMRM method"
SUPER_PATHWAY = f"{SUPER_PATHWAY_SRC})"

# Annotation values are padded with trailing spaces and non-breaking spaces, which splits a
# class into two look-alike strings for anyone filtering or pivoting the shipped sheet.
PAD_CHARS = " \t\xa0"

# Uncurated spellings of super-pathways that are otherwise consistent.
# metabolomics_visualization.Rmd folds exactly these before plotting; blank stays blank, since
# the Rmd's "Other" bucket is a plot-legend convenience rather than an annotation.
SUPER_PATHWAY_FIXES = {
    "Amino Acid?": "Amino Acid",
    "Nucleotide?": "Nucleotide",
    "Lipid?": "Lipid",
    "Lipids": "Lipid",
}

# Annotation columns (order) from the annotation CSV, joined Compound == NAME (NAME dropped).
ANNOTATION_COLS = [
    "HMDB",
    "KEGG",
    "Alias",
    "Pathway",
    "Metabolite Class (Vardhana lab)",
    SUPER_PATHWAY,
    "Chemical Taxonomy, Super Class (HMDB)",
    "Chemical Taxonomy, Sub Class (HMDB)",
]
# Lead identifiers: Compound + 4 database IDs up front (rest of annotation follows).
LEAD_COLS = [ID_COL, "HMDB", "KEGG", "Alias"]

# Spelling variants, as data name -> annotation NAME. Deoxyguanosine, dUMP and Oxalacetic acid
# are absent from the panel altogether and keep null annotation fields; `build_metabolomics`
# prints whatever is left unannotated so a changed panel doesn't pass silently.
COMPOUND_ALIASES = {
    "Indole-3-propionic acid": "3-Indolepropionic acid",
    "myo-Inositol": "Myoinositol",
}

# Per-sample value blocks (channel-ratio + raw signal), tagged with these suffixes.
NORMALIZED_SUFFIX = "_cell-volume-normalized-QRILC-imputation"
RAW_SUFFIX = "_raw-signal-intensity"

# The single comparison and its volcano metrics (the "index_*" helper column is dropped).
# The stimulus is named "TLR4", matching the stat columns of Data S1 and S2; the per-sample
# columns keep the LPS_* names the assay used.
COMPARISON = "TLR4 vs M0"
STAT_METRICS = ["log2_FC", "p_value", "-log10_pval", "-log10_pval_adj", "Regulation"]

SHEET_NAME = "S3-1 Polar metabolites"
TITLE = (
    "Abundance of polar metabolites in M0 and TLR4 stimulated macrophages "
    "(related to Figure 4 and Extended Data Fig. 3)"
)
DESCRIPTION = (
    "Polar metabolite abundance values in macrophages stimulated with the TLR4 agonist LPS "
    "compared to M0 macrophages as determined by LC-MS/MS analysis. The channel ratio values "
    "after quantile regression imputation of left censored data (QRILC) are shown alongside raw "
    "signal intensity values and differential abundance calculations. Data are from n = 4 donors."
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
    """Metabolite annotations keyed by Compound (annotation ``NAME`` -> ``Compound``).

    ``COMPOUND_ALIASES`` is applied to the annotation's names, not the spine's, so the sheet
    reports each metabolite under the name the analysis used.
    """
    df = pd.read_csv(require(ANNOTATION_CSV))
    df = df.rename(columns={"NAME": ID_COL, SUPER_PATHWAY_SRC: SUPER_PATHWAY})
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].str.strip(PAD_CHARS)
    df[SUPER_PATHWAY] = df[SUPER_PATHWAY].replace(SUPER_PATHWAY_FIXES)

    # A mapping that matched nothing would quietly reintroduce the blank rows it exists to fill.
    names = set(df[ID_COL])
    for data_name, annotation_name in COMPOUND_ALIASES.items():
        if annotation_name not in names:
            sys.exit(f"ERROR: annotation table has no '{annotation_name}' to alias")
        if data_name in names:
            sys.exit(f"ERROR: aliasing '{annotation_name}' would collide with '{data_name}'")
    df[ID_COL] = df[ID_COL].replace({v: k for k, v in COMPOUND_ALIASES.items()})

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

    unannotated = merged.loc[merged[ANNOTATION_COLS].isna().all(axis=1), ID_COL].tolist()
    if unannotated:
        print(f"  no annotation for {len(unannotated)}: {', '.join(unannotated)}")

    return merged[column_order]


def main():
    sheets = [SuppSheet(1, SHEET_NAME, TITLE, DESCRIPTION, build_metabolomics())]
    write_supplementary_workbook(OUTPUT_XLSX, sheets, freeze_cols=len(LEAD_COLS))
    print(f"Wrote {OUTPUT_XLSX}")
    for s in sheets:
        print(f"  '{s.sheet_name}': {s.df.shape[0]} rows x {s.df.shape[1]} columns")


if __name__ == "__main__":
    main()
