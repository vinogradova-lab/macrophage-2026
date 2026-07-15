"""Assemble Data S1 (bulk RNA-Seq): a contents index + one flat sheet, one row per
protein-coding gene (identifiers + per-stimulus-vs-M0 DESeq2 stats).

Reads the upstream DESeq2 tables (``DE_<COND>_rd.xlsx``, sheet ``DEseq2, coding``; Beatrice
Zhang) read-only from the Figure 1 panels folder, with the same TLR1->TLR1-2 relabeling as
``whole_proteome/multi_modal_significance.py``. ``DE_M0high_rd.xlsx`` is excluded (no
counterpart in the other modalities).

    conda run -n polars python rna_analysis/build_data_s1.py
"""

import sys
from pathlib import Path

import pandas as pd

# make `import src.*` work when run from anywhere in the repo
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
from src.macrophage import CONDITIONS, load_deg  # noqa: E402
from src.supp_data import SuppSheet, write_supplementary_workbook  # noqa: E402

OUTPUT_XLSX = REPO_ROOT / "Data S1.xlsx"

# ensembl.gene.id is unique per gene (unlike gene.symbol), so it is the join key; the other
# identifiers ride along from the first file a gene appears in.
JOIN_KEY = "ensembl.gene.id"
LEAD_COLS = ["gene.symbol", "uniprot.id", "ensembl.gene.id", "description"]
# DESeq2 statistics kept per comparison (native column names), renamed "<COND> vs M0 <metric>".
STAT_METRICS = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]

SHEET_NAME = "S1-1 Bulk RNA Sequencing Data"
TITLE = (
    "Differential gene expression analysis of bulk RNA-Seq data from M0 and TLR/STING stimulated "
    "macrophages (related to Figures 1 and Extended Data Fig. 1)."
)
DESCRIPTION = (
    "Differential gene expression analysis was performed with DESeq2, using its median-of-ratios "
    "method to normalize RNA counts. Genes with zero mapped reads were excluded from display. "
    "Data are from n = 3 donors."
)


def load_condition(condition):
    """Read one condition's DESeq2 coding sheet; return (identifier frame, renamed stat block)."""
    df = load_deg(condition)
    ids = df[LEAD_COLS]
    block = df[[JOIN_KEY, *STAT_METRICS]].rename(
        columns={m: f"{condition} vs M0 {m}" for m in STAT_METRICS}
    )
    return ids, block


def build_rna_seq():
    """S1-1: gene identifiers + per-comparison DESeq2 stat blocks, one row per gene.

    Rows are the union of genes across the seven comparisons (left merge on ``ensembl.gene.id``);
    identifiers come from the first comparison a gene appears in.
    """
    id_frames, blocks = [], []
    for condition in CONDITIONS:
        ids, block = load_condition(condition)
        id_frames.append(ids)
        blocks.append(block)

    spine = pd.concat(id_frames).drop_duplicates(subset=JOIN_KEY, keep="first")
    final = spine[LEAD_COLS]
    for block in blocks:
        final = final.merge(block, on=JOIN_KEY, how="left")
    return final


def main():
    sheets = [SuppSheet(1, SHEET_NAME, TITLE, DESCRIPTION, build_rna_seq())]
    write_supplementary_workbook(OUTPUT_XLSX, sheets, freeze_cols=len(LEAD_COLS))
    print(f"Wrote {OUTPUT_XLSX}")
    for s in sheets:
        print(f"  '{s.sheet_name}': {s.df.shape[0]} rows x {s.df.shape[1]} columns")


if __name__ == "__main__":
    main()
