"""Assemble Data S1 (bulk RNA-Seq): a contents index + one flat sheet, one row per
protein-coding gene on the primary assembly (identifiers + per-stimulus-vs-M0 DESeq2 stats +
per-sample normalized counts).

Reads the upstream DESeq2 tables (``DE_<COND>_rd.xlsx``, sheet ``DEseq2, coding``; Beatrice
Zhang) read-only from the Figure 1 panels folder, with the same TLR1->TLR1-2 relabeling as
``whole_proteome/multi_modal_significance.py``. ``DE_M0high_rd.xlsx`` is excluded (no
counterpart in the other modalities), and with it the ``M0 High`` samples.

    conda run -n polars python rna/build_data_s1.py
"""

import sys
from pathlib import Path

import pandas as pd

# make `import src.*` work when run from anywhere in the repo
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
from src.macrophage import CONDITIONS, load_deg, require  # noqa: E402
from src.supp_data import SuppSheet, write_supplementary_workbook  # noqa: E402

OUTPUT_XLSX = REPO_ROOT / "Data S1.xlsx"

# The count matrices GEO receives, written by src/data_deposition/build_processed_files.py.
NORMALIZED_COUNTS_TSV = REPO_ROOT / "rna/normalized_counts.tsv"
SAMPLE_MANIFEST = REPO_ROOT / "src/data_deposition/sample_manifest.tsv"

# ensembl.gene.id is unique per gene (unlike gene.symbol), so it is the join key; the other
# identifiers ride along from the first file a gene appears in.
JOIN_KEY = "ensembl.gene.id"
# The versioned GENCODE id, which is what the count matrices key on.
COUNTS_KEY = "gene.id"
LEAD_COLS = ["gene.symbol", "uniprot.id", "ensembl.gene.id", "gene.id", "chr", "description"]
# DESeq2 statistics kept per comparison (native column names), renamed "<COND> vs M0 <metric>".
STAT_METRICS = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]

# GENCODE annotates each pseudoautosomal locus twice, on X and as a `_PAR_Y` copy, and exactly
# one copy clears DESeq2's filtering per comparison - an arbitrary choice that varies by
# comparison, so 15 genes arrive under two `gene.id` values. The copies are interchangeable
# (each carries the whole locus, so they must not be summed); the X one is kept.
PAR_Y_SUFFIX = "_PAR_Y"

# Sequenced and deposited, but not a comparison here (see DE_M0high_rd.xlsx above).
EXCLUDED_STIMULUS = "M0 High"

# Genes on alt/fix/random contigs are dropped: they duplicate primary-assembly genes under a
# second ensembl.gene.id, each quantified in only some comparisons. No gene is annotated on
# both a canonical and a non-canonical contig, so this drops whole genes, never half of one.
PRIMARY_CONTIGS = frozenset([*(str(i) for i in range(1, 23)), "X", "Y", "M", "MT"])

SHEET_NAME = "S1-1 Bulk RNA Sequencing Data"
TITLE = (
    "Differential gene expression analysis of bulk RNA-Seq data from M0 and TLR/STING stimulated "
    "macrophages (related to Figure 1 and Extended Data Fig. 1)."
)
DESCRIPTION = (
    "Differential gene expression analysis was performed with DESeq2, using its median-of-ratios "
    "method to normalize RNA counts. Genes with zero mapped reads were excluded from display. "
    "Each stimulus was compared against the unstimulated low-vehicle control, written M0 in the "
    "statistics columns and M0 Low in the count columns.  Normalized counts for each stimulus and donor are shown to the right of the "
    "differential expression statistics. Genes are the union across the seven comparisons. Genes are restricted to the primary assembly. Data are from n = 3 donors."
)


def load_condition(condition):
    """Read one condition's DESeq2 coding sheet; return (identifier frame, renamed stat block)."""
    df = load_deg(condition)
    ids = df[LEAD_COLS]
    block = df[[JOIN_KEY, *STAT_METRICS]].rename(
        columns={m: f"{condition} vs M0 {m}" for m in STAT_METRICS}
    )
    return ids, block


def load_sample_labels():
    """``column_name -> "<stimulus>_<donor>"``, in manifest stimulus order then donor."""
    manifest = pd.read_csv(require(SAMPLE_MANIFEST), sep="\t")
    manifest = manifest[manifest["src_label"] != EXCLUDED_STIMULUS]
    manifest = manifest.sort_values(["cond_num", "donor"], kind="stable")
    labels = dict(
        zip(manifest["column_name"], manifest["src_label"] + "_" + manifest["donor"])
    )
    # Two samples collapsing onto one name would silently drop a column.
    assert len(set(labels.values())) == len(labels), "sample labels are not unique"
    return labels


def load_normalized_counts():
    """DESeq2-normalized counts keyed by ``gene.id``, columns relabeled and ordered."""
    labels = load_sample_labels()
    df = pd.read_csv(
        require(NORMALIZED_COUNTS_TSV), sep="\t", usecols=[COUNTS_KEY, *labels]
    )
    return df[[COUNTS_KEY, *labels]].rename(columns=labels)


def resolve_par_duplicates(ids):
    """One row per ``ensembl.gene.id``, preferring the non-``_PAR_Y`` copy. See PAR_Y_SUFFIX."""
    ids = ids.drop_duplicates()
    duplicated = ids.loc[ids[JOIN_KEY].duplicated(keep=False)]
    if len(duplicated):
        symbols = sorted(set(duplicated["gene.symbol"]))
        print(f"  PAR duplicates resolved to the X copy for {len(symbols)}: {', '.join(symbols)}")

    # Sorting the `_PAR_Y` rows last makes `keep="first"` pick the X copy; the positional index
    # survives, so sort_index restores first-appearance order.
    ranked = ids.assign(_par=ids[COUNTS_KEY].str.endswith(PAR_Y_SUFFIX))
    return (
        ranked.sort_values("_par", kind="stable")
        .drop_duplicates(subset=JOIN_KEY, keep="first")
        .sort_index()
        .drop(columns="_par")
    )


def drop_alt_contigs(spine):
    """Keep only genes on the primary assembly. See PRIMARY_CONTIGS."""
    primary = spine["chr"].astype(str).isin(PRIMARY_CONTIGS)
    dropped = spine.loc[~primary, "chr"].astype(str).nunique()
    print(
        f"  dropped {int((~primary).sum())} genes on {dropped} alt/fix/random contigs; "
        f"{int(primary.sum())} remain"
    )
    return spine[primary]


def build_rna_seq():
    """S1-1: gene identifiers + per-comparison DESeq2 stat blocks + per-sample counts.

    Rows are the union of genes across the seven comparisons (left merge on ``ensembl.gene.id``);
    identifiers come from the first comparison a gene appears in. Counts are merged last, so the
    block lands to the right of the statistics.
    """
    id_frames, blocks = [], []
    for condition in CONDITIONS:
        ids, block = load_condition(condition)
        id_frames.append(ids)
        blocks.append(block)

    spine = resolve_par_duplicates(pd.concat(id_frames, ignore_index=True))
    spine = drop_alt_contigs(spine)
    final = spine[LEAD_COLS]
    for block in blocks:
        final = final.merge(block, on=JOIN_KEY, how="left")

    counts = load_normalized_counts()
    # A left join on a stale key writes a block of blanks rather than failing, so assert first.
    missing = sorted(set(final[COUNTS_KEY]) - set(counts[COUNTS_KEY]))
    assert not missing, (
        f"{len(missing)} gene.id absent from {NORMALIZED_COUNTS_TSV.name}: {missing[:5]}"
    )
    return final.merge(counts, on=COUNTS_KEY, how="left")


def main():
    sheets = [SuppSheet(1, SHEET_NAME, TITLE, DESCRIPTION, build_rna_seq())]
    write_supplementary_workbook(OUTPUT_XLSX, sheets, freeze_cols=len(LEAD_COLS))
    print(f"Wrote {OUTPUT_XLSX}")
    for s in sheets:
        print(f"  '{s.sheet_name}': {s.df.shape[0]} rows x {s.df.shape[1]} columns")


if __name__ == "__main__":
    main()
