"""Stacked bar comparing FeNTA vs TiO2 phospho-peptide recovery.

Reads the tidy ``peptide_results_clean.csv`` written by ``aggregate_phospho.py``
and draws one bar per enrichment (FeNTA, TiO2), where each bar counts the peptide
backbones (uniprot + bare sequence, phospho-site markers stripped) that enrichment
detected, colored by the number of phospho-sites per peptide (1, 2, 3+). A backbone
seen in both enrichments contributes to both bars, so this is a per-enrichment
method comparison.

Styling mirrors the lab's guide-RNA figure: legend on the right, despined axes,
with outlines on the stacked boxes.

Conda environment (the lab's ``polars`` env already has these):

    conda activate polars   # needs polars, matplotlib

Usage:

    python num_mods_stacked_bar.py --input-dir /path/to/phospho_dir
    python num_mods_stacked_bar.py --peptides /path/to/peptide_results_clean.csv
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import polars as pl

PEPTIDES_FILENAME = "peptide_results_clean.csv"
# x-axis groups (enrichment), in plotting order, mapped to the "enr" key
GROUPS = [("FeNTA", "FeNTA"), ("TiO2", "TiO2")]
EDGE = "#3a3a3a"  # outline color on the stacked boxes

# color = number of phospho-sites per peptide, capped at 3+.
# stack order is bottom -> top, so "1" ends up on top (lightest).
CATEGORY_STACK = ["3+", "2", "1"]
CATEGORY_COLORS = {"1": "#dce6f0", "2": "#9dc0e0", "3+": "#4C72B0"}  # light -> dark, sequential
# legend order (top -> bottom) and labels
CATEGORY_LEGEND = [("1", "1"), ("2", "2"), ("3+", "3+")]


def counts_by_group(peptides: pl.DataFrame, min_donors: int = 2) -> pl.DataFrame:
    """Count peptide backbones per (enrichment, phospho-site count) for the bars.

    The unit is the bare backbone (uniprot + phospho-stripped sequence). A backbone
    is kept only if seen in ``min_donors`` or more distinct donors pooled across both
    enrichments. Each kept backbone contributes to every enrichment that detected it
    (so shared backbones count in both bars), bucketed by the max phospho-sites it
    reached in that enrichment: "1", "2", or "3+".
    """
    tagged = peptides.with_columns(
        pl.col("technical_replicate").str.extract(r"_(FeNTA|TiO2)$", 1).alias("enr"),
        pl.col("sequence").str.replace_all("*", "", literal=True).alias("bare_seq"),
    )
    # reproducibility gate: pooled donor count across enrichments
    kept = tagged.group_by("uniprot", "bare_seq").agg(
        pl.col("donor").n_unique().alias("n_donors")
    ).filter(pl.col("n_donors") >= min_donors).select("uniprot", "bare_seq")
    return (
        tagged.join(kept, on=["uniprot", "bare_seq"])
        .group_by("uniprot", "bare_seq", "enr")
        .agg(pl.col("num_modifications").max().alias("max_sites"))
        .with_columns(
            pl.when(pl.col("max_sites") >= 3).then(pl.lit("3+"))
            .when(pl.col("max_sites") == 2).then(pl.lit("2"))
            .otherwise(pl.lit("1")).alias("category")
        )
        .group_by("enr", "category").len()
    )


def make_stacked_bar(counts: pl.DataFrame, output_path: Path) -> None:
    """Draw and save the stacked bar in the sample style (2x2 in, 6 pt, thin lines)."""
    lw = 0.25
    plt.rcParams.update(
        {
            "font.size": 6,
            "font.family": "Arial",
            "text.color": "black",
            "axes.labelcolor": "black",
            "xtick.color": "black",
            "ytick.color": "black",
            "axes.linewidth": lw,
            "pdf.fonttype": 42,  # embed as editable TrueType, not paths
        }
    )

    # lookup: (enr, category) -> count
    lut = {(g, c): n for g, c, n in counts.select("enr", "category", "len").iter_rows()}

    fig, ax = plt.subplots(figsize=(1.55, 2))
    x = range(len(GROUPS))
    # stack bottom -> top per CATEGORY_STACK, so "mono" (lightest) lands on top
    for xi, (grp_key, _) in zip(x, GROUPS):
        bottom = 0
        for cat in CATEGORY_STACK:
            h = lut.get((grp_key, cat), 0)
            if h:
                ax.bar(xi, h, bottom=bottom, width=0.8, color=CATEGORY_COLORS[cat],
                       edgecolor=EDGE, linewidth=lw)
                bottom += h

    ax.set_xticks(list(x))
    ax.set_xticklabels([label for _, label in GROUPS])
    ax.set_ylabel("number of peptides")
    ax.set_ylim(bottom=0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out", length=2, width=lw)

    # legend on the right; fewest phospho-sites at top
    handles = [
        Patch(facecolor=CATEGORY_COLORS[cat], edgecolor=EDGE, linewidth=lw, label=label)
        for cat, label in CATEGORY_LEGEND
    ]
    ax.legend(
        handles=handles,
        title="phospho-sites\nper peptide",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        handlelength=1.0,
        handleheight=1.0,
        borderaxespad=0.2,
        labelspacing=0.3,
    )
    # narrow bar area (~0.55 in wide); reserve left for ylabel, right for the legend
    fig.subplots_adjust(left=0.27, right=0.62, bottom=0.18, top=0.95)
    fig.savefig(output_path)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Stacked bar of phospho-sites per peptide, one bar per enrichment (FeNTA/TiO2)."
    )
    parser.add_argument("--input-dir", type=Path, default=None,
                        help=f"Directory holding {PEPTIDES_FILENAME}.")
    parser.add_argument("--peptides", type=Path, default=None,
                        help=f"Path to the peptide table (default: <input-dir>/{PEPTIDES_FILENAME}).")
    parser.add_argument("--output-dir", type=Path, default=None,
                        help="Where to write the PNG (default: alongside the peptide table).")
    parser.add_argument("--min-donors", type=int, default=2,
                        help="Min distinct donors a peptide must appear in, pooled across "
                             "enrichments (default: 2).")
    args = parser.parse_args()

    if args.peptides is None and args.input_dir is None:
        parser.error("Pass --input-dir or --peptides.")

    peptides_path = args.peptides or (args.input_dir / PEPTIDES_FILENAME)
    if not peptides_path.exists():
        raise FileNotFoundError(f"Peptide table not found: {peptides_path}")

    output_dir = args.output_dir or peptides_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    peptides = pl.read_csv(peptides_path)
    counts = counts_by_group(peptides, min_donors=args.min_donors)

    print(f"Requiring >= {args.min_donors} donor(s) pooled across enrichments\n")
    print("backbones per enrichment x phospho-sites per peptide:")
    lut = {(g, c): n for g, c, n in counts.select("enr", "category", "len").iter_rows()}
    for enr_key, label in GROUPS:
        parts = ", ".join(
            f"{lbl}: {lut[(enr_key, cat)]}" for cat, lbl in CATEGORY_LEGEND
            if (enr_key, cat) in lut
        )
        total = sum(lut.get((enr_key, cat), 0) for cat, _ in CATEGORY_LEGEND)
        print(f"  {label:6} (n={total}): {parts}")

    output_path = output_dir / "num_mods_stacked_bar.pdf"
    make_stacked_bar(counts, output_path)
    print(f"\nWrote stacked bar: {output_path}")


if __name__ == "__main__":
    main()
