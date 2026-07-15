"""Venn diagram of peptides detected in the FeNTA vs TiO2 phospho enrichments.

Reads the tidy ``peptide_results_clean.csv`` written by ``aggregate_phospho.py``
and draws a two-set Venn of the peptides seen in each enrichment. Peptides are
keyed on (uniprot, bare sequence): the phosphosite-location markers (``*``) are
stripped from the sequence first, so the same backbone phosphorylated at
different residues counts as one peptide. The enrichment is carried in the
``technical_replicate`` column as a ``_TiO2`` / ``_FeNTA`` suffix
(aggregate_phospho appends it to the channel name), so we split it back out here.

A peptide only counts if it was detected in ``--min-donors`` distinct donors
(default 2) pooled *across enrichments* -- i.e. a donor counts whether it saw the
peptide in FeNTA or TiO2. A peptide passing this gate then lands in the FeNTA
and/or TiO2 set(s) according to where it was actually seen. Each donor is one TZ3
experiment (d2=TZ3-114, d3=TZ3-116, d1=TZ3-117, d5=TZ3-121); the ``donor`` column
carries this, while the trailing 1/2/3 of ``technical_replicate`` is within-donor
replication and is not counted here.

Conda environment (the lab's ``polars`` env already has these):

    conda activate polars   # needs polars, matplotlib, matplotlib-venn

Usage:

    python venn_enrichment.py --input-dir /path/to/phospho_dir
    python venn_enrichment.py --peptides /path/to/peptide_results_clean.csv
    python venn_enrichment.py --input-dir DIR --output-dir OUT
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import polars as pl
from matplotlib_venn import venn2

PEPTIDES_FILENAME = "peptide_results_clean.csv"
ENRICHMENTS = ("FeNTA", "TiO2")
# a peptide is identified by protein + bare sequence, ignoring phosphosite location
PEPTIDE_KEY = ["uniprot", "sequence"]


def peptides_by_enrichment(peptides: pl.DataFrame, min_donors: int = 2) -> dict[str, set]:
    """Map each enrichment to its detected (uniprot, bare-sequence) peptides.

    A peptide is kept only if it was seen in ``min_donors`` or more distinct
    donors pooled across both enrichments (the ``donor`` column, one per TZ3
    experiment). Kept peptides then populate the FeNTA and/or TiO2 set(s)
    according to the enrichment(s) in which they were detected.
    """
    tagged = peptides.with_columns(
        # enrichment is the suffix on technical_replicate, e.g. "1_TiO2"
        pl.col("technical_replicate").str.extract(r"_(FeNTA|TiO2)$", 1).alias("enrichment"),
        # drop phosphosite-location markers so we key on the bare backbone sequence
        pl.col("sequence").str.replace_all("*", "", literal=True).alias("sequence"),
    )
    # donor count pooled across enrichments; keep peptides meeting the threshold
    kept = (
        tagged.group_by(PEPTIDE_KEY)
        .agg(
            pl.col("donor").n_unique().alias("n_donors"),
            pl.col("enrichment").unique().alias("enrichments"),
        )
        .filter(pl.col("n_donors") >= min_donors)
    )
    sets = {}
    for enrichment in ENRICHMENTS:
        rows = (
            kept.filter(pl.col("enrichments").list.contains(enrichment))
            .select(PEPTIDE_KEY)
            .iter_rows()
        )
        sets[enrichment] = set(rows)
    return sets


def make_venn(sets: dict[str, set], output_path: Path) -> None:
    """Draw and save the two-set Venn (2x2 in, Arial 6 pt black, thin outlines)."""
    lw = 0.25
    plt.rcParams.update(
        {
            "font.size": 6,
            "font.family": "Arial",
            "text.color": "black",
            "axes.linewidth": lw,
            "pdf.fonttype": 42,  # embed as editable TrueType, not paths
        }
    )
    fentA, tio2 = sets["FeNTA"], sets["TiO2"]
    fig, ax = plt.subplots(figsize=(1.3, 1.3))
    v = venn2(
        [fentA, tio2],
        set_labels=("FeNTA", "TiO2"),
        set_colors=("#4C72B0", "#DD8452"),  # echoes the stacked-bar / source colors
        alpha=0.6,
        ax=ax,
    )
    # thin dark outlines on the circles, matching the bar's box outlines
    for patch in v.patches:
        if patch:
            patch.set_edgecolor("#3a3a3a")
            patch.set_linewidth(lw)
    for text in (v.set_labels or []):
        if text:
            text.set_fontsize(6)
            text.set_color("black")
    for text in (v.subset_labels or []):
        if text:
            text.set_fontsize(6)
            text.set_color("black")
    ax.set_title("Phosphopeptides detected\nby enrichment", fontsize=6)
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.04, top=0.82)
    fig.savefig(output_path)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Venn diagram of peptides detected in the FeNTA vs TiO2 enrichments, "
            "from aggregate_phospho.py's peptide_results_clean.csv."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=None,
        help=f"Directory holding {PEPTIDES_FILENAME} (default source of the peptide table).",
    )
    parser.add_argument(
        "--peptides",
        type=Path,
        default=None,
        help=f"Path to the peptide table (default: <input-dir>/{PEPTIDES_FILENAME}).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Where to write the Venn PNG (default: alongside the peptide table).",
    )
    parser.add_argument(
        "--min-donors",
        type=int,
        default=2,
        help="Min distinct donors a peptide must appear in, pooled across enrichments (default: 2).",
    )
    args = parser.parse_args()

    if args.peptides is None and args.input_dir is None:
        parser.error("Pass --input-dir or --peptides.")

    peptides_path = args.peptides or (args.input_dir / PEPTIDES_FILENAME)
    if not peptides_path.exists():
        raise FileNotFoundError(f"Peptide table not found: {peptides_path}")

    output_dir = args.output_dir or peptides_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    peptides = pl.read_csv(peptides_path)
    sets = peptides_by_enrichment(peptides, min_donors=args.min_donors)

    fentA, tio2 = sets["FeNTA"], sets["TiO2"]
    print(f"Requiring >= {args.min_donors} donor(s) pooled across enrichments\n")
    print(f"FeNTA peptides:      {len(fentA)}")
    print(f"TiO2 peptides:       {len(tio2)}")
    print(f"Shared (both):       {len(fentA & tio2)}")
    print(f"FeNTA only:          {len(fentA - tio2)}")
    print(f"TiO2 only:           {len(tio2 - fentA)}")
    print(f"Union (total):       {len(fentA | tio2)}")

    output_path = output_dir / "venn_enrichment_peptides.pdf"
    make_venn(sets, output_path)
    print(f"\nWrote Venn diagram:  {output_path}")


if __name__ == "__main__":
    main()
