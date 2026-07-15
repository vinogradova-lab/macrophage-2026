"""Reactivity-changes data: the canonical long-format reader and its long->wide pivot.

The run's ``output/`` folder holds two tables that look interchangeable but are not:

  reactivity_changes_long_format.csv  one row per uniprot/residue/condition (canonical, complete)
  reactivity_changes.csv             the same data pivoted per protein, but MISSING 63 proteins

The wide file is built in reactivity_polars.ipynb by inner-joining the six per-condition frames on
``identifier_cols``, which includes the pipe-joined ``residue`` and ``sequence`` strings. Those
strings carry their group's row order, so the join silently requires every condition to be present
*and* to spell its residue list byte-identically. 63 proteins fail that: 24 lack a condition, 38
have a genuinely condition-dependent residue set, and PFKFB3 (Q16875) merely lists ``C193,C206``
before ``C193`` under TLR7. 20 of the 63 have reactivity changes, which is the whole 355-vs-335 gap.

``pivot_rc_wide`` does the same pivot keyed on ``uniprot`` alone, over a complete
(uniprot x residue x condition) grid, so ragged proteins survive with blanks rather than vanishing.
Aside from those 63 proteins the two agree exactly, so this is a superset of the wide file.
"""

import re

import polars as pl

from src.macrophage import RC_CONDITIONS, RC_LONG_FORMAT, require

# Constant per protein; taken from a first-non-null spine, never used as a join key.
PROTEIN_COLS = [
    "protein",
    "description",
    "rc_n_peptides",
    "uniprot_function",
    "uniprot_goterms",
    "location",
]
# Constant within a (uniprot, condition) group -> one scalar "<COND>_<name>" column.
CONDITION_SCALAR_COLS = [
    "wp_percent_control",
    "passes_bio_replicate_variation_filter",
    "median_percent_control_of_residues",
    "rc_percent_control_min",
    "rc_percent_control_max",
    "rc_min_percent_control_log2_fc",
    "reference_percent_control",
]
# Booleans reported as the *name* of the residue they hold for, '|'-joined (the notebook's format):
# "C65|C201" means those two residues pass, and blank means none do.
BINARY_RESIDUE_COLS = [
    "passes_expression_filter",
    "passes_ratio_or_median_filter",
    "passes_ratio_in_two_replicates",
    "passes_two_rep_variability_filter",
    "reactivity_change",
]
# Per-residue values, '|'-joined and positionally aligned to the protein's ``residue`` list; a
# residue a condition lacks holds an empty slot so the columns stay aligned.
RESIDUE_VALUE_COLS = [
    "rc_percent_control",
    "rc_wp_ratio",
    "rc_max_percent_control",
    "rc_min_percent_control",
    "ratio_to_median_of_residues",
    "rc_ratio_min_max",
    "n_replicates",
    "direction_of_reactivity_change",
]
# Per-condition column order in the pivoted output.
CONDITION_BLOCK = [
    "wp_percent_control",
    "passes_bio_replicate_variation_filter",
    "median_percent_control_of_residues",
    "rc_percent_control",
    "rc_wp_ratio",
    "rc_percent_control_min",
    "rc_percent_control_max",
    "rc_max_percent_control",
    "rc_min_percent_control",
    "rc_min_percent_control_log2_fc",
    "reference_percent_control",
    "passes_expression_filter",
    "ratio_to_median_of_residues",
    "rc_ratio_min_max",
    "passes_ratio_or_median_filter",
    "passes_ratio_in_two_replicates",
    "passes_two_rep_variability_filter",
    "n_replicates",
    "direction_of_reactivity_change",
    "reactivity_change",
    "num_reactivity_changes",
]
# Funnel annotations already merged into the long format. They are dropped here: cell 38 built them
# behind a .drop_nulls() that discards any residue missing conservation OR pPSE OR AlphaMissense, so
# 293 of the 1419 reactivity-change residues carry no annotation. src/cross_ref.py recomputes them
# from the reference databases instead, keyed on every cysteine rather than an arbitrary first one.
_STALE_ANNOTATION_COLS = [
    "pos",
    "Depth of Conservation (# of Organisms)",
    "pPSE",
    "structure_group",
    "Complex",
    "am_class",
]

_RESIDUE_TOKEN = re.compile(r"C(\d+)")


def rc_positions(value):
    """Cysteine positions in any of this project's residue encodings, as a set of ints.

    Handles ``C349``, ``C351|C366``, ``C498, C501``, ``(C177;C367)`` and the long format's
    unsorted ``C10,C6`` alike. ``(C244;C429)`` is a *set of changed residues*, not one peptide --
    splitting on ``|`` alone is what makes the wide and long tables look like they disagree.
    """
    if value is None:
        return set()
    return {int(m) for m in _RESIDUE_TOKEN.findall(str(value))}


def load_rc_long(path=RC_LONG_FORMAT):
    """The canonical long-format reactivity table (one row per uniprot/residue/condition).

    ``reactivity_change`` is Boolean with nulls (a residue whose filters could not be evaluated);
    nulls are not changes, so never compare it with ``== True`` after a cast to string.
    """
    return pl.read_csv(require(path, "Re-run reactivity/reactivity_polars.ipynb."))


def _first_position(col):
    return pl.col(col).str.extract(r"C(\d+)", 1).cast(pl.Int64)


def pivot_rc_wide(long_df, conditions=RC_CONDITIONS):
    """One row per protein: identifiers + a "<COND>_<metric>" block per condition.

    Keyed on ``uniprot`` over a complete (uniprot x residue x condition) grid, so a protein missing
    a condition, or listing its residues differently between conditions, keeps its row with blanks
    instead of being dropped by an inner join.
    """
    long_df = long_df.filter(pl.col("condition").is_in(conditions)).drop(
        [c for c in _STALE_ANNOTATION_COLS if c in long_df.columns]
    )

    # Complete the grid so every (uniprot, condition) sees the protein's full residue union in one
    # deterministic order -- this is what lets ragged proteins survive.
    grid = long_df.select("uniprot", "residue").unique().join(
        pl.DataFrame({"condition": conditions}), how="cross"
    )
    full = (
        grid.join(
            long_df.with_columns(pl.lit(True).alias("_measured")),
            on=["uniprot", "residue", "condition"],
            how="left",
        )
        .with_columns(_first_position("residue").alias("_pos"))
        .sort(["uniprot", "_pos", "residue"])
    )

    spine = (
        long_df.sort(["uniprot", "condition"])
        .group_by("uniprot", maintain_order=True)
        .agg([pl.col(c).drop_nulls().first().alias(c) for c in PROTEIN_COLS])
    )
    residues = (
        full.filter(pl.col("condition") == conditions[0])
        .group_by("uniprot", maintain_order=True)
        .agg(
            pl.col("residue").str.join("|").alias("residue"),
            pl.col("sequence").fill_null("").str.join("|").alias("sequence"),
        )
    )

    out = spine.join(residues, on="uniprot", how="left")
    for condition in conditions:
        out = out.join(_condition_block(full, condition), on="uniprot", how="left")
    return out


def _condition_block(full, condition):
    """The "<COND>_<metric>" columns for one condition, one row per uniprot."""
    changed = pl.col("reactivity_change").fill_null(False)
    agg = (
        full.filter(pl.col("condition") == condition)
        .group_by("uniprot", maintain_order=True)
        .agg(
            pl.col("_measured").drop_nulls().len().alias("_n_measured"),
            # Residues that changed, in the protein's residue order.
            pl.col("residue").filter(changed).alias("_changed"),
            # The notebook counts an ambiguous protein's changes as one; see below.
            pl.col("wp_percent_control").drop_nulls().len().alias("_n_wp"),
            pl.col("rc_n_peptides").drop_nulls().first().alias("_n_peptides"),
            *[pl.col(c).drop_nulls().first().alias(c) for c in CONDITION_SCALAR_COLS],
            *[
                pl.col("residue").filter(pl.col(c).fill_null(False)).str.join("|").alias(c)
                for c in BINARY_RESIDUE_COLS
                if c != "reactivity_change"
            ],
            *[
                pl.col(c).cast(pl.String).fill_null("").str.join("|").alias(c)
                for c in RESIDUE_VALUE_COLS
            ],
        )
    )

    # A protein with no reference expression and few peptides cannot attribute the change to one
    # residue, so the notebook collapses those residues to "(C1;C2)" and counts them as ONE change.
    ambiguous = (
        pl.col("_n_wp").eq(0)
        & pl.col("_n_peptides").lt(5)
        & pl.col("_changed").list.len().gt(0)
    )
    agg = agg.with_columns(
        pl.when(ambiguous)
        .then(
            pl.concat_str(pl.lit("("), pl.col("_changed").list.join(";"), pl.lit(")"))
        )
        .otherwise(pl.col("_changed").list.join("|"))
        .alias("reactivity_change"),
        pl.when(ambiguous)
        .then(pl.lit(1))
        .otherwise(pl.col("_changed").list.len())
        .cast(pl.Int64)
        .alias("num_reactivity_changes"),
    )

    # A condition the protein was never measured in stays blank rather than reading as "0 changes".
    unmeasured = pl.col("_n_measured").eq(0)
    agg = agg.with_columns(
        [
            pl.when(unmeasured).then(None).otherwise(pl.col(c)).alias(c)
            for c in CONDITION_BLOCK
        ]
    )
    return agg.select(
        "uniprot", *[pl.col(c).alias(f"{condition}_{c}") for c in CONDITION_BLOCK]
    )


def rc_change_positions(wide_df, conditions=RC_CONDITIONS):
    """{uniprot: set of cysteine positions that change in at least one condition}."""
    cols = [f"{c}_reactivity_change" for c in conditions]
    out = {}
    for row in wide_df.select("uniprot", *cols).iter_rows(named=True):
        pos = set()
        for c in cols:
            pos |= rc_positions(row[c])
        if pos:
            out[row["uniprot"]] = pos
    return out


def has_reactivity_change(conditions=RC_CONDITIONS):
    """Expression selecting proteins with a reactivity change in at least one condition."""
    return pl.any_horizontal(
        [pl.col(f"{c}_num_reactivity_changes").fill_null(0).gt(0) for c in conditions]
    )
