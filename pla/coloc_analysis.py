#!/usr/bin/env python3
"""Superplots for the collaborator per-cell colocalization data.

Reads the per-cell CSVs (three antibody pairs x M0/LPS) from one or more deliveries, and
for each pair renders both directions of the split colocalization coefficient as an
M0-vs-LPS superplot: every small dot a cell, every large dot a donor mean/median, p from a
paired t-test across donors.

Two deliveries, two layouts
---------------------------
Donors arrive in batches (`--data-root`, repeatable) and are pooled into one figure: D1-D3
from `20260729_Yannis_D1-D3` and D4-D6 from `20260910_chromaticabberation` give the n=6 run.

D4-D6 is read from that 20260910 re-export rather than from the original
`20260908_yannis_D4-D6` delivery. The re-export re-runs the same cells with the chromatic
aberration of the Fy2 channel corrected: every Fy2-side count, size and intensity moves and
with it every coloc coefficient (by up to 0.29 on Fy2/Fy1), while the first marker's columns
come back byte-identical. It also restores the `Masked Area (Physical Units)` column that
the 20260908 export had dropped, which is what lets the density QC and the size-normalized
vesicle count reach all six donors. Note that only D4-D6 has been re-exported this way; the
D1-D3 tables are as delivered.

The exporter changed its metadata layout between the two deliveries and nothing in a file
announces which one it wrote, so `load_csv` sniffs it -- see DONOR_IN_FILENAME_RE. The
differences that matter here:

    D1-D3   GeneName = condition, siRNA_ID = donor (a within-batch 1/2/3).
    D4-D6   siRNA_ID = condition, GeneName = the marker pair, and no donor column at all --
            the donor is FileName's `bd26-NN` prefix.

Both now carry `Masked Area (Physical Units)`, the per-cell size measure that the density QC
and the normalized vesicle count both divide by. See CELL_SIZE_COL for why that column, and
not the `Mean: Volume[Total] : Cells3D` volume, is the one comparable across deliveries.

Donor ids are not comparable across deliveries -- D1-D3 numbers its donors 1/2/3 -- so
`load_pair` prefixes them with the batch and `assign_donor_labels` renumbers the union to
D1..Dn once for the whole run. `donor_key.csv` maps those back to the acquisition ids.

The two directions are not redundant. Coloc(X/Y) is the fraction of the X signal sitting
on Y ("how much of X is Y-positive"); Coloc(Y/X) is the reverse. They differ whenever the
two markers differ in abundance or spatial extent, and the pair is what distinguishes
"X is a subset of a broader Y pool" from "X and Y mark the same compartment".

Density QC
----------
Alongside the superplots the script checks the coloc coefficients against vesicle density,
a confound the superplots cannot show. Two per-cell quantities are correlated (Spearman):

    y  Coloc. (X/Y)Thesh=0.25 SubRand=TRUE by InegrIntense   [or `... by Number`]
       -- the exporter column itself, picked by `coloc_columns()`.
    x  Numb. Ves. (Y) / Masked Area (Physical Units)
       -- vesicle density of the *denominator* marker, derived by `density()`.

Y rather than X because the confound is directional: the denser the Y objects, the more
likely any X object overlaps one, so Coloc(X/Y) drifts up with density(Y) on chance-overlap
grounds alone. `SubRand=TRUE` is supposed to remove exactly that, and empirically does not
fully -- pooled within condition, rho reaches +0.76 for Fy4/Fy2 and +0.56 for Fy1/Fy2, and
the 72 donor x condition strata run from -0.24 to +0.91. Because density also varies by
donor and by condition, some of the M0/LPS coloc
difference may be segmentation rather than biology. `density_qc()` quantifies this at four
pooling scopes (see its docstring for the pooled-vs-per-donor trade-off) into
`coloc_density_qc.csv`, and `density_scatter()` renders it per direction.

Vesicle abundance and size
--------------------------
The same superplot treatment is also applied to the per-cell vesicle metrics of each
pair's non-Fy2 marker -- raw count, count normalized to cell size, and mean vesicle size
(see `vesicle_specs`). The normalized count divides by CELL_SIZE_COL, read inline from the
analysis table of every delivery. A pair whose donors are not *all* covered loses the
normalized panel rather than drawing it at a reduced n -- see `process_vesicle_metrics`.

Output layout
-------------
`--outdir` is split by figure type, each subfolder carrying its own stats CSV:

    colocalization/             the M0-vs-LPS superplots + coloc_stats.csv
    colocalization_percent_control/  the same six, each cell as a percentage of its own
                                donor's M0 summary, so every donor's M0 sits at 100 and the
                                donor batch effect stops dominating the axis. Annotated with
                                the same P as colocalization/ -- see
                                `process_direction_percent` for why.
    density_confound_qc/        the coloc-vs-density scatters + coloc_density_qc.csv
    vesicle_abundance_and_size/ the vesicle count/size superplots + vesicle_stats.csv

Colocalization panels are named by direction -- `coloc_Fy1-to-Fy2_intensity_mean.svg` is
Coloc(Fy1/Fy2), the fraction of the Fy1 signal sitting on Fy2, and is a different quantity
from its `Fy2-to-Fy1` counterpart.

Usage:
    python coloc_analysis.py --outdir plots_n6      # the n=6 run
    python coloc_analysis.py --metric both          # also the number-based coefficient
    python coloc_analysis.py --no-vesicle-metrics   # coloc superplots only
    python coloc_analysis.py --data-root DIR --data-root DIR2    # pick the deliveries
"""

import argparse
import glob
import os
import re

import pandas as pd

from pla_plotting import (
    AGGS,
    CONDITIONS,
    donor_palette,
    fisher_combine_rho,
    grand_summary,
    paired_ttest,
    plot_density_scatter,
    plot_swarm,
    spearman_rho,
    summarize_donors,
)

_IF_ROOT = (
    "/Users/henrysanford/Dropbox @RU Dropbox/Vinogradova Laboratory/"
    "Macrophage project/Manuscript 1/01_Data Analysis/08_IF"
)

# The two donor batches of the n=6 run, each the 0.15-um size-gated analysis tables.
#
# D4-D6 is taken from the 20260910 chromatic-aberration re-export, not from the 20260908
# original. Preferring the corrected coloc coefficients is the obvious reason; the reason it
# is *required* here is that the 20260908 export dropped `Masked Area (Physical Units)` and
# the 20260910 one restores it -- see CELL_SIZE_COL.
#
# D1-D3's ungated `no filter/` companions used to be a third root, as the only place a cell
# volume could be found for those donors. Nothing reads that volume any more, so they are
# not listed; `batch_label` still recognizes the folder for anyone who passes it explicitly.
DEFAULT_DATA_ROOTS = [
    os.path.join(_IF_ROOT, "20260729_Yannis_D1-D3", "0.15-filter"),
    os.path.join(_IF_ROOT, "20260910_chromaticabberation", "0.15 filter"),
]

# Non-numeric columns, exempt from the numeric coercion below.
META_COLS = ("FileName", "GeneName", "siRNA_ID", "Extension", "FullName")

# The exporter writes this string where a cell has no value. Left as-is it silently makes
# the whole column `object` dtype, and every later aggregation raises TypeError.
NA_VALUES = ["EMPTY", " EMPTY"]

# Metric selector -> the substring that identifies it in the column name. Both spellings
# are the exporter's own typos ("Thesh", "InegrIntense"); do not "fix" them.
METRIC_TOKENS = {"intensity": "by InegrIntense", "number": "by Number"}

# Coloc columns embed the marker names, so they differ per file:
#   " Coloc. (Rab5/Fy2)Thesh=0.25 SubRand=TRUE by InegrIntense"
COLOC_RE = re.compile(r"Coloc\. \((\w+)/(\w+)\)")

# The exporter changed its metadata layout between the two deliveries, and nothing in a
# file announces which layout it uses:
#
#   D1-D3   GeneName = condition (M0/LPS),  siRNA_ID = donor (1/2/3)
#   D4-D6   siRNA_ID = condition (M0/LPS),  GeneName = the marker pair ("Fy1_Fy2"), and
#           there is no donor column at all -- the donor is the `bd26-NN` prefix of FileName
#
# Read the D1-D3 way, a D4-D6 file yields a Condition of "Fy1_Fy2" that casts to all-NaN
# and a "Donor" of M0/LPS, so `load_csv` sniffs the layout rather than assuming one.
DONOR_IN_FILENAME_RE = re.compile(r"(bd26-\d+)", re.IGNORECASE)

# The collaborator's own name for a donor, as it appears in FileName: either a `bd26-NN`
# bleed code or an 8-digit acquisition date. Recorded alongside the D-numbers in
# donor_key.csv purely for traceability -- D1-D3's `siRNA_ID` is a within-batch 1/2/3 that
# means nothing outside its own delivery, so it is not enough to identify a donor.
SOURCE_ID_RE = re.compile(r"(bd26-\d+|\d{8})", re.IGNORECASE)

# `fy42nd-fy2` (D1-D3) and `fy4-fy2` (D4-D6) are the same antibody pair -- both carry
# `Coloc. (Fy4/Fy2)` columns. Without this they discover as two 3-donor pairs, and the n=6
# run quietly becomes two n=3 runs.
PAIR_ALIASES = {"fy42nd": "fy4"}

# Leaf folders that name the vesicle size gate rather than a delivery. Three spellings are
# in use across the deliveries -- `0.15-filter` (D1-D3, D4-D6), `0.15 filter` and
# `no filter` (the 20260910 chromatic-aberration re-export, and D1-D3's companions) -- and
# `batch_label` has to recognize all of them, or two roots of the same delivery come out
# under different batch names and their donors stop lining up.
GATE_DIR_RE = re.compile(r"(?:[\d.]+|no)[ _-]?filter", re.IGNORECASE)

# Re-downloaded copies land beside the originals as `...percell[18].csv`. They are
# byte-identical (verified by checksum), so letting them through would silently double
# every cell in the pair.
DUPLICATE_RE = re.compile(r"\[\d+\]\.csv$")

# Output subdirectories, one per kind of figure, each holding its own stats CSV. Grouped
# by figure type rather than by antibody pair because the panels of one type are what get
# assembled into a figure together, and because the three types answer different questions:
# COLOC_DIR is the result, DENSITY_QC_DIR is the check on it, VESICLE_DIR is a separate
# measurement that happens to come from the same files.
COLOC_DIR = "colocalization"
PERCENT_CONTROL_DIR = "colocalization_percent_control"
DENSITY_QC_DIR = "density_confound_qc"
VESICLE_DIR = "vesicle_abundance_and_size"

# The per-cell size measure. Two quantities are normalized by it -- the vesicle density the
# coloc QC correlates against, and the size-normalized vesicle count -- and both are compared
# across deliveries, so what matters about this column is not that it is the most physical
# measure available but that it means the same thing in both.
#
# The obvious alternative, `Mean: Volume[Total] : Cells3D`, does not. Take the per-cell ratio
# of that column to this one and it is 3.639 for every D1-D3 cell (SD 0.039, n=999) and 0.996
# for every D4-D6 cell (SD 0.020, n=755): D1-D3's "volume" is the masked area times a ~3.6 um
# cell height, D4-D6's is the masked area itself, un-scaled in z. That is a per-delivery
# calibration, not a per-cell measurement, so dividing counts by it put the two batches on
# denominators 3.6x apart and split the count-per-volume panel into a D1-D3 cluster and a
# D4-D6 cluster with no biology between them. Within D1-D3 the two denominators differ by a
# constant, so nothing is lost by using this one: the panel is the same figure rescaled.
#
# It is an *area*, not a volume, and the panels below say so. The 20260908 D4-D6 export had
# dropped it, which is why DEFAULT_DATA_ROOTS points at the 20260910 re-export instead.
CELL_SIZE_COL = "Masked Area (Physical Units)"


# ------------------------------------------------------------------
# Data loading
# ------------------------------------------------------------------
def subdir(outdir, name):
    """Create and return `outdir/name`, one of the per-figure-type output folders."""
    path = os.path.join(outdir, name)
    os.makedirs(path, exist_ok=True)
    return path


def load_csv(path):
    """Read one per-cell CSV, normalizing the exporter's formatting quirks."""
    df = pd.read_csv(path, na_values=NA_VALUES, skipinitialspace=False)
    df.columns = [c.strip() for c in df.columns]

    # Everything that isn't a metadata column should be numeric; anything left as object
    # dtype here means a stray sentinel we haven't listed in NA_VALUES.
    for col in df.columns:
        if col not in META_COLS and df[col].dtype == object:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Values are quoted with an inconsistent leading space (' "1"' in most files, '"3"' in
    # fy42nd-fy2_m0). Strip whitespace then quotes, or donor 3 splits into a 4th donor.
    def clean(s):
        return s.astype(str).str.strip().str.strip('"').str.strip()

    df["SourceID"] = (df["FileName"].astype(str)
                      .str.extract(SOURCE_ID_RE, expand=False).str.lower())

    gene, sirna = clean(df["GeneName"]), clean(df["siRNA_ID"])
    conditions = set(CONDITIONS)
    if set(gene.unique()) <= conditions:            # D1-D3 layout
        df["Condition"], df["Donor"] = gene, sirna
    elif set(sirna.unique()) <= conditions:         # D4-D6 layout
        df["Condition"] = sirna
        donor = df["FileName"].astype(str).str.extract(DONOR_IN_FILENAME_RE, expand=False)
        if donor.isna().any():
            raise SystemExit(
                f"{os.path.basename(path)} uses the D4-D6 layout, where the donor is the "
                f"'bd26-NN' prefix of FileName, but {int(donor.isna().sum())} of "
                f"{len(df)} rows have no such prefix")
        df["Donor"] = donor.str.lower()
    else:
        raise SystemExit(
            f"{os.path.basename(path)}: cannot tell which column holds the condition. "
            f"Expected {sorted(conditions)} in either GeneName (got "
            f"{sorted(gene.unique())[:4]}) or siRNA_ID (got {sorted(sirna.unique())[:4]})")
    return df


def batch_label(root):
    """A short name for the delivery a root belongs to, used to keep donors distinct.

    A leaf like `0.15-filter` or `no filter` names the size gate, not the batch, so for
    those roots the batch is the parent folder -- otherwise both deliveries would come out
    labelled "0.15-filter" and their donors would collide.

    Folding `no filter` in as well is what keeps D1-D3's two roots -- the gated analysis
    tables in `0.15-filter/` and the ungated companions in `no filter/` -- under one batch
    name. `load_pair` looks a companion up by `(condition, batch)`, so if they disagreed
    the lookup would return None and every D1-D3 cell would lose its volume with no error
    anywhere; the volume-normalized vesicle panel would just quietly drop to n=3.
    """
    root = os.path.normpath(root)
    base = os.path.basename(root)
    if GATE_DIR_RE.fullmatch(base):
        return os.path.basename(os.path.dirname(root))
    return base


def discover_pairs(roots):
    """Group the CSVs across all `roots` into {pair: {condition: [(path, batch), ...]}}.

    Extends the single-root version to several deliveries at once, which is what makes an
    n=6 run possible: one batch supplies donors 1-3 and another donors 4-6, and one antibody
    pair is assembled from both. Filenames are punctuated inconsistently (`size-gr` vs
    `size_gr`, `fy42nd` vs `fy4`), so the pair name is taken from the `<pair>-fy2` prefix and
    mapped through PAIR_ALIASES, and the condition from an `_m0`/`_lps` token.

    Each pair must have both conditions, or the M0-vs-LPS panel it feeds is half a
    comparison drawn as if it were whole.
    """
    pairs = {}
    for root in roots:
        batch = batch_label(root)
        paths = sorted(glob.glob(os.path.join(root, "*.csv")))
        if not paths:
            raise SystemExit(f"No .csv files found in {root!r}")
        for path in paths:
            base = os.path.basename(path).lower()
            if DUPLICATE_RE.search(base):
                print(f"  Skipping (re-downloaded duplicate): {os.path.basename(path)}")
                continue
            m = re.match(r"([a-z0-9]+)-fy2", base)
            cond = re.search(r"_(m0|lps)[_-]", base)
            if not (m and cond):
                print(f"  Skipping (couldn't parse pair/condition): {base}")
                continue
            pair = PAIR_ALIASES.get(m.group(1), m.group(1))
            (pairs.setdefault(pair, {})
                  .setdefault(cond.group(1).upper(), []).append((path, batch)))

    for pair, cells in sorted(pairs.items()):
        missing = set(CONDITIONS) - set(cells)
        if missing:
            raise SystemExit(f"Pair {pair!r} is missing condition(s) {sorted(missing)}; "
                             f"found {sorted(cells)}")
    return pairs


def load_pair(cells):
    """Load a pair's M0 and LPS analysis tables, across every batch, into one tidy frame.

    Donor labels are prefixed with the batch here. The raw ids are not unique across
    deliveries -- D1-D3 numbers its donors 1-3 -- and silently merging two donors would turn
    n=6 back into n=3 with no error anywhere.

    CELL_SIZE_COL is required rather than treated as optional. It is missing from exactly
    one export -- the superseded 20260908 D4-D6 delivery -- and its absence would otherwise
    surface only as the normalized vesicle panel quietly disappearing from the output.
    """
    frames = []
    for cond in CONDITIONS:
        for path, batch in cells[cond]:
            df = load_csv(path)
            if CELL_SIZE_COL not in df.columns:
                raise SystemExit(
                    f"{os.path.basename(path)} has no {CELL_SIZE_COL!r} column. The "
                    f"20260908 D4-D6 export dropped it; read D4-D6 from the 20260910 "
                    f"chromatic-aberration re-export instead.")
            df["Batch"] = batch
            df["Donor"] = batch + ":" + df["Donor"].astype(str)
            frames.append(df)

    data = pd.concat(frames, ignore_index=True)
    data["Condition"] = pd.Categorical(data["Condition"], categories=CONDITIONS,
                                       ordered=True)
    if data["Condition"].isna().any():
        raise ValueError(f"Unexpected condition labels (expected {CONDITIONS})")
    return data


def assign_donor_labels(frames):
    """{batch-prefixed donor -> "D1", "D2", ...} over every pair, in batch then id order.

    Computed once across all pairs rather than per pair, so that a donor carries the same
    D-number and the same color in every figure, and so that a pair which happened to be
    missing a donor could not silently shift the numbering of the others.

    The raw ids are not comparable between deliveries (D1-D3 numbers its donors 1-3, D4-D6
    identifies them as `bd26-46/52/53`), so the sort is on the batch-prefixed key; the
    batch folder names are date-prefixed, so this orders the deliveries chronologically.
    """
    donors = sorted({d for f in frames for d in f["Donor"].unique()})
    return {d: f"D{i}" for i, d in enumerate(donors, start=1)}


def coloc_columns(data, metric):
    """[(column, X, Y)] for the requested metric, in file order (X/Fy2 then Fy2/X)."""
    token = METRIC_TOKENS[metric]
    out = []
    for col in data.columns:
        m = COLOC_RE.search(col)
        if m and token in col:
            out.append((col, m.group(1), m.group(2)))
    if len(out) != 2:
        raise ValueError(f"Expected 2 {metric!r} coloc columns, found {len(out)}: "
                         f"{[c for c, _, _ in out]}")
    return out


# ------------------------------------------------------------------
# Density QC
# ------------------------------------------------------------------
def density(data, marker):
    """Vesicle density: segmented objects per unit masked area, one value per cell.

        `Numb. Ves. (<marker>)` / `Masked Area (Physical Units)`

    Note what the numerator is: `Numb. Ves.` is the object count the segmentation produced
    for this cell, not a physical quantity measured independently of it. So anything that
    makes one donor or one condition segment more objects -- staining efficiency, exposure,
    the `size_gr 0.15` size gate, focus -- moves this number without any change in biology.
    That is what makes it a confounder for the coloc coefficients rather than just a covariate.
    """
    return data[f"Numb. Ves. ({marker})"] / data[CELL_SIZE_COL]


def _scope_row(pair, X, Y, col, g, scope, condition, donor):
    """One QC row: densities, coloc and the coloc-vs-density Spearman rho for a cell subset.

    `condition` / `donor` are the string "all" wherever the scope pools over them, so that
    every scope lands on the same columns and the CSV stays rectangular. Returns None when
    fewer than 3 cells have a coloc value, since rho is undefined there.
    """
    ok = g[col].notna()
    if ok.sum() < 3:
        return None
    dens_y = density(g, Y)[ok]
    rho, pval = spearman_rho(dens_y, g.loc[ok, col])
    return {
        "pair": pair, "direction": f"{X}/{Y}", "numerator": X, "denominator": Y,
        "scope": scope, "condition": str(condition), "donor": str(donor),
        "n_cells": int(ok.sum()), "n_strata": 1,
        "mean_density_numerator": density(g, X)[ok].mean(),
        "mean_density_denominator": dens_y.mean(),
        "mean_cell_area": g.loc[ok, CELL_SIZE_COL].mean(),
        "mean_coloc": g.loc[ok, col].mean(),
        "rho_coloc_vs_density_denominator": rho,
        "p_rho": pval,
    }


def density_qc(pair, data, col, X, Y):
    """Coloc-vs-vesicle-density coupling for one direction, at five pooling scopes.

    Colocalization rises with the density of the *denominator* marker Y: the more Y objects
    per unit area, the more likely any X object overlaps one. `SubRand=TRUE` already
    subtracts an expected-random overlap, so this should be corrected for -- the point of
    the QC is that empirically it isn't, and since density also varies by donor and by
    condition, part of any M0/LPS coloc difference may be segmentation rather than biology.

    The `scope` column says how cells were grouped before correlating. Pooled is the primary
    scope -- the per-donor values exist to check it, not to replace it:

    - `pooled_within_condition` -- **the number to quote.** All donors of one condition
      together, n ~ 120-175. Pro: rho is estimated on enough cells to be stable, and M0 and
      LPS are pooled identically, so comparing them is like-for-like. Con: donor is a common
      cause of both variables (donor-level segmentation quality drives density, donor-level
      biology drives coloc), so the value absorbs between-donor structure and its magnitude
      is an upper bound. `p_rho` is anticonservative for the same reason -- pooled cells are
      not independent -- so read the rho, not the p.
    - `pooled_all` -- both conditions together. Context only; it mixes the M0/LPS contrast
      into the correlation.
    - `within_donor` -- one rho per donor x condition, n ~ 30-80. Too noisy to quote
      individually, and six loose values invite picking whichever suits the argument.
    - `combined_within_donor[_by_condition]` -- the `within_donor` rhos Fisher-averaged
      (`fisher_combine_rho`), over all strata and within each condition. This is the sign
      check on the pooled value. Where it agrees in sign, pooling is just amplifying a real
      within-donor effect and the pooled number is safe to use. Where it *disagrees*, the
      pooled correlation is a between-donor artifact -- Simpson's paradox -- and quoting it
      invites a reviewer who stratifies to get the opposite sign. In this dataset that
      happens for `rab5 Fy2/Rab5` (pooled -0.11, within-donor +0.23). `p_rho` is NaN for
      these scopes: an average of rhos is not a test on one sample.

    None of these scopes, pooled or not, establishes that an M0/LPS coloc difference is
    segmentation. They establish that the coefficient is density-sensitive. Whether the
    conditions actually differed in density is a donor-mean question -- see
    `print_donor_shift_check`.
    """
    rows = []

    def add(g, scope, condition="all", donor="all"):
        row = _scope_row(pair, X, Y, col, g, scope, condition, donor)
        if row is not None:
            rows.append(row)
        return row

    pooled = {"all": add(data, "pooled_all")}
    for cond, g in data.groupby("Condition", observed=True):
        pooled[str(cond)] = add(g, "pooled_within_condition", condition=cond)

    first_within = len(rows)
    for (cond, donor), g in data.groupby(["Condition", "Donor"], observed=True):
        add(g, "within_donor", condition=cond, donor=donor)
    within = rows[first_within:]

    # Fisher-averaged within-donor rho, over all strata and then within each condition.
    for scope, cond in [("combined_within_donor", "all"),
                        *[("combined_within_donor_by_condition", c) for c in CONDITIONS]]:
        strata = [r for r in within if cond in ("all", r["condition"])]
        context = pooled.get(cond)
        if context is None or not strata:
            continue
        rho, n_strata = fisher_combine_rho(
            [r["rho_coloc_vs_density_denominator"] for r in strata],
            [r["n_cells"] for r in strata],
        )
        # Context columns come from the pooled row over the same cells -- only the rho and
        # the stratum count are specific to this scope. p_rho is NaN: this is an average of
        # rhos, not a test on one sample.
        rows.append({**context, "scope": scope, "n_strata": n_strata,
                     "rho_coloc_vs_density_denominator": rho, "p_rho": float("nan")})
    return rows


def density_scatter(pair, data, col, X, Y, qc_rows, outdir, palette):
    """Save the per-cell coloc-vs-density scatter for one direction.

    The correlation annotation is read back out of `qc_rows` rather than recomputed, so the
    figure and `coloc_density_qc.csv` cannot disagree.
    """
    x_col = f"Vesicle density ({Y})"
    # assign() rather than a plain assignment: `data` is reused for the other direction and
    # for the superplots, and shouldn't collect derived columns as a side effect.
    plot_data = data.assign(**{x_col: density(data, Y)}).dropna(subset=[col, x_col])

    # One r over every cell in the panel, both conditions and all donors -- the `pooled_all`
    # scope, so it describes the same points the trendline is fitted to. Written "Spearman r"
    # rather than as the rho glyph: at 6 pt a rho reads as a p, and a "p" next to a
    # correlation gets read as a p-value. The per-condition values and the within-donor sign
    # check live in the console summary and coloc_density_qc.csv.
    pooled = [r["rho_coloc_vs_density_denominator"] for r in qc_rows
              if r["scope"] == "pooled_all"]

    out_svg = os.path.join(outdir, f"coloc_{X}-to-{Y}_vs_density_{Y}.svg")
    plot_density_scatter(
        plot_data, x_col, col, out_svg, palette,
        title=f"{X}/{Y} — {len(plot_data)} cells, "
              f"{plot_data['Donor'].nunique()} donors",
        xlabel=f"Vesicle density ({Y}) per cell",
        ylabel=f"Coloc. ({X}/{Y}) per cell",
        annotation=f"Spearman r = {pooled[0]:+.2f}" if pooled else None,
    )
    print(f"    saved {out_svg}")


# ------------------------------------------------------------------
# Orchestration
# ------------------------------------------------------------------
def percent_of_control(data, col, agg):
    """Per-cell values as a percentage of the same donor's own M0 summary.

    Normalizes by the donor's M0 value computed with the *same* `agg` the panel will
    summarize with. Dividing a donor's cells by a constant commutes with both the mean and
    the median, so each donor's M0 summary then lands on exactly 100 either way -- which is
    the property that makes the panel readable. A donor with no M0 cells maps to NaN and
    drops out; `paired_ttest` already intersects donors present in both conditions.
    """
    m0 = data[data["Condition"] == "M0"].groupby("Donor", observed=True)[col].agg(agg)
    return 100 * data[col] / data["Donor"].map(m0)


def superplot(data, col, aggs, outdir, palette, stem, title_prefix, ylabel, ident,
              test_override=None, hline=None, ylim=None):
    """Save one superplot per agg for a single per-cell metric; return the stats rows.

    Shared by the coloc directions and the vesicle metrics so both get the same donor
    summary, the same paired test across donors, and the same title conventions. `ident`
    is merged in ahead of every computed field, so it sets the leading columns of the
    stats CSV.

    `test_override` is an {agg: (t, p)} mapping supplying the statistic to annotate and
    report instead of the one computed from the plotted values -- keyed by agg because the
    test is re-run per agg, and carrying t alongside p so the two never come from different
    tests. When it is given, the test on the plotted values is still recorded, under
    `t_stat_as_plotted` / `p_as_plotted`. `hline` draws a horizontal reference line.

    Both exist for the percent-of-control panels, which draw a transform of the values their
    reportable test was run on. `ylim` pins the y-axis to a fixed range.
    """
    sub = data[data[col].notna()]
    dropped = len(data) - len(sub)
    if dropped:
        print(f"    ({dropped} of {len(data)} cell(s) dropped: no value for {col!r})")

    # Cell counts are given per condition, not as one total. Each donor dot is the summary
    # of that donor's cells *within one condition*, so a combined total names a group that
    # no mark on the plot corresponds to.
    n_by_cond = sub.groupby("Condition", observed=True).size()
    counts = " / ".join(f"{n_by_cond.get(c, 0)} {c}" for c in CONDITIONS)

    rows = []
    for agg in aggs:
        donor_values = summarize_donors(sub, col, agg)
        grand = grand_summary(donor_values)
        plotted_t, plotted_p, n_donors = paired_ttest(donor_values)
        tstat, pval = ((plotted_t, plotted_p) if test_override is None
                       else test_override[agg])

        out_svg = os.path.join(outdir, f"{stem}_{agg}.svg")
        # Two lines: what each dot is, then what the panel is made of. Split here rather
        # than after the pair name so the lines come out close to even -- the alternative
        # leaves a very short first line over a long second one. One line of this at the
        # narrowed panel width would shrink to near-illegible under `set_fitted_title`.
        plot_swarm(
            sub, donor_values, grand, col, out_svg, palette,
            title=f"{title_prefix} — {agg} per donor\n{counts} cells, {n_donors} donors",
            ylabel=ylabel, pvalue=pval, hline=hline, ylim=ylim,
        )
        print(f"    saved {out_svg}   P = {pval:.4f}")

        # Column names stay generic (summary/n_donors rather than the agg and marker
        # names) so every pair x agg lands on the same columns and the table is
        # rectangular.
        row = {**ident, "agg": agg, "n_donors": n_donors,
               "t_stat": tstat, "p_paired": pval}
        if test_override is not None:
            row["t_stat_as_plotted"] = plotted_t
            row["p_as_plotted"] = plotted_p
        for _, r in donor_values.iterrows():
            row[f"donor_{r['Donor']}_{r['Condition']}"] = r["donor_value"]
        for _, r in grand.iterrows():
            row[f"{r['Condition']}_summary"] = r["mean"]
            row[f"{r['Condition']}_sem"] = r["sem"]
        for cond in CONDITIONS:
            row[f"n_cells_{cond}"] = int((sub["Condition"] == cond).sum())
        rows.append(row)
    return rows


# The coloc coefficient is bounded in [0, 1] by construction, so the panels are drawn on
# that full range rather than fitted to each pair's own data. The pairs sit at very
# different levels -- Fy1/Fy2 averages ~0.18 where Fy4/Fy2 averages ~0.58 -- and a
# per-panel axis silently rescales each one, so a reader comparing panels side by side
# reads the same dot height as two different values. A shared range costs some vertical
# spread on the low pairs and buys comparability across the figure.
COLOC_YLIM = (0.0, 1.0)


def process_direction(pair, data, col, X, Y, metric, aggs, outdir, palette):
    """Save one superplot per agg for a single direction, and return the stats rows."""
    # Every column of coloc_stats.csv corresponds to something a reader can point at on the
    # superplot: the donor dots, the mean+/-SEM bars, the P, the cell counts in the title.
    # The per-condition vesicle densities used to be appended here too, but they are not on
    # the figure and are reproduced exactly by the `pooled_within_condition` rows of
    # coloc_density_qc.csv, which is now a sibling folder rather than a file to hunt for.
    return superplot(
        data, col, aggs, outdir, palette,
        stem=f"coloc_{X}-to-{Y}_{metric}",
        title_prefix=f"{X}/{Y}",
        ylabel=f"Coloc. ({X}/{Y}) per cell",
        ident={"pair": pair, "direction": f"{X}/{Y}", "numerator": X,
               "denominator": Y, "metric": metric},
        ylim=COLOC_YLIM,
    )


def process_direction_percent(pair, data, col, X, Y, metric, aggs, outdir, palette):
    """The same direction, with every cell as a percentage of its donor's own M0 summary.

    Why this panel exists: the un-normalized version is dominated by a donor effect. Across
    the six donors the M0 coloc spans up to 4.5x within a single direction (Fy1/Fy2: 0.07 to
    0.30), which spreads the donor dots across the axis and leaves the M0->LPS step within
    each donor hard to see. Normalizing puts every donor's M0 on exactly 100 so the shift
    reads directly.

    The annotated P is deliberately the paired t-test on the **raw** donor means -- the same
    number the un-normalized panel shows. Normalizing would otherwise convert a difference
    test into a ratio test and change every p-value (fy1 Fy1/Fy2: 0.037 -> 0.144), leaving
    two folders reporting different results for one comparison. The difference scale is also
    the better-specified model here: across donors the LPS effect is more consistent as an
    absolute difference than as a ratio in every direction that reaches significance. The
    ratio-scale test still reaches the CSV via `test_override`, as `p_as_plotted`.
    """
    sub = data[data[col].notna()]
    # Keyed by agg: the reportable test is re-run per agg, exactly as summarize_donors is.
    override = {agg: paired_ttest(summarize_donors(sub, col, agg))[:2] for agg in aggs}

    rows = []
    for agg in aggs:
        pct_col = f"{col} (% of donor M0)"
        # assign() rather than a plain assignment: `data` is shared with the other direction,
        # the density QC and the vesicle metrics, and must not collect derived columns.
        pct = sub.assign(**{pct_col: percent_of_control(sub, col, agg)})
        rows += superplot(
            pct, pct_col, [agg], outdir, palette,
            stem=f"coloc_{X}-to-{Y}_{metric}_percent_control",
            title_prefix=f"{X}/{Y} (% of M0)",
            ylabel=f"Coloc. ({X}/{Y}), % of donor M0",
            ident={"pair": pair, "direction": f"{X}/{Y}", "numerator": X,
                   "denominator": Y, "metric": metric},
            test_override={agg: override[agg]},
            hline=100,
        )
    return rows


# ------------------------------------------------------------------
# Per-cell vesicle abundance and size
# ------------------------------------------------------------------
def primary_marker(data):
    """The pair's non-Fy2 marker: the numerator of the first coloc column in file order."""
    for col in data.columns:
        m = COLOC_RE.search(col)
        if m:
            return m.group(1)
    raise ValueError("No 'Coloc. (X/Y)' column found; cannot identify the pair's marker")


def normalized_count_col(marker):
    """Name of the derived per-cell-size vesicle count for `marker`.

    Shared by `vesicle_specs` and `process_vesicle_metrics` so the column one of them
    derives is the one the other asks for.
    """
    return f"Numb. Ves. ({marker}) per cell area"


def vesicle_specs(data, marker):
    """[(column, slug, title, ylabel)] for the per-cell abundance/size metrics.

    Three panels:

    - `Numb. Ves. (<marker>)` -- the raw segmented-object count. Extensive: a larger cell
      holds more vesicles at the same size, so a shift here can be a shift in cell size.
    - the same count divided by CELL_SIZE_COL -- intensive, and the one to read for "did
      the condition change vesicle abundance".
    - `Mean: Size : <marker>` -- mean segmented vesicle size, already intensive.

    The middle panel is labelled per cell *area*, not per cell volume, because that is what
    the denominator is: `Masked Area (Physical Units)` is the cell's masked footprint, and
    the one column that could be called a volume is not comparable between the two
    deliveries (see CELL_SIZE_COL). Within either delivery on its own the two differ only by
    a constant factor, so this is the same panel a per-volume normalization would draw, on a
    denominator that does not also encode which batch the donor came from.

    The normalized column is derived rather than read: the exporter carries the count and
    the masked area, but not their ratio.
    """
    count = f"Numb. Ves. ({marker})"
    size = f"Mean: Size : {marker}"
    for col in (count, size):
        if col not in data.columns:
            raise KeyError(f"{col!r} not in columns: {list(data.columns)}")

    # Titles stay short -- the panel is pinned to a fixed width and the y-axis label
    # already carries the full quantity, so the title only has to name the panel.
    return [
        (count, "count", f"{marker} count", f"{marker} vesicles per cell"),
        (normalized_count_col(marker), "count_per_cell_area", f"{marker} count/area",
         f"{marker} vesicles per unit cell area"),
        (size, "size", f"{marker} size", f"Mean {marker} vesicle size"),
    ]


def process_vesicle_metrics(pair, data, aggs, outdir, palette):
    """Save the abundance/size superplots for a pair's non-Fy2 marker."""
    marker = primary_marker(data)

    norm_col = normalized_count_col(marker)
    # assign() rather than a plain assignment: `data` is reused by the coloc directions
    # and the density QC, and shouldn't collect derived columns.
    data = data.assign(
        **{norm_col: data[f"Numb. Ves. ({marker})"] / data[CELL_SIZE_COL]})

    specs = vesicle_specs(data, marker)

    # `load_pair` requires CELL_SIZE_COL, so this should no longer be reachable -- it stays
    # because the failure it guards against is silent: a donor whose cells all lacked a size
    # would put a 5-donor panel beside 6-donor ones with nothing on it to say so, which is
    # exactly the quiet n-reduction this run exists to rule out. Drop the panel, and say why.
    all_donors = set(data["Donor"].unique())
    with_size = set(data.loc[data[norm_col].notna(), "Donor"].unique())
    if with_size != all_donors:
        missing = sorted(all_donors - with_size)
        print(f"  (no cell size for donor(s) {missing}; skipping {marker} count/area "
              f"so every panel keeps all {len(all_donors)} donors)")
        specs = [spec for spec in specs if spec[1] != "count_per_cell_area"]

    rows = []
    for col, slug, title, ylabel in specs:
        print(f"  {col} [{slug}]")
        rows += superplot(
            data, col, aggs, outdir, palette,
            stem=f"{pair}_{marker}_vesicle_{slug}",
            title_prefix=title,
            ylabel=ylabel,
            ident={"pair": pair, "marker": marker, "metric": slug, "column": col},
        )
    return rows


RHO_COL = "rho_coloc_vs_density_denominator"


def print_density_summary(qc):
    """Print the density QC: pooled rho per condition, then the donor-mean shift check.

    The headline table pools cells across donors, per condition. That is the right scope for
    the question "is M0 or LPS more density-coupled": rho is estimated on n ~ 120-175 instead
    of ~40-60, and both conditions are pooled the same way, so the comparison between them is
    like-for-like.

    The one thing pooling cannot be trusted for is the *sign*, because donor is a common
    cause of both variables (donor-level segmentation quality drives density, donor-level
    biology drives coloc). The `within-donor` column exists only to flag that: it is the
    Fisher-averaged per-donor rho, and where it disagrees in sign with the pooled value the
    pooled number is a between-donor artifact and should not be quoted. In this dataset that
    fires for exactly one direction.
    """
    def wide(scope, value=RHO_COL):
        return (qc[qc["scope"] == scope]
                .pivot_table(index=["pair", "direction"], columns="condition",
                             values=value, sort=False))

    rho = wide("pooled_within_condition").reindex(columns=CONDITIONS)
    pval = wide("pooled_within_condition", "p_rho").reindex(columns=CONDITIONS)
    ncell = wide("pooled_within_condition", "n_cells").reindex(columns=CONDITIONS)
    within = wide("combined_within_donor")["all"]

    # p and n are formatted as strings: `float_format` below is for the rhos, and applying it
    # to a count renders "+174.000", to a p-value a signed "+0.000" that hides the magnitude.
    table = pd.concat({
        "rho": rho,
        "p": pval.map(lambda v: "<0.001" if v < 0.001 else f"{v:.3f}"),
        "n cells": ncell.astype("Int64"),
    }, axis=1)
    table[("", "higher")] = rho.apply(
        lambda r: "" if r.isna().any() else
        ("LPS" if abs(r["LPS"]) > abs(r["M0"]) else "M0"), axis=1)
    # Sign check against the pooled-over-both-conditions rho, which is the value a reviewer
    # re-deriving this from the CSV would land on.
    pooled_all = wide("pooled_all")["all"]
    table[("", "within-donor")] = [
        f"{w:+.3f} {'SIGN FLIP' if w * p < 0 else 'agrees'}"
        for w, p in zip(within, pooled_all)
    ]

    print("\n" + "=" * 88)
    print("density QC — Spearman rho: coloc vs denominator-marker vesicle density")
    print("cells pooled across donors, within condition")
    print("=" * 88)
    print(table.to_string(float_format=lambda v: f"{v:+.3f}"))
    n_lps = (table[("", "higher")] == "LPS").sum()
    print(f"\n  |rho| larger in LPS for {n_lps} of {len(table)} directions.")

    print_donor_shift_check(qc)


def print_donor_shift_check(qc):
    """Print the LPS-minus-M0 shift in coloc and in density, per donor.

    This is the level the paired t-test actually runs on, so it -- not the cell-level rho at
    any pooling scope -- is what decides whether an M0/LPS coloc difference could be
    segmentation. A cell-level rho says the coefficient is density-sensitive; it does not say
    the conditions differed in density. Only this does.

    Read the concordance count: all 6 donors shifting coloc and density the same way means
    the coloc difference and the segmentation difference are indistinguishable in this data.
    Density moving the *opposite* way is the good case -- the coloc shift happened against
    the density gradient, so it cannot be explained by it.
    """
    w = qc[qc["scope"] == "within_donor"]
    rows = []
    for (pair, direction), g in w.groupby(["pair", "direction"], sort=False):
        p = g.pivot_table(index="donor", columns="condition",
                          values=["mean_coloc", "mean_density_denominator"])
        d_coloc = p[("mean_coloc", "LPS")] - p[("mean_coloc", "M0")]
        d_dens = (p[("mean_density_denominator", "LPS")]
                  - p[("mean_density_denominator", "M0")])
        same = int((d_coloc * d_dens > 0).sum())
        rows.append({
            "pair": pair, "direction": direction,
            "mean d_coloc": d_coloc.mean(), "mean d_density": d_dens.mean(),
            "concordant": f"{same}/{len(d_coloc)}",
            "reading": ("density could explain it" if same == len(d_coloc)
                        else "coloc shifts against density" if same == 0
                        else "mixed"),
        })

    print("\n" + "=" * 88)
    print("could the M0/LPS coloc difference be segmentation? donor-mean shifts (LPS - M0)")
    print("=" * 88)
    print(pd.DataFrame(rows).to_string(index=False,
                                       float_format=lambda v: f"{v:+.3f}"))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-root", action="append", default=None, metavar="DIR",
                        help="directory holding per-cell CSVs; repeat to combine several "
                             "deliveries into one run (donors are kept distinct per "
                             "directory). Defaults to the D1-D3 and D4-D6 batches.")
    parser.add_argument("--outdir", default="plots", help="directory for plots & CSVs")
    # Mean, not median. The per-cell distributions are near-symmetric (median skewness
    # -0.03 across the 36 donor x condition x direction strata, excess kurtosis ~0) and
    # 1.5-IQR outliers are only ~2% of cells, so there is nothing for the median to be
    # robust against -- while bootstrapping the donor summary at each donor's actual n puts
    # SE(median)/SE(mean) at 1.28, close to the sqrt(pi/2) penalty expected for normal data.
    # The median just costs precision here. `--agg both` remains available.
    parser.add_argument("--agg", choices=AGGS + ["both"], default="mean",
                        help="how to collapse cells -> one value per donor")
    parser.add_argument("--metric", choices=list(METRIC_TOKENS) + ["both"],
                        default="intensity", help="which colocalization coefficient")
    parser.add_argument("--no-vesicle-metrics", dest="vesicle_metrics",
                        action="store_false",
                        help="skip the per-cell vesicle abundance/size superplots")
    parser.add_argument("--no-percent-control", dest="percent_control",
                        action="store_false",
                        help="skip the percent-of-donor-M0 colocalization superplots")
    parser.add_argument("--no-density-qc", dest="density_qc", action="store_false",
                        help="skip the coloc-vs-vesicle-density QC")
    args = parser.parse_args()

    aggs = AGGS if args.agg == "both" else [args.agg]
    metrics = list(METRIC_TOKENS) if args.metric == "both" else [args.metric]

    roots = args.data_root or DEFAULT_DATA_ROOTS

    coloc_dir = subdir(args.outdir, COLOC_DIR)
    qc_dir = subdir(args.outdir, DENSITY_QC_DIR) if args.density_qc else None
    vesicle_dir = subdir(args.outdir, VESICLE_DIR)
    percent_dir = subdir(args.outdir, PERCENT_CONTROL_DIR) if args.percent_control else None
    pairs = discover_pairs(roots)

    # Every pair is loaded before anything is plotted, so donor labels can be assigned once
    # across the whole run rather than per pair -- see `assign_donor_labels`.
    loaded = {pair: load_pair(cells) for pair, cells in sorted(pairs.items())}
    donor_labels = assign_donor_labels(loaded.values())
    for data in loaded.values():
        data["Donor"] = data["Donor"].map(donor_labels)

    # The D-numbers are what appear on the figures, so the mapping back to the delivery and
    # the collaborator's own id is written out beside them rather than left in the console.
    key = pd.DataFrame(sorted(donor_labels.items()), columns=["source", "donor"])
    key[["batch", "raw_id"]] = key["source"].str.split(":", n=1, expand=True)
    # The acquisition id out of FileName, which is what actually names the donor outside
    # its own delivery. Taken over every pair, so a disagreement would show up as a list.
    src = (pd.concat([d[["Donor", "SourceID"]] for d in loaded.values()])
           .dropna().drop_duplicates()
           .groupby("Donor")["SourceID"].apply(lambda v: "/".join(sorted(set(v)))))
    key["file_id"] = key["donor"].map(src)
    key = key[["donor", "batch", "raw_id", "file_id"]].sort_values("donor")
    key.to_csv(os.path.join(args.outdir, "donor_key.csv"), index=False)
    print(f"\n{len(key)} donors (written to donor_key.csv):")
    print(key.to_string(index=False))

    # One palette for the whole run, so a donor keeps its color across every figure.
    palette = donor_palette(pd.DataFrame({"Donor": list(donor_labels.values())}))

    stat_rows, qc_rows, vesicle_rows, percent_rows = [], [], [], []
    for pair, data in loaded.items():
        print(f"\n=== {pair} ===")
        print(f"  {len(data)} cells, "
              f"{data.groupby('Condition', observed=True).size().to_dict()}, "
              f"donors {sorted(data['Donor'].unique())}")

        for metric in metrics:
            for col, X, Y in coloc_columns(data, metric):
                print(f"  {X}/{Y} [{metric}]")
                stat_rows += process_direction(pair, data, col, X, Y, metric, aggs,
                                               coloc_dir, palette)
                if percent_dir is not None:
                    percent_rows += process_direction_percent(
                        pair, data, col, X, Y, metric, aggs, percent_dir, palette)
                if args.density_qc and metric == "intensity":
                    rows = density_qc(pair, data, col, X, Y)
                    qc_rows += rows
                    density_scatter(pair, data, col, X, Y, rows, qc_dir, palette)

        if args.vesicle_metrics:
            vesicle_rows += process_vesicle_metrics(pair, data, aggs, vesicle_dir,
                                                    palette)

    # Each stats CSV goes in the folder holding the plots it describes, so a subdirectory
    # is self-contained: the figures and the numbers behind them travel together.
    stats_csv = os.path.join(coloc_dir, "coloc_stats.csv")
    pd.DataFrame(stat_rows).to_csv(stats_csv, index=False)
    print(f"\nsaved {stats_csv}")

    if percent_rows:
        percent_csv = os.path.join(percent_dir, "coloc_percent_control_stats.csv")
        pd.DataFrame(percent_rows).to_csv(percent_csv, index=False)
        print(f"saved {percent_csv}")

    if vesicle_rows:
        vesicle_csv = os.path.join(vesicle_dir, "vesicle_stats.csv")
        pd.DataFrame(vesicle_rows).to_csv(vesicle_csv, index=False)
        print(f"saved {vesicle_csv}")

    if qc_rows:
        qc_csv = os.path.join(qc_dir, "coloc_density_qc.csv")
        qc = pd.DataFrame(qc_rows)
        qc.to_csv(qc_csv, index=False)
        print(f"saved {qc_csv}")
        print_density_summary(qc)

    summary = pd.DataFrame(stat_rows)[
        ["pair", "direction", "metric", "agg", "n_donors", "p_paired"]
    ]
    print("\n" + "=" * 60)
    print("paired t-test across donors")
    print("=" * 60)
    print(summary.to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
