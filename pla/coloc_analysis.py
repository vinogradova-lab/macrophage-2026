#!/usr/bin/env python3
"""Superplots for the collaborator per-cell colocalization data.

Reads the six per-cell CSVs (three antibody pairs x M0/LPS), and for each pair renders
both directions of the split colocalization coefficient as an M0-vs-LPS superplot: every
small dot a cell, every large dot a donor mean/median, p from a paired t-test across the
three donors.

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
fully -- rho reaches ~+0.6 to +0.7 in some donor x condition strata for the Fy1 and Fy4
pairs. Because density also varies by donor and by condition, some of the M0/LPS coloc
difference may be segmentation rather than biology. `density_qc()` quantifies this at four
pooling scopes (see its docstring for the pooled-vs-per-donor trade-off) into
`coloc_density_qc.csv`, and `density_scatter()` renders it per direction.

Vesicle abundance and size
--------------------------
The same superplot treatment is also applied to the per-cell vesicle metrics of each
pair's non-Fy2 marker -- raw count, count normalized to cell volume, and mean vesicle
size (see `vesicle_specs`). The normalized count needs a cell volume, which the analysis
tables do not carry; it comes from the companion `-per-cell.csv` export, joined per cell
in `load_pair`. Only rab5 has that companion file, so only rab5 gets all three panels.

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

Usage:
    python coloc_analysis.py                        # intensity metric + vesicle metrics
    python coloc_analysis.py --metric both          # also the number-based coefficient
    python coloc_analysis.py --no-vesicle-metrics   # coloc superplots only
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

DEFAULT_DATA_ROOT = (
    "/Users/henrysanford/Dropbox @RU Dropbox/Vinogradova Laboratory/"
    "Macrophage project/Manuscript 1/01_Data Analysis/08_IF/if_collaborator_data"
)

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

# The export ships two kinds of file per pair x condition, distinguished by a single
# hyphen. Read the suffixes carefully:
#   *-percell.csv    the size-gated analysis table -- one row per cell, carrying the coloc
#                    coefficients, the vesicle counts and the vesicle sizes.
#   *-per-cell.csv   the ungated companion -- seven rows per cell (a `total` row plus one
#                    per vesicle-size stratum), and the only file carrying cell volume.
# Only rab5 has a companion file; fy1 and fy4 were delivered as analysis tables alone.
VOLUME_SUFFIX = "-per-cell.csv"

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

# Cell volume lives only in the companion file, under this column (its column H).
CELL_VOLUME_COL = "Mean: Volume[Total] : Cells3D"
CELL_VOLUME = "Cell volume"

# The companion file repeats each cell once per vesicle-size stratum; `total` is the
# whole-cell row. Volume is identical across a cell's strata (verified), so selecting one
# stratum picks a representative row rather than discarding information.
TOTAL_EXTENSION = "total"


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

    df["Donor"] = clean(df["siRNA_ID"])
    df["Condition"] = clean(df["GeneName"])
    return df


def discover_pairs(data_root):
    """Group the CSVs into {pair: {"cells": {cond: path}, "volume": {cond: path}}}.

    Filenames are punctuated inconsistently (`size-gr` vs `size_gr`), so the pair name is
    taken from the `<pair>-fy2` prefix and the condition from an `_m0`/`_lps` token. Each
    pair must have both conditions of the analysis table; the `volume` companion is
    optional and present only for rab5.
    """
    paths = sorted(glob.glob(os.path.join(data_root, "*.csv")))
    if not paths:
        raise SystemExit(f"No .csv files found in {data_root!r}")

    pairs = {}
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
        role = "volume" if base.endswith(VOLUME_SUFFIX) else "cells"
        pairs.setdefault(m.group(1), {}).setdefault(role, {})[cond.group(1).upper()] = path

    for pair, files in sorted(pairs.items()):
        missing = set(CONDITIONS) - set(files.get("cells", {}))
        if missing:
            raise SystemExit(f"Pair {pair!r} is missing condition(s) {sorted(missing)}; "
                             f"found {sorted(files.get('cells', {}))}")
        # A companion for only one condition would silently give a volume-normalized plot
        # with a whole condition missing, so require both or neither.
        vol = files.get("volume", {})
        if vol and set(vol) != set(CONDITIONS):
            raise SystemExit(f"Pair {pair!r} has a companion volume file for "
                             f"{sorted(vol)} but not {sorted(set(CONDITIONS) - set(vol))}")
    return pairs


def cell_ids(series):
    """The per-cell UUID out of a FullName like ` M0:1:cel_<uuid>:ext_total`.

    The two files' FullName cannot be joined directly: they differ in the trailing `ext_`
    tag -- `ext_total` in the companion file against `ext_Size > 0.15` in the analysis
    table -- so the strings never match even for the same cell. The `cel_` UUID is the
    shared part, and is unique per cell within one condition's export.
    """
    return series.astype(str).str.extract(r"(cel_[0-9a-f-]+)", expand=False)


def load_cell_volume(path):
    """{cell UUID -> cell volume} from one companion `-per-cell.csv`."""
    df = load_csv(path)
    if CELL_VOLUME_COL not in df.columns:
        raise SystemExit(f"{os.path.basename(path)} has no {CELL_VOLUME_COL!r} column; "
                         f"it does not look like a companion volume file")
    total = df[df["Extension"].astype(str).str.strip() == TOTAL_EXTENSION]
    ids = cell_ids(total["FullName"])
    if ids.isna().any() or ids.duplicated().any():
        raise SystemExit(f"{os.path.basename(path)}: cell IDs in the {TOTAL_EXTENSION!r} "
                         f"rows are missing or not unique")
    return pd.Series(total[CELL_VOLUME_COL].values, index=ids.values)


def load_pair(files):
    """Load a pair's M0 and LPS analysis tables into one tidy frame.

    When a companion volume file is present, each cell's volume is joined on in the same
    pass. The join is done per condition and before the concat, so an M0 and an LPS cell
    can never collide on a shared UUID.
    """
    frames = []
    for cond in CONDITIONS:
        df = load_csv(files["cells"][cond])
        vol_path = files.get("volume", {}).get(cond)
        if vol_path:
            volumes = load_cell_volume(vol_path)
            ids = cell_ids(df["FullName"])
            unmatched = int((~ids.isin(volumes.index)).sum())
            if unmatched:
                raise SystemExit(
                    f"{unmatched} of {len(df)} {cond} cells have no matching row in "
                    f"{os.path.basename(vol_path)}; the two files are not the same cells")
            df[CELL_VOLUME] = ids.map(volumes).values
            # Distinct from an unmatched cell above: these cells are in both files, but
            # the 3D cell segmentation failed on them, so they carry no volume. They drop
            # out of the volume-normalized metric only.
            no_vol = int(df[CELL_VOLUME].isna().sum())
            if no_vol:
                print(f"    ({no_vol} of {len(df)} {cond} cells have no 3D cell volume)")
        frames.append(df)

    data = pd.concat(frames, ignore_index=True)
    data["Condition"] = pd.Categorical(data["Condition"], categories=CONDITIONS,
                                       ordered=True)
    if data["Condition"].isna().any():
        bad = load_csv(files["cells"]["M0"])["Condition"].unique()
        raise ValueError(f"Unexpected condition labels (expected {CONDITIONS}): {bad}")
    return data


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
    return data[f"Numb. Ves. ({marker})"] / data["Masked Area (Physical Units)"]


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
        "mean_cell_area": g.loc[ok, "Masked Area (Physical Units)"].mean(),
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
      happens for `rab5 Fy2/Rab5` (pooled -0.28, within-donor +0.14). `p_rho` is NaN for
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

    out_svg = os.path.join(outdir, f"{pair}_{X}-{Y}_coloc_vs_density_{Y}.svg")
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
        stem=f"{pair}_{X}-{Y}_coloc_{metric}",
        title_prefix=f"{X}/{Y}",
        ylabel=f"Coloc. ({X}/{Y}) per cell",
        ident={"pair": pair, "direction": f"{X}/{Y}", "numerator": X,
               "denominator": Y, "metric": metric},
        ylim=COLOC_YLIM,
    )


def process_direction_percent(pair, data, col, X, Y, metric, aggs, outdir, palette):
    """The same direction, with every cell as a percentage of its donor's own M0 summary.

    Why this panel exists: the un-normalized version is dominated by a donor batch effect
    (donor 1 > 2 > 3 in all six directions, up to 4x), which spreads the three donor dots
    across the axis and leaves the M0->LPS step within each donor hard to see. Normalizing
    puts every donor's M0 on exactly 100 so the shift reads directly.

    The annotated P is deliberately the paired t-test on the **raw** donor means -- the same
    number the un-normalized panel shows. Normalizing would otherwise convert a difference
    test into a ratio test and change every p-value (fy1 Fy1/Fy2: 0.054 -> 0.183), leaving
    two folders reporting different results for one comparison. The difference scale is also
    the better-specified model here: across donors the LPS effect is more consistent as an
    absolute difference than as a ratio in 5 of the 6 directions. The ratio-scale test still
    reaches the CSV via `test_override`, as `p_as_plotted`.
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
            stem=f"{pair}_{X}-{Y}_coloc_{metric}_percent_control",
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


def vesicle_specs(data, marker):
    """[(column, slug, title, ylabel)] for the per-cell abundance/size metrics.

    Three panels, only reached when a companion volume file was found:

    - `Numb. Ves. (<marker>)` -- the raw segmented-object count. Extensive: a larger cell
      holds more vesicles at the same density, so a shift here can be a shift in cell size.
    - the same count divided by `Cell volume` -- intensive, and the one to read for "did
      the condition change vesicle abundance". Only available with the companion file.
    - `Mean: Size : <marker>` -- mean segmented vesicle size, already intensive.

    The normalized column is derived rather than read: the companion file carries cell
    volume but not the counts, and the analysis table the reverse, so the ratio only
    exists after the join in `load_pair`.
    """
    count = f"Numb. Ves. ({marker})"
    size = f"Mean: Size : {marker}"
    for col in (count, size):
        if col not in data.columns:
            raise KeyError(f"{col!r} not in columns: {list(data.columns)}")

    # Titles stay short -- the panel is pinned to a fixed width and the y-axis label
    # already carries the full quantity, so the title only has to name the panel.
    norm = f"{count} per cell volume"
    return [
        (count, "count", f"{marker} count", f"{marker} vesicles per cell"),
        (norm, "count_per_volume", f"{marker} count/vol",
         f"{marker} vesicles per unit cell volume"),
        (size, "size", f"{marker} size", f"Mean {marker} vesicle size"),
    ]


def process_vesicle_metrics(pair, data, aggs, outdir, palette):
    """Save the abundance/size superplots for a pair's non-Fy2 marker.

    Runs only for pairs that came with a companion volume file. The raw count and size
    columns exist for every pair, but the point of this set of panels is the comparison
    between the raw and the size-normalized count, and without a cell volume that
    comparison is missing its middle term. Delivering the companion file for another pair
    is therefore what opts it in -- nothing here is keyed to a pair by name.
    """
    if CELL_VOLUME not in data.columns:
        print(f"  (no companion volume file for {pair}; skipping the vesicle metrics)")
        return []

    marker = primary_marker(data)
    specs = vesicle_specs(data, marker)

    norm_col = f"Numb. Ves. ({marker}) per cell volume"
    # assign() rather than a plain assignment: `data` is reused by the coloc directions
    # and the density QC, and shouldn't collect derived columns.
    data = data.assign(**{norm_col: data[f"Numb. Ves. ({marker})"] / data[CELL_VOLUME]})

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

    Read the concordance count: 3/3 donors shifting coloc and density the same way means the
    coloc difference and the segmentation difference are indistinguishable in this data.
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
    parser.add_argument("--data-root", default=DEFAULT_DATA_ROOT,
                        help="directory holding the six per-cell CSVs")
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
    args = parser.parse_args()

    aggs = AGGS if args.agg == "both" else [args.agg]
    metrics = list(METRIC_TOKENS) if args.metric == "both" else [args.metric]

    coloc_dir = subdir(args.outdir, COLOC_DIR)
    qc_dir = subdir(args.outdir, DENSITY_QC_DIR)
    vesicle_dir = subdir(args.outdir, VESICLE_DIR)
    percent_dir = subdir(args.outdir, PERCENT_CONTROL_DIR) if args.percent_control else None
    pairs = discover_pairs(args.data_root)

    stat_rows, qc_rows, vesicle_rows, percent_rows = [], [], [], []
    for pair, files in sorted(pairs.items()):
        print(f"\n=== {pair} ===")
        data = load_pair(files)
        print(f"  {len(data)} cells, "
              f"{data.groupby('Condition', observed=True).size().to_dict()}, "
              f"donors {sorted(data['Donor'].unique())}")
        palette = donor_palette(data)

        for metric in metrics:
            for col, X, Y in coloc_columns(data, metric):
                print(f"  {X}/{Y} [{metric}]")
                stat_rows += process_direction(pair, data, col, X, Y, metric, aggs,
                                               coloc_dir, palette)
                if percent_dir is not None:
                    percent_rows += process_direction_percent(
                        pair, data, col, X, Y, metric, aggs, percent_dir, palette)
                if metric == "intensity":
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
