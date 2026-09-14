#!/usr/bin/env python3
"""How much of the colocalization result depends on the image-processing settings?

The D4-D6 donors were delivered three times, quantified three ways:

    20260908_yannis_D4-D6/0.15-filter         no chromatic-aberration correction, objects
                                              below 0.15 um excluded
    20260910_chromaticabberation/0.15 filter  corrected, objects below 0.15 um excluded
    20260910_chromaticabberation/no filter    corrected, no size gate

Holding the correction on, does the 0.15 um gate change the result? Holding the gate on,
does the correction? There is no fourth delivery -- uncorrected and ungated -- so the design
is three corners of a 2x2, not four.

"The result" is what `coloc_analysis.py` already reports: a paired t-test of M0 against LPS
across donors, on each donor's mean per-cell colocalization coefficient. This script re-runs
that test once per delivery and draws it two ways -- p as a heatmap, and the LPS effect size
as bars. The effect size is the more robust of the two: p also moves with n and with donor
spread, so it wobbles under settings that barely touch the quantity being estimated.

Every number comes from `coloc_analysis` and `pla_plotting`; nothing statistical is
reimplemented here.

Only D4-D6 can be compared this way. D1-D3 has a `no filter/` folder too, but those are the
ungated *companions* -- seven rows per cell, one per size stratum, carrying cell volume and
no coloc coefficients -- not a second analysis table. So every column is the same three
donors and every p is n=3, which is not comparable to the n=6 p-values in `plots_n6/` and is
not meant to be. What is comparable is the per-donor values; `cross_check` asserts those.

Donors are named D4-D6 throughout, matching the n=6 run and the within-state comparison. The
collaborator's `bd26-NN` ids are the internal join key and are recorded in `donor_key.csv`;
`apply_donor_numbers` explains why the D-numbers cannot come from `assign_donor_labels`.

Verified before writing: for every pair x condition the three deliveries carry exactly the
same cells -- identical `cel_` UUID sets, not merely equal counts. Counts are pair-specific
(fy1 135 M0 / 103 LPS, fy4 132/132, rab5 134/119), so a cell total only means anything
within a pair. `check_columns` re-asserts the count half of that on every run.

Usage:
    python filter_settings_comparison.py
    python filter_settings_comparison.py --metric both --agg both
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle

import coloc_analysis as ca
# Imported for the module-level rcParams as much as the functions: these figures sit beside
# the superplots in a panel, so they inherit their 6pt type and 0.25pt line work.
from pla_plotting import AGGS, CONDITIONS, LINE_WIDTH, paired_ttest, summarize_donors

# (label, correction level, gate level, root), reference delivery first then one setting
# changed at a time. The two levels are what the heatmap header rows show: a reader can see
# which factor moved without diffing two phrases.
VARIANTS = [
    ("no CA correction / 0.15 µm filter", "no", "0.15 µm",
     os.path.join(ca._IF_ROOT, "20260908_yannis_D4-D6", "0.15-filter")),
    ("CA corrected / 0.15 µm filter", "yes", "0.15 µm",
     os.path.join(ca._IF_ROOT, "20260910_chromaticabberation", "0.15 filter")),
    ("CA corrected / no filter", "yes", "none",
     os.path.join(ca._IF_ROOT, "20260910_chromaticabberation", "no filter")),
]
FACTOR_ROWS = [("chromatic aberration", 1), ("vesicle size gate", 2)]
SLOT_LABELS = ["no CA, 0.15 µm", "CA, 0.15 µm", "CA, no gate"]

REFERENCE_VARIANT = VARIANTS[0][0]
REFERENCE_RUN = "plots_n6"

# Fixed rather than taken from the data, so a pair failing to load leaves a visible hole
# instead of silently renumbering the rows.
ROW_ORDER = [("fy1", "Fy1/Fy2"), ("fy1", "Fy2/Fy1"),
             ("fy4", "Fy4/Fy2"), ("fy4", "Fy2/Fy4"),
             ("rab5", "Rab5/Fy2"), ("rab5", "Fy2/Rab5")]

INK, INK_ON_FILL, SURFACE, MUTED = "#0b0b0b", "#ffffff", "#ffffff", "#52514e"

# --- heatmap ----------------------------------------------------------------------------
# Geometry is pinned and the figures are saved without bbox_inches="tight", the way
# pla_plotting pins the superplots, so a panel lands at the same size every run.
# "0.0559" at 6pt is ~0.28in, so 0.46in is the narrowest a cell can be. Height follows:
# square-ish cells stop the grid reading as horizontal bands, which drags the eye along rows
# when the across-column comparison is the point.
CELL_W_IN, CELL_H_IN = 0.46, 0.42
MARGIN_L_IN, MARGIN_T_IN, MARGIN_B_IN = 0.92, 0.74, 0.22
CBAR_GAP_IN, CBAR_W_IN, CBAR_LABEL_IN = 0.08, 0.10, 0.30
HEADER_ROW_IN = 0.13
CELL_GAP_PT = 1.0  # surface-colored edge, so cells read as separate marks

# Sequential blue, steps 100->700. One hue, monotonic in lightness: p is a magnitude, and a
# multi-hue ramp would invent categories. Dark = small p.
BLUE_RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
             "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
CMAP = LinearSegmentedColormap.from_list("coloc_p", BLUE_RAMP)
INK_FLIP = 0.55  # above this ramp fraction the fill is too dark for black type

# Color is on -log10(p), which spreads 0.2-0.01 instead of packing it into a tenth of a
# linear axis. Bounds fixed, not fitted, so a cell means the same thing in every panel.
P_MIN, P_MAX = 0.01, 1.0
CBAR_TICKS = [1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01]

# --- effect-size bars --------------------------------------------------------------------
# Same mark vocabulary as the within-state summary in `settings_within_state.py` -- bars from
# zero, donor marks dodged up the bar, a P column outside the axes -- so the two halves of
# the settings question read side by side. Written out rather than imported so the two
# scripts stay independent.
#
# No palette is introduced for the deliveries, which is the same call the within-state panels
# make: a setting is not something the manuscript has a color for, so it is carried by
# POSITION (a labelled slot per direction) and hue is left to mean what it means elsewhere.
# The bar is the LPS response, so it takes the manuscript's LPS color -- LPS is a TLR4
# agonist, and #FFCC31 is the TLR4 yellow in macrophage_figure_setting.R.
BAR_COLOR = "#FFCC31"
BAR_H_IN, BAR_GAP_IN, GROUP_GAP_IN = 0.115, 0.022, 0.075
EFF_GRID_W_IN = 1.30
EFF_MARGIN_L_IN, EFF_MARGIN_R_IN = 1.06, 0.52
EFF_MARGIN_T_IN, EFF_MARGIN_B_IN = 0.44, 0.44
EFF_X_PAD, EFF_P_COL_X = 1.12, 1.03
# Donors are told apart by SHAPE, not color: the bar already carries the LPS yellow, and D4's
# donor color on the superplots *is* that yellow. White fill with a black edge reads both on
# the bar and off it, since a donor mark can land either side of zero.
DONOR_MARKERS = ["o", "^", "s", "D", "v", "P"]
DONOR_FILL, DONOR_SIZE = "#ffffff", 3.2
STACK_FRAC = 0.68  # fraction of bar height the donor marks are dodged within


def donor_names(df):
    """Donor ids carried as `donor_<id>_<condition>` columns."""
    suffixes = tuple(f"_{c}" for c in CONDITIONS)
    return sorted({c.split("_")[1] for c in df.columns
                   if c.startswith("donor_") and c.endswith(suffixes)})


def p_to_norm(p):
    """p-value -> 0..1 position on the ramp. NaN passes through for a blank cell."""
    if not np.isfinite(p):
        return np.nan
    lo, hi = -np.log10(P_MAX), -np.log10(P_MIN)
    return float(np.clip((-np.log10(max(p, 1e-12)) - lo) / (hi - lo), 0.0, 1.0))


def format_p(p):
    """Three significant figures, not three decimals.

    Three decimals rounds 0.0497 and 0.0503 to the same "0.050", collapsing exactly the
    distinction a reader checking against 0.05 is looking for.
    """
    if not np.isfinite(p):
        return "n/a"
    return "<0.001" if p < 0.001 else f"{p:#.3g}"


# ------------------------------------------------------------------
# Statistics: one delivery at a time, through the existing pipeline
# ------------------------------------------------------------------
def variant_stats(label, root, metrics, aggs):
    """Run one delivery through `coloc_analysis` and return its paired-t rows.

    Entered at the same point `coloc_analysis.main` enters it, so a column here and the
    corresponding panel there are the same computation.

    `assign_donor_labels` is called for its check, not its labels: one delivery must yield
    exactly three donors, and anything else means a root was mis-parsed.
    """
    pairs = ca.discover_pairs([root])
    loaded = {pair: ca.load_pair(files) for pair, files in sorted(pairs.items())}

    donor_labels = ca.assign_donor_labels(loaded.values())
    if len(donor_labels) != 3:
        raise SystemExit(f"{label!r}: expected 3 donors in {root!r}, found "
                         f"{len(donor_labels)}: {sorted(donor_labels)}")
    # "<batch>:bd26-46" -> "bd26-46"; unique within one delivery, which is all a single-root
    # run holds. See `apply_donor_numbers` for why the D-numbering is deferred.
    raw_ids = {d: d.split(":", 1)[1] for d in donor_labels}
    for data in loaded.values():
        data["Donor"] = data["Donor"].map(raw_ids)

    rows = []
    for pair, data in loaded.items():
        for metric in metrics:
            for col, X, Y in ca.coloc_columns(data, metric):
                sub = data[data[col].notna()]
                for agg in aggs:
                    donor_values = summarize_donors(sub, col, agg)
                    tstat, pval, n_donors = paired_ttest(donor_values)
                    row = {"variant": label, "root": root, "pair": pair,
                           "direction": f"{X}/{Y}", "numerator": X, "denominator": Y,
                           "metric": metric, "agg": agg, "n_donors": n_donors,
                           "t_stat": tstat, "p_paired": pval}
                    for _, r in donor_values.iterrows():
                        row[f"donor_{r['Donor']}_{r['Condition']}"] = r["donor_value"]
                    for cond in CONDITIONS:
                        row[f"n_cells_{cond}"] = int((sub["Condition"] == cond).sum())
                    rows.append(row)
    return rows


def apply_donor_numbers(stats_df, reference_run):
    """Rename `donor_bd26-46_M0` -> `donor_D4_M0`; return (renamed, key frame).

    The load path keys on `bd26-NN` because that is the only id stable across runs:
    `assign_donor_labels` renumbers from 1 within whatever roots it was given, so a
    single-delivery run of D4-D6 calls its donors D1-D3. Emitting that would put a
    `donor_D1_M0` column beside `plots_n6`'s `donor_D1_M0`, which is a different person.

    So the D-numbers come from the published run's key, after which every donor column here
    matches the n=6 table by name. Without that key the raw ids are kept and the caller is
    told, rather than the CSV schema changing silently with what happens to be on disk.
    """
    raw = donor_names(stats_df)
    key_csv = os.path.join(reference_run, "donor_key.csv")
    if not os.path.exists(key_csv):
        print(f"  ! {key_csv} not found; donor columns keep the collaborator's ids "
              f"({', '.join(raw)}) instead of D-numbers")
        return stats_df, pd.DataFrame({"donor": raw, "raw_id": raw})

    published = pd.read_csv(key_csv)
    mapping = dict(zip(published["raw_id"].astype(str), published["donor"].astype(str)))
    missing = [r for r in raw if r not in mapping]
    if missing:
        raise SystemExit(f"{key_csv} has no D-number for {missing}; it does not describe "
                         f"the same donors as {REFERENCE_VARIANT!r}")

    renamed = stats_df.rename(columns={f"donor_{r}_{cond}": f"donor_{mapping[r]}_{cond}"
                                       for r in raw for cond in CONDITIONS})
    key = pd.DataFrame({"donor": [mapping[r] for r in raw], "raw_id": raw})
    return renamed, key.sort_values("donor").reset_index(drop=True)


# ------------------------------------------------------------------
# The p-value heatmap
# ------------------------------------------------------------------
def heatmap(stats, metric, agg, out_svg):
    """One p-value grid: rows are pair x direction, columns are deliveries."""
    sub = stats[(stats["metric"] == metric) & (stats["agg"] == agg)]
    donors = donor_names(stats)
    variants = [v[0] for v in VARIANTS]
    grid = (sub.pivot_table(index=["pair", "direction"], columns="variant",
                            values="p_paired", sort=False)
            .reindex(index=ROW_ORDER, columns=variants))

    n_rows, n_cols = len(ROW_ORDER), len(variants)
    grid_w, grid_h = CELL_W_IN * n_cols, CELL_H_IN * n_rows
    fig_w = MARGIN_L_IN + grid_w + CBAR_GAP_IN + CBAR_W_IN + CBAR_LABEL_IN
    fig_h = MARGIN_T_IN + grid_h + MARGIN_B_IN

    fig = plt.figure(figsize=(fig_w, fig_h), facecolor=SURFACE)
    ax = fig.add_axes([MARGIN_L_IN / fig_w, MARGIN_B_IN / fig_h,
                       grid_w / fig_w, grid_h / fig_h])
    ax.set_facecolor(SURFACE)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xlim(0, n_cols)
    ax.set_ylim(n_rows, 0)  # row 0 at the top, reading order
    ax.set_xticks([])
    ax.set_yticks([])

    for r, key in enumerate(ROW_ORDER):
        for c, variant in enumerate(variants):
            p = grid.loc[key, variant] if key in grid.index else np.nan
            t = p_to_norm(p)
            # The gap is an edge on the fill rather than a shrunken rectangle, so cells stay
            # exactly on their grid coordinates.
            ax.add_patch(Rectangle((c, r), 1, 1, edgecolor=SURFACE, linewidth=CELL_GAP_PT,
                                   facecolor=SURFACE if not np.isfinite(t) else CMAP(t),
                                   zorder=1))
            ax.text(c + 0.5, r + 0.5, format_p(p), ha="center", va="center", zorder=4,
                    color=INK if (not np.isfinite(t) or t < INK_FLIP) else INK_ON_FILL)

    for r, (_, direction) in enumerate(ROW_ORDER):
        ax.text(-0.06, r + 0.5, direction, ha="right", va="center",
                transform=ax.get_yaxis_transform(), color=INK)

    # Header: one row per factor, name in the left margin, level over each column. Placed in
    # axes coordinates off the top edge so it tracks the grid at any cell size.
    row_h = HEADER_ROW_IN / grid_h
    for i, (factor, idx) in enumerate(FACTOR_ROWS):
        y = 1 + (len(FACTOR_ROWS) - i - 0.5) * row_h
        ax.text(-0.06, y, factor, ha="right", va="center", transform=ax.transAxes,
                color=MUTED)
        for c, variant in enumerate(VARIANTS):
            ax.text((c + 0.5) / n_cols, y, variant[idx], ha="center", va="center",
                    transform=ax.transAxes, color=INK)

    # Neither fact is inferable from the grid: the test is M0 vs LPS (the columns are
    # deliveries, not conditions), and n is 3, not the 6 of the published run.
    fig.text(0.02, 1 - 0.08 / fig_h,
             f"M0 vs LPS paired t-test\n"
             f"colocalization ({metric}, {agg} per donor)\n"
             f"n = {len(donors)} donors ({donors[0]}–{donors[-1]}) in every column",
             ha="left", va="top", color=INK, linespacing=1.35)
    fig.text(0.02, 0.04 / fig_h, "cell = p-value", ha="left", va="bottom", color=MUTED)

    cax = fig.add_axes([(MARGIN_L_IN + grid_w + CBAR_GAP_IN) / fig_w, MARGIN_B_IN / fig_h,
                        CBAR_W_IN / fig_w, grid_h / fig_h])
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0, 1), cmap=CMAP), cax=cax)
    bar.outline.set_visible(False)
    bar.set_ticks([p_to_norm(p) for p in CBAR_TICKS])
    bar.set_ticklabels([f"{p:g}" for p in CBAR_TICKS])
    cax.tick_params(length=1.5, width=LINE_WIDTH, colors=INK, pad=1.5)
    # Above the bar, not beside it: as a y-label it lands mid-axis right of the tick numbers
    # and reads as a stray letter in the margin.
    cax.set_title("p", color=INK, pad=3)

    fig.savefig(out_svg)
    plt.close(fig)
    return grid


# ------------------------------------------------------------------
# The effect-size figure
# ------------------------------------------------------------------
def donor_effects(row, donors):
    """[(donor, LPS - M0)] for the donors this row actually carries.

    Returns the donor alongside its value, never a bare list: the donor decides both the
    marker shape and the dodge slot, so a row missing one donor must not shift the others
    onto the wrong shape.

    LPS - M0, not the M0 - LPS that `paired_ttest` hands to ttest_rel: the sign of a p-value
    is not on the figure, but the sign of an effect is, and "what LPS did" is the direction
    a reader expects.
    """
    out = []
    for d in donors:
        v = row.get(f"donor_{d}_LPS", np.nan) - row.get(f"donor_{d}_M0", np.nan)
        if np.isfinite(v):
            out.append((d, float(v)))
    return out


def effect_plot(stats_df, metric, agg, out_svg):
    """The LPS effect, three deliveries per direction, as bars from zero.

    The companion to the heatmap and the more robust of the two, for the reason in the module
    docstring. Every donor is drawn on its bar, per the superplot convention the rest of the
    figure set uses: at n=3 a bar alone hides whether the donors agree, and on the Rab5 rows
    they emphatically do not.
    """
    sub = stats_df[(stats_df["metric"] == metric) & (stats_df["agg"] == agg)]
    donors = donor_names(stats_df)
    slot_of = {d: i for i, d in enumerate(donors)}

    n_series = len(VARIANTS)
    group_h = n_series * BAR_H_IN + (n_series - 1) * BAR_GAP_IN + GROUP_GAP_IN
    axes_h = len(ROW_ORDER) * group_h
    fig_w = EFF_MARGIN_L_IN + EFF_GRID_W_IN + EFF_MARGIN_R_IN
    fig_h = EFF_MARGIN_T_IN + axes_h + EFF_MARGIN_B_IN

    # Collected before drawing so the axis can be fitted to the donor marks as well as the
    # bars -- a dissenting donor otherwise lands outside the panel -- and so the tick labels
    # come from the bars that exist rather than from an assumed full grid.
    bars, extremes = [], [0.0]
    for gi, key in enumerate(ROW_ORDER):
        for si, variant in enumerate(VARIANTS):
            match = sub[(sub["pair"] == key[0]) & (sub["direction"] == key[1])
                        & (sub["variant"] == variant[0])]
            if len(match) != 1:
                continue
            row = match.iloc[0]
            effects = donor_effects(row, donors)
            if not effects:
                continue
            values = [v for _, v in effects]
            y = gi * group_h + BAR_H_IN / 2 + si * (BAR_H_IN + BAR_GAP_IN)
            bars.append((y, float(np.mean(values)), effects, row["p_paired"],
                         SLOT_LABELS[si]))
            extremes += [float(np.mean(values)), *values]
    xmax = max(abs(v) for v in extremes) * EFF_X_PAD
    if not xmax:
        # Nothing matched -- a metric/agg the table does not carry. Keep the axis
        # non-degenerate so the empty panel is legible rather than a matplotlib warning.
        print(f"  ! no {metric}/{agg} rows to plot; {out_svg} will be empty")
        xmax = 1.0

    fig = plt.figure(figsize=(fig_w, fig_h), facecolor=SURFACE)
    ax = fig.add_axes([EFF_MARGIN_L_IN / fig_w, EFF_MARGIN_B_IN / fig_h,
                       EFF_GRID_W_IN / fig_w, axes_h / fig_h])
    ax.set_facecolor(SURFACE)
    ax.spines["left"].set_visible(False)
    ax.set_ylim(axes_h, 0)
    ax.set_xlim(min(-xmax, 0), xmax)
    ax.axvline(0, color=INK, lw=LINE_WIDTH, zorder=5)

    dodge = BAR_H_IN * STACK_FRAC / max(len(donors) - 1, 1)
    for y, mean, effects, pval, _ in bars:
        ax.barh(y, mean, height=BAR_H_IN, color=BAR_COLOR, edgecolor=INK,
                linewidth=LINE_WIDTH, zorder=2)
        for donor, v in effects:
            i = slot_of[donor]
            ax.plot(v, y - BAR_H_IN * STACK_FRAC / 2 + i * dodge,
                    marker=DONOR_MARKERS[i % len(DONOR_MARKERS)], markersize=DONOR_SIZE,
                    markerfacecolor=DONOR_FILL, markeredgecolor=INK,
                    markeredgewidth=LINE_WIDTH, linestyle="none", zorder=6)
        ax.text(EFF_P_COL_X, y, f"P = {format_p(pval)}", ha="left", va="center",
                transform=ax.get_yaxis_transform(), color=INK)

    # Two label levels: the delivery slot on the axis, the direction outside it.
    ax.set_yticks([b[0] for b in bars])
    ax.set_yticklabels([b[4] for b in bars], color=MUTED)
    ax.tick_params(axis="y", length=0, pad=1.5)
    for gi, (_, direction) in enumerate(ROW_ORDER):
        centre = gi * group_h + (n_series * BAR_H_IN + (n_series - 1) * BAR_GAP_IN) / 2
        ax.text(-(EFF_MARGIN_L_IN - 0.06) / EFF_GRID_W_IN, centre, direction,
                ha="left", va="center", transform=ax.get_yaxis_transform(), color=INK)

    ax.tick_params(axis="x", length=1.5, width=LINE_WIDTH, pad=1.5)
    ax.set_xlabel("Δ colocalization coefficient (LPS − M0)")
    fig.text(0.02, 1 - 0.07 / fig_h,
             f"Does the processing setting change the LPS effect?\n"
             f"{metric}, {agg} per donor, paired t-test, "
             f"n = {len(donors)} donors ({', '.join(donors)})",
             ha="left", va="top", color=INK, linespacing=1.35)

    # Donors only: the bars carry no categorical hue to explain, and the delivery is named on
    # the axis rather than in a key.
    for donor in donors:
        i = slot_of[donor]
        ax.plot([], [], marker=DONOR_MARKERS[i % len(DONOR_MARKERS)],
                markersize=DONOR_SIZE, markerfacecolor=DONOR_FILL, markeredgecolor=INK,
                markeredgewidth=LINE_WIDTH, linestyle="none", label=donor)
    leg = ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.005), ncol=len(donors),
                    frameon=False, handletextpad=0.3, columnspacing=0.9, borderpad=0.0)
    for text in leg.get_texts():
        text.set_color(INK)

    fig.savefig(out_svg)
    plt.close(fig)


# ------------------------------------------------------------------
# Checks
# ------------------------------------------------------------------
def cross_check(stats, reference_run):
    """Assert the reference column reproduces the D4-D6 donors of the n=6 run.

    The p-values cannot match -- that run is n=6, this one n=3 -- but the per-donor values
    behind them must, because the reference column reads the very files the n=6 run read. A
    mismatch means the loading path diverged, which the p-values alone would not reveal.
    """
    stats_csv = os.path.join(reference_run, "colocalization", "coloc_stats.csv")
    if not os.path.exists(stats_csv):
        print(f"\ncross-check skipped: {stats_csv} not present")
        return True

    published = pd.read_csv(stats_csv)
    # Both tables are on D-numbers (see `apply_donor_numbers`), so same-named columns are
    # the same donor by construction.
    donors = donor_names(stats)
    mine = stats[stats["variant"] == REFERENCE_VARIANT]

    checked, failures = 0, []
    for _, row in mine.iterrows():
        match = published[(published["pair"] == row["pair"])
                          & (published["direction"] == row["direction"])
                          & (published["metric"] == row["metric"])
                          & (published["agg"] == row["agg"])]
        if len(match) != 1:
            continue
        ref = match.iloc[0]
        for donor in donors:
            for cond in CONDITIONS:
                col = f"donor_{donor}_{cond}"
                if col not in row or col not in ref:
                    continue
                a, b = row[col], ref[col]
                if pd.isna(a) and pd.isna(b):
                    continue
                checked += 1
                if not (pd.notna(a) and pd.notna(b) and abs(a - b) < 1e-9):
                    failures.append(f"{row['pair']} {row['direction']} "
                                    f"{row['metric']}/{row['agg']} {donor} {cond}: "
                                    f"{a} vs {b}")

    print(f"\ncross-check vs {stats_csv}")
    if failures:
        print(f"  FAIL  {len(failures)} of {checked} per-donor values differ")
        for f in failures[:10]:
            print(f"        {f}")
        return False
    print(f"  PASS  {checked} per-donor values identical to the n=6 run's "
          f"{'/'.join(donors)}")
    return True


def check_columns(stats):
    """Every column is the same cells, so the cell counts must agree across deliveries."""
    failures = [f"{r['variant']} {r['pair']} {r['direction']}: n_donors = {r['n_donors']}, "
                f"expected 3"
                for _, r in stats[stats["n_donors"] != 3].iterrows()]
    for (pair, metric, agg, direction), g in stats.groupby(["pair", "metric", "agg",
                                                            "direction"]):
        for cond in CONDITIONS:
            counts = set(g[f"n_cells_{cond}"])
            if len(counts) > 1:
                failures.append(f"{pair} {direction} {metric}/{agg} {cond}: cell counts "
                                f"differ across deliveries ({sorted(counts)}) -- the "
                                f"deliveries should hold the same cells")

    print("\ncolumn sanity")
    if failures:
        print(f"  FAIL  {len(failures)} problem(s)")
        for f in failures[:10]:
            print(f"        {f}")
        return False
    print("  PASS  every column n=3 donors, cell counts identical across deliveries")
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", default="plots_filter_settings",
                        help="directory for the figures & stats CSV")
    parser.add_argument("--metric", choices=list(ca.METRIC_TOKENS) + ["both"],
                        default="intensity", help="which colocalization coefficient")
    parser.add_argument("--agg", choices=AGGS + ["both"], default="mean",
                        help="how to collapse cells -> one value per donor")
    parser.add_argument("--reference-run", default=REFERENCE_RUN,
                        help="run directory supplying the donor key and the cross-check")
    parser.add_argument("--no-cross-check", dest="cross_check", action="store_false",
                        help="skip the comparison against the published run")
    args = parser.parse_args()

    metrics = list(ca.METRIC_TOKENS) if args.metric == "both" else [args.metric]
    aggs = AGGS if args.agg == "both" else [args.agg]
    os.makedirs(args.outdir, exist_ok=True)

    rows = []
    for label, _corrected, _gate, root in VARIANTS:
        print(f"\n=== {label} ===\n  {root}")
        rows += variant_stats(label, root, metrics, aggs)
    stats, donor_key = apply_donor_numbers(pd.DataFrame(rows), args.reference_run)
    donors = donor_names(stats)

    key_csv = os.path.join(args.outdir, "donor_key.csv")
    donor_key.to_csv(key_csv, index=False)
    print(f"\nsaved {key_csv}")
    print(donor_key.to_string(index=False))

    stats_csv = os.path.join(args.outdir, "coloc_filter_settings_stats.csv")
    stats.to_csv(stats_csv, index=False)
    print(f"saved {stats_csv}")

    for metric in metrics:
        for agg in aggs:
            stem = f"by_filter_setting_{metric}_{agg}.svg"
            grid = heatmap(stats, metric, agg,
                           os.path.join(args.outdir, f"coloc_p_{stem}"))
            effect_plot(stats, metric, agg,
                        os.path.join(args.outdir, f"coloc_effect_size_{stem}"))
            print(f"saved coloc_p_{stem} and coloc_effect_size_{stem}")
            print(f"\np-values, {metric} / {agg} "
                  f"(M0 vs LPS, n={len(donors)} donors {donors[0]}–{donors[-1]}):")
            print(grid.round(4).to_string())

    ok = check_columns(stats)
    if args.cross_check:
        ok = cross_check(stats, args.reference_run) and ok
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
