#!/usr/bin/env python3
"""Does the image-processing setting itself move the colocalization estimate?

Output layout
-------------
`--outdir` holds one subfolder per comparison, each self-contained (its panels and the
numbers behind them travel together), matching `coloc_analysis.py`:

    size_gate_effect/       12 panels + coloc_settings_stats.csv
    ca_correction_effect/   12 panels + coloc_settings_stats.csv
    donor_key.csv

A panel is one direction at one state: `coloc_Fy1-to-Fy2_M0_intensity_mean.svg` is
Coloc(Fy1/Fy2) in M0 cells under the two settings. Cells are the small dots, donors the
large ones, and each donor's two summaries are joined by a line -- that line is the paired
unit the P is computed on.

Usage:
    python settings_within_state.py
    python settings_within_state.py --metric both --agg both
    python settings_within_state.py --outdir plots_settings_within_state
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import coloc_analysis as ca
from pla_plotting import (
    AGGS,
    ANNOTATION_HEADROOM,
    CONDITIONS,
    DONOR_COLORS,
    LINE_WIDTH,
    AXIS_PAD_FRAC,
    CELL_JITTER,
    CELL_MARKER_ALPHA,
    CELL_MARKER_AREA,
    CONDITION_XLIM,
    DONOR_MARKER_DIAM,
    DONOR_MARKER_EDGE,
    MEAN_BAR_HALF,
    SWARM_FIGSIZE,
    SWARM_RECT,
    grand_summary,
    paired_ttest,
    set_fitted_title,
    summarize_donors,
    TITLE_MIN_FONT_SIZE,
    _PER_UNIT_IN,  # inches per category unit; see DONOR_SPREAD below
)

# The three deliveries. Keyed rather than listed because both comparisons share the
# corrected + gated delivery as their pivot, and it must be loaded once: loading it twice
# would be wasted work, but more importantly the two comparisons would then be free to
# disagree about what the shared column contains.
DELIVERIES = {
    "uncorrected_0.15": ("no CA correction, 0.15 µm filter",
                         os.path.join(ca._IF_ROOT, "20260908_yannis_D4-D6", "0.15-filter")),
    "corrected_0.15": ("CA corrected, 0.15 µm filter",
                       os.path.join(ca._IF_ROOT, "20260910_chromaticabberation",
                                    "0.15 filter")),
    "corrected_nofilter": ("CA corrected, no filter",
                           os.path.join(ca._IF_ROOT, "20260910_chromaticabberation",
                                        "no filter")),
}

# Reference is always the setting *without* the treatment under test, so
# `delta = test - reference` reads as the effect of turning that treatment on.
#
# `axis` and `ticks` describe the panel's x-axis, and they belong to the comparison rather
# than to the deliveries either side of it. A delivery differs from its partner in exactly
# one dimension, but *which* dimension depends on which partner it is paired against:
# `corrected_0.15` is the "CA on" level in one comparison and the "gate on" level in the
# other. Labelling the axis from the delivery therefore mislabels one of the two panels --
# the size-gate panels came out reading "no filter" vs "CA corr.", naming the dimension
# being held constant instead of the one being varied.
#
# `held` names what is pinned, and goes on the axis under the variable, so a panel states
# both halves of its design rather than leaving the constant to be inferred.
COMPARISONS = [
    {"folder": "size_gate_effect",
     "reference": "corrected_nofilter", "test": "corrected_0.15",
     "short": "the 0.15 µm gate",
     "isolates": "0.15 µm size gate (CA correction on throughout)",
     "axis": "Object size gate", "ticks": ("none", "0.15 µm"),
     "held": "chromatic-aberration correction on"},
    {"folder": "ca_correction_effect",
     "reference": "uncorrected_0.15", "test": "corrected_0.15",
     "short": "CA correction",
     "isolates": "chromatic-aberration correction (0.15 µm gate on throughout)",
     "axis": "Chromatic-aberration correction", "ticks": ("off", "on"),
     "held": "0.15 µm size gate on"},
]

# The run whose donor numbering and colors this script adopts, and whose per-donor values
# the uncorrected column is checked against.
REFERENCE_RUN = "plots_n6"
# The delivery that `plots_n6/` was built from -- the one whose numbers must reproduce.
REFERENCE_DELIVERY = "uncorrected_0.15"

# `paired_ttest` is hardcoded to the two biological conditions, so the two settings are
# renamed onto them for the call. Test takes CONDITIONS[0] because the function computes
# ttest_rel(first, second): binding test to the first level makes `t_stat` carry the sign
# of `delta` (test - reference) rather than its opposite, which is checked in `main`.
DUMMY_TEST, DUMMY_REF = CONDITIONS

# The biological states, each analyzed on its own. Same strings as CONDITIONS, but a
# separate name: here they are the variable held *fixed*, not the variable under test, and
# conflating the two is the one way to get this analysis silently wrong.
STATES = CONDITIONS

# Donor dots are fanned across this much of each setting's x position so that three donors
# at similar values stay distinguishable *and* land on known coordinates -- the paired line
# has to start and end exactly on its own dots, which `sns.swarmplot`'s automatic packing
# (used by the superplots) cannot guarantee. Given in inches and converted, like every
# other horizontal dimension in `pla_plotting`.
DONOR_SPREAD_IN = 0.11
DONOR_SPREAD = DONOR_SPREAD_IN / _PER_UNIT_IN

PAIR_LINE_WIDTH = 0.4   # the donor-pairing lines; heavier than the 0.25pt axis line work
PAIR_LINE_ALPHA = 0.55

# The x-label here is two lines -- the variable under test over the setting held constant --
# where the superplots' is one, so the panel needs a taller bottom margin. The extra height
# is added to the *figure*, not taken from the axes: the axes box stays 0.899 x 1.32in,
# byte-identical to a `plots_n6/` panel, so the two figure sets still sit side by side at
# the same scale. Only the strip of white below the axes grows.
_EXTRA_BOTTOM_IN = 0.12
SETTING_FIGSIZE = (SWARM_FIGSIZE[0], SWARM_FIGSIZE[1] + _EXTRA_BOTTOM_IN)
_AXES_H_IN = SWARM_FIGSIZE[1] * (SWARM_RECT["top"] - SWARM_RECT["bottom"])
_BOTTOM_IN = SWARM_FIGSIZE[1] * SWARM_RECT["bottom"] + _EXTRA_BOTTOM_IN
SETTING_RECT = dict(
    left=SWARM_RECT["left"], right=SWARM_RECT["right"],
    bottom=_BOTTOM_IN / SETTING_FIGSIZE[1],
    top=(_BOTTOM_IN + _AXES_H_IN) / SETTING_FIGSIZE[1],
)

# ---- the summary panel ----
# Row order, pinned rather than taken from the data so a pair that failed to load leaves a
# gap instead of shifting every other row. Pairs in the manuscript's order (fy1, fy4,
# rab5), and within a pair X/Fy2 before Fy2/X, as the exporter writes them.
SUMMARY_ROW_ORDER = [("fy1", "Fy1/Fy2"), ("fy1", "Fy2/Fy1"),
                     ("fy4", "Fy4/Fy2"), ("fy4", "Fy2/Fy4"),
                     ("rab5", "Rab5/Fy2"), ("rab5", "Fy2/Rab5")]

# The manuscript's own stimulus palette (macrophage_figure_setting.R): M0 is its grey, and
# LPS is a TLR4 agonist, so it takes the TLR4 yellow. Reused rather than re-picked so a
# state is the same color here as everywhere else in the figure set.
STATE_COLORS = {"M0": "#AAAAAA", "LPS": "#FFCC31"}

# Summary geometry, in inches and pinned like every other panel: two bars per direction,
# a gap between directions, and margins sized by the text that goes in them.
BAR_H_IN, GROUP_GAP_IN = 0.115, 0.055
# Gap between the two state bars *within* one direction. Without it the bars abut, and a
# donor mark sitting near the edge of one bar reads as easily as belonging to the other --
# the marks are what carry the n, so ambiguous ownership defeats the point of drawing them.
# Kept well below GROUP_GAP_IN so the two bars still group as one direction.
STATE_GAP_IN = 0.022
# The per-donor marks drawn on each bar. Donors are told apart by SHAPE here, not by the
# color they carry on the superplots: the bar fill already encodes state with the
# manuscript's stimulus palette, and D4's donor color *is* the LPS yellow, so a colored dot
# on a colored bar would put two different meanings on one hue. Shape is free.
#
# Marks are white-filled with a black edge so they read on both the grey M0 bar and the
# yellow LPS one, and they are dodged into non-overlapping slots up the bar -- the slot
# spacing is derived from the marker size rather than guessed, so the dodge cannot silently
# start overlapping if the marker or the bar height is retuned.
SUMMARY_DOT_AREA = 2.4        # points^2
SUMMARY_DOT_GAP = 1.15        # slot pitch, as a multiple of the marker diameter
# The stack of donor marks is confined to this fraction of the bar height and centred on
# it, rather than being fanned across the whole bar. Clear space above and below the stack
# is what makes a mark unambiguously the property of one bar; filling the bar edge to edge
# made the marks read as floating between the two.
SUMMARY_STACK_FRAC = 0.68
DONOR_MARKERS = ["o", "^", "s", "D", "v", "P"]
DONOR_MARK_FILL = "#ffffff"
SUMMARY_GRID_W_IN = 1.30
SUMMARY_MARGIN_L_IN, SUMMARY_MARGIN_R_IN = 0.54, 0.50
SUMMARY_MARGIN_T_IN, SUMMARY_MARGIN_B_IN = 0.46, 0.48
# Headroom beyond the largest bar, so the longest bar does not run into the axis edge.
SUMMARY_X_PAD = 1.12
# Where the P column sits, in axes-width fractions to the right of the plot. The labels are
# a column outside the axes rather than text at each bar tip: at bar tips they collide with
# the legend on the long bars and cross the zero line on the short ones, and which of those
# happens changes with the data, so the panel could not be laid out once and trusted.
SUMMARY_P_COL_X = 1.03


# ------------------------------------------------------------------
# Loading
# ------------------------------------------------------------------
def load_delivery(root):
    """{pair: tidy per-cell frame} for one delivery, with a `Cell` column added.

    Enters the pipeline exactly where `coloc_analysis.main` does, so a value here is the
    same value the published panels plot. Two things are changed afterwards:

    `Donor` loses the batch prefix `load_pair` adds. That prefix exists to keep donors of
    different *deliveries* apart, which is the opposite of what is needed here -- the whole
    point is that `bd26-46` in one delivery and in another are the same donor, and they
    have to collide so they can be paired.

    `Cell` is the `cel_` UUID, the pairing key for the cell-level test.
    """
    pairs = ca.discover_pairs([root])
    loaded = {}
    for pair, files in sorted(pairs.items()):
        data = ca.load_pair(files)
        data["Donor"] = data["Donor"].str.split(":", n=1).str[1]
        data["Cell"] = ca.cell_ids(data["FullName"])
        if data["Cell"].isna().any() or data["Cell"].duplicated().any():
            raise SystemExit(f"{root!r} pair {pair!r}: cell UUIDs missing or not unique; "
                             f"they are the pairing key and must identify a cell exactly")
        loaded[pair] = data
    return loaded


def donor_labels_from_reference(reference_run):
    """{raw id -> D-number} from the published run's donor key, or {} if it is absent.

    Read rather than recomputed: `assign_donor_labels` numbers the donors of whatever run
    it is given, so on this three-delivery run it would return D1-D3 for the donors the
    manuscript calls D4-D6. A figure captioned D1 that is D4 elsewhere is worse than no
    D-number at all, so without the key this falls back to the collaborator's `bd26-NN`.
    """
    key_csv = os.path.join(reference_run, "donor_key.csv")
    if not os.path.exists(key_csv):
        print(f"  ({key_csv} absent -- labelling donors by their acquisition id)")
        return {}
    key = pd.read_csv(key_csv)
    return dict(zip(key["raw_id"].astype(str), key["donor"].astype(str)))


def donor_palette_from_labels(donors, labeled):
    """Donor -> color, matching the published run when the D-numbers came from its key.

    `pla_plotting.donor_palette` hands out DONOR_COLORS in sorted-donor order, which for a
    D4-D6 run means D4 takes D1's color. When the labels are D-numbers the color is picked
    by that number instead, so a donor is the same color here as on the n=6 panels.
    """
    if not labeled:
        return dict(zip(sorted(donors), DONOR_COLORS))
    return {d: DONOR_COLORS[int(d[1:]) - 1] for d in sorted(donors)}


# ------------------------------------------------------------------
# Statistics
# ------------------------------------------------------------------
def paired_setting_ttest(frame, unit, value):
    """Paired t-test of the test setting against the reference, paired over `unit`.

    Delegates to `pla_plotting.paired_ttest` -- the function that produces every P in
    `plots_n6/` -- by renaming the two settings onto its two hardcoded levels. The test,
    the sort onto matched pairs and the small-n guard are therefore identical to the
    published ones; only the meaning of the two levels differs. `unit` is "Donor" for the
    reportable test and "Cell" for the technical one, which is the only difference between
    them.

    Returns (t, p, n_pairs, delta), delta being the mean of test - reference over the
    matched pairs.
    """
    wide = frame.pivot_table(index=unit, columns="Setting", values=value)
    if not {"reference", "test"} <= set(wide.columns):
        return float("nan"), float("nan"), 0, float("nan")
    wide = wide.dropna(subset=["reference", "test"])

    t, p, n = paired_ttest(pd.DataFrame({
        "Donor": np.concatenate([wide.index.values, wide.index.values]),
        "Condition": [DUMMY_TEST] * len(wide) + [DUMMY_REF] * len(wide),
        "donor_value": np.concatenate([wide["test"].values, wide["reference"].values]),
    }))
    delta = float((wide["test"] - wide["reference"]).mean()) if len(wide) else float("nan")
    return t, p, n, delta


def direction_stats(tidy, col, state, agg, ident):
    """One stats row: both pairing levels of one direction at one state.

    The donor-level numbers come from `summarize_donors` -- the same collapse of cells to
    one value per donor that the superplots draw -- so the donor dots on the panel and the
    donor row here are the same quantity.
    """
    sub = tidy[(tidy["State"] == state) & tidy[col].notna()]

    donor_values = summarize_donors(
        sub.rename(columns={"Setting": "Condition"}), col, agg
    ).rename(columns={"Condition": "Setting"})
    d_t, d_p, d_n, d_delta = paired_setting_ttest(donor_values, "Donor", "donor_value")
    c_t, c_p, c_n, c_delta = paired_setting_ttest(sub, "Cell", col)

    ref_level = donor_values.loc[donor_values["Setting"] == "reference", "donor_value"]
    row = {**ident, "state": state, "agg": agg,
           "n_donors": d_n, "t_stat": d_t, "p_paired": d_p, "delta": d_delta,
           "delta_percent_of_reference": (100 * d_delta / ref_level.mean()
                                          if len(ref_level) and ref_level.mean() else
                                          float("nan")),
           "n_cells": c_n, "t_stat_cell": c_t, "p_paired_cell": c_p,
           "delta_cell": c_delta}

    # Per-donor values, so a reader can see the three pairs the P rests on without
    # re-running anything -- the same role the donor columns play in coloc_stats.csv.
    for _, r in donor_values.iterrows():
        row[f"donor_{r['Donor']}_{r['Setting']}"] = r["donor_value"]
    for _, r in grand_summary(
            donor_values.rename(columns={"Setting": "Condition"})).iterrows():
        row[f"{r['Condition']}_summary"] = r["mean"]
        row[f"{r['Condition']}_sem"] = r["sem"]
    return row, donor_values


# ------------------------------------------------------------------
# The panel
# ------------------------------------------------------------------
def plot_setting_pair(sub, donor_values, col, out_svg, palette, tick_labels, xlabel,
                      title, ylabel, pvalue, ylim):
    """Superplot of one direction at one state, with the two settings as the categories.

    Deliberately the same mark vocabulary as `pla_plotting.plot_swarm` -- small dots are
    cells, large outlined dots are donors, black bars are mean +/- SEM of the donor
    summaries, one P above -- at the same figsize and axes rect, so a panel from here sits
    beside a panel from `plots_n6/` without rescaling.

    Two things differ, both because this comparison is paired within a donor rather than
    between two groups of cells. Each donor's two summary dots are joined by a line: that
    line *is* the paired unit, and its slope is what the t-test is on. And the dots are
    fanned to fixed offsets instead of beeswarm-packed, because a line has to land on its
    own dots -- `plot_swarm` cannot be reused for that reason, not merely because its
    x-axis is hardcoded to M0/LPS.
    """
    order = ["reference", "test"]
    x_pos = {s: i for i, s in enumerate(order)}
    donors = sorted(donor_values["Donor"].unique())
    # Fixed, symmetric offsets: donor k of n sits at the same place at both settings, so
    # the pairing line is vertical when the setting changes nothing.
    offset = {d: (i - (len(donors) - 1) / 2) * (DONOR_SPREAD / max(len(donors) - 1, 1))
              for i, d in enumerate(donors)}

    fig, ax = plt.subplots(figsize=SETTING_FIGSIZE)

    rng = np.random.default_rng(0)
    for setting in order:
        at = sub[sub["Setting"] == setting]
        for donor in donors:
            vals = at.loc[at["Donor"] == donor, col]
            if len(vals) == 0:
                continue
            jitter = rng.uniform(-CELL_JITTER, CELL_JITTER, size=len(vals))
            ax.scatter(x_pos[setting] + jitter, vals, color=palette[donor],
                       alpha=CELL_MARKER_ALPHA, s=CELL_MARKER_AREA, zorder=2,
                       edgecolor="none")

    wide = donor_values.pivot_table(index="Donor", columns="Setting", values="donor_value")
    for donor in donors:
        if donor not in wide.index:
            continue
        pair = wide.loc[donor]
        if pair.isna().any():
            continue
        ax.plot([x_pos["reference"] + offset[donor], x_pos["test"] + offset[donor]],
                [pair["reference"], pair["test"]], color=palette[donor],
                lw=PAIR_LINE_WIDTH, alpha=PAIR_LINE_ALPHA, zorder=3)
        for setting in order:
            ax.scatter(x_pos[setting] + offset[donor], pair[setting],
                       color=palette[donor], s=DONOR_MARKER_DIAM ** 2, zorder=4,
                       edgecolor="black", linewidth=DONOR_MARKER_EDGE)

    for _, r in grand_summary(
            donor_values.rename(columns={"Setting": "Condition"})).iterrows():
        xc = x_pos[r["Condition"]]
        ax.plot([xc - MEAN_BAR_HALF, xc + MEAN_BAR_HALF], [r["mean"], r["mean"]],
                color="black", lw=LINE_WIDTH, zorder=5)
        ax.errorbar(xc, r["mean"], yerr=r["sem"], color="black", capsize=5,
                    lw=LINE_WIDTH, capthick=LINE_WIDTH, zorder=5)

    # The P annotation goes in headroom *above* the coefficient's range, not inside the top
    # of it as `plot_swarm` does. `plot_swarm` reserves the top 7% of a fixed (0, 1) axis
    # for the bracket, which works only while no cell reaches there -- in the n=6 run none
    # does, topping out at 0.898.
    #
    # The CA-corrected + 0.15 um delivery does reach there: Coloc(Fy2/Fy4) hits 0.948 and
    # Coloc(Fy2/Rab5) 0.930. Note it is the *gated* corrected delivery, not the ungated one
    # -- the ungated tops out at 0.869 and clears the bracket. So this is not a quirk of a
    # delivery nobody plots: `coloc_analysis.py --data-root <CA/0.15 filter>` dies partway
    # through with ValueError, and that is the delivery that would replace the manuscript's
    # figures if the panels move to corrected data.
    #
    # Reserving the band in-band would therefore either collide with real cells or force a
    # per-panel axis, destroying the comparability the fixed range exists for. Since the
    # coefficient is bounded in [0, 1] by construction, 0..1 stays the data range on every
    # panel, the ticks stop at 1.0, and the bracket sits in space no datum can occupy.
    #
    # `plot_swarm` now does exactly this for its own fixed-range panels, and shares the
    # constant, so the two figure sets cannot drift apart on where the bracket sits. The
    # geometry is repeated here rather than delegated only because `plot_swarm` hardcodes
    # its x-axis to M0/LPS; see this function's docstring.
    axis_span = ylim[1] - ylim[0]
    y_bar = ylim[1] + axis_span * ANNOTATION_HEADROOM * 0.30
    hi = sub[col].max()
    if hi > ylim[1]:
        raise ValueError(f"{col!r} reaches {hi:.3g}, outside the range {ylim} the "
                         f"colocalization coefficient is bounded to by construction")
    if np.isfinite(pvalue):
        ax.plot([0, 1], [y_bar, y_bar], color="black", lw=LINE_WIDTH)
        ax.text(0.5, y_bar + axis_span * 0.015, format_p(pvalue), ha="center",
                va="bottom")

    ax.set_ylim(ylim[0] - axis_span * AXIS_PAD_FRAC,
                ylim[1] + axis_span * ANNOTATION_HEADROOM)
    ax.set_yticks(np.linspace(ylim[0], ylim[1], 5))
    ax.set_xlim(*CONDITION_XLIM)

    handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=palette[d],
                          markersize=DONOR_MARKER_DIAM, label=d) for d in donors]
    ax.legend(handles=handles, title="Donor", loc="upper left", bbox_to_anchor=(1.02, 1),
              borderaxespad=0, frameon=False, handletextpad=0.4, labelspacing=0.4)

    ax.set_xticks(list(x_pos.values()))
    ax.set_xticklabels([tick_labels[s] for s in order])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.subplots_adjust(**SETTING_RECT)
    set_fitted_title(fig, ax, title)
    fig.savefig(out_svg)
    plt.close(fig)


def format_p(p):
    """The P label. `f"{p:.3f}"` alone renders a 1e-4 p-value as the false "P = 0.000"."""
    if not np.isfinite(p):
        return "n/a"
    return "P < 0.001" if p < 0.001 else f"P = {p:.3f}"


def fit_figure_text(fig, text):
    """Shrink `text` until it fits the figure width, mirroring `set_fitted_title`.

    Same failure it guards against: matplotlib does not clip overlong text, it draws past
    the canvas edge and the ends are simply lost from the saved SVG with no warning. The
    figure width is pinned so panels stay droppable into a layout side by side, so the
    text has to give rather than the canvas. `set_fitted_title` cannot be used here --
    it sets an *axes* title, and this one has to clear the legend above the axes.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    size = text.get_fontsize()
    while size > TITLE_MIN_FONT_SIZE:
        ext = text.get_window_extent(renderer)
        if ext.x0 >= 0 and ext.x1 <= fig.bbox.width:
            break
        size -= 0.25
        text.set_fontsize(size)
    return text


def plot_effect_summary(stats, short, metric, agg, palette, out_svg):
    """All twelve tests of one comparison on one axis: signed effect, not just a p-value.

    This is the panel that answers the question the analysis is actually for -- "does CA
    correction *increase* Coloc(Fy2/Fy4) in M0?" -- in one look. A p-value grid cannot: it
    says whether something moved, never which way, and the two are genuinely independent
    here, because the size gate pushes the X/Fy2 directions up and the Fy2/X directions
    down. So the encoding is the signed shift itself, `test - reference`, as a bar
    diverging from zero: right of the line the setting raises the coefficient, left of it
    lowers it, and bar length is how much.

    Direction on the y-axis and state as the bar fill, rather than the reverse, because the
    question is always asked about one direction ("Fy2/Fy4") and the two states are the
    thing being compared within it -- the M0 and LPS bars of a row sit adjacent so a shift
    that hits one state and not the other is visible as a mismatched pair.

    The units are colocalization coefficient, which the panel shares with every superplot
    in the run, so a bar can be read against the coefficient's own scale rather than as an
    abstract effect size. P (the paired t-test across donors, n = 3) rides in a column at
    the right: the direction and the magnitude come first, the significance second.

    Each bar is the mean of three donor deltas, and those three are drawn on it as dodged
    white marks, one shape per donor (see DONOR_MARKERS for why shape and not color). At
    n = 3 a bar alone is a summary of too few
    numbers to stand on its own -- it cannot distinguish three donors shifting by the same
    small amount (which is what a processing artifact looks like) from one donor carrying
    the whole effect. The dots also make the paired design visible: each is one donor's
    own test-minus-reference, the quantity the t-test is actually run on, so a bar whose
    dots all sit to one side of zero is the picture of a significant P.
    """
    sub = stats[(stats["metric"] == metric) & (stats["agg"] == agg)]
    rows = [key for key in SUMMARY_ROW_ORDER
            if ((sub["pair"] == key[0]) & (sub["direction"] == key[1])).any()]

    n_rows = len(rows)
    row_h = len(STATES) * BAR_H_IN + (len(STATES) - 1) * STATE_GAP_IN
    grid_h = n_rows * row_h + max(n_rows - 1, 0) * GROUP_GAP_IN
    grid_w = SUMMARY_GRID_W_IN
    fig_w = SUMMARY_MARGIN_L_IN + grid_w + SUMMARY_MARGIN_R_IN
    fig_h = SUMMARY_MARGIN_T_IN + grid_h + SUMMARY_MARGIN_B_IN

    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([SUMMARY_MARGIN_L_IN / fig_w, SUMMARY_MARGIN_B_IN / fig_h,
                       grid_w / fig_w, grid_h / fig_h])

    # Row 0 at the top: the rows are an ordered list to be read down the page, not a
    # quantity, so they run in reading order rather than matplotlib's bottom-up default.
    step = row_h + GROUP_GAP_IN
    # Donors, taken from the columns rather than assumed, so a run with a different donor
    # set still plots. Their per-donor deltas set the axis alongside the means: an
    # individual donor reaches further than the mean it belongs to, and would otherwise be
    # drawn outside the axes.
    donors = [c[len("donor_"):-len("_reference")] for c in sub.columns
              if c.startswith("donor_") and c.endswith("_reference")]
    donor_delta = {d: sub[f"donor_{d}_test"] - sub[f"donor_{d}_reference"] for d in donors}
    reach = [abs(sub["delta"].dropna()).max()] if sub["delta"].notna().any() else [0]
    reach += [abs(v.dropna()).max() for v in donor_delta.values() if v.notna().any()]
    limit = max(max(reach), 1e-6) * SUMMARY_X_PAD

    # Donor slots, sized so the marks cannot overlap: pitch is the marker diameter times
    # SUMMARY_DOT_GAP, and the stack of them has to fit inside one bar. y data units are
    # inches scaled by the axis range, which carries one extra GROUP_GAP_IN beyond grid_h.
    marker_diam_in = np.sqrt(SUMMARY_DOT_AREA) / 72
    slot_in = marker_diam_in * SUMMARY_DOT_GAP
    stack_in = (len(donors) - 1) * slot_in + marker_diam_in
    budget_in = BAR_H_IN * SUMMARY_STACK_FRAC
    if stack_in > budget_in:
        raise ValueError(
            f"{len(donors)} donor marks need {stack_in:.4f}in but only {budget_in:.4f}in "
            f"is allowed inside a bar (BAR_H_IN={BAR_H_IN} x "
            f"SUMMARY_STACK_FRAC={SUMMARY_STACK_FRAC}). Shrink SUMMARY_DOT_AREA, tighten "
            f"SUMMARY_DOT_GAP, or raise BAR_H_IN -- do not widen the stack into the bar "
            f"edges, which is what makes a mark's bar ambiguous.")
    y_range = (n_rows - 1) * step + row_h + GROUP_GAP_IN
    per_in = y_range / grid_h if grid_h else 1.0
    slots = (np.arange(len(donors)) - (len(donors) - 1) / 2) * slot_in * per_in

    for r, key in enumerate(rows):
        for s, state in enumerate(STATES):
            match = sub[(sub["pair"] == key[0]) & (sub["direction"] == key[1])
                        & (sub["state"] == state)]
            if len(match) != 1:
                continue
            row = match.iloc[0]
            delta, p = row["delta"], row["p_paired"]
            # Bar centre. STATE_GAP_IN separates the two state bars, so the bar itself is
            # drawn nearly full height -- the whitespace between bars comes from the gap.
            y_c = -(r * step + s * (BAR_H_IN + STATE_GAP_IN) + BAR_H_IN / 2)
            ax.barh(y_c, delta, height=BAR_H_IN * 0.94,
                    color=STATE_COLORS[state], edgecolor="black", linewidth=LINE_WIDTH,
                    zorder=3)
            # The three donor deltas behind that bar, one per dodge slot. The vertical
            # position is dodge only and carries no meaning; the horizontal one is the
            # donor's own test-minus-reference.
            idx = match.index[0]
            for donor, off in zip(donors, slots):
                value = donor_delta[donor].loc[idx]
                if not np.isfinite(value):
                    continue
                ax.scatter(value, y_c + off, s=SUMMARY_DOT_AREA,
                           marker=DONOR_MARKERS[donors.index(donor) % len(DONOR_MARKERS)],
                           facecolor=DONOR_MARK_FILL, edgecolor="black",
                           linewidth=DONOR_MARKER_EDGE * 0.6, zorder=5)
            # x in axes fractions, y in data coordinates: the label is pinned to its own
            # bar's row but to a column that does not move with the bar's length.
            ax.text(SUMMARY_P_COL_X, y_c,
                    format_p(p),
                    transform=ax.get_yaxis_transform(), ha="left", va="center", zorder=4)

    ax.axvline(0, color="black", lw=LINE_WIDTH, zorder=2)
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-((n_rows - 1) * step + row_h) - GROUP_GAP_IN / 2, GROUP_GAP_IN / 2)
    ax.set_yticks([-(r * step + row_h / 2) for r in range(n_rows)])
    ax.set_yticklabels([d for _, d in rows])
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel(f"Δ colocalization coefficient\n(with − without {short})")

    # Laid out along the top rather than in the right margin, which now belongs to the P
    # column. Two entries side by side cost one line of height and no width.
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=STATE_COLORS[s], edgecolor="black",
                             linewidth=LINE_WIDTH, label=s) for s in STATES]
    handles += [plt.Line2D([0], [0], color="w", label=d,
                           marker=DONOR_MARKERS[i % len(DONOR_MARKERS)],
                           markerfacecolor=DONOR_MARK_FILL, markeredgecolor="black",
                           markeredgewidth=DONOR_MARKER_EDGE * 0.6,
                           markersize=np.sqrt(SUMMARY_DOT_AREA) * 1.6)
                for i, d in enumerate(donors)]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0, 1.01),
              ncol=len(STATES) + len(donors), frameon=False, handlelength=1.0,
              handletextpad=0.3, columnspacing=0.8, borderaxespad=0)

    # A figure-level title, not `set_fitted_title`: that writes an axes title, which would
    # land on top of the legend now sitting above the axes. Left-aligned to the axes' left
    # edge so the title, the direction labels and the bars share one margin.
    n_donors = int(sub["n_donors"].max()) if sub["n_donors"].notna().any() else 0
    fit_figure_text(fig, fig.text(
        SUMMARY_MARGIN_L_IN / fig_w, 1 - 0.06 / fig_h,
        f"Does {short} change colocalization?\n"
        f"{agg} per donor, paired t-test, n = {n_donors} donors",
        ha="left", va="top", linespacing=1.35))
    fig.savefig(out_svg)
    plt.close(fig)


# ------------------------------------------------------------------
# Checks
# ------------------------------------------------------------------
def check_paired_cells(pair, ref, test, ref_name, test_name):
    """The two deliveries must hold exactly the same cells, or nothing here is paired.

    Equal row counts are not enough: two deliveries could carry the same number of cells
    from a different segmentation. The identity is checked on (state, donor, cell UUID),
    which is what the cell-level test pairs on, and a mismatch is fatal rather than a
    warning -- a partial overlap would still produce a p-value, just one computed on
    whichever cells happened to survive an inner join.
    """
    key = ["Condition", "Donor", "Cell"]
    a = set(map(tuple, ref[key].astype(str).values))
    b = set(map(tuple, test[key].astype(str).values))
    if a != b:
        raise SystemExit(
            f"pair {pair!r}: {ref_name} and {test_name} do not hold the same cells "
            f"({len(a)} vs {len(b)}, {len(a & b)} shared). The settings comparison is "
            f"paired cell by cell and donor by donor; it is not valid across different "
            f"segmentations.")
    return len(a)


def cross_check_reference(stats, reference_run):
    """The uncorrected column must reproduce the published run's D4-D6 donor values.

    This script reads the very files `plots_n6/` was built from, through the same loader,
    so its per-donor values for that delivery must match to the last digit. They are the
    reference side of one whole comparison, so if the loading path has drifted -- a donor
    mis-parsed, a duplicate file let through -- every `ca_correction_effect` p-value is
    wrong, and the p-values alone would not show it.

    Only the M0 and LPS donor values are checkable this way; the published run has no
    counterpart for the corrected deliveries, which is the point of this script.
    """
    stats_csv = os.path.join(reference_run, "colocalization", "coloc_stats.csv")
    if not os.path.exists(stats_csv):
        print(f"\ncross-check skipped: {stats_csv} not present")
        return True
    published = pd.read_csv(stats_csv)

    mine = stats[stats["reference_delivery"] == REFERENCE_DELIVERY]
    checked, failures = 0, []
    for _, row in mine.iterrows():
        match = published[(published["pair"] == row["pair"])
                          & (published["direction"] == row["direction"])
                          & (published["metric"] == row["metric"])
                          & (published["agg"] == row["agg"])]
        if len(match) != 1:
            continue
        ref = match.iloc[0]
        for donor in [c.split("_")[1] for c in row.index
                      if c.startswith("donor_") and c.endswith("_reference")]:
            mine_col = f"donor_{donor}_reference"
            ref_col = f"donor_{donor}_{row['state']}"
            if ref_col not in ref.index:
                continue
            a, b = row[mine_col], ref[ref_col]
            if pd.isna(a) and pd.isna(b):
                continue
            checked += 1
            if not (pd.notna(a) and pd.notna(b) and abs(a - b) < 1e-9):
                failures.append(f"{row['pair']} {row['direction']} {row['state']} "
                                f"{row['metric']}/{row['agg']} {donor}: {a} vs {b}")

    print(f"\ncross-check vs {stats_csv}")
    if failures:
        print(f"  FAIL  {len(failures)} of {checked} per-donor values differ")
        for f in failures[:10]:
            print(f"        {f}")
        return False
    print(f"  PASS  {checked} per-donor values identical to the n=6 run's D4–D6 donors")
    return True


def check_sign_convention(stats):
    """`t_stat` and `delta` must agree in sign, at both pairing levels.

    The two come from different places -- `delta` is computed here, `t_stat` inside
    `paired_ttest`, which fixes which level it subtracts from which -- so a disagreement
    means the settings were bound to the dummy levels the wrong way round and every effect
    in the CSV points backwards. Cheap to check, and silent if it ever regresses.
    """
    failures = []
    for t_col, d_col in [("t_stat", "delta"), ("t_stat_cell", "delta_cell")]:
        bad = stats[stats[t_col].notna() & stats[d_col].notna()
                    & (np.sign(stats[t_col]) != np.sign(stats[d_col]))]
        for _, r in bad.iterrows():
            failures.append(f"{r['comparison']} {r['pair']} {r['direction']} {r['state']}: "
                            f"{t_col}={r[t_col]:.3f} but {d_col}={r[d_col]:+.5f}")
    print("\nsign convention (t_stat and delta both test - reference)")
    if failures:
        print(f"  FAIL  {len(failures)} row(s) disagree")
        for f in failures[:10]:
            print(f"        {f}")
        return False
    print(f"  PASS  {len(stats)} rows")
    return True


# ------------------------------------------------------------------
# Driver
# ------------------------------------------------------------------
def run_comparison(spec, deliveries, metrics, aggs, palette, outdir):
    """Every direction x state x metric x agg of one comparison. Returns the stats rows."""
    folder, short, isolates = spec["folder"], spec["short"], spec["isolates"]
    ref_key, test_key = spec["reference"], spec["test"]
    ref_label, _ = DELIVERIES[ref_key]
    test_label, _ = DELIVERIES[test_key]
    tick_labels = dict(zip(("reference", "test"), spec["ticks"]))
    xlabel = f"{spec['axis']}\n({spec['held']})"
    out = ca.subdir(outdir, folder)

    print(f"\n=== {folder}: effect of the {isolates} ===")
    print(f"  reference  {ref_label}   [x-axis: {spec['axis']} = {spec['ticks'][0]}]")
    print(f"  test       {test_label}   [x-axis: {spec['axis']} = {spec['ticks'][1]}]")

    rows = []
    for pair in sorted(deliveries[ref_key]):
        ref, test = deliveries[ref_key][pair], deliveries[test_key][pair]
        n_cells = check_paired_cells(pair, ref, test, ref_label, test_label)
        # `Setting` carries the role, not the delivery name, so everything downstream --
        # the pivot, the panel, the sign of delta -- is written once and is the same for
        # both comparisons. The delivery names are recorded on the row instead.
        tidy = pd.concat([ref.assign(Setting="reference"), test.assign(Setting="test")],
                         ignore_index=True)
        # `Condition` is renamed to `State` rather than copied. The reused helpers
        # (`summarize_donors`, `grand_summary`, `paired_ttest`) all key on a column called
        # "Condition", and this analysis feeds them the *setting* under that name -- so a
        # surviving `Condition` holding the state would be a second column of the same name
        # after the rename, and pandas raises on the ambiguity rather than picking one.
        # Leaving exactly one of the two in the frame is what keeps "which variable is
        # under test" unambiguous for every function downstream.
        tidy = tidy.rename(columns={"Condition": "State"})
        tidy["State"] = tidy["State"].astype(str)
        print(f"  {pair}: {n_cells} cells per setting, "
              f"donors {sorted(tidy['Donor'].unique())}")

        for metric in metrics:
            for col, X, Y in ca.coloc_columns(ref, metric):
                ident = {"comparison": folder, "isolates": isolates, "short_name": short,
                         "reference_delivery": ref_key, "test_delivery": test_key,
                         "reference_setting": ref_label,
                         "test_setting": test_label,
                         "pair": pair, "direction": f"{X}/{Y}", "numerator": X,
                         "denominator": Y, "metric": metric}
                for state in STATES:
                    for agg in aggs:
                        row, donor_values = direction_stats(tidy, col, state, agg, ident)
                        rows.append(row)

                        svg = os.path.join(
                            out, f"coloc_{X}-to-{Y}_{state}_{metric}_{agg}.svg")
                        sub = tidy[(tidy["State"] == state) & tidy[col].notna()]
                        plot_setting_pair(
                            sub, donor_values, col, svg, palette, tick_labels,
                            xlabel,
                            # Names the direction, the state and which setting is under
                            # test, then the signed effect -- so the panel answers "does
                            # <setting> raise Coloc(X/Y) in <state>?" without its folder.
                            title=(f"{X}/{Y} in {state} — effect of {short}\n"
                                   f"Δ = {row['delta']:+.4f} "
                                   f"({row['delta_percent_of_reference']:+.1f}%), "
                                   f"{row['n_cells']} cells, {row['n_donors']} donors"),
                            ylabel=f"Coloc. ({X}/{Y}) per cell",
                            pvalue=row["p_paired"], ylim=ca.COLOC_YLIM)
                        print(f"    {X}/{Y} {state} [{metric}/{agg}]  "
                              f"Δ = {row['delta']:+.4f}  P = {row['p_paired']:.4f}  "
                              f"(cells P = {row['p_paired_cell']:.2e})")

    stats = pd.DataFrame(rows)
    for metric in metrics:
        for agg in aggs:
            svg = os.path.join(out, f"settings_effect_summary_{metric}_{agg}.svg")
            plot_effect_summary(stats, short, metric, agg, palette, svg)
            print(f"  saved {svg}")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", default="plots_settings_within_state",
                        help="directory for the panels & per-comparison stats CSVs")
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

    # Only the deliveries the comparisons actually name, so an unused entry in DELIVERIES
    # cannot fail a run.
    needed = {spec[k] for spec in COMPARISONS for k in ("reference", "test")}
    deliveries = {}
    for key in sorted(needed):
        label, root = DELIVERIES[key]
        print(f"loading {label}\n  {root}")
        deliveries[key] = load_delivery(root)

    labels = donor_labels_from_reference(args.reference_run)
    donors = sorted({d for pairs in deliveries.values() for f in pairs.values()
                     for d in f["Donor"].unique()})
    unmapped = [d for d in donors if d not in labels]
    if labels and unmapped:
        raise SystemExit(f"{args.reference_run}/donor_key.csv does not name {unmapped}; "
                         f"it maps {sorted(labels)}")
    if labels:
        for pairs in deliveries.values():
            for f in pairs.values():
                f["Donor"] = f["Donor"].map(labels)
        donors = sorted(labels[d] for d in donors)
    palette = donor_palette_from_labels(donors, bool(labels))

    key = pd.DataFrame({"donor": donors})
    key["raw_id"] = key["donor"].map({v: k for k, v in labels.items()}) if labels else donors
    key.to_csv(os.path.join(args.outdir, "donor_key.csv"), index=False)
    print(f"\n{len(key)} donors (written to {args.outdir}/donor_key.csv):")
    print(key.to_string(index=False))

    all_rows = []
    for spec in COMPARISONS:
        rows = run_comparison(spec, deliveries, metrics, aggs, palette, args.outdir)
        csv = os.path.join(args.outdir, spec["folder"], "coloc_settings_stats.csv")
        pd.DataFrame(rows).to_csv(csv, index=False)
        print(f"  saved {csv}")
        all_rows += rows

    stats = pd.DataFrame(all_rows)
    ok = check_sign_convention(stats)
    if args.cross_check:
        ok = cross_check_reference(stats, args.reference_run) and ok

    for metric in metrics:
        for agg in aggs:
            print("\n" + "=" * 78)
            print(f"paired t-test across donors, setting vs setting at a fixed state "
                  f"({metric}, {agg})")
            print("=" * 78)
            sub = stats[(stats["metric"] == metric) & (stats["agg"] == agg)]
            print(sub[["comparison", "pair", "direction", "state", "n_donors", "delta",
                       "delta_percent_of_reference", "p_paired", "n_cells",
                       "p_paired_cell"]]
                  .to_string(index=False,
                             float_format=lambda v: f"{v:.4f}" if abs(v) >= 1e-4 else
                             f"{v:.2e}"))
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
