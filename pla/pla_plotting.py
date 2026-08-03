#!/usr/bin/env python3
"""Shared style, statistics and superplot rendering for the PLA / colocalization scripts.

Both `pla_analysis.py` (per-donor .xls region data) and `coloc_analysis.py` (per-cell
colocalization CSVs) import from here so the two figure sets share one visual style and
one statistical treatment.

The "superplot" convention used throughout: every small dot is a cell, colored by donor;
every large outlined dot is one donor's summary; the black bar is mean +/- SEM *of the
donor summaries*; and the p-value is a paired t-test across donors, never across cells.
"""

import re

import matplotlib

matplotlib.use("Agg")  # headless: render to files, never open a window

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

CONDITIONS = ["M0", "LPS"]
AGGS = ["mean", "median"]

FONT_SIZE = 6
LINE_WIDTH = 0.25

# Panel-scale geometry. Figure size and axes rect are both pinned, and figures are saved
# without bbox_inches="tight", so every plot lands with a byte-identical axes box and can
# be dropped into a figure panel side by side. The rect leaves room on the right for the
# legend, which sits outside the axes.
# Widths are chosen in inches, then expressed as fractions of the figure. The axes holds
# only two categories, so a wide box just pushes M0 and LPS apart with nothing between them:
# 0.90in (was 1.30in) keeps the two clouds close enough to compare at a glance. The left
# margin stays at 0.54in for the y-label and ticks, and the right margin at 0.86in for the
# legend, which sits outside the axes -- both were already at their minimum, so the figure
# narrows by exactly what the axes gives up. `top` leaves room for a two-line title.
# The vertical margins are sized by their text -- two title lines above, tick labels and the
# x-label below -- so they are held at a fixed 0.34in each and any added height goes entirely
# to the axes, which is 1.32in tall here.
SWARM_FIGSIZE = (2.3, 2.0)
SWARM_RECT = dict(left=0.235, right=0.626, bottom=0.17, top=0.83)
VIOLIN_FIGSIZE = (4.0, 1.7)
VIOLIN_RECT = dict(left=0.13, right=0.83, bottom=0.20, top=0.88)
CELL_MARKER_AREA = 1.5  # per-cell dots, scatter `s` (points^2)
# Per-cell dots are the background of a superplot, not its subject: the donor summaries and
# the mean+/-SEM bars are what the reader is meant to compare, and at 120-175 cells per
# condition an opaque cloud competes with them. Kept high enough that a lone cell in a sparse
# region is still visible -- below ~0.2 the tails start dropping out of the figure.
CELL_MARKER_ALPHA = 0.3
# Horizontal layout. The two conditions sit at x = 0 and x = 1 by definition, 1.0 category
# unit apart no matter what, so what sets their separation on the page is how much of the
# fixed axes box that 0..1 span occupies: matplotlib's categorical default of (-0.5, 1.5)
# gives it half the box, and a wider range packs the conditions toward the middle with the
# slack going to the outer edges.
#
# Both dimensions below are given in inches and converted to category units, rather than
# written as category units directly. Widening the range squeezes a fixed number of units
# into fewer inches, so a hand-written jitter silently thins the cell clouds every time the
# tick gap is retuned -- and the same for the mean+/-SEM bar. Deriving them keeps the panel
# looking identical apart from the one dimension being changed.
CONDITION_TICK_GAP_IN = 0.39  # M0 to LPS, centre to centre
CELL_CLOUD_WIDTH_IN = 0.135   # full width of one condition's jittered cell cloud
MEAN_BAR_WIDTH_IN = 0.18      # full width of the mean +/- SEM bar

_AXES_W_IN = SWARM_FIGSIZE[0] * (SWARM_RECT["right"] - SWARM_RECT["left"])
_X_UNITS = _AXES_W_IN / CONDITION_TICK_GAP_IN  # category units spanning the axes box
_PER_UNIT_IN = _AXES_W_IN / _X_UNITS

CONDITION_XLIM = (0.5 - _X_UNITS / 2, 0.5 + _X_UNITS / 2)
CELL_JITTER = CELL_CLOUD_WIDTH_IN / _PER_UNIT_IN / 2  # half-width, in category units
MEAN_BAR_HALF = MEAN_BAR_WIDTH_IN / _PER_UNIT_IN / 2
# Breathing room below the bottom of the y-axis, as a fraction of the axis range. Markers
# are drawn centered on their value, so a point sitting at the very bottom is half cut off
# by the axis line -- the coloc coefficients reach 0.003, which is visually zero. Applied
# identically whether the axis is pinned by `ylim` or fitted to the data, so every panel
# gets the same proportional gap and they stay comparable side by side.
AXIS_PAD_FRAC = 0.05
DONOR_MARKER_DIAM = 4.5  # donor-summary dots, swarmplot `size` (points)
# Outline width for those dots. Deliberately not LINE_WIDTH: that is the axis-spine weight,
# and at 0.25pt the black ring is too fine to separate a donor dot from the cell cloud behind
# it. The ring is what makes the summary read as a different kind of mark rather than just a
# larger cell, so it is set independently of the line work.
DONOR_MARKER_EDGE = 0.5
# Density-scatter dots. Larger than CELL_MARKER_AREA because these encode condition in the
# marker fill, and filled-vs-open is not distinguishable at 1.5 points^2.
SCATTER_MARKER_AREA = 3.0

# Donors are assigned these in sorted-donor order. These are the manuscript's stimulus
# palette (macrophage_figure_setting.R), reordered: #AAAAAA is the manuscript's M0 color,
# so leading with it painted donor 1 in the color that means "M0" everywhere else -- on a
# plot whose x-axis is M0 vs LPS. Grey is kept last, as a fallback for an 8th donor.
DONOR_COLORS = [
    "#4869B2",
    "#EC2427",
    "#F58420",
    "#FFCC31",
    "#5CBED7",
    "#237D41",
    "#742D16",
    "#AAAAAA",
]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": FONT_SIZE,
    "axes.titlesize": FONT_SIZE,
    "axes.labelsize": FONT_SIZE,
    "xtick.labelsize": FONT_SIZE,
    "ytick.labelsize": FONT_SIZE,
    "legend.fontsize": FONT_SIZE,
    "legend.title_fontsize": FONT_SIZE,
    "axes.spines.top": False,  # x and y axis lines only, no bounding box
    "axes.spines.right": False,
    "axes.linewidth": LINE_WIDTH,
    "lines.linewidth": LINE_WIDTH,
    "patch.linewidth": LINE_WIDTH,
    "grid.linewidth": LINE_WIDTH,
    "xtick.major.width": LINE_WIDTH,
    "ytick.major.width": LINE_WIDTH,
    "xtick.minor.width": LINE_WIDTH,
    "ytick.minor.width": LINE_WIDTH,
})


def slugify(name):
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")


# ------------------------------------------------------------------
# Summaries & statistics
# ------------------------------------------------------------------
def summarize_donors(data, metric, agg):
    """Collapse per-cell rows to one value per Donor x Condition using `agg`."""
    return (
        data.groupby(["Condition", "Donor"], observed=True)[metric]
        .agg(agg)
        .reset_index()
        .rename(columns={metric: "donor_value"})
    )


def grand_summary(donor_means):
    """Per-condition mean / SEM / n across the donor summary values (the black bars)."""
    return (
        donor_means.groupby("Condition", observed=True)["donor_value"]
        .agg(mean="mean", sem=lambda x: x.std(ddof=1) / np.sqrt(len(x)), n="count")
        .reset_index()
    )


def paired_ttest(donor_means):
    """Paired t-test across donors present in both conditions.

    Returns (tstat, pval, n_donors). If fewer than 2 paired donors, returns
    (nan, nan, n_donors).
    """
    common = set(donor_means.loc[donor_means["Condition"] == "M0", "Donor"]) & set(
        donor_means.loc[donor_means["Condition"] == "LPS", "Donor"]
    )
    n = len(common)
    if n < 2:
        return float("nan"), float("nan"), n

    m0 = (
        donor_means[(donor_means.Condition == "M0") & (donor_means.Donor.isin(common))]
        .sort_values("Donor")["donor_value"]
        .values
    )
    lps = (
        donor_means[(donor_means.Condition == "LPS") & (donor_means.Donor.isin(common))]
        .sort_values("Donor")["donor_value"]
        .values
    )
    tstat, pval = stats.ttest_rel(m0, lps)
    return tstat, pval, n


def spearman_rho(x, y):
    """Spearman rho and its two-sided p-value, or (nan, nan) when n < 3.

    Wrapped rather than called inline so every scope of the density QC gets the same
    small-n guard: below n = 3 scipy returns nan behind a warning, and the QC pools cells
    at several group sizes.
    """
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if len(x) < 3:
        return float("nan"), float("nan")
    res = stats.spearmanr(x, y)
    return float(res.statistic), float(res.pvalue)


def fisher_combine_rho(rhos, ns):
    """Average Spearman rhos across independent strata. Returns (rho, n_strata_used).

    Averaging rho directly under-weights the strong strata, because rho is compressed as it
    approaches +/-1. So the strata are averaged after the Fisher transform z = arctanh(rho),
    which is variance-stabilizing, with weights n - 3 (the variance of z for a Spearman rho
    is approximately 1/(n - 3)), then transformed back with tanh.

    Strata with n <= 3 or |rho| == 1 are dropped: the first has no usable weight, the second
    sends z to infinity.
    """
    rhos, ns = np.asarray(rhos, dtype=float), np.asarray(ns, dtype=float)
    ok = np.isfinite(rhos) & (ns > 3) & (np.abs(rhos) < 1)
    if not ok.any():
        return float("nan"), 0
    z = np.average(np.arctanh(rhos[ok]), weights=ns[ok] - 3)
    return float(np.tanh(z)), int(ok.sum())


# ------------------------------------------------------------------
# Plots (saved to disk)
# ------------------------------------------------------------------
def donor_palette(data):
    """Donor -> color. Built once from the full dataset so colors are stable across variants."""
    donors = sorted(data["Donor"].unique())
    if len(donors) > len(DONOR_COLORS):
        raise ValueError(
            f"{len(donors)} donors but only {len(DONOR_COLORS)} colors in DONOR_COLORS"
        )
    return dict(zip(donors, DONOR_COLORS))


def plot_density_scatter(data, x_col, y_col, out_svg, palette, title, xlabel, ylabel,
                         annotation=None):
    """One dot per cell: vesicle density on x, colocalization on y, colored by donor.

    Condition is encoded as the marker fill -- M0 filled, LPS open -- rather than split into
    two panels, so that a condition shift along x (segmentation density) and a shift along y
    (coloc) are visible against each other in one panel. That comparison is the point of the
    figure: if the LPS cloud sits to the right *and* higher, the coloc difference is
    confounded with how densely the two conditions segmented.

    A single least-squares trend line is drawn over all points -- both conditions, all donors
    -- so it describes the same cells the pooled r in `annotation` is computed on. Note the
    mismatch: the line is a linear fit while Spearman's r is rank-based, so the line shows the
    direction and rough steepness of the trend, not the quantity being reported.
    """
    donors = sorted(data["Donor"].unique())
    fig, ax = plt.subplots(figsize=SWARM_FIGSIZE)

    for cond in CONDITIONS:
        sub = data[data["Condition"] == cond]
        filled = cond == "M0"
        for donor in donors:
            g = sub[sub["Donor"] == donor]
            if g.empty:
                continue
            ax.scatter(
                g[x_col], g[y_col], s=SCATTER_MARKER_AREA, alpha=0.7, zorder=2,
                facecolor=palette[donor] if filled else "none",
                edgecolor="none" if filled else palette[donor],
                linewidth=0 if filled else LINE_WIDTH,
            )

    # Least-squares trend line over all points, spanning only the observed x-range so it
    # doesn't imply the fit extends past the data. Drawn under the annotation (zorder 6) and
    # over the points, so it stays visible inside the dense part of the cloud.
    xy = data[[x_col, y_col]].to_numpy(dtype=float)
    xy = xy[np.isfinite(xy).all(axis=1)]
    if len(xy) >= 2 and xy[:, 0].min() < xy[:, 0].max():
        slope, intercept = np.polyfit(xy[:, 0], xy[:, 1], 1)
        xs = np.array([xy[:, 0].min(), xy[:, 0].max()])
        ax.plot(xs, slope * xs + intercept, color="#333333", lw=LINE_WIDTH * 2, zorder=5)

    if annotation:
        # Translucent backing: which corner is empty differs per panel, so the text is
        # always placed top-left and made to survive whatever points sit under it.
        ax.text(0.02, 0.98, annotation, transform=ax.transAxes, ha="left", va="top",
                linespacing=1.4, zorder=6,
                bbox=dict(facecolor="white", alpha=0.75, edgecolor="none",
                          boxstyle="square,pad=0.25"))

    # Two legends stacked outside the axes on the right: donor colors at the top, the
    # fill convention at the bottom. The first has to be re-added as an artist, or the
    # second ax.legend() call replaces it.
    donor_handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=palette[d],
                   markersize=DONOR_MARKER_DIAM, label=d)
        for d in donors
    ]
    cond_handles = [
        plt.Line2D([0], [0], marker="o", color="w", label=cond,
                   markerfacecolor="#444444" if cond == "M0" else "none",
                   markeredgecolor="#444444", markeredgewidth=LINE_WIDTH * 2,
                   markersize=DONOR_MARKER_DIAM)
        for cond in CONDITIONS
    ]
    leg_kw = dict(borderaxespad=0, frameon=False, handletextpad=0.4, labelspacing=0.4)
    donor_leg = ax.legend(handles=donor_handles, title="Donor", loc="upper left",
                          bbox_to_anchor=(1.02, 1), **leg_kw)
    ax.add_artist(donor_leg)
    ax.legend(handles=cond_handles, title="Condition", loc="lower left",
              bbox_to_anchor=(1.02, 0), **leg_kw)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.subplots_adjust(**SWARM_RECT)
    fig.savefig(out_svg)
    plt.close(fig)


# Floor for the shrink-to-fit below. Under ~4 pt a title is illegible in print, so a
# title that still doesn't fit is left overflowing rather than silently made unreadable.
TITLE_MIN_FONT_SIZE = 4.0


def set_fitted_title(fig, ax, title):
    """Set `title`, stepping its font size down until it fits inside the figure.

    matplotlib does not clip an overlong title -- it draws past the canvas edge and the
    ends are simply lost from the saved SVG, with no warning. Panels are pinned to a fixed
    figsize so they can be dropped into a figure side by side, so widening the figure or
    wrapping to a second line are both out: either would leave one panel a different size
    from its neighbours. Measuring the rendered extent is the only reliable check, since
    character count is a poor proxy once the string mixes digits, letters and em dashes.
    """
    text = ax.set_title(title)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    size = text.get_fontsize()
    while size > TITLE_MIN_FONT_SIZE:
        ext = text.get_window_extent(renderer)
        # The title is centered on the axes, which is inset asymmetrically (the legend
        # takes the right margin), so check both edges rather than the width alone.
        if ext.x0 >= 0 and ext.x1 <= fig.bbox.width:
            break
        size -= 0.25
        text.set_fontsize(size)
    return text


def plot_swarm(data, donor_means, grand, metric, out_svg, palette, title, ylabel=None,
               pvalue=None, hline=None, ylim=None):
    """Per-cell dots + donor summary dots + mean±SEM bars + paired-t P annotation.

    `ylabel` defaults to "<metric> per cell"; pass it explicitly when the column name is
    too unwieldy for an axis label (e.g. the colocalization columns).

    `pvalue` overrides the annotation, which otherwise comes from `paired_ttest` on
    `donor_means`. Needed when the plotted values are a transform of the ones the reportable
    test was run on -- the percent-of-control panels plot ratios but annotate the paired test
    on the raw donor means, so that they and the un-normalized panels cannot disagree.

    `hline` draws a horizontal reference line at that y value, under everything else. Used
    for the 100% baseline on those same panels.

    `ylim` pins the axis to a fixed (lo, hi) instead of fitting it to the data, so panels of
    different magnitude stay comparable side by side. Note this changes where the P
    annotation goes: by default it is drawn *above* the data and the axis is grown to fit it,
    which a fixed range cannot do, so the top of the range is reserved for it instead. Only
    pass `ylim` for a quantity with a meaningful fixed range -- the caller is asserting the
    data cannot leave it.
    """
    donors = sorted(data["Donor"].unique())
    x_pos = {c: i for i, c in enumerate(CONDITIONS)}

    fig, ax = plt.subplots(figsize=SWARM_FIGSIZE)

    if hline is not None:
        ax.axhline(hline, color="black", lw=LINE_WIDTH, alpha=0.35, zorder=1)

    # small individual-cell dots, jittered, colored by donor
    rng = np.random.default_rng(0)
    for cond in CONDITIONS:
        sub = data[data["Condition"] == cond]
        for donor in donors:
            vals = sub.loc[sub["Donor"] == donor, metric]
            if len(vals) == 0:
                continue
            jitter = rng.uniform(-CELL_JITTER, CELL_JITTER, size=len(vals))
            ax.scatter(
                x_pos[cond] + jitter, vals,
                color=palette[donor], alpha=CELL_MARKER_ALPHA, s=CELL_MARKER_AREA, zorder=2,
               edgecolor="none",
            )

    # big donor-summary dots, beeswarm-packed so donors at the same value stay visible
    sns.swarmplot(
        data=donor_means, x="Condition", y="donor_value",
        order=CONDITIONS, hue="Donor", hue_order=donors, palette=palette,
        size=DONOR_MARKER_DIAM, edgecolor="black", linewidth=DONOR_MARKER_EDGE,
        ax=ax, legend=False, zorder=4,
    )

    # black mean +/- SEM bars (of the donor summary values)
    for _, r in grand.iterrows():
        xc = x_pos[r["Condition"]]
        ax.plot([xc - MEAN_BAR_HALF, xc + MEAN_BAR_HALF], [r["mean"], r["mean"]],
                color="black", lw=LINE_WIDTH, zorder=5)
        ax.errorbar(xc, r["mean"], yerr=r["sem"], color="black", capsize=5,
                    lw=LINE_WIDTH, capthick=LINE_WIDTH, zorder=5)

    # paired t-test annotation
    _, computed, n = paired_ttest(donor_means)
    pval = computed if pvalue is None else pvalue
    lo, hi = data[metric].min(), data[metric].max()
    span = (hi - lo) or 1.0
    if ylim is None:
        axis_lo, axis_hi = lo - span * AXIS_PAD_FRAC, None  # top grown to fit the P below
        y_bar, text_gap = hi + span * 0.06, span * 0.02
    else:
        # Reserve the top of the fixed range for the annotation rather than growing the axis.
        axis_span = ylim[1] - ylim[0]
        axis_lo, axis_hi = ylim[0] - axis_span * AXIS_PAD_FRAC, ylim[1]
        y_bar, text_gap = ylim[0] + axis_span * 0.93, axis_span * 0.015
        if hi > y_bar:
            raise ValueError(
                f"{metric!r} reaches {hi:.3g}, which collides with the P annotation at "
                f"{y_bar:.3g}; ylim={ylim} is too tight for this data")

    y_top = hi
    if not np.isnan(pval):
        ax.plot([0, 1], [y_bar, y_bar], color="black", lw=LINE_WIDTH)
        ax.text(0.5, y_bar + text_gap, f"P = {pval:.3f}", ha="center", va="bottom")
        y_top = y_bar + span * 0.12

    # swarmplot resets the axes data limits to its own points (the donor summaries), which
    # drops every per-cell dot outside that range off the plot -- so set the range by hand.
    ax.set_ylim(axis_lo, y_top if axis_hi is None else axis_hi)
    # After swarmplot, which resets the x range to its own categorical default.
    ax.set_xlim(*CONDITION_XLIM)

    # legend (donor colors), outside plot to the right
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=palette[d],
                   markersize=DONOR_MARKER_DIAM, label=d)
        for d in donors
    ]
    ax.legend(handles=handles, title="Donor", loc="upper left",
              bbox_to_anchor=(1.02, 1), borderaxespad=0, frameon=False,
              handletextpad=0.4, labelspacing=0.4)

    ax.set_xticks(list(x_pos.values()))
    ax.set_xticklabels(list(x_pos.keys()))
    ax.set_xlabel("Condition")
    ax.set_ylabel(ylabel if ylabel is not None else f"{metric} per cell")
    # subplots_adjust first: the fit check measures against the final axes position.
    fig.subplots_adjust(**SWARM_RECT)
    set_fitted_title(fig, ax, title)
    fig.savefig(out_svg)
    plt.close(fig)
