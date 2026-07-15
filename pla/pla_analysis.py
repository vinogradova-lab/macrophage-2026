#!/usr/bin/env python3
"""Parameterized PLA region analysis.

Parses Donor/Condition (M0/LPS) from the filename, collapses per-cell rows
to one value per donor (mean or median), runs a paired t-test across donors, and saves
swarm + violin plots


Usage:
    python pla_analysis.py                                  # all folders under the data root
"""

import argparse
import glob
import os
import re

import matplotlib

matplotlib.use("Agg")  # headless: render to files, never open a window

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

CONDITIONS = ["M0", "LPS"]
DEFAULT_DATA_ROOT = os.path.join(os.path.dirname(__file__), "20260713_PLA_region_analysis")
DEFAULT_METRIC = "Ch1_Particle_Count"
EXCLUDE_DONORS = ["20260616", "20260513"]  # date-named samples dropped for the robustness variants
EXCLUDE_SLUG = "_".join(EXCLUDE_DONORS)
AGGS = ["mean", "median"]

FONT_SIZE = 6
LINE_WIDTH = 0.25

# Panel-scale geometry. Figure size and axes rect are both pinned, and figures are saved
# without bbox_inches="tight", so every plot lands with a byte-identical axes box and can
# be dropped into a figure panel side by side. The rect leaves room on the right for the
# legend, which sits outside the axes.
SWARM_FIGSIZE = (2.7, 1.7)
SWARM_RECT = dict(left=0.20, right=0.68, bottom=0.20, top=0.88)
VIOLIN_FIGSIZE = (4.0, 1.7)
VIOLIN_RECT = dict(left=0.13, right=0.83, bottom=0.20, top=0.88)
CELL_MARKER_AREA = 1.5  # per-cell dots, scatter `s` (points^2)
DONOR_MARKER_DIAM = 4.5  # donor-summary dots, swarmplot `size` (points)

# Donors are assigned these in sorted-donor order.
DONOR_COLORS = [
    "#AAAAAA",
    "#EC2427",
    "#F58420",
    "#FFCC31",
    "#4869B2",
    "#5CBED7",
    "#237D41",
    "#742D16",
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


# ------------------------------------------------------------------
# Data loading
# ------------------------------------------------------------------
def load_folder(folder_path, metric):
    """Load & parse every .xls in `folder_path` into one tidy DataFrame.

    Files are tab-delimited despite the .xls extension. Donor and Condition are parsed
    from the filename via the `<donor>_(M0|LPS)-` pattern.
    """
    files = glob.glob(os.path.join(folder_path, "*.xls"))
    assert files, f"No .xls files found in {folder_path!r}"

    rows = []
    for f in files:
        df = pd.read_csv(f, sep="\t")
        m = re.search(r"([^/\\]+?)_(M0|LPS)-", f)
        if not m:
            print(f"  Skipping (couldn't parse donor/condition): {f}")
            continue
        df["Donor"], df["Condition"] = m.group(1), m.group(2)
        rows.append(df)

    if not rows:
        raise ValueError(f"No parseable files in {folder_path!r}")

    data = pd.concat(rows, ignore_index=True)
    if metric not in data.columns:
        raise KeyError(f"Metric {metric!r} not in columns: {list(data.columns)}")
    data["Condition"] = pd.Categorical(data["Condition"], categories=CONDITIONS, ordered=True)
    return data


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


# ------------------------------------------------------------------
# Plots (saved to disk)
# ------------------------------------------------------------------
def donor_palette(data):
    """Donor -> color. Built once from the full folder so colors are stable across variants."""
    donors = sorted(data["Donor"].unique())
    if len(donors) > len(DONOR_COLORS):
        raise ValueError(
            f"{len(donors)} donors but only {len(DONOR_COLORS)} colors in DONOR_COLORS"
        )
    return dict(zip(donors, DONOR_COLORS))


def plot_swarm(data, donor_means, grand, metric, out_svg, palette, title):
    """Per-cell dots + donor summary dots + mean±SEM bars + paired-t P annotation."""
    donors = sorted(data["Donor"].unique())
    x_pos = {c: i for i, c in enumerate(CONDITIONS)}

    fig, ax = plt.subplots(figsize=SWARM_FIGSIZE)

    # small individual-cell dots, jittered, colored by donor
    rng = np.random.default_rng(0)
    for cond in CONDITIONS:
        sub = data[data["Condition"] == cond]
        for donor in donors:
            vals = sub.loc[sub["Donor"] == donor, metric]
            if len(vals) == 0:
                continue
            jitter = rng.uniform(-0.15, 0.15, size=len(vals))
            ax.scatter(
                x_pos[cond] + jitter, vals,
                color=palette[donor], alpha=0.6, s=CELL_MARKER_AREA, zorder=2,
               edgecolor="none",
            )

    # big donor-summary dots, beeswarm-packed so donors at the same value stay visible
    sns.swarmplot(
        data=donor_means, x="Condition", y="donor_value",
        order=CONDITIONS, hue="Donor", hue_order=donors, palette=palette,
        size=DONOR_MARKER_DIAM, edgecolor="black", linewidth=LINE_WIDTH,
        ax=ax, legend=False, zorder=4,
    )

    # black mean +/- SEM bars (of the donor summary values)
    for _, r in grand.iterrows():
        xc = x_pos[r["Condition"]]
        ax.plot([xc - 0.2, xc + 0.2], [r["mean"], r["mean"]], color="black",
                lw=LINE_WIDTH, zorder=5)
        ax.errorbar(xc, r["mean"], yerr=r["sem"], color="black", capsize=5,
                    lw=LINE_WIDTH, capthick=LINE_WIDTH, zorder=5)

    # paired t-test annotation
    _, pval, n = paired_ttest(donor_means)
    lo, hi = data[metric].min(), data[metric].max()
    span = (hi - lo) or 1.0
    y_top = hi
    if not np.isnan(pval):
        y_bar = hi + span * 0.06
        ax.plot([0, 1], [y_bar, y_bar], color="black", lw=LINE_WIDTH)
        ax.text(0.5, y_bar + span * 0.02, f"P = {pval:.3f}", ha="center", va="bottom")
        y_top = y_bar + span * 0.12

    # swarmplot resets the axes data limits to its own points (the donor summaries), which
    # drops every per-cell dot outside that range off the plot -- so set the range by hand.
    ax.set_ylim(lo - span * 0.05, y_top)

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
    ax.set_ylabel(f"{metric} per cell")
    ax.set_title(title)
    fig.subplots_adjust(**SWARM_RECT)
    fig.savefig(out_svg)
    plt.close(fig)


def plot_violin(data, metric, out_svg, title):
    """M0 vs LPS violins side-by-side, per donor, with jittered per-cell points."""
    donors = sorted(data["Donor"].unique())
    fig, ax = plt.subplots(figsize=VIOLIN_FIGSIZE)

    sns.violinplot(
        data=data, x="Donor", y=metric, hue="Condition",
        order=donors, hue_order=CONDITIONS,
        split=False, inner="quartile", ax=ax,
        linewidth=LINE_WIDTH,
        cut=0,  # a particle count can't go below 0; don't draw KDE tails past the data
        palette={"M0": "#4C72B0", "LPS": "#DD8452"},
    )
    sns.stripplot(
        data=data, x="Donor", y=metric, hue="Condition",
        order=donors, hue_order=CONDITIONS,
        dodge=True, jitter=True, size=0.8, color="black", alpha=0.4,
        linewidth=0, ax=ax, legend=False,
    )

    ax.set_xlabel("Donor")
    ax.set_ylabel(f"{metric} per cell")
    ax.set_title(title)
    # outside to the right: at panel size there is no empty space left inside the axes
    ax.legend(title="Condition", loc="upper left", bbox_to_anchor=(1.02, 1),
              borderaxespad=0, frameon=False, handletextpad=0.4, labelspacing=0.4)
    fig.subplots_adjust(**VIOLIN_RECT)
    fig.savefig(out_svg)
    plt.close(fig)


# ------------------------------------------------------------------
# Orchestration
# ------------------------------------------------------------------
def slugify(name):
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")


def discover_folders(data_root, requested):
    """Return list of (name, path) folders that contain .xls files."""
    if requested:
        names = requested
    else:
        names = sorted(
            d for d in os.listdir(data_root)
            if os.path.isdir(os.path.join(data_root, d))
        )
    found = []
    for name in names:
        path = os.path.join(data_root, name)
        if glob.glob(os.path.join(path, "*.xls")):
            found.append((name, path))
        else:
            print(f"  (no .xls) skipping folder: {name}")
    return found


def donor_set_variants(data):
    """(suffix, label, subset) for each donor set: all donors, then minus EXCLUDE_DONORS.

    The exclusion variant is skipped when the folder has none of EXCLUDE_DONORS, since it
    would duplicate the all-donors plots.
    """
    n_all = data["Donor"].nunique()
    variants = [("", f"all {n_all} donors", data)]

    kept = data[~data["Donor"].isin(EXCLUDE_DONORS)]
    if kept["Donor"].nunique() == n_all:
        print(f"  (no donors {', '.join(EXCLUDE_DONORS)}) skipping the exclusion plots")
    else:
        variants.append((
            f"_without_{EXCLUDE_SLUG}",
            f"without {', '.join(EXCLUDE_DONORS)} ({kept['Donor'].nunique()} donors)",
            kept,
        ))
    return variants


def process_folder(name, path, metric, aggs, outdir):
    """Load one folder and save a swarm plot per (donor set x agg) plus a violin per donor set."""
    print(f"\n=== {name} ===")
    data = load_folder(path, metric)
    print(data.groupby(["Condition", "Donor"])[metric].agg(["mean", "median", "count"]))

    slug = slugify(name)
    palette = donor_palette(data)  # from the full folder, so colors match across variants

    for suffix, label, subset in donor_set_variants(data):
        violin_svg = os.path.join(outdir, f"{slug}_violin{suffix}.svg")
        plot_violin(subset, metric, violin_svg, f"{name} — {label}")
        print(f"  saved {violin_svg}")

        for agg in aggs:
            donor_means = summarize_donors(subset, metric, agg)
            grand = grand_summary(donor_means)
            print(f"\n  [{label}, agg={agg}]")
            print(grand.to_string(index=False))

            swarm_svg = os.path.join(outdir, f"{slug}_swarm_{agg}{suffix}.svg")
            plot_swarm(subset, donor_means, grand, metric, swarm_svg,
                       palette, f"{name} — agg={agg}, {label}")
            print(f"  saved {swarm_svg}")


def _p(donor_means):
    return paired_ttest(donor_means)[1]


def compare_agg(folders, metric, outdir):
    """Part 2: paired-t p-value with mean- vs median-summarized donors, per folder.

    Also reports each p-value recomputed with donors `EXCLUDE_DONORS` removed, as a
    robustness check against those (date-named) samples.
    """
    records = []
    for name, path in folders:
        data = load_folder(path, metric)
        kept = data[~data["Donor"].isin(EXCLUDE_DONORS)]
        n = paired_ttest(summarize_donors(data, metric, "mean"))[2]
        records.append({
            "folder": name,
            "n_donors": n,
            "p_mean": _p(summarize_donors(data, metric, "mean")),
            "p_median": _p(summarize_donors(data, metric, "median")),
            f"p_mean_without_{EXCLUDE_SLUG}": _p(summarize_donors(kept, metric, "mean")),
            f"p_median_without_{EXCLUDE_SLUG}": _p(summarize_donors(kept, metric, "median")),
        })

    table = pd.DataFrame.from_records(records)
    csv_path = os.path.join(outdir, "pvalue_mean_vs_median.csv")
    table.to_csv(csv_path, index=False)

    print("\n" + "=" * 60)
    print("Part 2 — paired-t p-value: mean vs median donor summary")
    print("=" * 60)
    print(table.to_string(index=False,
                          float_format=lambda v: f"{v:.4f}" if pd.notna(v) else "nan"))
    print(f"\n  saved {csv_path}")
    return table


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-root", default=DEFAULT_DATA_ROOT,
                        help="root dir containing the per-target subfolders")
    parser.add_argument("--folder", action="append", dest="folders",
                        help="subfolder name to process (repeatable; default: all)")
    parser.add_argument("--metric", default=DEFAULT_METRIC, help="column to analyze")
    parser.add_argument("--agg", choices=AGGS + ["both"], default="both",
                        help="how to collapse cells -> one value per donor (for plots)")
    parser.add_argument("--outdir", default="plots", help="directory for saved plots & CSV")
    args = parser.parse_args()

    aggs = AGGS if args.agg == "both" else [args.agg]

    os.makedirs(args.outdir, exist_ok=True)
    folders = discover_folders(args.data_root, args.folders)
    if not folders:
        raise SystemExit(f"No folders with .xls files under {args.data_root!r}")

    for name, path in folders:
        process_folder(name, path, args.metric, aggs, args.outdir)

    compare_agg(folders, args.metric, args.outdir)


if __name__ == "__main__":
    main()
