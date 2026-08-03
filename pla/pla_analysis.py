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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from pla_plotting import (
    AGGS,
    CONDITIONS,
    LINE_WIDTH,
    VIOLIN_FIGSIZE,
    VIOLIN_RECT,
    donor_palette,
    grand_summary,
    paired_ttest,
    plot_swarm,
    slugify,
    summarize_donors,
)

DEFAULT_DATA_ROOT = os.path.join(os.path.dirname(__file__), "20260713_PLA_region_analysis")
DEFAULT_METRIC = "Ch1_Particle_Count"
EXCLUDE_DONORS = ["20260616", "20260513"]  # date-named samples dropped for the robustness variants
EXCLUDE_SLUG = "_".join(EXCLUDE_DONORS)


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
# Plots (saved to disk)
# ------------------------------------------------------------------
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
    # stripplot's jitter draws from numpy's *global* RNG, so without this the point
    # positions differ on every run and the figure isn't reproducible. plot_swarm seeds
    # its own generator for the same reason.
    np.random.seed(0)
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
            # Two lines, matching the coloc panels, split to keep them close to even.
            plot_swarm(subset, donor_means, grand, metric, swarm_svg,
                       palette, f"{name} — agg={agg}\n{label}")
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
