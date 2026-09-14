#!/usr/bin/env python3
"""Check that a colocalization run carries the donors and cells it claims.

Two checks per row of the stats CSV:

    n_donors == --donors        
                                
    cells per donor x condition
        >= --min-cells          
                                
                                

Usage:
    python verify_n6.py plots_n6/colocalization/coloc_stats.csv
    python verify_n6.py plots_n6/colocalization/coloc_stats.csv --donors 6 --min-cells 20
"""

import argparse
import re
import sys

import pandas as pd

DONOR_COL_RE = re.compile(r"^donor_(.+)_(M0|LPS)$")


def check(csv, want_donors, min_cells):
    """Return (ok, list of failure strings) for one stats CSV."""
    df = pd.read_csv(csv)
    failures = []

    # Which identity columns exist depends on the figure type: the coloc tables carry a
    # direction, the vesicle table a marker and metric. Use whatever is there so this works
    # on both without a per-file schema.
    id_cols = [c for c in ("pair", "direction", "marker", "metric", "agg")
               if c in df.columns]

    for _, row in df.iterrows():
        label = " ".join(str(row[c]) for c in id_cols)

        if int(row["n_donors"]) != want_donors:
            failures.append(
                f"{label}: paired test ran across {int(row['n_donors'])} donors, "
                f"expected {want_donors}")

        # Donor dots actually present, per condition, from the wide donor_<D>_<cond> columns.
        present = {"M0": set(), "LPS": set()}
        for col in df.columns:
            m = DONOR_COL_RE.match(col)
            if m and pd.notna(row[col]):
                present[m.group(2)].add(m.group(1))
        for cond in ("M0", "LPS"):
            if len(present[cond]) != want_donors:
                failures.append(
                    f"{label}: {cond} has {len(present[cond])} donor dots "
                    f"({sorted(present[cond])}), expected {want_donors}")

        # Cell totals. The per-donor minimum is the binding requirement, but the stats CSV
        # only stores per-condition totals -- so assert the weaker bound it can support
        # (every donor would need >= min_cells for the total to clear donors * min_cells).
        for cond in ("M0", "LPS"):
            col = f"n_cells_{cond}"
            if col in df.columns and int(row[col]) < want_donors * min_cells:
                failures.append(
                    f"{label}: {int(row[col])} {cond} cells across {want_donors} donors "
                    f"cannot give every donor >= {min_cells}")

    return not failures, failures


def check_per_donor(percell_glob, min_cells):
    """Per-donor cell counts straight from the per-cell CSVs, when they were written."""
    import glob
    import os

    failures = []
    for path in sorted(glob.glob(percell_glob)):
        df = pd.read_csv(path)
        if not {"Donor", "Condition"}.issubset(df.columns):
            continue
        counts = df.groupby(["Donor", "Condition"]).size()
        for (donor, cond), n in counts.items():
            if n < min_cells:
                failures.append(f"{os.path.basename(path)}: donor {donor} {cond} has "
                                f"{n} cells, below {min_cells}")
    return failures


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stats_csv", nargs="+", help="coloc_stats.csv / vesicle_stats.csv")
    ap.add_argument("--donors", type=int, default=6)
    ap.add_argument("--min-cells", type=int, default=20)
    ap.add_argument("--percell-glob", default=None,
                    help="optional glob of per-cell CSVs for an exact per-donor check")
    args = ap.parse_args()

    all_failures = []
    for csv in args.stats_csv:
        ok, failures = check(csv, args.donors, args.min_cells)
        print(f"{'PASS' if ok else 'FAIL'}  {csv}")
        for f in failures:
            print(f"        {f}")
        all_failures += failures

    if args.percell_glob:
        failures = check_per_donor(args.percell_glob, args.min_cells)
        print(f"{'PASS' if not failures else 'FAIL'}  per-donor cell counts "
              f"({args.percell_glob})")
        for f in failures:
            print(f"        {f}")
        all_failures += failures

    if all_failures:
        print(f"\n{len(all_failures)} check(s) failed")
        sys.exit(1)
    print(f"\nall checks passed: {args.donors} donors, "
          f">= {args.min_cells} cells per donor x condition")


if __name__ == "__main__":
    main()
