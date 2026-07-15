# goal: copy quantification files for every experiment + enrichment + diff-mod
# search to a dest directory, sorted into output_<N>diffmods/ subdirs.
#
# The same raw runs are searched several times with different max_num_diffmod
# settings (1, 2, 3, or 10). There is one upload dir per experiment+enrichment,
# and inside it one quant run per diff-mod search. Each quant run carries its
# OWN copy of search.xml right next to census-out.txt, so the diff-mod count for
# a census-out is read from the search.xml in the same quant dir:
#
#     SOURCE_ROOT/TZ3_114_FeNTA_<ts>/
#       quant/<ts1>/search.xml       # <max_num_diffmod>2</max_num_diffmod>
#       quant/<ts1>/census-out.txt
#       quant/<ts2>/search.xml       # <max_num_diffmod>1</max_num_diffmod>
#       quant/<ts2>/census-out.txt
#       ...
#
# Files land in output_<N>diffmods/ so the downstream aggregate_phospho.py
# (extract_n_diff_mods) can recover N from the path.
#
# Conda environment setup (only the standard library is used, so any recent
# Python works; these steps just pin an isolated interpreter):
#
#     conda create -n census python=3.11
#     conda activate census
#     python retrieve_census_out.py

import glob
import os
import re
import shutil
from datetime import datetime

EXPERIMENTS = [
    "TZ3_113",
    "TZ3_114",
    "TZ3_116",
    "TZ3_117",
    "TZ3_121",
]

ENRICHMENTS = [
    "TiO2",
    "FeNTA",
]

# Every kept donor + enrichment pair is searched at each of these
# max_num_diffmod settings; used to report which variants are missing.
EXPECTED_DIFFMODS = [1, 2, 3, 10]

# Base directory holding one sub-directory per experiment + enrichment.
SOURCE_ROOT = "/lustre/fs4/vino_lab/store/ip2/ip2_data/tzhang01/Macrophages/"

DESTINATION = "/ru-auth/local/home/hsanford/scratch/phospho_quants/"

LOG_FILE = os.path.join(DESTINATION, "retrieve_census_out.log")


def find_upload_dirs(experiment, enrichment):
    """Return the upload dirs for a combo, oldest first.

    Normally one per combo; the timestamp suffix is matched with a glob.
    """
    pattern = os.path.join(SOURCE_ROOT, f"{experiment}_{enrichment}_*")
    return sorted(d for d in glob.glob(pattern) if os.path.isdir(d))


def find_quant_dirs(upload_dir):
    """Return the quant run dirs in an upload dir (one per diff-mod), oldest first.

    Quant dir names carry a zero-padded timestamp, so a string sort is
    chronological and puts the newest run last.
    """
    pattern = os.path.join(upload_dir, "quant", "*")
    return sorted(d for d in glob.glob(pattern) if os.path.isdir(d))


def parse_max_num_diffmod(quant_dir):
    """Return max_num_diffmod (int) from the search.xml inside a quant dir, or None."""
    path = os.path.join(quant_dir, "search.xml")
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as handle:
        text = handle.read()
    m = re.search(r"<max_num_diffmod>\s*(\d+)\s*</max_num_diffmod>", text)
    return int(m.group(1)) if m else None


def main():
    os.makedirs(DESTINATION, exist_ok=True)

    log_lines = []

    def record(line):
        print(line)
        log_lines.append(line)

    copied = 0
    # (experiment, enrichment) -> set of diff-mod values successfully copied
    found = {}
    problems = []

    for experiment in EXPERIMENTS:
        for enrichment in ENRICHMENTS:
            combo = f"{experiment}_{enrichment}"
            found[(experiment, enrichment)] = set()

            for upload_dir in find_upload_dirs(experiment, enrichment):
                for quant_dir in find_quant_dirs(upload_dir):
                    source = os.path.join(quant_dir, "census-out.txt")
                    if not os.path.exists(source):
                        record(f"NO CENSUS: {quant_dir}")
                        problems.append(quant_dir)
                        continue

                    n_diffmod = parse_max_num_diffmod(quant_dir)
                    if n_diffmod is None:
                        record(f"NO DIFFMOD: {quant_dir}")
                        problems.append(quant_dir)
                        continue

                    dest_dir = os.path.join(DESTINATION, f"output_{n_diffmod}diffmods")
                    os.makedirs(dest_dir, exist_ok=True)
                    dest = os.path.join(dest_dir, f"{combo}-census-out.txt")
                    if os.path.exists(dest):
                        record(f"WARNING: overwriting {dest} (duplicate diffmod {n_diffmod})")
                    shutil.copy2(source, dest)
                    record(f"COPIED:  {source} -> {dest}")
                    copied += 1
                    found[(experiment, enrichment)].add(n_diffmod)

    # Report per-pair coverage and every missing (pair x diffmod) triple.
    record("")
    missing = 0
    for experiment in EXPERIMENTS:
        for enrichment in ENRICHMENTS:
            combo = f"{experiment}_{enrichment}"
            got = found[(experiment, enrichment)]
            gaps = [n for n in EXPECTED_DIFFMODS if n not in got]
            found_str = ", ".join(str(n) for n in sorted(got)) or "none"
            gap_str = f" (missing {', '.join(str(n) for n in gaps)})" if gaps else ""
            record(f"{combo}: {{{found_str}}}{gap_str}")
            missing += len(gaps)

    for experiment in EXPERIMENTS:
        for enrichment in ENRICHMENTS:
            combo = f"{experiment}_{enrichment}"
            for n in EXPECTED_DIFFMODS:
                if n not in found[(experiment, enrichment)]:
                    log_lines.append(f"MISSING: {combo} diffmod {n}")

    total_expected = len(EXPERIMENTS) * len(ENRICHMENTS) * len(EXPECTED_DIFFMODS)
    summary = (
        f"\nDone. Copied {copied}, missing {missing} of {total_expected} "
        f"(pairs x {{{', '.join(str(n) for n in EXPECTED_DIFFMODS)}}}). "
        f"Problems: {len(problems)}. Log: {LOG_FILE}"
    )
    print(summary)

    with open(LOG_FILE, "w", encoding="utf-8") as log:
        log.write(f"retrieve_census_out.py run at {datetime.now().isoformat()}\n\n")
        log.write("\n".join(log_lines))
        log.write(summary + "\n")


if __name__ == "__main__":
    main()
