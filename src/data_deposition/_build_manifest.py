"""Emit sample_manifest.tsv: one row per submitted sample.

Regenerate with ``python _build_manifest.py``.  The md5 columns are filled later, by
``load_md5sums`` in build_geo_metadata.py.
"""

from pathlib import Path

# Condition suffix -> (source label, GEO stimulus, agonist, dose, vendor/cat, DMSO carrier).
CONDITIONS = {
    1: ("M0 Low", "none", "none", "none", "", "0.025% DMSO"),
    2: ("TLR1-2", "TLR1-2", "CU-T12-09", "10 µM", "Millipore Sigma 5325830001", "0.025% DMSO"),
    3: ("TLR3", "TLR3", "Poly(I:C)", "20 µg/mL", "Tocris Bioscience 42-871-0", "0.025% DMSO"),
    4: ("TLR4", "TLR4", "LPS", "100 ng/mL", "Millipore Sigma L2630-10", "0.025% DMSO"),
    5: ("TLR7", "TLR7", "Imiquimod", "10 µg/mL", "Fisher Scientific AAJ63990MC", "1% DMSO"),
    6: ("TLR8", "TLR8", "CL075", "0.5 µg/mL", "Millipore Sigma SML2235", "0.025% DMSO"),
    7: ("TLR9", "TLR9", "CpG ODN 2006", "5 µg/mL", "Integrated DNA Technologies, custom", "0.025% DMSO"),
    8: ("STING", "STING", "cGAMP", "10 µg/mL", "Invivogen tlrl-nacga23-5", "0.025% DMSO"),
    9: ("M0 High", "none", "none", "none", "", "1% DMSO"),
}
DONORS = ("TZ10", "TZ11", "TZ12")  # TZ3 excluded: only 8 of 9 conditions (no M0 High)

COLUMNS = [
    "column_name", "donor", "cond_num", "src_label", "stimulus", "agonist", "dose",
    "vendor_catalog", "vehicle", "fastq_r1", "fastq_r2", "md5_r1", "md5_r2",
]


def rows():
    """27 rows in counts.tsv column order; S-numbers run 1..27 over that same order."""
    n = 0
    for donor in DONORS:
        for suffix, (src, stim, agonist, dose, cat, vehicle) in CONDITIONS.items():
            n += 1
            name = f"{donor}_{suffix}"
            yield {
                "column_name": name,
                "donor": donor,
                "cond_num": str(suffix),
                "src_label": src,
                "stimulus": stim,
                "agonist": agonist,
                "dose": dose,
                "vendor_catalog": cat,
                "vehicle": vehicle,
                "fastq_r1": f"{name}_S{n}_R1_001.fastq.gz",
                "fastq_r2": f"{name}_S{n}_R2_001.fastq.gz",
                "md5_r1": "",
                "md5_r2": "",
            }


if __name__ == "__main__":
    out = Path(__file__).parent / "sample_manifest.tsv"
    all_rows = list(rows())
    assert len(all_rows) == 27, len(all_rows)
    with out.open("w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for row in all_rows:
            fh.write("\t".join(row[c] for c in COLUMNS) + "\n")
    print(f"wrote {out} ({len(all_rows)} samples)")
