"""Fill GEO's seq_template.xlsx for the macrophage bulk RNA-seq series.

    python build_geo_metadata.py [--template PATH] [--md5 PATH] [--verify-protocols]

Runbook: rna/GEO_SUBMISSION.md (this is its Phase 3).
"""

import argparse
import csv
import sys
from pathlib import Path

import openpyxl
from openpyxl.formatting.formatting import ConditionalFormattingList
from openpyxl.formatting.rule import Rule
from openpyxl.styles import Font, PatternFill
from openpyxl.styles.differential import DifferentialStyle

import protocols_verbatim as P

HERE = Path(__file__).parent
DEFAULT_TEMPLATE = Path.home() / "Downloads" / "seq_template.xlsx"
MANIFEST = HERE / "sample_manifest.tsv"
OUTPUT = HERE / "seq_template_filled.xlsx"

PROCESSED_FILES = ["counts_27samples.tsv.gz", "normalized_counts.tsv.gz"]

SUPPLEMENTARY_FILES = PROCESSED_FILES

TITLE = (
    "Bulk RNA-seq of primary human monocyte-derived macrophages following "
    "Toll-like receptor and STING activation"
)

SUMMARY = P.ABSTRACT

DESIGN = (
    "Primary human monocyte-derived macrophages from three healthy donors (TZ10, TZ11, TZ12) "
    "were profiled across nine conditions each, for 27 samples in total: seven agonist "
    "conditions -- TLR1/2 (CU-T12-09), TLR3 (Poly(I:C)), TLR4 (LPS), TLR7 (imiquimod), TLR8 "
    "(CL075), TLR9 (CpG ODN 2006) and STING (cGAMP) -- plus two unstimulated DMSO vehicle "
    "controls, at 0.025% and 1% DMSO. No appreciable difference was observed between them; only "
    "the 0.025% control is shown in the associated manuscript, but both are deposited here. In "
    "the count matrices they carry the source labels 'M0 Low' and 'M0 High' respectively."
)
# Filled in at publication.  GEO format: Firstname,MI,Lastname -- one per row.
CONTRIBUTORS = []

# PROTOCOLS -- verbatim from the STAR Methods, plus the EDIT_* additions below.

# (a) agonist table inlined
AGONIST_TABLE = (
    "TLR1-2, CU-T12-09, 10 µM (Millipore Sigma 5325830001); "
    "TLR3, Poly(I:C), 20 µg/mL (Tocris Bioscience 42-871-0); "
    "TLR4, LPS, 100 ng/mL (Millipore Sigma L2630-10); "
    "TLR7, Imiquimod, 10 µg/mL (Fisher Scientific AAJ63990MC); "
    "TLR8, CL075, 0.5 µg/mL (Millipore Sigma SML2235); "
    "TLR9, CpG ODN 2006, 5 µg/mL (Integrated DNA Technologies, custom); "
    "STING, cGAMP, 10 µg/mL (Invivogen tlrl-nacga23-5)"
)
EDIT_A = ("see table below", AGONIST_TABLE)

# (e) DMSO vehicle
EDIT_E = (
    " Agonists were delivered in DMSO; unstimulated controls received vehicle alone, at 0.025% "
    "DMSO matching the carrier load of all conditions other than TLR7, or at 1% DMSO matching "
    "the TLR7 (imiquimod) condition."
)

# (b) strandedness
EDIT_B = " Libraries are reverse-stranded (featureCounts -s 2)."

# (c) TE annotation
EDIT_C = "Gene-level quantification used GENCODE release 44."

# (d) DESeq2 normalization
EDIT_D = (
    "Differential gene expression analysis was performed with DESeq2, using its "
    "median-of-ratios method to normalize RNA counts."
)

GROWTH = f"{P.PBMC_MONOCYTE} {P.DIFFERENTIATION}"
TREATMENT = P.STIMULATION.replace(*EDIT_A) + EDIT_E
EXTRACT = P.RNA_ISOLATION
LIBRARY = P.LIBRARY_PREP + EDIT_B
DATA_PROCESSING_STEPS = [P.DATA_PROCESSING, EDIT_C, EDIT_D]
GENOME_BUILD = "GRCh38"
PROCESSED_FORMATS = [
    "counts_27samples.tsv.gz: tab-delimited text, raw featureCounts gene-level counts. Rows are "
    "69,233 features (67,928 version-suffixed Ensembl/GENCODE gene IDs plus 1,305 Dfam "
    "transposable-element family accessions); columns are the 27 samples, named as in the "
    "'description' field of each Sample record.",
    "normalized_counts.tsv.gz: tab-delimited text, the same matrix after DESeq2 median-of-ratios "
    "size-factor normalization (estimateSizeFactors), computed across all 27 samples.",
]

def load_manifest():
    with MANIFEST.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if len(rows) != 27:
        raise SystemExit(f"expected 27 samples in {MANIFEST}, found {len(rows)}")
    return rows


def load_md5sums(path):
    """Parse the core's .md5 file -> {basename: checksum}, either field order."""
    sums = {}
    for line in Path(path).read_text().splitlines():
        parts = line.split()
        if len(parts) != 2:
            continue
        a, b = parts
        digest, name = (a, b) if len(a) == 32 else (b, a)
        sums[Path(name).name] = digest
    return sums


def sample_title(row):
    if row["stimulus"] == "none":
        return (
            f"Primary human macrophages, donor {row['donor']}, "
            f"unstimulated control ({row['vehicle']} vehicle)"
        )
    return (
        f"Primary human macrophages, donor {row['donor']}, "
        f"{row['stimulus']} ({row['agonist']})"
    )


def treatment_cell(row):
    if row["stimulus"] == "none":
        return f"vehicle only ({row['vehicle']}), 24 h"
    return f"{row['agonist']} ({row['stimulus']} agonist), {row['dose']}, 24 h"


def stimulus_cell(row):
    if row["stimulus"] == "none":
        return f"none (vehicle control, {row['vehicle']})"
    return row["stimulus"]


def find_label_row(ws, label, limit=200):
    for r in range(1, limit):
        if ws.cell(row=r, column=1).value == label:
            return r
    raise SystemExit(f"could not find {label!r} in column A -- template layout changed?")


def fill_label_block(ws, label, values):
    """Write `values` down column B against repeated `label` rows, adding rows if short."""
    first = find_label_row(ws, label)
    slots = 0
    while ws.cell(row=first + slots, column=1).value == label:
        slots += 1
    if len(values) > slots:
        ws.insert_rows(first + slots, len(values) - slots)
        for r in range(first + slots, first + len(values)):
            ws.cell(row=r, column=1, value=label)
    for i, value in enumerate(values):
        ws.cell(row=first + i, column=2, value=value)


def write_sheet(ws, samples):
    ws["B12"] = TITLE
    ws["B13"] = SUMMARY
    ws["B14"] = DESIGN
    fill_label_block(ws, "contributor", CONTRIBUTORS)
    fill_label_block(ws, "supplementary file", SUPPLEMENTARY_FILES)

    header = find_label_row(ws, "*library name")
    first_row = header + 1
    ws.cell(row=header, column=8, value="donor")
    ws.cell(row=header, column=10, value="stimulus")

    protocols = find_label_row(ws, "PROTOCOLS")
    ws.insert_rows(protocols, len(samples) - (protocols - first_row))

    for i, row in enumerate(samples):
        r = first_row + i
        values = [
            row["column_name"],
            sample_title(row),
            "RNA-Seq",
            "Homo sapiens",
            None,
            None,
            "monocyte-derived macrophage",
            row["donor"],
            treatment_cell(row),
            stimulus_cell(row),
            "polyA RNA",
            "paired-end",
            "Illumina NovaSeq 6000",
            f"Source condition label: {row['src_label']}",
            PROCESSED_FILES[0],
            PROCESSED_FILES[1],
            row["fastq_r1"],
            row["fastq_r2"],
        ]
        for c, value in enumerate(values, start=1):
            if value is not None:
                ws.cell(row=r, column=c, value=value)

    protocols = find_label_row(ws, "growth protocol")
    ws.cell(row=protocols, column=2, value=GROWTH)
    ws.cell(row=find_label_row(ws, "treatment protocol"), column=2, value=TREATMENT)
    ws.cell(row=find_label_row(ws, "*extract protocol"), column=2, value=EXTRACT)
    ws.cell(row=find_label_row(ws, "*library construction protocol"), column=2, value=LIBRARY)

    step = find_label_row(ws, "*data processing step")
    for i, text in enumerate(DATA_PROCESSING_STEPS):
        ws.cell(row=step + i, column=2, value=text)
    ws.cell(row=find_label_row(ws, "*genome build/assembly"), column=2, value=GENOME_BUILD)

    fmt = find_label_row(ws, "*processed data files format and content")
    for i, text in enumerate(PROCESSED_FORMATS):
        ws.cell(row=fmt + i, column=2, value=text)

    paired = find_label_row(ws, "file name 1", limit=400)
    for i, row in enumerate(samples):
        ws.cell(row=paired + 1 + i, column=1, value=row["fastq_r1"])
        ws.cell(row=paired + 1 + i, column=2, value=row["fastq_r2"])

    realign_duplicate_highlighting(ws, len(samples), first_row, paired + 1)


# ``insert_rows`` does not move the template's duplicate-value rules, so they are rebuilt
# below against the rows the samples actually occupy.
DUPLICATE_STYLE = DifferentialStyle(
    font=Font(color="9C0006"), fill=PatternFill(bgColor="FFC7CE")
)


def realign_duplicate_highlighting(ws, n_samples, first_row, paired_first_row):
    last = first_row + n_samples - 1
    ws.conditional_formatting = ConditionalFormattingList()
    for ref in (
        f"A{first_row}:A{last}",
        f"B{first_row}:B{last}",
        f"Q{first_row}:T{last}",
        f"A{paired_first_row}:B{paired_first_row + n_samples - 1}",
    ):
        ws.conditional_formatting.add(
            ref, Rule(type="duplicateValues", dxf=DUPLICATE_STYLE, stopIfTrue=False)
        )


def write_md5_sheet(ws, samples, md5sums):
    missing = []
    for i, row in enumerate(samples):
        for j, key in enumerate(("fastq_r1", "fastq_r2")):
            r = 9 + i * 2 + j
            name = row[key]
            ws.cell(row=r, column=1, value=name)
            digest = md5sums.get(name, "")
            if not digest:
                missing.append(name)
            ws.cell(row=r, column=2, value=digest or "TODO: checksum")
    for i, name in enumerate(PROCESSED_FILES):
        ws.cell(row=9 + i, column=6, value=name)
        ws.cell(row=9 + i, column=7, value=md5sums.get(name, "") or "TODO: checksum")
    return missing


def verify_protocols():
    """Re-extract from the .docx and confirm the committed constants still match."""
    from _extract_protocols import extract

    fresh = extract()
    bad = [name for name, text in fresh.items() if getattr(P, name) != text]
    if bad:
        raise SystemExit(f"protocols_verbatim.py has drifted from the .docx: {bad}")
    print(f"abstract and protocols match the .docx sources verbatim ({len(fresh)} blocks)")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--template", type=Path, default=DEFAULT_TEMPLATE)
    ap.add_argument("--md5", type=Path, help="sequencing core's .md5 file")
    ap.add_argument("--verify-protocols", action="store_true")
    args = ap.parse_args()

    if args.verify_protocols:
        verify_protocols()

    samples = load_manifest()
    md5sums = load_md5sums(args.md5) if args.md5 else {}

    wb = openpyxl.load_workbook(args.template)
    write_sheet(wb["Metadata"], samples)
    missing = write_md5_sheet(wb["MD5 Checksums"], samples, md5sums)
    wb.save(OUTPUT)

    print(f"wrote {OUTPUT}")
    print(f"  samples          {len(samples)}")
    print(f"  contributors     {len(CONTRIBUTORS) or '(blank -- add at publication)'}")
    print(f"  paired-end rows  {len(samples)}")
    print(f"  md5 raw rows     {len(samples) * 2} ({len(missing)} still missing)")
    print("  protocol edits   (a) agonist table inlined  (b) strandedness  (c) TE annotation"
          "  (d) DESeq2 normalization  (e) DMSO vehicle")
    if missing:
        print(f"  NOTE: no checksums yet -- rerun with --md5 <core .md5 file>", file=sys.stderr)


if __name__ == "__main__":
    main()
