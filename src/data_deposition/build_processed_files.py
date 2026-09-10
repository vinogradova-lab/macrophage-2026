"""Build the two processed data matrices GEO receives, then prove both are faithful.

  counts_27samples.tsv   counts.tsv minus its eight ``TZ3_*`` columns.
  normalized_counts.tsv  the DESeq2 matrix from Dropbox, with the ``gene.id`` header field
                         R's ``write.table(row.names=TRUE)`` omits prepended.

    python build_processed_files.py                 # build both, then verify
    python build_processed_files.py --verify-only   # re-check existing outputs

A non-zero exit means the matrices must not be uploaded.

Runbook: rna/GEO_SUBMISSION.md (project note 2, Phase 1).
"""

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).parent
RNA = HERE.parent.parent / "rna"
DEFAULT_SOURCE = RNA / "counts.tsv"
DEFAULT_OUT = RNA / "counts_27samples.tsv"
MANIFEST = HERE / "sample_manifest.tsv"

# Written to Dropbox by DESEQ2_normalization.Rmd:54, not to the repo.
DROPBOX_RNA = Path(
    "/Users/henrysanford/Dropbox @RU Dropbox/Vinogradova Laboratory/Macrophage project/"
    "Manuscript 1/01_Data Analysis/02_RNA-seq"
)
DEFAULT_NORM_SOURCE = DROPBOX_RNA / "normalized_counts.tsv"
DEFAULT_NORM_OUT = RNA / "normalized_counts.tsv"

EXCLUDED_DONOR_PREFIX = "TZ3_"
ID_COLUMN = "gene.id"

fails = 0


def check(ok, description, detail=""):
    global fails
    if ok:
        print(f"  PASS  {description}")
    else:
        fails += 1
        print(f"  FAIL  {description}{': ' + detail if detail else ''}")
    return ok


def rows(path):
    """Yield (line_number, fields) with no quoting or line-ending translation applied."""
    with open(path, newline="") as fh:
        for n, line in enumerate(fh, start=1):
            yield n, line.rstrip("\n").split("\t")


def manifest_columns():
    """The submitted sample names, in the order sample_manifest.tsv lists them."""
    with MANIFEST.open(newline="") as fh:
        header = next(fh).rstrip("\n").split("\t")
        idx = header.index("column_name")
        return [line.rstrip("\n").split("\t")[idx] for line in fh if line.strip()]


def split_header(header):
    """(keep indices, kept names, dropped names) -- selection is by name, not position."""
    keep = [i for i, name in enumerate(header) if not name.startswith(EXCLUDED_DONOR_PREFIX)]
    dropped = [name for name in header if name.startswith(EXCLUDED_DONOR_PREFIX)]
    return keep, [header[i] for i in keep], dropped


def build(source, out):
    with open(source, newline="") as src:
        header = next(src).rstrip("\n").split("\t")
        keep, kept_names, dropped = split_header(header)
        with open(out, "w", newline="") as dst:
            dst.write("\t".join(kept_names) + "\n")
            for line in src:
                fields = line.rstrip("\n").split("\t")
                dst.write("\t".join(fields[i] for i in keep) + "\n")
    print(f"wrote {out}")
    print(f"  kept {len(kept_names) - 1} sample columns, dropped {len(dropped)}: {', '.join(dropped)}")


def verify(source, out):
    """Re-read both files from disk and compare; nothing in memory from build() is trusted."""
    src_bytes = source.read_bytes()
    out_bytes = out.read_bytes()

    check(b"\r" not in src_bytes and b'"' not in src_bytes and src_bytes.endswith(b"\n"),
          "source is LF-only, unquoted, newline-terminated")
    check(b"\r" not in out_bytes and b'"' not in out_bytes and out_bytes.endswith(b"\n"),
          "output is LF-only, unquoted, newline-terminated")

    src = rows(source)
    _, header = next(src)
    keep, kept_names, dropped = split_header(header)
    n_src_fields = len(header)

    check(header[0] == ID_COLUMN, f"source first column is {ID_COLUMN!r}", f"got {header[0]!r}")
    expected_dropped = sorted(n for n in header if n.startswith(EXCLUDED_DONOR_PREFIX))
    check(len(dropped) == 8 and sorted(dropped) == expected_dropped,
          "dropped columns are exactly the excluded donor's",
          f"dropped {len(dropped)}: {', '.join(dropped)}")

    expected_kept = manifest_columns()
    check(kept_names[1:] == expected_kept,
          f"kept columns match sample_manifest.tsv order ({len(expected_kept)} samples)",
          f"{kept_names[1:]} != {expected_kept}")

    out_rows = rows(out)
    _, out_header = next(out_rows)
    check(out_header == kept_names, "output header is the kept source columns, in source order",
          f"{out_header} != {kept_names}")

    n_src = src_bytes.count(b"\n")
    n_out = out_bytes.count(b"\n")
    check(n_out == n_src, "output row count equals source row count",
          f"output {n_out} != source {n_src}")

    n_out_fields = len(kept_names)
    ragged_src = ragged_out = None
    mismatch = None
    for (line_no, src_fields), (_, out_fields) in zip(src, out_rows):
        if len(src_fields) != n_src_fields and ragged_src is None:
            ragged_src = f"line {line_no} has {len(src_fields)} fields, expected {n_src_fields}"
        if len(out_fields) != n_out_fields and ragged_out is None:
            ragged_out = f"line {line_no} has {len(out_fields)} fields, expected {n_out_fields}"
        if mismatch is None:
            for j, i in enumerate(keep):
                got = out_fields[j] if j < len(out_fields) else "<missing>"
                want = src_fields[i] if i < len(src_fields) else "<missing>"
                if got != want:
                    mismatch = (f"line {line_no}, column {kept_names[j]!r}: "
                                f"output {got!r} != source {want!r}")
                    break

    check(ragged_src is None, f"every source row has {n_src_fields} fields", ragged_src or "")
    check(ragged_out is None, f"every output row has {n_out_fields} fields", ragged_out or "")
    check(mismatch is None,
          f"every retained cell is byte-identical to the source "
          f"({min(n_src, n_out) - 1:,} rows x {n_out_fields} columns)",
          mismatch or "")

    return n_src - 1


def build_normalized(source, out):
    """Copy the DESeq2 matrix in, prepending the ID field R's write.table omits."""
    with open(source, newline="") as src:
        header = src.readline().rstrip("\n").split("\t")
        # n names = headerless-rowname form; n+1 = the column is already labelled.
        n_samples = len(manifest_columns())
        header = [ID_COLUMN] + (header[1:] if len(header) > n_samples else header)
        with open(out, "w", newline="") as dst:
            dst.write("\t".join(header) + "\n")
            while chunk := src.read(1 << 20):  # body passes through untouched
                dst.write(chunk)
    print(f"wrote {out}")
    print(f"  {len(header) - 1} sample columns, ID column {ID_COLUMN!r} prepended to the header")


def verify_normalized(source, out, counts_out):
    """Check the copy against its source *and* against the raw-counts matrix beside it."""
    out_bytes = out.read_bytes()
    check(b"\r" not in out_bytes and b'"' not in out_bytes and out_bytes.endswith(b"\n"),
          "normalized output is LF-only, unquoted, newline-terminated")

    with open(source, newline="") as fh:
        src_header = fh.readline().rstrip("\n").split("\t")
        src_body = fh.read()
    out_text = out_bytes.decode()
    out_header = out_text.split("\n", 1)[0].split("\t")
    out_body = out_text.split("\n", 1)[1]

    check(out_body == src_body, "normalized body is byte-identical to the Dropbox matrix",
          f"{len(out_body)} bytes out vs {len(src_body)} in")

    samples = manifest_columns()
    check(out_header == [ID_COLUMN] + samples,
          f"normalized header is {ID_COLUMN!r} + the {len(samples)} submitted samples, in order",
          f"got {out_header}")
    check(len(src_header) in (len(samples), len(samples) + 1),
          f"source header carried {len(samples)} sample names",
          f"got {len(src_header)}: {src_header}")

    if counts_out.exists():
        norm_ids = [line.split("\t", 1)[0] for line in out_body.split("\n") if line]
        counts_rows = rows(counts_out)
        next(counts_rows)
        counts_ids = [f[0] for _, f in counts_rows]
        first_diff = next((i for i, (a, b) in enumerate(zip(norm_ids, counts_ids)) if a != b), None)
        check(norm_ids == counts_ids,
              f"gene IDs match {counts_out.name} row for row ({len(norm_ids):,} features)",
              f"row {first_diff + 2}: {norm_ids[first_diff]!r} != {counts_ids[first_diff]!r}"
              if first_diff is not None
              else f"{len(norm_ids):,} rows vs {len(counts_ids):,}")

    n_fields = {len(line.split("\t")) for line in out_body.split("\n") if line}
    check(n_fields == {len(out_header)},
          f"every normalized row has {len(out_header)} fields", f"saw widths {sorted(n_fields)}")


def feature_breakdown(out):
    ensembl = dfam = other = 0
    src = rows(out)
    next(src)
    for _, fields in src:
        if fields[0].startswith("ENSG"):
            ensembl += 1
        elif fields[0].startswith("DF"):
            dfam += 1
        else:
            other += 1
    return ensembl, dfam, other


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--source", type=Path, default=DEFAULT_SOURCE, help="full counts matrix")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT, help="subset matrix to write")
    ap.add_argument("--norm-source", type=Path, default=DEFAULT_NORM_SOURCE,
                    help="DESeq2 matrix written by DESEQ2_normalization.Rmd")
    ap.add_argument("--norm-out", type=Path, default=DEFAULT_NORM_OUT,
                    help="header-repaired normalized matrix to write")
    ap.add_argument("--verify-only", action="store_true",
                    help="check existing outputs against their sources without rewriting them")
    args = ap.parse_args()

    if not args.source.exists():
        sys.exit(f"source not found: {args.source}\n"
                 "counts.tsv is gitignored; copy it from the Dropbox path in "
                 "DESEQ2_normalization.Rmd:25-28")
    if not args.norm_source.exists():
        sys.exit(f"normalized source not found: {args.norm_source}\n"
                 "normalized_counts.tsv is written to Dropbox by DESEQ2_normalization.Rmd:54")
    for out in (args.out, args.norm_out):
        if args.verify_only and not out.exists():
            sys.exit(f"--verify-only but output not found: {out}")

    if not args.verify_only:
        build(args.source, args.out)
        build_normalized(args.norm_source, args.norm_out)

    print(f"\nverifying {args.out.name} against {args.source.name}")
    n_rows = verify(args.source, args.out)

    ensembl, dfam, other = feature_breakdown(args.out)
    print(f"{n_rows:,} features: {ensembl:,} Ensembl/GENCODE (ENSG*), "
          f"{dfam:,} Dfam TE families (DF*), {other:,} other")

    print(f"\nverifying {args.norm_out.name} against the Dropbox matrix")
    verify_normalized(args.norm_source, args.norm_out, args.out)

    if fails:
        print(f"{fails} check(s) FAILED -- do not upload these matrices", file=sys.stderr)
        return 1
    print("all checks PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
