"""Pull phosphoproteomics data for all donors and aggregate S/T/Y residues.

This is a portable lift of the first two preprocessing sections of
``phosphoproteomics.ipynb`` ("Pull data for all donors" and "Aggregating S, T,
and Y across donors"). Point it at a raw-data directory and it writes the tidy
per-peptide table and the aggregated residue x channel matrix that the rest of
the notebook (median normalization, PCA, volcano) consumes.

Conda environment setup (polars is the only third-party dependency):

    conda create -n phospho python=3.11 polars
    conda activate phospho

Usage:

    python aggregate_phospho.py --input-dir /path/to/phospho_dir
    python aggregate_phospho.py --input-dir DIR --output-dir OUT --no-step2
    python aggregate_phospho.py --input-dir DIR --metadata /path/to/sheet.csv

Two directory layouts are supported; both are discovered by globbing for
``**/01_processed_files/processed_*census-out*.csv`` and both key the sample-name
mapping on the TZ3 number:
  - multi-diffmod (combined enrichments): one or more ``output_<N>diffmods/`` subdirs,
    each containing a ``.../01_processed_files/`` dir of ``processed_*census-out*.csv``
    files (the intermediate processing-run subdir may be any number, not just ``2``);
  - per-enrichment: a single ``.../01_processed_files/processed_<experiment>-census-out.csv``
    where ``<experiment>`` embeds the enrichment (e.g. ``TZ3_114_TiO2``). Enrichments are
    kept as distinct samples (the enrichment is appended to the channel name).

The metadata sheet (``Macrophage phosphoproteomics experiments - Sheet1.csv`` by
default, or ``--metadata``) maps ``tmt_channel_*`` to sample names; its ``experiment``
column is matched to each file by TZ3 number, not by exact string.
"""

import argparse
import re
from collections import defaultdict
from pathlib import Path

import polars as pl
import polars.selectors as cs

METADATA_FILENAME = "Macrophage phosphoproteomics experiments - Sheet1.csv"
# processed CSVs live directly under a 01_processed_files dir in both layouts
PROCESSED_GLOB = "**/01_processed_files/processed_*census-out*.csv"
# known experiment-label typos -> corrected TZ3 number
TZ_TYPO_MAP = {"113": "114"}

# ID columns carried through cross-donor aggregation (notebook's ``all_id_cols``).
ALL_ID_COLS = [
    "uniprot",
    "protein",
    "cell_state",
    "donor",
    "technical_replicate",
    "description",
]


def unnest_channel_name(df):
    """Split the ``cell_state_donor_technicalreplicate`` channel name into columns."""
    return df.with_columns(
        pl.col("channel_name")
        .str.splitn("_", 3)
        .struct.rename_fields(["cell_state", "donor", "technical_replicate"]),
    ).unnest("channel_name")


def parse_experiment(filename: str) -> str:
    """Pull the experiment label out of a processed census-out filename.

    Handles both ``processed_census-out_<exp>.csv`` (experiment after the marker)
    and ``processed_<exp>-census-out.csv`` (experiment before it).
    """
    stem = filename[: -len(".csv")] if filename.endswith(".csv") else filename
    if stem.startswith("processed_"):
        stem = stem[len("processed_") :]
    if "census-out_" in stem:
        return stem.split("census-out_", 1)[1]
    if "-census-out" in stem:
        return stem.split("-census-out", 1)[0]
    return stem


def extract_tz(experiment: str):
    """Return the TZ3 number embedded in an experiment label (typo-corrected)."""
    m = re.search(r"TZ3[-_]?(\d+)", experiment)
    if not m:
        return None
    return TZ_TYPO_MAP.get(m.group(1), m.group(1))


def extract_enrichment(experiment: str):
    """Return the enrichment (TiO2/FeNTA) named in the label, or None if absent."""
    m = re.search(r"TiO2|FeNTA", experiment, re.IGNORECASE)
    return m.group(0) if m else None


def extract_n_diff_mods(path: Path) -> str:
    """Return N from an ``output_<N>diffmods`` path component, else 'NA'."""
    for part in path.parts:
        m = re.fullmatch(r"output_(\d+)diffmods", part)
        if m:
            return m.group(1)
    return "NA"


def load_peptide_data(input_dir: Path, metadata_path: Path) -> pl.DataFrame:
    """Section 1: collect processed data from every donor/search into one tidy table.

    Discovers all ``processed_*census-out*`` CSVs under any ``01_processed_files``
    dir, attaches the TMT channel/sample names from the metadata sheet (matched by
    TZ3 number), unpivots to long format, and cleans into ``peptide_results_clean``.
    Enrichment (when present in the filename) is appended to the sample name so
    FeNTA and TiO2 stay as distinct samples.
    """
    # Collect all processed dfs
    dfs = []
    processed_files = sorted(input_dir.glob(PROCESSED_GLOB))
    if not processed_files:
        raise FileNotFoundError(
            f"No processed census-out CSVs found under {input_dir} "
            f"(looked for {PROCESSED_GLOB})"
        )
    for path in processed_files:
        experiment = parse_experiment(path.name)
        tz = extract_tz(experiment)
        enrichment = extract_enrichment(experiment)
        n_diffs = extract_n_diff_mods(path)
        print(
            f"Reading {path.name}  "
            f"(tz={tz}, enrichment={enrichment}, {n_diffs} diff mods)"
        )
        df = pl.read_csv(path).with_columns(
            pl.lit(n_diffs).alias("n_diff_mods"),
            pl.lit(experiment).alias("experiment"),
            pl.lit(tz, dtype=pl.Utf8).alias("tz"),
            pl.lit(enrichment, dtype=pl.Utf8).alias("enrichment"),
        )
        dfs.append(df)
    df = pl.concat(dfs, how="vertical_relaxed")

    # Add channel names, keyed on TZ3 number so combined/per-enrichment files and
    # the sheet's dated experiment labels all line up.
    metadata = (
        pl.read_csv(metadata_path)
        .drop(cs.contains("peptide ids"), "abbreviation")
        .unpivot(
            on=cs.contains("tmt_channel_"),
            index=~cs.contains("tmt_channel_"),
            value_name="channel_name",
        )
        .with_columns(pl.col("experiment").str.extract(r"TZ3[-_]?(\d+)", 1).alias("tz"))
        .select("tz", "variable", "channel_name")
    )
    tag_columns = [col for col in df.columns if col.startswith("tag_")]
    rename_mapping = {col: f"tmt_channel_{i}" for i, col in enumerate(tag_columns, 1)}

    clean_df = (
        df.rename(mapping=rename_mapping)
        .unpivot(
            on=cs.contains("tmt_channel_"),
            index=~cs.contains("tmt_channel_"),
            value_name="channel_ratio",
        )
        .join(other=metadata, on=["tz", "variable"], how="left")
        # keep enrichments as distinct samples
        .with_columns(
            pl.when(pl.col("enrichment").is_not_null())
            .then(pl.concat_str("channel_name", pl.lit("_"), "enrichment"))
            .otherwise(pl.col("channel_name"))
            .alias("channel_name")
        )
    )

    unmatched = clean_df.filter(pl.col("channel_name").is_null())
    if unmatched.height:
        missing_tz = sorted(set(unmatched["tz"].to_list()))
        print(
            f"WARNING: {unmatched.height} rows had no sample-name match "
            f"(TZ3 numbers not in metadata sheet: {missing_tz})"
        )

    clean_df = clean_df.drop(["experiment", "variable", "tz", "enrichment"])

    peptide_results_clean = unnest_channel_name(
        clean_df.select(
            [
                "uniprot",
                "residue",
                "channel_name",
                "channel_ratio",
                "description",
                "sequence",
                "n_diff_mods",
            ]
        )
    )

    peptide_results_clean = (
        peptide_results_clean.filter(~pl.col("uniprot").str.contains("contaminant"))
        .with_columns(
            # we are not removing mods from the sequence
            # so sequence aggregation won't average multiple mods
            pl.col("sequence").str.split(".").list.get(1),
            pl.col("sequence")
            .str.count_matches(pattern="*", literal=True)
            .alias("num_modifications"),
            pl.when(pl.col("description").str.contains("GN="))
            .then(pl.col("description").str.extract(r"GN=(\S+)"))
            .otherwise(pl.col("uniprot"))
            .alias("protein"),
        )
        .sort(by="num_modifications")
    )

    return peptide_results_clean


def get_shortest_string(sequences, aas_to_aggregate):
    """Find shortest common string in list containing amino acids to aggregate."""
    # polars series to python list
    sequences = list(sequences)

    # handle empty or single item cases
    if not sequences:
        return None
    if len(sequences) == 1:
        return sequences[0]

    # get sequences with minimum length
    min_len = min(len(x) for x in sequences if x is not None)
    sequence_list = [x for x in sequences if x is not None and len(x) == min_len]

    if len(sequence_list) == 1:
        return sequence_list[0]

    if len(set(sequence_list)) == 1:
        return sequence_list[0]  # all sequences are equal

    # find sequence with amino acid closest to N-terminus
    tmp_list = []
    seq_indices = []

    for i, seq in enumerate(sequence_list):
        for aa in aas_to_aggregate:
            finder = seq.find(aa)
            tmp_list.append(finder)
            seq_indices.append(i)

    if tmp_list:
        tmp_list = [float("inf") if x == -1 else x for x in tmp_list]
        min_index = tmp_list.index(min(tmp_list))
        sequence = sequence_list[seq_indices[min_index]]
    else:
        sequence = sequence_list[0]

    return sequence


def average_seqs(original_df, aggregated_seqs_df):
    """Join with aggregation mappings and average peptides that share a representative.

    Takes the average (mean) channel ratio for aggregated residues.
    """
    return (
        original_df.join(aggregated_seqs_df, on=["uniprot", "sequence"], how="left")
        .with_columns(
            pl.col("sequence_representative").fill_null(
                pl.col("sequence")
            )  # fallback to original if no mapping found
        )
        .drop("sequence")
        .rename({"sequence_representative": "sequence"})
        .group_by(ALL_ID_COLS + ["sequence"])
        .agg(
            pl.col("channel_ratio").mean(),
            pl.col("residue").unique().sort().str.join(", "),
        )
    )


def solveContainedSequences(sequence_string, list_of_unique_seq, aa_to_aggregate):
    """Collapse a peptide onto the shortest peptide that contains it, when safe."""
    # Find all sequences that contain or are contained by sequence_string
    contained = [
        item
        for item in list_of_unique_seq
        if sequence_string in item or item in sequence_string
    ]

    if len(contained) <= 1:
        return sequence_string

    # Get shortest sequence (parent)
    parent_string = min(contained, key=len)

    # Check if leftover contains target amino acids
    leftover = sequence_string.replace(parent_string, "")
    aa_set = set(aa_to_aggregate)

    return sequence_string if aa_set & set(leftover) else parent_string


def aggregate_residues(
    peptide_results_clean: pl.DataFrame, perform_step_2: bool = True
) -> pl.DataFrame:
    """Section 2: aggregate S, T, Y across donors and pivot to a residue x channel matrix.

    Step 1 averages multiple peptides that map to the same shortest representative
    sequence. Step 2 (optional) additionally collapses a peptide fully contained in
    another peptide.
    """
    print(
        f"Starting peptides: "
        f"{len(peptide_results_clean.select(['uniprot', 'sequence']).unique())}"
    )

    # aggregate residues so multiple peptides containing the same residue are averaged
    # build mapping keyed on (uniprot, sequence, num_modifications)
    shortest_seq_rows = []
    for uniprot, residue, num_modifications, seqs in (
        peptide_results_clean.group_by("uniprot", "residue", "num_modifications")
        .agg(pl.col("sequence").unique())
        .rows()
    ):
        for seq in seqs:
            shortest_seq_rows.append(
                {
                    "uniprot": uniprot,
                    "sequence": seq,
                    "sequence_representative": get_shortest_string(seqs, ["S", "T", "Y"]),
                }
            )

    shortest_seqs_df = pl.DataFrame(shortest_seq_rows).unique(
        subset=["uniprot", "sequence"]
    )  # guard against any duplicates

    sequence_means = average_seqs(
        original_df=peptide_results_clean, aggregated_seqs_df=shortest_seqs_df
    )

    print(
        f"after step 1: "
        f"{len(sequence_means.select(['uniprot', 'sequence']).unique())}"
    )

    if perform_step_2:
        # Create aggregation mapping across all donors for each uniprot
        aggregated_seq_rows = []
        for uniprot, seqs in (
            sequence_means.group_by("uniprot").agg(pl.col("sequence").unique()).rows()
        ):
            for seq in seqs:
                aggregated_seq_rows.append(
                    {
                        "uniprot": uniprot,
                        "sequence": seq,
                        "sequence_representative": solveContainedSequences(
                            seq, seqs, ["S", "T", "Y"]
                        ),
                    }
                )

        aggregated_seqs_df = pl.DataFrame(aggregated_seq_rows).unique(
            subset=["uniprot", "sequence"]
        )

        # Apply the sequence aggregation mapping and then re-aggregate by donor
        sequence_means = average_seqs(
            original_df=sequence_means, aggregated_seqs_df=aggregated_seqs_df
        )

        print(
            f"after step 2: "
            f"{len(sequence_means.select(['uniprot', 'sequence']).unique())}"
        )

    sequence_to_residue = defaultdict(set)
    for row in (
        sequence_means.select(["sequence", "residue"]).unique().iter_rows(named=True)
    ):
        # Split on "," and strip whitespace
        residues = row["residue"].replace(";", ",").split(",")
        sequence_to_residue[row["sequence"]].update(r.strip() for r in residues)
    clean_sequence_to_residue = {}
    for sequence in sequence_to_residue.keys():
        clean_sequence_to_residue[sequence] = ",".join(sequence_to_residue[sequence])

    # reformat to residue x channel matrix
    df = (
        (
            sequence_means.with_columns(
                pl.concat_str(
                    ["cell_state", "donor", "technical_replicate"], separator="_"
                ).alias("channel_name")
            )
            .drop(["cell_state", "donor", "technical_replicate"])
            .pivot(
                on="channel_name",
                values="channel_ratio",
                index=["uniprot", "protein", "description", "sequence"],
            )
        )
        .with_columns(
            pl.col("sequence").replace(clean_sequence_to_residue).alias("residue")
        )
        .with_columns(
            pl.concat_str(
                ["uniprot", "protein", "description", "residue", "sequence"],
                separator="_",
            ).alias("identifier")
        )
    )

    return df


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Pull phosphoproteomics data for all donors and aggregate S/T/Y residues "
            "into a residue x channel matrix."
        )
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        type=Path,
        help="Directory holding the metadata sheet and output_<N>diffmods subdirs.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Where to write the output CSVs (default: --input-dir).",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=None,
        help=(
            "Path to the sample-name mapping sheet "
            f"(default: <input-dir>/{METADATA_FILENAME})."
        ),
    )
    parser.add_argument(
        "--no-step2",
        dest="perform_step_2",
        action="store_false",
        help="Skip step-2 aggregation of peptides fully contained in another peptide.",
    )
    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir or input_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata_path = args.metadata or (input_dir / METADATA_FILENAME)
    if not metadata_path.exists():
        raise FileNotFoundError(
            f"Metadata sheet not found: {metadata_path}. "
            f"Pass --metadata to point at the sample-name mapping sheet."
        )

    peptide_results_clean = load_peptide_data(input_dir, metadata_path)
    aggregated = aggregate_residues(
        peptide_results_clean, perform_step_2=args.perform_step_2
    )

    peptides_path = output_dir / "peptide_results_clean.csv"
    matrix_path = output_dir / "aggregated_phospho_matrix.csv"
    peptide_results_clean.write_csv(peptides_path)
    aggregated.write_csv(matrix_path)

    print(f"\nWrote tidy peptide table:   {peptides_path}")
    print(f"Wrote aggregated matrix:    {matrix_path}")


if __name__ == "__main__":
    main()
