"""Cysteine cross-reference: the funnel schematic's reference databases, in one place.

The funnel narrows reactivity changes through complex membership (CORUM / ComplexPortal),
conservation depth, and AlphaMissense pathogenicity; ClinVar proximity is analysed alongside it.
Everything here is a fact about a ``(uniprot, position)`` pair (or a whole protein), so it is
independent of which reactivity run is being annotated.

The residue-level databases are large (AlphaMissense 1.4 GB, ClinVar 3.7 GB), so
``load_cys_cross_ref`` caches a small per-study table:

    conda run -n polars python -m src.cross_ref --rebuild

The cache is keyed to the study's cysteines -- rebuild it whenever the reactivity run changes.

Cutoffs are the funnel's own (reactivity_polars.ipynb cells 38/46) and are the single definition
used by the notebook, the alluvial plot, and Data S2-3.
"""

import argparse
from pathlib import Path

import polars as pl

from src.macrophage import CLINVAR_TXT, REFERENCE_DBS, REPO_ROOT, require

# AlphaMissense supplement cutoffs.
AM_PATHOGENIC, AM_BENIGN = 0.564, 0.34
# Conservation is a count out of 102 organisms (doi.org/10.1038/s41586-025-08756-y table 10).
CONSERVATION_ORGANISMS = 102
CONSERVATION_BINS = ((49, "51+"), (10, "11-50"))
CONSERVATION_FLOOR = "1-10"
# pPSE: solvent accessibility. Used by the alluvial plot, not by Data S2-3.
PPSE_BURIED, PPSE_INTERMEDIATE = 9, 5
# ClinVar variants within this many residues of a cysteine are considered to implicate it.
CLINVAR_RESIDUE_TOLERANCE = 15
CLINVAR_PATHOGENIC = ["Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic"]

CORUM_XLSX = REFERENCE_DBS / "CORUM download 2022_09_12.xlsx"
COMPLEXPORTAL_TSV = REFERENCE_DBS / "9606.tsv"
CONSERVATION_CSV = REFERENCE_DBS / "cysteine_conservation.csv"
# Comma-separated despite the .tsv suffix (it was written with write_csv; see cell 38).
ALPHAMISSENSE_TSV = REFERENCE_DBS / "AlphaMissense_aa_substitutions_filtered.tsv"
PPSE_CSV = REFERENCE_DBS / "20230315_Supplementary_AlphaFold_pPSE.csv"

CACHE_CSV = REPO_ROOT / "reactivity" / "cys_cross_ref.csv"


def conservation_bin(depth):
    """Conservation depth -> the funnel's bin label. None passes through."""
    if depth is None:
        return None
    for floor, label in CONSERVATION_BINS:
        if depth > floor:
            return label
    return CONSERVATION_FLOOR


def am_class(score):
    """AlphaMissense pathogenicity -> Pathogenic / Ambiguous / Benign. None passes through."""
    if score is None:
        return None
    if score > AM_PATHOGENIC:
        return "Pathogenic"
    if score < AM_BENIGN:
        return "Benign"
    return "Ambiguous"


def ppse_class(ppse):
    """pPSE -> buried / intermediate / accessible. None passes through."""
    if ppse is None:
        return None
    if ppse > PPSE_BURIED:
        return "buried"
    if ppse > PPSE_INTERMEDIATE:
        return "intermediate"
    return "accessible"


def load_corum():
    """Human CORUM complex names per subunit accession."""
    return (
        pl.read_excel(require(CORUM_XLSX))
        .filter(pl.col("Organism").eq("Human"))
        .select("subunits(UniProt IDs)", "ComplexName")
        .with_columns(pl.col("subunits(UniProt IDs)").str.split(";").alias("uniprot"))
        .explode("uniprot")
        .group_by("uniprot")
        .agg(pl.col("ComplexName").unique().sort().str.join("|").alias("corum_complexes"))
    )


def load_complexportal():
    """ComplexPortal complex names per participating accession (isoform suffixes stripped)."""
    return (
        pl.read_csv(require(COMPLEXPORTAL_TSV), separator="\t")
        .with_columns(
            pl.col("Identifiers (and stoichiometry) of molecules in complex")
            .str.split("|")
            .alias("uniprot")
        )
        .explode("uniprot")
        .with_columns(
            pl.col("uniprot").str.split("(").list.get(0).str.split("-").list.get(0)
        )
        .group_by("uniprot")
        .agg(
            pl.col("Recommended name")
            .unique()
            .sort()
            .str.join("|")
            .alias("complexportal_complexes")
        )
    )


def load_complex_membership(uniprots):
    """Complex membership per protein: the two source name lists + the funnel's binary call.

    Returns one row per requested accession, so a protein in neither database is reported as
    "Not complex" rather than going missing. Small enough to read on demand -- not cached.
    """
    return (
        pl.DataFrame({"uniprot": sorted(set(uniprots))})
        .join(load_corum(), on="uniprot", how="left")
        .join(load_complexportal(), on="uniprot", how="left")
        .with_columns(
            pl.when(
                pl.col("corum_complexes").is_not_null()
                | pl.col("complexportal_complexes").is_not_null()
            )
            .then(pl.lit("Complex"))
            .otherwise(pl.lit("Not complex"))
            .alias("complex_membership")
        )
    )


def load_conservation(uniprots):
    """Conservation depth per (uniprot, cysteine position), out of 102 organisms.

    The file carries a title row above the header, and its Accession/Gene/Cysteine columns hold
    space-delimited parallel lists (one entry per isoform), so they must be exploded together.
    Depth is averaged across the UCSC transcripts reporting the same cysteine, as the funnel does.
    """
    id_cols = ["Accession", "Gene", "Cysteine"]
    return (
        pl.scan_csv(require(CONSERVATION_CSV), skip_rows=1)
        .select(*id_cols, "Depth of Conservation (# of Organisms)")
        .drop_nulls()
        .with_columns(pl.col(id_cols).str.split(" "))
        .explode(id_cols)
        .filter(pl.col("Accession").is_in(set(uniprots)))
        .with_columns(pl.col("Cysteine").cast(pl.Int64, strict=False).alias("pos"))
        .drop_nulls("pos")
        .group_by(["Accession", "pos"])
        .agg(
            pl.col("Depth of Conservation (# of Organisms)")
            .mean()
            .alias("conservation_depth")
        )
        .rename({"Accession": "uniprot"})
        .collect()
    )


def load_alphamissense(uniprots):
    """Mean AlphaMissense pathogenicity per (uniprot, position), across substitutions.

    Averaging every substitution at a position *before* thresholding is what the funnel does, so
    the shipped per-variant ``am_class`` column is deliberately not reused -- it would not match.
    """
    return (
        pl.scan_csv(require(ALPHAMISSENSE_TSV))
        .filter(pl.col("uniprot").is_in(set(uniprots)))
        .with_columns(
            pl.col("protein_variant")
            .str.slice(1, pl.col("protein_variant").str.len_chars() - 2)
            .cast(pl.Int64, strict=False)
            .alias("pos")
        )
        .drop_nulls("pos")
        .group_by(["uniprot", "pos"])
        .agg(pl.col("am_pathogenicity").mean())
        .collect()
    )


def load_ppse(uniprots):
    """AlphaFold pPSE (solvent accessibility) and secondary structure per (uniprot, position)."""
    return (
        pl.scan_csv(require(PPSE_CSV), ignore_errors=True)
        .rename({"protein_id": "uniprot", "position": "pos"})
        .filter(pl.col("uniprot").is_in(set(uniprots)))
        .select("uniprot", "pos", "pPSE", "structure_group")
        .collect()
    )


def load_clinvar(genes):
    """Pathogenic ClinVar variants for the given gene symbols, with their protein position.

    ClinVar keys on gene symbol rather than accession. The protein position is parsed out of the
    HGVS ``p.`` term in ``Name``; rows without one are dropped.
    """
    position = (
        pl.when(pl.col("Name").str.contains("p.", literal=True))
        .then(
            pl.col("Name")
            .str.split("p.")
            .list.last()
            .str.split(")")
            .list.first()
            .str.extract(r"(\d+)", 1)
            .cast(pl.Int64, strict=False)
        )
        .otherwise(None)
        .alias("clinvar_position")
    )
    return (
        pl.scan_csv(
            require(CLINVAR_TXT),
            separator="\t",
            truncate_ragged_lines=True,
            ignore_errors=True,
        )
        .filter(
            pl.col("GeneSymbol").is_in(set(genes)),
            pl.col("Assembly") == "GRCh37",
            pl.col("ClinicalSignificance").is_in(CLINVAR_PATHOGENIC),
        )
        .with_columns(position)
        .drop_nulls("clinvar_position")
        .select(
            pl.col("GeneSymbol").alias("protein"),
            "clinvar_position",
            pl.col("ClinicalSignificance").alias("clinical_significance"),
            pl.col("PhenotypeList").alias("phenotype_list"),
        )
        .collect()
    )


def _clinvar_per_residue(residues):
    """Pathogenic ClinVar variants within CLINVAR_RESIDUE_TOLERANCE of each cysteine."""
    clinvar = load_clinvar(residues["protein"].unique())
    return (
        residues.join(clinvar, on="protein", how="inner")
        .filter(
            (pl.col("pos") - pl.col("clinvar_position")).abs()
            < CLINVAR_RESIDUE_TOLERANCE
        )
        .group_by(["uniprot", "pos"])
        .agg(
            pl.len().alias("clinvar_pathogenic_variants"),
            pl.col("phenotype_list").unique().sort().str.join("|").alias("clinvar_phenotypes"),
        )
    )


def build_cys_cross_ref(residues, path=CACHE_CSV):
    """Annotate ``residues`` (uniprot, protein, pos) from the reference DBs and cache the result.

    Scans the large sources once. Returns one row per (uniprot, pos) with conservation depth,
    AlphaMissense pathogenicity, and ClinVar proximity counts.
    """
    residues = residues.select("uniprot", "protein", "pos").unique()
    uniprots = residues["uniprot"].unique()
    out = (
        residues.select("uniprot", "pos")
        .join(load_conservation(uniprots), on=["uniprot", "pos"], how="left")
        .join(load_alphamissense(uniprots), on=["uniprot", "pos"], how="left")
        .join(_clinvar_per_residue(residues), on=["uniprot", "pos"], how="left")
        .with_columns(pl.col("clinvar_pathogenic_variants").fill_null(0))
        .sort(["uniprot", "pos"])
    )
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    out.write_csv(path)
    return out


def load_cys_cross_ref(residues=None, path=CACHE_CSV, rebuild=False):
    """Cached per-cysteine annotations, building them on first use.

    ``residues`` defaults to every cysteine in the current reactivity run. The cache covers the
    study's cysteines only; pass ``rebuild=True`` (or run ``python -m src.cross_ref --rebuild``)
    after the reactivity run changes.
    """
    if rebuild or not Path(path).exists():
        return build_cys_cross_ref(_study_residues() if residues is None else residues, path)
    return pl.read_csv(path)


def _study_residues():
    """Every cysteine of every protein in the current reactivity run."""
    from src.reactivity import load_rc_long, rc_positions

    long_df = load_rc_long()
    rows = {}
    for u, p, r in zip(long_df["uniprot"], long_df["protein"], long_df["residue"]):
        for pos in rc_positions(r):
            rows[(u, p, pos)] = None
    return pl.DataFrame(
        {
            "uniprot": [k[0] for k in rows],
            "protein": [k[1] for k in rows],
            "pos": [k[2] for k in rows],
        }
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--rebuild", action="store_true", help="rescan the reference databases"
    )
    args = parser.parse_args()

    out = load_cys_cross_ref(rebuild=args.rebuild)
    print(f"Wrote {CACHE_CSV}: {out.shape[0]} rows x {out.shape[1]} columns")
    for col in ["conservation_depth", "am_pathogenicity", "clinvar_pathogenic_variants"]:
        print(f"  {col:30} non-null: {out[col].drop_nulls().len()}")


if __name__ == "__main__":
    main()
