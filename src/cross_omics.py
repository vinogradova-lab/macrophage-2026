"""The cross-omics GO-term Venn diagrams (Extended Data Fig. 5h/5j), in one place.


A change set is a set of gene symbols, matching the ``protein`` column all three sources key on.
"""

from functools import lru_cache

import polars as pl

from src.macrophage import PROTEOMICS, RC_LONG_FORMAT, REPO_ROOT, require

# MSigDB C5 GO 2026.1.Hs, exported by whole_proteome/export_msigdb_go.R. The same file the R
# go_enrich() draws term membership from, so venn and enrichment gene sets are identical.
MSIGDB_GO_CSV = REPO_ROOT / "reference_dbs" / "msigdb_go_2026.csv"

# The three omics, all TLR4 vs. M0. Each entry is (label, source CSV, filter expressions); the
# labels are the Venn's set labels and become the membership columns in Data S2-8.
#
# NB phospho is both directions ("Significant Up"/"Significant Down"). The separate up-only set
# belongs to the phospho volcano's own GO:BP enrichment (Figure 3g / Extended Data Fig. 5d), not
# here -- panels 5e/5i plot both directions, so the Venns and the cross-omics enrichment do too.
_DIRECTIONS = ["Significant Up", "Significant Down"]
_PHOSPHO_REGULATION = "Regulation_phospho - TLR4 vs. M0 (9591 Residues)"

OMIC_SOURCES = [
    (
        "Reactivity changes (TLR4)",
        RC_LONG_FORMAT,
        [pl.col("reactivity_change"), pl.col("condition").eq("TLR4")],
    ),
    (
        "Phosphorylated proteins in TLR4 (FC > 2, p < 0.05)",
        PROTEOMICS
        / "05_phosphoproteomics/processed_results"
        / "expression-normalized_phosphorylation_table.csv",
        [pl.col(_PHOSPHO_REGULATION).is_in(_DIRECTIONS)],
    ),
    (
        "Expression changes in TLR4 (FC > 1.5, p < 0.05)",
        PROTEOMICS
        / "01_unenriched proteomics/03_results/20250219_3reps/04_results/volcano_plots"
        / "volcano_data_long_format.csv",
        [pl.col("condition").eq("TLR4"), pl.col("Regulation").is_in(_DIRECTIONS)],
    ),
]

OMIC_NAMES = [name for name, _, _ in OMIC_SOURCES]


def omic_change_sets(restrict=None):
    """``{omic label: set of gene symbols with a significant TLR4-vs-M0 change}``.

    ``restrict`` (a gene set) additionally intersects each omic with that GO term, which is exactly
    what the Venn does; omitting it gives the whole-omic denominators.
    """
    out = {}
    for name, path, filters in OMIC_SOURCES:
        conds = list(filters)
        if restrict is not None:
            conds.append(pl.col("protein").is_in(restrict))
        out[name] = set(pl.read_csv(require(path)).filter(*conds)["protein"])
    return out


ANNOTATION_COLS = ["uniprot", "description", "uniprot_function"]


@lru_cache(maxsize=1)
def protein_annotation():
    """``protein -> (uniprot, description, uniprot_function)``, from the omic sources themselves.

    All three carry uniprot and description; only reactivity and phospho carry the UniProt
    function blurb, so a protein seen exclusively in the expression data has that field blank.
    First non-null wins, so no UniProt lookup (and no entry cache) is needed at build time.
    """
    frames = []
    for _, path, _ in OMIC_SOURCES:
        present = pl.read_csv(require(path), n_rows=0).columns
        wanted = [c for c in ANNOTATION_COLS if c in present]
        df = pl.read_csv(path, columns=["protein"] + wanted).with_columns(
            pl.lit(None, dtype=pl.String).alias(c) for c in ANNOTATION_COLS if c not in wanted
        )
        # read_csv returns file order, which differs between the sources; pin it for the concat.
        frames.append(df.select("protein", *ANNOTATION_COLS))
    return (
        pl.concat(frames)
        .drop_nulls("protein")
        .group_by("protein")
        .agg(pl.col(c).drop_nulls().first() for c in ANNOTATION_COLS)
    )


@lru_cache(maxsize=1)
def _msigdb_go():
    """The MSigDB GO membership table (GO_id, term, gene), read once per process."""
    return pl.read_csv(
        require(MSIGDB_GO_CSV, "Run whole_proteome/export_msigdb_go.R to regenerate it."),
        columns=["GO_id", "term", "gene"],
    )


def go_term_genes(go_id):
    """Gene symbols of one MSigDB C5 GO set, by GO id."""
    return set(_msigdb_go().filter(pl.col("GO_id") == go_id)["gene"])


def go_term_genes_by_name(terms):
    """Union of gene symbols across MSigDB GO sets matched by (case-insensitive) term name.

    Mirrors go_genes_terms() in phosphoproteomics/visualization.Rmd so the Venn uses the same
    multi-term gene set as the phospho volcano's highlighting.
    """
    wanted = [t.lower() for t in terms]
    return set(_msigdb_go().filter(pl.col("term").str.to_lowercase().is_in(wanted))["gene"])


# "gtpase regulator activity" has no MSigDB term and is silently dropped, here and in
# go_genes_terms(); the union resolves to the other five.
GTPASE_TERMS = [
    "gtpase activity",
    "gtpase activator activity",
    "gtpase regulator activity",
    "gtpase activating protein binding",
    "gtpase inhibitor activity",
    "gtpase binding",
]

# The terms with a Venn panel, in figure order. ``source`` is the human-readable gene-set
# provenance carried into Data S2-8: GTPase activity is a six-term union rather than GO:0003924
# alone, so its circles (39/73/126) must not be reconciled against that single term's enrichment.
VENN_TERMS = {
    "Endosomal transport": {
        "genes": lambda: go_term_genes("GO:0016197"),
        "source": "GO:0016197 (GOBP endosomal transport)",
        "figure": "Extended Data Fig. 5f",
    },
    "GTPase activity": {
        "genes": lambda: go_term_genes_by_name(GTPASE_TERMS),
        "source": (
            "Union of 5 MSigDB GTPase GO terms (GO:0003924 and related), "
            "not GO:0003924 alone"
        ),
        "figure": "Extended Data Fig. 5h",
    },
    "Vesicle-mediated transport": {
        "genes": lambda: go_term_genes("GO:0016192"),
        "source": "GO:0016192 (GOBP vesicle-mediated transport)",
        "figure": "Extended Data Fig. 5k",
    },
    # Drawn by the notebook for context but not placed in a panel, so it is excluded from S2-8.
    "Endocytosis": {
        "genes": lambda: go_term_genes("GO:0006897"),
        "source": "GO:0006897 (GOBP endocytosis)",
        "figure": None,
    },
}


def venn_membership(term):
    """Binary membership table for one VENN_TERMS entry: one row per protein in the union of the
    three omic sets, with a 0/1 column per omic. Column order matches OMIC_NAMES, so the per-column
    sums are the Venn's circle totals."""
    sets = omic_change_sets(restrict=VENN_TERMS[term]["genes"]())
    proteins = sorted(set().union(*sets.values()))
    return pl.DataFrame({"protein": proteins}).with_columns(
        pl.col("protein").is_in(sets[name]).cast(pl.Int8).alias(name) for name in OMIC_NAMES
    )
