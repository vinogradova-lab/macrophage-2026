"""Pull gene -> disease/phenotype associations from the Open Targets Platform.

Why: HPO has no DISC1 annotations (its gene->phenotype links are Mendelian
OMIM/Orphanet-derived, whereas DISC1's disease links are GWAS/literature-based).
Open Targets integrates 20+ evidence sources and returns a quantitative 0-1
overall association score per disease, which maps naturally onto the reactivity
dot-plot style (x = score, bubble size = evidence breadth).

Usage:
    Run once to (re)generate the committed table the R figure reads:
        python src/opentargets_utils.py
    -> reference_dbs/phenotype_associations.csv
       (columns: gene, disease, score, n_evidence)
"""

from pathlib import Path

import pandas as pd
import requests

OPEN_TARGETS_ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"

# Genes to profile: PPP1R21 + TBCK (FERRY complex, carry the reactive cysteines)
# and DISC1 (the reason we left HPO -- see module docstring).
DEFAULT_GENES = ["PPP1R21", "TBCK", "DISC1"]

# Open Targets therapeutic-area ids for a nervous-system / psychiatric focus. Filtering
# to these before taking the top-N keeps the panel on-topic: unfiltered, DISC1's highest
# raw scores are pleiotropic GWAS hits (scoliosis, hypercholesterolemia, ...) and its
# schizophrenia link (tagged with both areas below) is buried at ~rank 21.
NEURO_PSYCH_TA = {
    "MONDO_0005071": "nervous system disorder",
    "MONDO_0002025": "psychiatric disorder",
}

# Open Targets datasources under the genetic-association data type. The evidences
# endpoint filters by datasource (not data type), so we enumerate them to count a
# gene-disease pair's genetic-evidence records (-> dot size). NB DISC1 has ~none of
# these -- its disease links are literature-mined -- so its dots stay small.
GENETIC_DATASOURCES = [
    "gwas_credible_sets", "ot_genetics_portal", "gene_burden", "eva",
    "genomics_england", "gene2phenotype", "clingen", "uniprot_literature",
    "uniprot_variants", "orphanet",
]

_SEARCH_QUERY = """
query($q: String!){
  search(queryString: $q, entityNames: ["target"]){
    hits { id name entity }
  }
}
"""

# enableIndirect:false keeps only diseases directly associated with the gene, so
# parent ontology terms don't flood the list.
_ASSOC_QUERY = """
query($id: String!, $size: Int!){
  target(ensemblId: $id){
    approvedSymbol
    associatedDiseases(page:{index:0, size:$size}, enableIndirect:false){
      rows{
        disease { id name therapeuticAreas { id } }
        score
      }
    }
  }
}
"""


def _post(query, variables):
    """POST a GraphQL query to Open Targets and return the `data` payload."""
    resp = requests.post(
        OPEN_TARGETS_ENDPOINT,
        json={"query": query, "variables": variables},
        timeout=60,
    )
    resp.raise_for_status()
    payload = resp.json()
    if "errors" in payload:
        raise RuntimeError(f"Open Targets GraphQL error: {payload['errors']}")
    return payload["data"]


def symbol_to_ensembl(symbol):
    """Resolve a gene symbol to its Ensembl gene id via the Open Targets search.

    Prefers an exact (case-insensitive) symbol match among the target hits;
    falls back to the first target hit if none matches exactly.
    """
    hits = _post(_SEARCH_QUERY, {"q": symbol})["search"]["hits"]
    if not hits:
        raise ValueError(f"No Open Targets target hit for '{symbol}'")
    for hit in hits:
        if hit["name"].upper() == symbol.upper():
            return hit["id"]
    print(f"  (no exact match for {symbol}; using first hit {hits[0]['name']})")
    return hits[0]["id"]


def fetch_associated_diseases(ensembl_id, size=50):
    """Return the associatedDiseases `rows` for an Ensembl gene id (direct only)."""
    data = _post(_ASSOC_QUERY, {"id": ensembl_id, "size": size})
    return data["target"]["associatedDiseases"]["rows"]


_EVIDENCE_COUNT_QUERY = """
query($efo: String!, $ens: [String!]!, $ds: [String!]!){
  disease(efoId: $efo){
    evidences(ensemblIds: $ens, datasourceIds: $ds, size: 1){ count }
  }
}
"""


def count_genetic_evidence(ensembl_id, efo_id):
    """Number of genetic-association evidence records for a gene-disease pair."""
    data = _post(
        _EVIDENCE_COUNT_QUERY,
        {"efo": efo_id, "ens": [ensembl_id], "ds": GENETIC_DATASOURCES},
    )
    disease = data["disease"]
    return disease["evidences"]["count"] if disease else 0


def build_phenotype_table(symbols=DEFAULT_GENES, output_path=None, top_n=8,
                          therapeutic_area_ids=None):
    """Build a tidy gene -> top disease-association table and write it to CSV.

    For each symbol: resolve to Ensembl, fetch its direct disease associations,
    optionally keep only diseases whose therapeutic areas intersect
    `therapeutic_area_ids` (e.g. NEURO_PSYCH_TA), then keep the `top_n` by overall
    score, emitting one row per (gene, disease) with:
        gene       -- input symbol
        disease    -- Open Targets disease name
        score      -- overall association score (0-1), -> dot x-position
        n_genetic  -- number of genetic-association evidence records, -> dot size

    When filtering, a larger page is fetched first so enough on-topic diseases
    remain to fill top_n.
    """
    ta_filter = set(therapeutic_area_ids) if therapeutic_area_ids else None
    page_size = 250 if ta_filter else max(top_n * 3, 50)

    records = []
    for symbol in symbols:
        print(f"Querying Open Targets for {symbol}...")
        ensembl_id = symbol_to_ensembl(symbol)
        rows = fetch_associated_diseases(ensembl_id, size=page_size)
        if ta_filter is not None:
            rows = [
                r for r in rows
                if ta_filter & {ta["id"] for ta in r["disease"]["therapeuticAreas"]}
            ]
        rows = sorted(rows, key=lambda r: r["score"], reverse=True)[:top_n]
        for r in rows:
            n_genetic = count_genetic_evidence(ensembl_id, r["disease"]["id"])
            records.append(
                {
                    "gene": symbol,
                    "disease": r["disease"]["name"],
                    "score": r["score"],
                    "n_genetic": n_genetic,
                }
            )

    table = pd.DataFrame.from_records(records)
    if output_path is not None:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(output_path, index=False)
        print(f"Wrote {len(table)} rows ({table['gene'].nunique()} genes) to {output_path}")
    return table


if __name__ == "__main__":
    ref_dir = Path(__file__).resolve().parents[1] / "reference_dbs"
    # Unfiltered: top diseases by raw score (all therapeutic areas).
    build_phenotype_table(DEFAULT_GENES, ref_dir / "phenotype_associations.csv")
    # Nervous-system / psychiatric focus (keeps DISC1's schizophrenia link on-panel).
    build_phenotype_table(DEFAULT_GENES, ref_dir / "phenotype_associations_neuro.csv",
                          therapeutic_area_ids=NEURO_PSYCH_TA)
