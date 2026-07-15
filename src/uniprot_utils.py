import json
from pathlib import Path

from unipressed import UniprotkbClient


def create_entry_cache(df, uniprot_col="uniprot", cache_path=None, chunk_size=500):
    """Fetch UniProtKB entries for every unique accession in df[uniprot_col],
    returned as a dict keyed by primaryAccession.

    If cache_path is given and the file exists, the cache is loaded from it instead
    of querying UniProt. If cache_path is given and the file does NOT exist, the
    queried result is saved there. cache_path=None preserves the old in-memory-only
    behavior.
    """
    if cache_path is not None and Path(cache_path).exists():
        with open(cache_path) as f:
            return json.load(f)

    uniprots = sorted(
        u for u in set(df[uniprot_col]) if u is not None and u != "None"
    )
    entry_dict = {}
    print("Querying UniProt...")
    for start in range(0, len(uniprots), chunk_size):
        print(f"Retrieved {start} entries out of {len(uniprots)}")
        for entry in UniprotkbClient.fetch_many(uniprots[start:start + chunk_size]):
            entry_dict[entry["primaryAccession"]] = entry
    print("Done querying UniProt.")

    if cache_path is not None:
        with open(cache_path, "w") as f:
            json.dump(entry_dict, f)
    return entry_dict


def get_function(acc, cache):
    """UniProt FUNCTION comment text(s) for an accession, '|'-joined."""
    ret = set()
    if acc in cache and "comments" in cache[acc]:
        for c in cache[acc]["comments"]:
            if "texts" in c and c["commentType"] == "FUNCTION":
                for t in c["texts"]:
                    ret.add(t["value"])
    return "|".join(ret)


def get_go_terms(acc, cache):
    """GO terms for an accession as 'Term (GO:id)', '|'-joined."""
    ret = []
    if acc in cache and "uniProtKBCrossReferences" in cache[acc]:
        for db in cache[acc]["uniProtKBCrossReferences"]:
            if db["database"] == "GO":
                for p in db["properties"]:
                    if p["key"] == "GoTerm":
                        ret.append(f'{p["value"].split(":")[1]} ({db["id"]})')
    return "|".join(ret)


def build_gene_to_uniprot(cache):
    """Map gene symbol -> primary accession from a UniProt entry cache."""
    return {
        cache[acc]["genes"][0]["geneName"]["value"]: acc
        for acc in cache
        if cache[acc].get("genes") and cache[acc]["genes"][0].get("geneName")
    }
