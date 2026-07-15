"""Build the reference-list -> GO-terms methods table for the whole_proteome analysis.

Materializes, as a flat CSV, which GO terms define each functional "reference protein
list" used to annotate the whole-proteome data. The list definitions are the
``search_strings`` dict in wp_downstream_analysis.ipynb (cell 13), copied verbatim
below, plus the HALLMARK-based "Glycolysis" list. GO term names come from go-basic.obo
(the source named in protein_lists/README.md); the 4 mitochondrial CC ids that
MSigDB omits fall back to the goatools search log.

Output: reference_list_go_terms.csv
    columns: Reference list name, GO terms used, GO term names   (alphabetized)
"""

import re
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
MSIGDB_CSV = REPO / "reference_dbs" / "msigdb_go_2026.csv"
SEARCH_LOG = HERE / "search_logs" / "mitochondri_CC_gos_human.log"
OUT_CSV = HERE / "reference_list_go_terms.csv"

# The GTPase list is a name/definition substring search; a ("HALLMARK", set_id)
# tuple marks an MSigDB hallmark set (not GO-defined).
GTPASE = "__gtpase_substring__"

# GO-based lists copied verbatim from wp_downstream_analysis.ipynb cell 13
# (search_strings); the MSigDB HALLMARK sets are the ones fed to plot_gsea() in
# rna_analysis/plot_GSEA.Rmd (Glycolysis is shared with the cell-13 annotation).
REFERENCE_LISTS = {
    "Electron transport chain": ["GO:0022900"],
    "Mitochondrial": [
        "GO:0005739",
        "GO:0031966",
        "GO:0005759",
        "GO:0098800",
        "GO:0005758",
        "GO:0044233",
        "GO:0005740",
        "GO:0042645",
        "GO:0098799",
        "GO:0005746",
    ],
    "Mitochondrial translation": ["GO:0032543"],
    "Vesicle-mediated transport": [
        "GO:0060627",
        "GO:0048193",
        "GO:0006888",
        "GO:0042147",
    ],
    "GTPase activity": GTPASE,
    "Glycolysis": ("HALLMARK", "HALLMARK_GLYCOLYSIS"),
    "Interferon gamma response": ("HALLMARK", "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
    "TNFα signaling via NF-κB": ("HALLMARK", "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
}
# domains only matter for GTPase (BP, MF), per cell 13.
GTPASE_DOMAINS = {"BP", "MF"}


def load_msigdb():
    """The exact table query_gene_ontology_msigdb() reads."""
    return pd.read_csv(MSIGDB_CSV, usecols=["ontology", "GO_id", "term", "description"])


def resolve_gtpase_ids(msigdb):
    """Reproduce query_gene_ontology_msigdb's match for the token 'GTPase':
    ontology in {BP, MF}, term OR description contains 'gtpase' (case-insensitive)."""
    df = msigdb[msigdb["ontology"].isin(GTPASE_DOMAINS)]
    hit = df["term"].str.contains("gtpase", case=False, regex=False, na=False)
    hit = hit | df["description"].str.contains("gtpase", case=False, regex=False, na=False)
    return sorted(df.loc[hit, "GO_id"].unique().tolist())


def build_name_lookup(msigdb):
    """GO id -> name. Prefer authoritative go-basic.obo names; fall back to the
    goatools search log (covers the 4 mito ids MSigDB lacks) and finally MSigDB term."""
    lookup = {}

    # 1) go-basic.obo (best: hyphenated, complete). Download+cache; skip on failure.
    try:
        from goatools.base import download_go_basic_obo
        from goatools.obo_parser import GODag

        obo = HERE / "go-basic.obo"
        if not obo.exists():
            download_go_basic_obo(str(obo))
        dag = GODag(str(obo))
        for gid, term in dag.items():
            lookup[gid] = term.name
    except Exception as exc:  # network unavailable etc. -> offline fallbacks below
        print(f"go-basic.obo unavailable ({exc}); using offline name sources")

    # 2) goatools search log: 'MATCH GO:0031966(mitochondrial membrane) name: ...'
    if SEARCH_LOG.exists():
        for gid, name in re.findall(r"MATCH (GO:\d+)\(([^)]+)\)", SEARCH_LOG.read_text()):
            lookup.setdefault(gid, name)

    # 3) MSigDB term (last resort).
    for gid, term in msigdb.drop_duplicates("GO_id").set_index("GO_id")["term"].items():
        lookup.setdefault(gid, term)

    return lookup


def main():
    msigdb = load_msigdb()
    names = build_name_lookup(msigdb)

    rows = []
    for list_name, definition in REFERENCE_LISTS.items():
        if isinstance(definition, tuple) and definition[0] == "HALLMARK":
            rows.append(
                {
                    "Reference list name": list_name,
                    "GO terms used": "—",
                    "GO term names": f"MSigDB {definition[1]} (not a GO-defined list)",
                }
            )
            continue

        go_ids = resolve_gtpase_ids(msigdb) if definition == GTPASE else definition
        missing = [g for g in go_ids if g not in names]
        if missing:
            print(f"WARNING: no name resolved for {missing} in '{list_name}'")
        rows.append(
            {
                "Reference list name": list_name,
                "GO terms used": ", ".join(go_ids),
                "GO term names": ", ".join(names.get(g, g) for g in go_ids),
            }
        )

    out = pd.DataFrame(rows).sort_values("Reference list name").reset_index(drop=True)
    out.to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV}  ({len(out)} rows)")
    for _, r in out.iterrows():
        n = 0 if r["GO terms used"] == "—" else r["GO terms used"].count(",") + 1
        print(f"  {r['Reference list name']}: {n} GO terms")


if __name__ == "__main__":
    main()
