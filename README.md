# Macrophage multi-omics analysis

Analysis code for the macrophage TLR/STING stimulation study. Each modality is a
self-contained folder pairing a Python processing notebook with an R (`.Rmd`) figure step.

## 1. Environment

```bash
conda env create -f environment.yml   # env "polars" (Python 3.13)
conda activate polars
nbstripout --install --attributes .gitattributes   # strip notebook outputs on commit
```

`environment.yml` pins deps and installs `src/` via `-e .` (see `setup.py`), so
`import src.uniprot_utils` etc. work from any notebook. Register the kernel if needed:
`python -m ipykernel install --user --name polars`.

R is used for the `.Rmd` figure steps and Fisher GO enrichment. Install the packages sourced by
`macrophage_figure_setting.R` (notably `msigdbr` 26.1.0, `clusterProfiler`, `ggplot2`, `tidyverse`).

## 2. Repository structure

Shared code:

- `src/`: installable package (`uniprot_utils`, `peptide_funcs`, `pca_utils`, `opentargets_utils`)
- `macrophage_figure_setting.R`: shared R theme + GO enrichment (`go_enrich()` / `.go_msigdb()`)
- `reference_dbs/`: shared annotations (`msigdb_go_2026.csv`, GO semdata, phospho sites)

Per-modality folders, each with its own `README.md`:

- `whole_proteome/`: unenriched proteome expression, volcano/GO downstream
- `reactivity/`: cysteine reactivity (isoTOP-ABPP)
- `ipms/`: endogenous IP-MS (interactomes, CRAPome filtering, Venn)
- `phosphoproteomics/`: phospho-site quantification, kinome, enrichment
- `rna/`: DESeq2 normalization, RNA/GSEA
- `metabolomics/`: metabolite abundance

## 3. GO annotation source

Functional gene-list curation and Fisher GO enrichment share one ontology source: the
**MSigDB C5 GO collection, 2026.1.Hs** (`msigdbr` 26.1.0).

- `whole_proteome/export_msigdb_go.R` writes `reference_dbs/msigdb_go_2026.csv`
  (`ontology, GO_id, term, gene`; 10,490 terms, 19,591 genes), the single source of truth,
  also consumed by `go_enrich()`.
- `wp_funcs.query_gene_ontology_msigdb()` curates gene lists: GO-id tokens resolve to the
  DAG-propagated MSigDB set; other tokens match term name/definition.
- Category flags in `volcano_data_long_format.csv` (Autophagy, Mitochondrial, etc.) derive
  from these sets.

Caveats: `GO:0006914` (autophagy) and `GO:0006096` (glycolysis) aren't MSigDB gene sets, so they
resolve via the term-name tokens `"autophagy"` / `"glycolytic"`. The `simplify()` step in
`go_enrich()` still uses 2021 `org.Hs.eg.db`/`GO.db` IC for semantic tie-breaking only.
