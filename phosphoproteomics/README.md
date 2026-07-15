# Phosphoproteomics analysis details:

## Data filtering

Database search: Experiment HK3_14_phospho was searched with UniProt 2024 canonical database, with differential modifications on S,T,Y (one modification search)

Filtering:  Phosphoproteomic data was filtered using standard parameters:
- 2 missed cleavages
- Unique peptides only
- Average signal intensity > 5000
- CV < 0.5
Similar parameters used in Kemper 2022 analysis

PTM aggregation: take mean of all peptides per modified residue

## Principal component analysis

Channel ratio values were normalized to whole proteome and then log2 transformed before PCA. Top 5 loadings in each direction for PC1 and 2 are shown.

## Motif enrichment analysis 

Differential Serine/Threonine (303 kinases) motif enrichment analysis using the Kinase Library
Default parameters were used for enrichment analysis:
Fold change of 2 (**after normalization to whole proteome**) was used to select differentially regulated phosphosites (will filter on p-value as well once we have biological replicates) 
Enrichment Threshold of 15 was used to determine kinases-substrate interactions
