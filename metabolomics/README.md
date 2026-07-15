# Metabolomics analysis

## Data processing:

Four independent biological replicates were processed by the metabolomics core.

Within each biological replicate, metabolites were required to have at least one condition with an average signal intensity above 5000. Out of metabolites provided by the metabolomics core, 179 passed signal intensity filtering. After filtering and log transformation, missing data was imputed through quantile regression imputation of left-censored data (QRILC) performed using imputeLCMD version 2.1. Channel ratios (signal intensity / sum of signal intensities per metabolite) were calculated and used for all analyses.

## Heatmap:

Only proteins quantified in all replicates are displayed. A z-score was calculated for each observation by z = (x-μ)/σ where x  is the channel ratio, μ is the mean channel ratio of the metabolite across all replicates, and σ  is the standard deviation of the channel ratio of the metabolite across all replicates. Hierarchical clustering was performed using ComplexHeatmap (Gu et al., 2022) library version 2.10.0 for R version 4.1.1. The distance matrix was calculated using the “euclidean” method and clustering was performed using the “complete” method. Colors represent Z-score of metabolite channel ratio. 

Metabolite class annotations were curated by the Vardhana lab through literature review.

## Principal component analysis:

Metabolite ratio values were log2 transformed and only proteins with a channel ratio above 0 in all channels were used. PC1 and PC2 were plotted. Principal component analysis was performed with scikit-learn version 1.4.2 for Python version 3.1.1. Loadings for plotting were selected by taking the 5 metabolites with highest and lowest loadings for both PC1 and PC2. 

## Volcano plots: 

P-values were calculated with T-test for the means of two independent samples. Volcano plots show p-value on y axis and fold change of the condition versus control on x axis. Significant proteins (< 0.05, -log10(0.05) = 1.3) which show a > 1.5 fold change are highlighted.
 
