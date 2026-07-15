# Reactivity analysis

## Data processing

Each experiment was processed individually using following filters:

- Peptides with more than 2 internal missed cleavage sites were removed.
- Half-tryptic peptides were removed.
- Peptides with low average of reporter ion intensities (<5,000 for all technical replicates) were removed ("floating control")
- Peptides with high variation between all technical replicate channels (CV >0.5) were removed ("floating control")


## Reactivity change analysis

The following criteria were used to identify cysteine reactivity changes:

Data quality filters:
- Only peptides quantified in at least two replicates were considered
- Only proteins with at least two peptides were considered
- Percent control values for channels with a mean signal intensity < 5000 were dropped if the corresponding M0 value was also < 5000
- Residues that were quantified in two replicates where one of the replicates had < 1.5 fold change and the ratio between the two replicates > 2 were not considered reactivity changes (**passes_two_rep_variability_filter**).
- Proteins manually annotated as low quality ("red") were not considered reactivity changes

Reactivity change criteria:
- At least one cysteine within 1.5-fold of the expression value (**passes_bio_replicate_variation_filter**) or max/min ratio greater than 3
- At least one cysteines is changing less than 2-fold or max/min ratio greater than 3
- The cysteine must have a ratio to another cysteine on the same protein greater than 2-fold in at least 2 replicates (**passes_ratio_in_two_replicates**). If a peptide was not detected in one replicate, the average value for that peptide from other replicates is used when applying this filter.

- If the cysteine has whole proteome expression data, the cysteine must have a ratio to the whole proteome value greater than 2. (**passes_expression_filter**)
- Proteins with fewer than five peptides and no expression data are listed as (Cys1;Cys2) and are counted as 1 reactivity change
- For proteins with two-four peptides, the max/min ratio must be greater than 2. For proteins with five+ peptides, the cysteine must have >2 fold change from the median cysteine value for that protein.  (**passes_ratio_or_median_filter**). 

