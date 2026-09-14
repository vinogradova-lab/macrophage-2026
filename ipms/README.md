# Endogenous IP-MS analysis

Endogenous pulldowns of the FERRY subunits Fy2 (PPP1R21) and Tbck (Fy1) were performed in M0 and
LPS-stimulated macrophages from three donors, each with an IgG control and technical replicates.

## Data processing

- Common contaminants and keratins were removed.
- Proteins were required to have at least two peptides.
- Proteins were required to be quantified in at least two biological replicates (donors).
- Signal intensities were normalized within each donor to the bait protein, scaling every IP
  sample so that its bait signal equals the donor's median bait signal.

## Enrichment analysis

P-values were calculated with a t-test for the means of two independent samples. Fold change was
calculated from the mean normalized signal of each condition. Three comparisons were made per bait:
IP vs. IgG in M0, IP vs. IgG in LPS, and LPS IP vs. M0 IP. Proteins with p < 0.05 (unadjusted) and
a > 1.5 fold change were called significant. Proteins significantly up over IgG in either state
were considered IP-enriched.

Proteins were annotated as FERRY subunits, ribosomal (UniProt GO terms containing "ribosom"), or
RNA-binding proteins (Census of RBPs, Gerstberger et al., 2014).

## CRAPome filtering

Each quantified protein was queried against the CRAPome v2.0 database (716 human control
experiments). Proteins found in more than 200 experiments were flagged as contaminants.

## GO enrichment

CRAPome-filtered IP-enriched proteins were tested for GO enrichment (BP/MF/CC) with a Fisher exact
test against the full human background, using the MSigDB C5 GO collection (see root README).

## Comparison to published FERRY interactomes

IP-enriched proteins were compared with the FERRY-complex enriched proteins from the GST-FERRY pulldown
(Schuhmacher et al., 2023) and with the Tbck, Fy2 and C12orf4 AP-MS interactomes (Wang et al.,
Tables S11–S13). Overlaps are shown as three-way Venn diagrams.
