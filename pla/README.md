# PLA and colocalization analysis

Proximity ligation assay (PLA) and immunofluorescence colocalization of Fy2 with Fy1, Fy4 and Rab5
in M0 and LPS-stimulated macrophages.

## Statistics

All comparisons use superplots (Lord et al., 2020). Per-cell values were collapsed to one value
per donor and condition (mean, or median as a sensitivity check). P-values were calculated with a
paired t-test across donors, never across cells. Bars show the mean ± SEM of donor values.

## PLA

The per-cell PLA particle count (Ch1) was compared between M0 and LPS. As a robustness check, the
test was repeated without the two date-named samples (20260513, 20260616).

## Colocalization

Six donors were pooled from two image deliveries (D1–D3, D4–D6). The D4–D6 images were corrected
for chromatic aberration in the Fy2 channel. Objects smaller than 0.15 µm were excluded. The split
colocalization coefficient was used: intensity-weighted, threshold 0.25, random overlap subtracted.
Both directions were analyzed for each pair: Coloc(X/Y) is the fraction of X signal on Y. Values
were also shown as a percentage of each donor's M0 value, annotated with the raw-value p-value.
Vesicle count, count per unit cell area, and mean vesicle size were compared the same way.

## Sensitivity to image-processing settings

D4–D6 were quantified three ways: uncorrected with the 0.15 µm gate, corrected with the gate, and
corrected without it. The effect of each setting was tested within M0 and within LPS by a paired
t-test across donors, with a cell-level paired test as a technical check. The M0 vs. LPS test was
also repeated under each setting and compared by p-value and effect size.
