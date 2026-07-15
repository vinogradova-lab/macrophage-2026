# RNA analysis

RNA count table was shared by Beatrice Zhang.

Counts were normalized by DESeq2's normalization procedure. 

## DEG

DEG was performed by Beatrice Zhang.

## PCA

Normalized count values were log2 transformed before PCA. Only protein coding genes were used for PCA. The top 5,000 genes based on row wise variance were used for PCA. PCA was performed using scikit-learn (Pedregosa et al., 2011) library for Python.