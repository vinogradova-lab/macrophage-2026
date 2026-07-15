"""Shared PCA core for the whole-proteome and RNA-seq pipelines.

Dataset-specific preprocessing, plotting, and CSV schemas stay in the callers.
"""

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def run_pca(matrix, n_components, log2=True):
    """Run PCA on a features (rows) by samples (columns) matrix.

    With ``log2=True`` the matrix is log2-transformed and features with any
    non-finite value are dropped (complete cases only).

    Returns ``(pcs_df, loadings_df, explained_variance_ratio)``: samples scored
    on ``PC1..PCn`` plus ``sample_name``, features loaded on ``PC1..PCn`` plus
    ``variable``, and the per-component variance fractions.
    """
    df = matrix
    if log2:
        df = np.log2(df)
    df = df.replace(-np.inf, np.nan).dropna()

    # samples as rows, features as columns
    df_t = df.transpose()
    X = df_t.values
    sample_names = df_t.index.tolist()

    pca = PCA(n_components=n_components)
    principal_components = pca.fit_transform(X)

    pc_cols = [f"PC{i + 1}" for i in range(n_components)]
    pcs_df = pd.DataFrame(data=principal_components, columns=pc_cols)
    pcs_df["sample_name"] = sample_names

    loadings_df = pd.DataFrame(data=np.transpose(pca.components_), columns=pc_cols)
    loadings_df["variable"] = df.index.tolist()

    return pcs_df, loadings_df, pca.explained_variance_ratio_


def percent_explained_df(explained_variance_ratio):
    """Variance-explained table (``principal_component``, ``percent_explained``).

    Written as ``percent_explained.csv``; the R figure code reads this schema to
    build PCA axis labels.
    """
    return pd.DataFrame(
        [
            {"principal_component": f"PC{i + 1}", "percent_explained": round(r * 100, 2)}
            for i, r in enumerate(explained_variance_ratio)
        ]
    )
