"""
ANALYSIS UTILS

This module contains analysis methods used across the project.
Methods are called from metabolomics_analysis.ipynb.
"""

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import scipy.stats as stat
import plotly.express as px
import math


color_discrete_map = {
    "Significant Up": "#ff8080",
    "Not Significant Up": "#ffcccc",
    "Significant but <1.5 FC": "darkgrey",
    "Not Significant": "lightgrey",
    "Significant Down": "#71A0C6",
    "Not Significant Down": "#CDDEFA",
}

def get_pca_plot(df, index_cols, title_string, out_dir, with_log = False):
    """Principal component analysis.
    Scales data with log. Runs and plots PCA with scikit learn
    """
    out_dir.mkdir(exist_ok=True)
    # drop na
    df = df.dropna()
    df = df.set_index(index_cols)
    # log2 transform
    if not with_log:
        df_log = np.log2(df)
    else:
        df_log = df

    df_log = df_log.replace(-np.inf, np.nan)
    df_log = df_log.dropna()

    number_of_proteins_in_common = df_log.shape[0]

    # transpose table
    df_log_t = df_log.transpose()
    df_log_t.reset_index(inplace=True)
    df_log_t = df_log_t.rename(columns={"index": "channel_name"})

    # get X and y
    X = df_log_t.drop("channel_name", axis=1)

    # get PCA 3, fit transform, get df
    pca = PCA(n_components=3)
    principalComponents = pca.fit_transform(X)
    principalDf = pd.DataFrame(
        data=principalComponents,
        columns=[
            "PC1",
            "PC2",
            "PC3",
        ],
    )
    # add channel name
    finalDf = pd.concat([principalDf, df_log_t[["channel_name"]]], axis=1)

    # get PCA %
    pca_1_percent = round((pca.explained_variance_ratio_[0] * 100), 2)
    pca_2_percent = round((pca.explained_variance_ratio_[1] * 100), 2)
    pca_3_percent = round((pca.explained_variance_ratio_[2] * 100), 2)

    percent_df = pd.DataFrame(
        [
            {"principal_component": "PC1", "percent_explained": pca_1_percent},
            {"principal_component": "PC2", "percent_explained": pca_2_percent},
            {"principal_component": "PC3", "percent_explained": pca_3_percent},
        ]
    )

    # add condition
    finalDf["condition"] = finalDf["channel_name"]

    finalDf["condition"] = finalDf["condition"].str.split("_").str[0]
    color_discrete_sequence_list = [
        "salmon",
        "#A3A533",
        "#56BC82",
        "#4EADF0",
        "#D772EC",
    ]
    title_info = title_string + " (" + str(number_of_proteins_in_common) + " Proteins)"

    fig = px.scatter(
        finalDf,
        x="PC1",
        y="PC2",
        hover_data=["channel_name"],
        color=finalDf["condition"],
        color_discrete_sequence=color_discrete_sequence_list,
        labels={
            "principal component 1": "PC1 ({}%)".format(pca_1_percent),
            "principal component 2": "PC2 ({}%)".format(pca_2_percent),
            "condition": "Condition",
        },
        title=title_info,
        template="plotly_white",
    )

    fig.update_xaxes(
        dtick=5,
        showline=True,
        linewidth=2,
        linecolor="black",
        mirror=True,
        showgrid=True,
        gridwidth=2,
        gridcolor="#E8E8E8",
        zerolinewidth=2,
        zerolinecolor="#E8E8E8",
    )

    fig.update_yaxes(
        dtick=5,
        showline=True,
        linewidth=2,
        linecolor="black",
        mirror=True,
        showgrid=True,
        gridwidth=2,
        gridcolor="#E8E8E8",
        zerolinewidth=2,
        zerolinecolor="#E8E8E8",
    )

    fig.update_traces(marker=dict(size=12), selector=dict(mode="markers"))

    fig.update_layout(height=600, width=600, showlegend=True, legend_title_text="")

    loadings_df = pd.DataFrame(
        data=np.transpose(pca.components_), columns=["PC1", "PC2", "PC3"]
    )
    loadings_df["variable"] = df_log.index.tolist()

    loadings_df = loadings_df.set_index("variable")

    finalDf.to_csv(out_dir / "pca_results.csv")
    loadings_df.to_csv(out_dir / "loadings_results.csv")
    percent_df.to_csv(out_dir / "percent_explained.csv")

    return (fig, loadings_df, finalDf, percent_df)


def get_p_value(row, cond_1, cond_2):
    ttest_result = stat.ttest_ind(row[cond_1], row[cond_2], nan_policy="omit")
    return ttest_result[1]


def get_expr(row):
    """define regulation of row based on fold change and p value"""
    p_value = row["-log10_pval"]
    p_value_cutoff = -1 * math.log10(0.05)
    fc_cutoff = math.log2(1.5)
    if (row["log2_FC"] > fc_cutoff) & (p_value > p_value_cutoff):
        return "Significant Up"
    if (row["log2_FC"] > fc_cutoff) & (p_value < p_value_cutoff):
        return "Not Significant Up"
    if (row["log2_FC"] < -fc_cutoff) & (p_value > p_value_cutoff):
        return "Significant Down"
    if (row["log2_FC"] < -fc_cutoff) & (p_value < p_value_cutoff):
        return "Not Significant Down"
    if (
        (row["log2_FC"] > -fc_cutoff)
        & (row["log2_FC"] < fc_cutoff)
        & (p_value > p_value_cutoff)
    ):
        return "Significant but <1.5 FC"
    else:
        return "Not Significant"
    
def get_volcano_plot_treatment_vs_control(
    conditions_list,
    control_labelling,
    df,
    file_name,
    folder_path,
    index_cols,
    target_name,
):
    volcano_plot_list = []
    volcano_df_list = []
    df = df.set_index(index_cols)

    for condition in conditions_list:
        copy_df = df.copy()
        list_cond_1 = copy_df.filter(like=condition).columns.tolist()
        list_cond_2 = copy_df.filter(like=control_labelling).columns.tolist()
        if len(list_cond_1) <= 1 or len(list_cond_2) <= 1:
            print(
                f"Condition {condition} does not have enough replicates to be shown in volcano plot!"
            )
            print(list_cond_1)
            print(list_cond_2)
            continue

        for labelling in [control_labelling, condition]:
            list_columns_labelling = copy_df.filter(like=labelling).columns.tolist()
            labelling_df = copy_df[list_columns_labelling]
            copy_df["mean_" + labelling] = labelling_df.mean(axis=1)

        copy_df["FC"] = (
            copy_df["mean_" + condition] / copy_df["mean_" + control_labelling]
        )
        x_axis_name = "log2(" + condition + "/" + control_labelling + ")"

        idx_cond_1 = copy_df.columns.get_indexer(list_cond_1)
        idx_cond_2 = copy_df.columns.get_indexer(list_cond_2)
        copy_df["p_value"] = copy_df.apply(
            get_p_value,
            axis=1,
            args=(idx_cond_1, idx_cond_2),
        )

        copy_df["log2_FC"] = np.log2(copy_df["FC"])
        copy_df = copy_df.replace(-np.inf, -7)
        copy_df = copy_df.replace(np.inf, 7)
        volcano_df = copy_df[["p_value", "log2_FC"]]

        volcano_df = volcano_df.dropna()
        volcano_df["-log10_pval"] = -1 * np.log10(volcano_df["p_value"])
        # Using BH correction here (previous versions used Boneferroni)
        volcano_df["-log10_pval_adj"] = -1 * np.log10(
            stat.false_discovery_control(volcano_df["p_value"])
        )
        volcano_df["Regulation"] = volcano_df.apply(get_expr, axis=1)

        volcano_df["Regulation"] = volcano_df["Regulation"].astype("category")

        volcano_df = volcano_df.reset_index()

        title_name = (
            file_name
            + " - "
            + control_labelling
            + " vs. "
            + condition
            + " ("
            + str(len(volcano_df))
            + " "
            + target_name
            + ")"
        )

        fig = px.scatter(
            volcano_df,
            x="log2_FC",
            y="-log10_pval",
            color="Regulation",
            color_discrete_map=color_discrete_map,
            hover_data=index_cols,
            title=title_name.replace("processed_census-out_", ""),
            labels={"log2_FC": x_axis_name},
            template="simple_white",
            category_orders={"Regulation": np.sort(volcano_df["Regulation"].unique())},
        )
        fig.add_vline(x=0.58, line_width=2, line_dash="dash", line_color="grey")
        fig.add_vline(x=-0.58, line_width=2, line_dash="dash", line_color="grey")
        fig.add_hline(y=1.3, line_width=2, line_dash="dash", line_color="grey")
        fig.update_layout(legend=dict(title=""), title_x=0.5, font_family="Arial")
        fig.write_image(
            folder_path
            / ("volcano_plot_" + title_name.split(" - ")[1].split(" (")[0] + ".svg"),
            engine="kaleido",
        )

        sign_up_df = (
            volcano_df.loc[volcano_df["Regulation"] == "Significant Up"]
            .sort_values(by="-log10_pval", ascending=False)
            .head(30)
        )
        sign_up_df_fc = (
            volcano_df.loc[volcano_df["Regulation"] == "Significant Up"]
            .sort_values(by="log2_FC", ascending=False)
            .head(30)
        )
        sign_down_df = (
            volcano_df.loc[volcano_df["Regulation"] == "Significant Down"]
            .sort_values(by="-log10_pval", ascending=False)
            .head(30)
        )
        sign_down_df_fc = (
            volcano_df.loc[volcano_df["Regulation"] == "Significant Down"]
            .sort_values(by="log2_FC", ascending=True)
            .head(30)
        )
        labels_df = pd.concat(
            [sign_up_df, sign_down_df, sign_up_df_fc, sign_down_df_fc], axis=0
        )
        labels_df = labels_df.drop_duplicates()
        labels_df = labels_df.loc[~np.isinf(labels_df["log2_FC"])]
        for i, r in labels_df.iterrows():
            if r["Regulation"] == "Significant Down":
                color = color_discrete_map["Significant Down"]
            elif r["Regulation"] == "Significant Up":
                color = color_discrete_map["Significant Up"]

            fig.add_annotation(
                x=r["log2_FC"],
                y=r["-log10_pval"],
                text=r[index_cols[0]],
                showarrow=False,
                xanchor="center",
                yanchor="bottom",
                font=dict(size=10, color=color),
            )

        volcano_plot_list.append(fig)
        volcano_df = volcano_df.reset_index().set_index(index_cols)
        volcano_df = volcano_df.add_suffix("_" + title_name)
        volcano_df_list.append(volcano_df)

    with open(folder_path / "volcano_plots.html", "w") as f:
        for fig in volcano_plot_list:
            f.write(fig.to_html(full_html=False, include_plotlyjs="cdn"))

    volcano_df = pd.concat(volcano_df_list, axis=1)
    volcano_df.to_csv(folder_path / "volcano_data.csv")
    return volcano_df, copy_df
