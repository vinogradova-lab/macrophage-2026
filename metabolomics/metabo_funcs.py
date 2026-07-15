import pandas as pd
from pathlib import Path
import math
import plotly.express as px
import scipy.stats as stat
import numpy as np
from functools import reduce
import subprocess
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import colors
from analysis_utils import *
sns.set_theme(font="Arial", style="white")

AVG_INTENSITY_CUTOFF = 5000


def add_metabolite_annotation(volcano_df, index_col="Compound"):
    """merge metabolite class annotations with quantification
    returns annotated data frame
    """
    annotation = pd.read_csv(
        "../../data/metabolomics/20240924_list of polar metabolites for the targeted assay_annotations_LSP_v3.csv"
    )
    annotation_index = "NAME"
    merged_df = volcano_df.merge(
        right=annotation,
        left_on=index_col,
        right_on=annotation_index,
        how="left",
    )
    merged_df[index_col] = volcano_df.reset_index()[index_col]
    return merged_df


def passes_intensity_cutoff(row, conditions):
    return any(
        row.filter(like=condition).mean() > AVG_INTENSITY_CUTOFF
        for condition in conditions
    )


def filter_on_intensity_cutoff(metabo_df, donors, conditions, index_cols, output_path):
    donor_datasets = []
    for donor in donors:
        labelling_df = metabo_df.copy().filter(like=donor)
        metabolites_passing = labelling_df.apply(
            passes_intensity_cutoff, axis=1, conditions=conditions
        )
        keep_df = labelling_df[metabolites_passing]
        discard_df = labelling_df[~metabolites_passing]
        discard_df.index.names = index_cols
        discard_df.to_csv(output_path / "{}_filtered_out_metabolites.csv".format(donor))
        keep_df.index.names = index_cols
        donor_datasets.append(keep_df)
    return pd.concat(donor_datasets, axis=1)


def calc_channel_ratio(metabo_df, donors, index_cols):
    donor_datasets = []
    for donor in donors:
        labelling_df = metabo_df.filter(like=donor)
        row_sums = labelling_df.sum(axis=1)
        normalized_df = labelling_df.div(row_sums, axis=0)
        normalized_df.index.names = index_cols
        donor_datasets.append(normalized_df)

    # Concatenate all donor datasets
    return pd.concat(donor_datasets, axis=1)


def make_result_dirs(parent_dir):
    formatted_data = parent_dir / "01_formatted_data"
    filtered_files = parent_dir / "02_filtered_out_metabolites"
    combined_files = parent_dir / "03_combined_files"
    results_files = parent_dir / "04_results"

    for dir in [combined_files, filtered_files, results_files, formatted_data]:
        dir.mkdir(exist_ok=True, parents=True)
    return (combined_files, filtered_files, results_files, formatted_data)


def imputed_data_histogram(
    index_cols, filtered_df, impute_lcmd_results, results_files, combined_files
):
    """ "Create a histogram showing imputed vs original data"""
    long_df = (
        filtered_df.reset_index()
        .melt(
            id_vars=index_cols,
            var_name="channel_name",
            value_name="signal_intensity",
        )
        .dropna()
        .set_index(index_cols + ["channel_name"])
    )

    long_df["log_2_signal_intensity"] = np.log2(long_df["signal_intensity"])
    long_df = long_df.replace(-np.inf, 0)
    long_df["Imputed data"] = long_df["log_2_signal_intensity"] == 0
    stat_df = long_df.replace(
        0, np.nan
    ).dropna()  # subset to non-zero data to calculate median/stddev
    missing_data = long_df[long_df["Imputed data"]].copy()
    not_missing_data = long_df[~long_df["Imputed data"]].copy()
    imputed_df = pd.concat([not_missing_data, missing_data])
    plt.rcParams["font.family"] = "Arial"
    fig, [ax1, ax2] = plt.subplots(ncols=2, figsize=(10, 5))
    sns.histplot(
        stat_df,
        x="log_2_signal_intensity",
        ax=ax1,
        hue="Imputed data",
        multiple="stack",
        palette="viridis",
        binrange=(10, 32),
        hue_order=[False, True],
        binwidth=0.5,
    )

    merge_cols = [index_cols[0]] + ["channel_name"]
    merged_df = imputed_df.reset_index().merge(
        (
            impute_lcmd_results.melt(
                id_vars=index_cols[0],
                var_name="channel_name",
                value_name="imputed_log2_signal_intensity",
            )
            .dropna()
            .set_index(merge_cols)
        ).reset_index(),
        on=merge_cols,
    )
    sns.histplot(
        merged_df,
        x="imputed_log2_signal_intensity",
        hue="Imputed data",
        multiple="stack",
        palette="viridis",
        ax=ax2,
        binrange=(10, 32),
        hue_order=[False, True],
        binwidth=0.5,
    )
    for ax in [ax1, ax2]:
        ax.set_xlabel("log2(signal intensity)")
        ax.set_xlim(9, 31)
    fig.suptitle("Distribution of metabolite signal intensities")
    ax1.set_title("Before imputation")
    ax2.set_title("After imputation through QRILC")
    (results_files / "data_imputation").mkdir(exist_ok=True)
    fig.savefig(results_files / "data_imputation" / "data_imputation_histogram.svg")
    plt.show()

    imputed_df.to_csv(results_files / "data_imputation" / "data_imputation_data.csv")

    print(
        "Median value of original data",
        np.exp2(
            np.median(
                merged_df[~merged_df["Imputed data"]]["imputed_log2_signal_intensity"]
            )
        ),
        "\n" "Median value of imputed data",
        np.exp2(
            np.median(
                merged_df[merged_df["Imputed data"]]["imputed_log2_signal_intensity"]
            )
        ),
    )

    filtered_df.to_csv(combined_files / "raw_data_with_imputed_values.csv")


def differential_expression_analysis(
    index_cols, results_files, file_name, target_name, channel_ratio_df
):
    volcano_dir = results_files / "volcano_plots"
    volcano_dir.mkdir(exist_ok=True)
    # define the comparisons we want to make, instead of every possible combination
    comparisons = {"M0": ["LPS"]}
    de_results = []
    for control, conditions in comparisons.items():
        comparison_dir = volcano_dir / "{}_comparison".format(control)
        comparison_dir.mkdir(exist_ok=True)
        volcano_df, copy_df = get_volcano_plot_treatment_vs_control(
            conditions_list=conditions,
            control_labelling=control,
            file_name=file_name,
            target_name=target_name,
            folder_path=comparison_dir,
            df=channel_ratio_df.reset_index(),
            index_cols=index_cols,
        )
        de_results.append(volcano_df)

    volcano_df = de_results[0]

    val_counts_list = []
    reg_sub_df = volcano_df.filter(regex="Regulation_")
    for col in reg_sub_df.columns:
        val_counts = reg_sub_df[col].value_counts().to_frame()
        val_counts.columns = [col]
        val_counts_list.append(val_counts)

    val_count_df = pd.concat(val_counts_list, axis=1)
    cols = []
    for col in val_count_df.columns:
        col = col.split(f"Regulation_{file_name} - ")[1].rsplit(" (")[0]
        cols.append(col)
    val_count_df.columns = cols

    new_index = [
        "Not Significant",
        "Not Significant Down",
        "Not Significant Up",
        "Significant but <1.5 FC",
        "Significant Down",
        "Significant Up",
    ]

    val_count_df = val_count_df.T
    val_count_df = val_count_df.reindex(new_index, axis=1)

    val_count_df.plot(
        kind="bar",
        stacked=True,
        color=["lightgrey", "#CDDEFA", "#ffcccc", "darkgrey", "#71A0C6", "#ff8080"],
    )
    plt.xlabel("Comparisons")
    plt.ylabel(target_name)
    plt.title("Differential Abundance Overview")
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.savefig(
        results_files / "volcano_plots" / "diff_expr_overview.svg", bbox_inches="tight"
    )

    val_count_df[["Significant Up", "Significant Down"]].plot(
        kind="bar", color=["#ff8080", "#71A0C6"]
    )
    plt.xlabel("Comparisons")
    plt.ylabel(target_name)
    plt.title("Significant Up/Down regulation Overview")
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.savefig(
        results_files / "volcano_plots" / "signupdown_overview.svg", bbox_inches="tight"
    )

    return (volcano_df, de_results)
