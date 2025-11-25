"""
24/11/2025
RNA-seq Processing Pipeline
Authors:
- Dulce Alejandra Carrillo Carlos
- Andrea Villarruel Garcia
- Antonio Platas Renteral

This script: 
1. Load raw count data and map samples to experimental conditions
2. Preprocess counts with log2(x + 1) transformation
3. Generate QC plots:
    - Boxplot of sample distributions.
    - Density plots per sample and condition.
    - PCA for sample clustering visualization.
4. Compute Counts Per Million (CPM) and filter lowly expressed genes
5. Perform differential expression analysis using DESeq2 (via pyDESeq2):
    - Save full results, significant genes, and separate UP/DOWN gene lists
6. Create a volcano plot with UP, DOWN, and non-significant genes

Usage:
    python3 differential_expression.py \
        --counts ../data/raw_counts.tsv \
        --output ../results/DESeq2_output
"""

import os
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from matplotlib.patches import Patch
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

def load_and_prepare_counts(path: str, sample_to_condition: dict):
    """ Load raw counts, apply log2(x+1) to numeric columns and add a condition row.
     - loads and preprocesses a raw counts table,
     - produces QC plots (boxplot, density, PCA),
     - applies CPM filtering,
     - runs differential expression with pyDESeq2,
     - creates a volcano plot and saves results.

    Args:
        path (str): Path to raw counts table (tab-separated, gene samples)
        sample_to_condition (dict): Mapping sample_name -> condition string

    """

    df_raw = pd.read_csv(path, sep="\t", index_col=0)

    # keep non-numeric first column as in original script
    non_num = df_raw.iloc[:, :1]
    numeric = df_raw.iloc[:, 1:].astype(float)

    # log2(x + 1) transform for visualization
    log_num = np.log2(numeric + 1)
    df_log = pd.concat([non_num, log_num], axis=1)

    # build condition row in the exact format used previously
    condition_row = pd.DataFrame(
        [["Condition"] + [sample_to_condition.get(c, "Unknown") for c in df_log.columns if c != "Name"]],
        columns=df_log.columns,
        index=["ID"],
    )

    df_log_with_condition = pd.concat([condition_row, df_log])

    return df_raw, df_log_with_condition, numeric


def plot_boxplot(df_log_with_condition: pd.DataFrame, figures_dir: str, colors_map: dict):
    """ Create and save the samples boxplot.

    Args:
        df_log_with_condition (pd.DataFrame): DataFrame created by load_and_prepare_counts
        figures_dir (str): Directory to save figure
        colors_map (dict): Mapping conditions
    """
    os.makedirs(figures_dir, exist_ok=True)

    # extract conditions and data
    conditions = df_log_with_condition.iloc[0, 1:].to_dict()
    df_data = df_log_with_condition.iloc[1:, :].copy()

    samples = df_data.columns[1:]
    data = [df_data[s].astype(float).values for s in samples]
    box_colors = [colors_map.get(conditions[s], "#C0C0C0") for s in samples]

    fig, ax = plt.subplots(figsize=(12, 6))

    box = ax.boxplot(
        data,
        patch_artist=True,
        widths=0.6,
        showfliers=False,
        medianprops=dict(color="black", linewidth=1.5),
        boxprops=dict(linewidth=1.5),
        whiskerprops=dict(linewidth=1.5),
        capprops=dict(linewidth=1.5),
    )

    for patch, color in zip(box["boxes"], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_xticks(range(1, len(samples) + 1))
    ax.set_xticklabels(samples, rotation=0, fontsize=12)
    ax.set_ylabel("Log2(x + 1)", fontsize=14, fontweight="bold")
    ax.tick_params(axis="y", labelsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.5)
    ax.spines["bottom"].set_linewidth(1.5)

    legend_elements = [
        Patch(facecolor=colors_map.get(k), alpha=0.7, label=k) for k in sorted(colors_map.keys())
    ]
    ax.legend(handles=legend_elements, loc="upper right", title="Condition", fontsize=10, title_fontsize=12)
    plt.tight_layout()

    out = os.path.join(figures_dir, "samples_boxplot.png")
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Boxplot saved to: {out}")


def plot_density(df_log_with_condition: pd.DataFrame, figures_dir: str, sample_colors: dict):
    """ Create and save density plot per sample type. 

    Args:
        df_log_with_condition (pd.DataFrame): DataFrame from load_and_prepare_counts
        figures_dir (str): Directory to save figure
        sample_colors (dict): Mapping sample -> color hex
    """
    os.makedirs(figures_dir, exist_ok=True)

    df_melt = df_log_with_condition.iloc[1:, :].melt(id_vars=["Name"], var_name="Sample", value_name="Value")
    conditions = df_log_with_condition.iloc[0, 1:].to_dict()
    df_melt["Cond"] = df_melt["Sample"].map(conditions)
    df_melt["Value"] = df_melt["Value"].astype(float)

    graf = sns.displot(
        data=df_melt,
        x="Value",
        hue="Sample",
        col="Cond",
        kind="kde",
        fill=True,
        alpha=0.5,
        palette=sample_colors,
        height=5,
        aspect=1.2,
    )

    graf.set_axis_labels("Log2(x + 1)", "Density")
    graf.fig.suptitle("Log2 coverage distribution by condition and sample", fontsize=16)
    graf.tight_layout()
    graf.fig.subplots_adjust(top=0.88)

    out = os.path.join(figures_dir, "samples_densityplot.png")
    graf.savefig(out, dpi=300)
    plt.close()
    print(f"Density plot saved to: {out}")


def plot_pca(df_log_with_condition: pd.DataFrame, figures_dir: str, colors_map: dict):
    """ Compute PCA on numeric data and save a PC1 vs PC2 scatter. 

    Args:
        df_log_with_condition (pd.DataFrame): DataFrame from load_and_prepare_counts
        figures_dir (str): Directory to save plot
        colors_map (dict): Mapping condition -> color hex

    Returns:
        sklearn.decomposition.PCA: fitted PCA object (returned for availability of explained variance)
    """
    os.makedirs(figures_dir, exist_ok=True)

    df_data = df_log_with_condition.iloc[1:, 1:].astype(float)
    expr = df_data.T
    expr.index.name = "Sample"

    pca = PCA(n_components=2)
    pcs = pca.fit_transform(expr)

    conditions = df_log_with_condition.iloc[0, 1:].to_dict()
    pca_df = pd.DataFrame(pcs, columns=["PC1", "PC2"], index=expr.index)
    pca_df["Condition"] = pca_df.index.map(conditions)

    plt.figure(figsize=(8, 6))
    markers_map = {"Control": "o", "Tratamiento": "*"},

    sns.scatterplot(
        x="PC1",
        y="PC2",
        data=pca_df,
        hue="Condition",
        palette=colors_map,
        s=100,
        style="Condition",
        markers = {"Control": "o", "Osteosarcoma": "*"},
        edgecolor="black",
    )

    plt.title("PCA", fontsize=16, fontweight="bold")
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f}%)", fontsize=12)
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f}%)", fontsize=12)

    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0, fontsize=10, title="Condition")
    plt.grid(True)
    plt.tight_layout()

    out = os.path.join(figures_dir, "PCA_samples.png")
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"PCA plot saved to: {out}")

    return pca

def compute_and_filter_cpm(df_numeric: pd.DataFrame):
    """ Compute CPM and filter genes.
    - Keep genes with CPM >= MIN_CPM in at least MIN_SAMPLES samples.

    Returns:
        tuple: (df_cpm, df_filtered, passed)
            df_cpm: CPM-normalized DataFrame
            df_filtered: CPM-filtered count matrix (gene samples)
            passed: boolean Series indicating which genes passed the filter
    
    """

    MIN_CPM = 2
    MIN_SAMPLES = 3

    # Compute CPM
    df_cpm = (df_numeric / df_numeric.sum()) * 1_000_000

    # Boolean passed
    passed = (df_cpm >= MIN_CPM).sum(axis=1) >= MIN_SAMPLES

    # Filtered matrix
    df_filtered = df_numeric[passed].copy()

    return df_cpm, df_filtered, passed

def run_deseq2(df_counts: pd.DataFrame, sample_to_condition: dict, outdir: str) -> pd.DataFrame:
    """ Run DESeq2 differential expression using pyDESeq2.

    Args:
        df_counts (pd.DataFrame): Filtered integer count matrix (genes samples).
        sample_to_condition (dict): Sample → condition mapping.
        outdir (str): Output directory.

    Returns:
        pd.DataFrame: DESeq2 results dataframe.

    """
    metadata = pd.DataFrame({"condition": [sample_to_condition.get(s, "Unknown") for s in df_counts.columns]})
    metadata.index = df_counts.columns

    dds = DeseqDataSet(counts=df_counts.T, metadata=metadata, design_factors="condition")
    dds.deseq2()

    ds = DeseqStats(dds, contrast=["condition", "Osteosarcoma", "Control"])
    ds.summary()
    results = ds.results_df

    results.to_csv(os.path.join(outdir, "DESeq2_all_results.tsv"), sep="\t")
    return results

def volcano_plot(df: pd.DataFrame, outpath: str):
    """ Volcano plot

    Args:
        df (pd.DataFrame): DE results including log2FoldChange and padj.
        outpath (str): Path to save the figure.

    """
    df = df.copy()
    df["log10Neg"] = -np.log10(df["padj"].clip(lower=1e-300))

    palette = {"UP": "#FF6B6B", "DOWN": "#4ECDC4", "Non-DE": "#C0C0C0"}

    plt.figure(figsize=(7, 8))
    sns.scatterplot(data=df, x="log2FoldChange", y="log10Neg", hue="Expresion", palette=palette)
    plt.xlabel("log2(Fold Change)")
    plt.ylabel("-log10(padj)")
    plt.title("Volcano Plot")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Volcano plot saved to: {outpath}")

def plot_sig_heatmap(df_log_with_condition, sig_genes, figures_dir):

    os.makedirs(figures_dir, exist_ok=True)

    df_log = df_log_with_condition.iloc[1:, :].copy()
    df_sig = df_log[df_log["Name"].isin(sig_genes)].copy()

    if df_sig.empty:
        print("No significant genes found for heatmap.")
        return

    df_sig = df_sig.set_index("Name")
    mat = df_sig.astype(float)

    row_means = mat.mean(axis=1).to_numpy()[:, None]
    row_stds  = mat.std(axis=1).to_numpy()[:, None]
    mat_z = (mat.to_numpy() - row_means) / row_stds
    mat_z = pd.DataFrame(mat_z, index=mat.index, columns=mat.columns)

    g = sns.clustermap(
        mat_z,
        cmap="RdYlBu_r",
        xticklabels=True,
        yticklabels=False,
        figsize=(28, 25),
        dendrogram_ratio=(0.15, 0.1),
        cbar_pos=(0.92, 0.2, 0.02, 0.6)
    )

    # Column labels
    g.ax_heatmap.set_xticklabels(
        g.ax_heatmap.get_xticklabels(),
        fontsize=28,
        rotation=0
    )

    # Colorbar
    g.ax_cbar.tick_params(labelsize=26)

    # Title
    g.fig.suptitle(
        "Heatmap of Significant Genes",
        fontsize=40,
        fontweight="bold",
        y=1.02
    )

    # Layout Adjustments
    g.fig.subplots_adjust(right=0.9)

    # Cbar position 
    g.ax_cbar.set_position((0.93, 0.25, 0.02, 0.5))

    #Save
    out = os.path.join(figures_dir, "SignificantGenes_Heatmap.png")
    g.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Heatmap saved to: {out}")


def main():
    parser = argparse.ArgumentParser(description="Modular RNA-seq pipeline with QC plots and DESeq2")
    parser.add_argument("--counts", required=True, help="Path to raw count table")
    parser.add_argument("--output", required=True, help="Output directory")

    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    figures_dir = os.path.join(args.output, "Figures")
    os.makedirs(figures_dir, exist_ok=True)

    sample_to_condition = {
        "ERR13451008": "Osteosarcoma",
        "ERR13448321": "Control",
        "ERR13451021": "Osteosarcoma",
        "ERR13448329": "Osteosarcoma",
        "ERR13451054": "Osteosarcoma",
        "ERR13447502": "Control",
    }

    colors_map = {"Osteosarcoma": "#FC4E92", "Control": "#02A4D1"}

    sample_colors = {
        "ERR13447502": "#02A4D1",
        "ERR13448321": "#FFD93D",
        "ERR13448329": "#845EC2",
        "ERR13451008": "#FC4E92",
        "ERR13451021": "#4D96FF",
        "ERR13451054": "#FFA94D",
    }

    df_raw, df_log_with_condition, numeric_part = load_and_prepare_counts(args.counts, sample_to_condition)

    # QC plots (boxplot, density, PCA)
    plot_boxplot(df_log_with_condition, figures_dir, colors_map)
    plot_density(df_log_with_condition, figures_dir, sample_colors)
    pca = plot_pca(df_log_with_condition, figures_dir, colors_map)

    df_cpm, df_filtered, passed = compute_and_filter_cpm(numeric_part)

    passed.to_csv(os.path.join(args.output, "FilterCount.csv"))
    print("Genes after filter:", df_filtered.shape[0])

    # DESeq2
    df_int = df_filtered.astype(int)
    results = run_deseq2(df_int, sample_to_condition, args.output)

    # add expression label used for volcano
    results["Expresion"] = np.select(
        [
            (results["padj"] <= 0.01) & (results["log2FoldChange"] >= 2),
            (results["padj"] <= 0.01) & (results["log2FoldChange"] <= -2),
        ],
        ["UP", "DOWN"],
        default="Non-DE",
    )

    # Volcano
    volcano_plot(results, os.path.join(figures_dir, "VolcanoPlot.png"))

    # Save DE results and significant subsets
    results.to_csv(os.path.join(args.output, "DESeq2_all_results.tsv"), sep="\t")
    sig = results[results["padj"] < 0.05]
    sig_genes = sig.index.tolist()
    plot_sig_heatmap(df_log_with_condition, sig_genes, figures_dir)
    sig.to_csv(os.path.join(args.output, "DESeq2_significant_genes.tsv"), sep="\t")
    up = sig[sig["log2FoldChange"] > 0]
    down = sig[sig["log2FoldChange"] < 0]
    up.to_csv(os.path.join(args.output, "DESeq2_UP_genes.tsv"), sep="\t")
    down.to_csv(os.path.join(args.output, "DESeq2_DOWN_genes.tsv"), sep="\t")

    print(f"Saved results in: {args.output}")
    print(f"Total significant DE genes: {len(sig)} (UP: {len(up)}, DOWN: {len(down)})")

    # add expression label used for volcano
    results["Expresion"] = np.select(
        [
            (results["padj"] <= 0.01) & (results["log2FoldChange"] > 1),
            (results["padj"] <= 0.01) & (results["log2FoldChange"] < -1),
        ],
        ["UP", "DOWN"],
        default="Non-DE"
    )

    # Save volcano categories
    results.to_csv(os.path.join(args.output, "DESeq2_all_results.tsv"), sep="\t")

    # Volcano plot
    volcano_plot(results, os.path.join(figures_dir, "VolcanoPlot.png"))

    # Load significant gene list to use in heatmap
    sig_path = "/export/space3/users/carrillo/LCG/BIOPYTHON/Metastasis-in-Osteosarcomas/results/DESeq2_output/DESeq2_significant_genes.tsv"

    if os.path.exists(sig_path):
        sig_df = pd.read_csv(sig_path, sep="\t")
        sig_genes = sig_df.iloc[:, 0].astype(str).tolist()

        # Remove version suffix
        sig_genes = [g.split(".")[0] for g in sig_genes]

        print(f"Loaded {len(sig_genes)} significant genes for heatmap.")
    else:
        print(f"WARNING: No significant file found at: {sig_path}")
        sig_genes = []

    # Ensure gene IDs in matrix also have version removed
    df_log_with_condition.index = df_log_with_condition.index.str.split(".").str[0]

    # Generate heatmap
    plot_sig_heatmap(df_log_with_condition, sig_genes, figures_dir)

    print("Pipeline completed successfully.")

 
if __name__ == "__main__":
    main()