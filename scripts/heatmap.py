#!/usr/bin/env python3
"""
Heatmap Generation for Metastasis in Osteosarcomas Analysis

This script creates heatmaps to visualize expression patterns of differentially
expressed genes between normal osteoblastic cells and metastatic osteosarcoma cell lines.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
import argparse
import os

def create_heatmap(expression_file, deseq_file=None, output_dir="figures", 
                   top_genes=50, figsize=(12, 10)):
    """
    Create heatmaps from expression data.
    
    Parameters:
    - expression_file: Path to CSV file with expression matrix
    - deseq_file: Optional path to DESeq2 results to select top genes
    - output_dir: Directory to save the figure
    - top_genes: Number of top differentially expressed genes to include
    - figsize: Figure size tuple
    """
    
    # Read expression data
    expr_df = pd.read_csv(expression_file)
    
    # Select top differentially expressed genes if DESeq2 results provided
    if deseq_file:
        deseq_df = pd.read_csv(deseq_file)
        # Sort by adjusted p-value and select top genes
        deseq_df = deseq_df.dropna(subset=['padj'])
        deseq_df = deseq_df.sort_values('padj')
        top_gene_ids = deseq_df.head(top_genes)['gene_id'].tolist()
        
        # Filter expression data to top genes
        expr_df = expr_df[expr_df['gene_id'].isin(top_gene_ids)]
    else:
        # Use first N genes if no DESeq2 results
        expr_df = expr_df.head(top_genes)
    
    # Prepare data for heatmap
    gene_names = expr_df['gene_name'].tolist() if 'gene_name' in expr_df.columns else expr_df['gene_id'].tolist()
    expression_cols = [col for col in expr_df.columns if col not in ['gene_id', 'gene_name']]
    expression_matrix = expr_df[expression_cols].values
    
    # Z-score normalization across samples for each gene
    expression_matrix_zscore = (expression_matrix - expression_matrix.mean(axis=1, keepdims=True)) / expression_matrix.std(axis=1, keepdims=True)
    
    # Create sample annotations
    sample_annotations = []
    for col in expression_cols:
        if 'Normal' in col:
            sample_annotations.append('Normal')
        else:
            sample_annotations.append('Metastatic')
    
    # Create the figure with subplots
    fig = plt.figure(figsize=figsize)
    
    # Main heatmap
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 8], width_ratios=[8, 1], 
                         hspace=0.05, wspace=0.05)
    
    # Sample annotation bar
    ax_sample_annot = fig.add_subplot(gs[0, 0])
    sample_colors = [0 if x == 'Normal' else 1 for x in sample_annotations]
    ax_sample_annot.imshow([sample_colors], aspect='auto', cmap='coolwarm')
    ax_sample_annot.set_xlim(0, len(expression_cols))
    ax_sample_annot.set_xticks([])
    ax_sample_annot.set_yticks([])
    ax_sample_annot.set_ylabel('Condition', rotation=0, ha='right', va='center')
    
    # Main heatmap
    ax_heatmap = fig.add_subplot(gs[1, 0])
    
    # Create heatmap
    im = ax_heatmap.imshow(expression_matrix_zscore, cmap='RdBu_r', aspect='auto', 
                          vmin=-3, vmax=3)
    
    # Set labels
    ax_heatmap.set_xticks(range(len(expression_cols)))
    ax_heatmap.set_xticklabels(expression_cols, rotation=90, ha='center')
    ax_heatmap.set_yticks(range(len(gene_names)))
    ax_heatmap.set_yticklabels(gene_names, fontsize=8)
    ax_heatmap.set_xlabel('Samples')
    ax_heatmap.set_ylabel('Genes')
    
    # Colorbar
    ax_colorbar = fig.add_subplot(gs[1, 1])
    cbar = plt.colorbar(im, cax=ax_colorbar)
    cbar.set_label('Z-score', rotation=270, labelpad=15)
    
    # Add title
    fig.suptitle(f'Expression Heatmap: Top {len(gene_names)} Differentially Expressed Genes\nNormal vs Metastatic Osteosarcoma', 
                fontsize=14, fontweight='bold')
    
    # Create legend for sample annotations
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='blue', label='Normal'),
                      Patch(facecolor='red', label='Metastatic')]
    ax_sample_annot.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5))
    
    # Save the plot
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, 'heatmap.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'heatmap.pdf'), bbox_inches='tight')
    
    print(f"Heatmap saved to: {output_file}")
    print(f"Showing expression patterns for {len(gene_names)} genes across {len(expression_cols)} samples")
    
    return plt.gcf()

def create_clustered_heatmap(expression_file, deseq_file=None, output_dir="figures", 
                           top_genes=50, figsize=(14, 10)):
    """
    Create a clustered heatmap with dendrograms.
    """
    
    # Read and prepare data (same as above)
    expr_df = pd.read_csv(expression_file)
    
    if deseq_file:
        deseq_df = pd.read_csv(deseq_file)
        deseq_df = deseq_df.dropna(subset=['padj'])
        deseq_df = deseq_df.sort_values('padj')
        top_gene_ids = deseq_df.head(top_genes)['gene_id'].tolist()
        expr_df = expr_df[expr_df['gene_id'].isin(top_gene_ids)]
    else:
        expr_df = expr_df.head(top_genes)
    
    gene_names = expr_df['gene_name'].tolist() if 'gene_name' in expr_df.columns else expr_df['gene_id'].tolist()
    expression_cols = [col for col in expr_df.columns if col not in ['gene_id', 'gene_name']]
    expression_matrix = expr_df[expression_cols].values
    
    # Z-score normalization
    expression_matrix_zscore = (expression_matrix - expression_matrix.mean(axis=1, keepdims=True)) / expression_matrix.std(axis=1, keepdims=True)
    
    # Create sample annotations for colors
    sample_annotations = ['Normal' if 'Normal' in col else 'Metastatic' for col in expression_cols]
    sample_colors = ['blue' if x == 'Normal' else 'red' for x in sample_annotations]
    
    # Create clustered heatmap using seaborn
    plt.figure(figsize=figsize)
    
    # Create a dataframe for seaborn
    heatmap_df = pd.DataFrame(expression_matrix_zscore, 
                             index=gene_names, 
                             columns=expression_cols)
    
    # Create the clustered heatmap
    g = sns.clustermap(heatmap_df, 
                       cmap='RdBu_r', 
                       center=0,
                       vmin=-3, vmax=3,
                       col_colors=sample_colors,
                       figsize=figsize,
                       cbar_pos=(0.02, 0.8, 0.03, 0.18))
    
    # Add title
    g.fig.suptitle(f'Clustered Heatmap: Top {len(gene_names)} DE Genes\nNormal vs Metastatic Osteosarcoma', 
                   fontsize=14, fontweight='bold', y=0.98)
    
    # Add legend for sample colors
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='blue', label='Normal'),
                      Patch(facecolor='red', label='Metastatic')]
    g.ax_heatmap.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1))
    
    # Save the plot
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, 'clustered_heatmap.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'clustered_heatmap.pdf'), bbox_inches='tight')
    
    print(f"Clustered heatmap saved to: {output_file}")
    
    return g.fig

def main():
    parser = argparse.ArgumentParser(description='Generate heatmaps from expression data')
    parser.add_argument('--expression', '-e', required=True, 
                       help='Input CSV file with expression matrix')
    parser.add_argument('--deseq', '-d', help='Optional DESeq2 results file to select top genes')
    parser.add_argument('--output', '-o', default='figures', help='Output directory for figures')
    parser.add_argument('--top-genes', '-n', type=int, default=50, 
                       help='Number of top genes to include')
    parser.add_argument('--clustered', '-c', action='store_true', 
                       help='Create clustered heatmap with dendrograms')
    
    args = parser.parse_args()
    
    if args.clustered:
        create_clustered_heatmap(args.expression, args.deseq, args.output, args.top_genes)
    else:
        create_heatmap(args.expression, args.deseq, args.output, args.top_genes)

if __name__ == "__main__":
    main()