#!/usr/bin/env python3
"""
Volcano Plot Generation for Metastasis in Osteosarcomas Analysis

This script creates volcano plots to visualize differential expression results
from DESeq2 analysis comparing normal osteoblastic cells vs metastatic osteosarcoma cell lines.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import argparse
import os

def create_volcano_plot(data_file, output_dir="figures", significance_threshold=0.05, 
                       fold_change_threshold=2.0, figsize=(10, 8)):
    """
    Create a volcano plot from DESeq2 results.
    
    Parameters:
    - data_file: Path to CSV file with DESeq2 results
    - output_dir: Directory to save the figure
    - significance_threshold: Adjusted p-value threshold for significance
    - fold_change_threshold: Log2 fold change threshold
    - figsize: Figure size tuple
    """
    
    # Read the data
    df = pd.read_csv(data_file)
    
    # Remove rows with NaN values
    df = df.dropna(subset=['padj', 'log2FoldChange'])
    
    # Calculate -log10(padj)
    df['neg_log10_padj'] = -np.log10(df['padj'])
    
    # Define significance categories
    df['significance'] = 'Not significant'
    df.loc[(df['padj'] < significance_threshold) & 
           (df['log2FoldChange'] > np.log2(fold_change_threshold)), 'significance'] = 'Upregulated'
    df.loc[(df['padj'] < significance_threshold) & 
           (df['log2FoldChange'] < -np.log2(fold_change_threshold)), 'significance'] = 'Downregulated'
    
    # Create the plot
    plt.figure(figsize=figsize)
    
    # Define colors for each category
    colors = {'Not significant': 'lightgray', 'Upregulated': 'red', 'Downregulated': 'blue'}
    
    # Plot points by significance category
    for category in ['Not significant', 'Upregulated', 'Downregulated']:
        subset = df[df['significance'] == category]
        plt.scatter(subset['log2FoldChange'], subset['neg_log10_padj'], 
                   c=colors[category], alpha=0.6, s=30, label=category)
    
    # Add significance threshold lines
    plt.axhline(y=-np.log10(significance_threshold), color='black', linestyle='--', alpha=0.5)
    plt.axvline(x=np.log2(fold_change_threshold), color='black', linestyle='--', alpha=0.5)
    plt.axvline(x=-np.log2(fold_change_threshold), color='black', linestyle='--', alpha=0.5)
    
    # Customize the plot
    plt.xlabel('Log2 Fold Change', fontsize=12)
    plt.ylabel('-Log10(Adjusted P-value)', fontsize=12)
    plt.title('Volcano Plot: Normal vs Metastatic Osteosarcoma', fontsize=14, fontweight='bold')
    plt.legend(frameon=True, fancybox=True, shadow=True)
    plt.grid(True, alpha=0.3)
    
    # Add summary statistics to the plot
    n_up = len(df[df['significance'] == 'Upregulated'])
    n_down = len(df[df['significance'] == 'Downregulated'])
    n_total = len(df)
    
    textstr = f'Total genes: {n_total}\nUpregulated: {n_up}\nDownregulated: {n_down}'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    plt.text(0.02, 0.98, textstr, transform=plt.gca().transAxes, fontsize=10,
             verticalalignment='top', bbox=props)
    
    # Save the plot
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, 'volcano_plot.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'volcano_plot.pdf'), bbox_inches='tight')
    
    print(f"Volcano plot saved to: {output_file}")
    print(f"Summary: {n_up} upregulated, {n_down} downregulated genes (padj < {significance_threshold}, |log2FC| > {np.log2(fold_change_threshold):.2f})")
    
    return plt.gcf()

def main():
    parser = argparse.ArgumentParser(description='Generate volcano plot from DESeq2 results')
    parser.add_argument('--input', '-i', required=True, help='Input CSV file with DESeq2 results')
    parser.add_argument('--output', '-o', default='figures', help='Output directory for figures')
    parser.add_argument('--pvalue-threshold', '-p', type=float, default=0.05, 
                       help='Adjusted p-value threshold')
    parser.add_argument('--fold-change-threshold', '-f', type=float, default=2.0, 
                       help='Fold change threshold')
    
    args = parser.parse_args()
    
    create_volcano_plot(args.input, args.output, args.pvalue_threshold, args.fold_change_threshold)

if __name__ == "__main__":
    main()