"""
24/11/2025
KEGG Pathway Enrichment Visualization
Authors:
- Dulce Alejandra Carrillo Carlos
- Andrea Villarruel Garcia
- Antonio Platas Renteral

This script:
1. Define KEGG pathway lists for UP and DOWN differentially expressed genes:
    - Pathway name
    - Gene count
    - Functional category

    Example to obtain principal pathways:  tail -n +2 ../results/KEGG_output/DOWN_genes_pathways.tsv | cut -f3 | sort | uniq -c | sort -nr | head -n 25

2. Assign ordered categorical labels to control plotting order
3. Generate two types of visualizations for both UP and DOWN pathways:
    - Horizontal barplot of gene counts
    - Bubble plot (scatter plot) with bubble size proportional to gene count
4. Save resulting figures as:
    - top_pathways_UP.png
    - top_pathways_DOWN.png

Usage:
    python3 pathways_plot.py
    
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os


#Up regulated pathways, their frequency and category
up_counts = [
    ("Metabolic pathways", 14, "Metabolism"),
    ("Pathways in cancer", 9, "Cancer"),
    ("Rap1 signaling pathway", 8, "Signaling"),
    ("Wnt signaling pathway", 7, "Signaling"),
    ("Regulation of actin cytoskeleton", 7, "Cytoskeleton"),
    ("MAPK signaling pathway", 7, "Signaling"),
    ("PI3K-Akt signaling pathway", 6, "Signaling"),
    ("Cytoskeleton in muscle cells", 6, "Cytoskeleton"),
    ("Cushing syndrome", 6, "Hormone"),
    ("Calcium signaling pathway", 6, "Signaling"),
    ("TGF-beta signaling pathway", 5, "Signaling"),
    ("Ras signaling pathway", 5, "Signaling"),
    ("Phospholipase D signaling pathway", 5, "Signaling"),
    ("Pathways of neurodegeneration", 5, "Neurodegeneration"),
    ("Motor proteins", 5, "Cytoskeleton"),
    ("Transcriptional misregulation in cancer", 4, "Cancer"),
    ("Signaling pathways regulating pluripotency of stem cells", 4, "Signaling"),
    ("Proteoglycans in cancer", 4, "Cancer"),
    ("Prostate cancer", 4, "Cancer"),
    ("MicroRNAs in cancer", 4, "Cancer"),
    ("Hypertrophic cardiomyopathy", 4, "Cardiovascular"),
    ("Hormone signaling", 4, "Hormone"),
    ("Circadian rhythm", 4, "Hormone"),
    ("Circadian entrainment", 4, "Hormone"),
    ("Basal cell carcinoma", 4, "Cancer"),
]

#Down regulated pathways, their frequency and category
down_counts = [
    ("Coronavirus disease - COVID-19", 9, "Virus/Infection"),
    ("Influenza A", 6, "Virus/Infection"),
    ("Epstein-Barr virus infection", 6, "Virus/Infection"),
    ("Cytosolic DNA-sensing pathway", 6, "Signaling"),
    ("AGE-RAGE signaling pathway in diabetic complications", 6, "Signaling"),
    ("TNF signaling pathway", 5, "Signaling"),
    ("Salmonella infection", 5, "Virus/Infection"),
    ("RIG-I-like receptor signaling pathway", 5, "Signaling"),
    ("Pathways of neurodegeneration", 5, "Neurodegeneration"),
    ("Pathogenic Escherichia coli infection", 5, "Virus/Infection"),
    ("NOD-like receptor signaling pathway", 5, "Signaling"),
    ("Measles", 5, "Virus/Infection"),
    ("Lipid and atherosclerosis", 5, "Metabolism"),
    ("Human cytomegalovirus infection", 5, "Virus/Infection"),
    ("Herpes simplex virus 1 infection", 5, "Virus/Infection"),
    ("Hepatitis B", 5, "Virus/Infection"),
    ("Complement and coagulation cascades", 5, "Signaling"),
    ("Yersinia infection", 4, "Virus/Infection"),
    ("Toll-like receptor signaling pathway", 4, "Signaling"),
    ("Shigellosis", 4, "Virus/Infection"),
    ("Pertussis", 4, "Virus/Infection"),
    ("Pathways in cancer", 4, "Cancer"),
    ("NF-kappa B signaling pathway", 4, "Signaling"),
    ("Metabolic pathways", 4, "Metabolism"),
    ("MAPK signaling pathway", 4, "Signaling"),
]


def make_df(count_list):
    df = pd.DataFrame(count_list, columns=["Pathway_Name", "Gene_Count", "Category"])
    df["Occurrences"] = df["Gene_Count"]
    return df

up_df = make_df(up_counts)
down_df = make_df(down_counts)

#Color Palette
palette = {
    "Virus/Infection": "#3CC0F0", 
    "Signaling": "#FFD52E", 
    "Metabolism": "#88DB5E",
    "Cancer": "#F45B96", 
    "Cytoskeleton": "#57A1F7",
    "Hormone": "#FF7B2E", 
    "Neurodegeneration": "#B309E8",
    "Cardiovascular": "#00B4D8",
}

#Category Order
category_order = ["Virus/Infection", "Signaling", "Metabolism", "Cancer",
                  "Cytoskeleton", "Hormone", "Neurodegeneration",
                  "Cardiovascular"]

def order_df_by_category(df):
    """
    Sort pathways within a DataFrame according to a predefined biological
    category order.

    Steps:
    1. Convert the 'Category' column into an ordered pandas Categorical,
       using the global 'category_order' list.
    2. Sort pathways first by Category (following that order) and then
       alphabetically by Pathway_Name.

    """
    df["Category"] = pd.Categorical(df["Category"], categories=category_order, ordered=True)
    df = df.sort_values(["Category", "Pathway_Name"])
    return df

up_df = order_df_by_category(up_df) 
down_df = order_df_by_category(down_df)

#Plot function
def plot_bar_and_bubble(df, title, filename, output_dir="figures"):
    os.makedirs(output_dir, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(19, 5))

    # Barplot
    sns.barplot(
        data=df,
        y="Pathway_Name",
        x="Gene_Count",
        hue="Category",
        dodge=False,
        palette=palette,
        ax=axes[0]
    )
    axes[0].invert_yaxis()
    axes[0].set_title(f"{title}", fontsize=16)
    axes[0].set_xlabel("Gene number")
    axes[0].legend(title="Category", bbox_to_anchor=(1.05, 1), loc="upper left")

    # Bubble plot
    sns.scatterplot(
        data=df,
        x="Gene_Count",
        y="Pathway_Name",
        size="Occurrences",
        hue="Category",
        palette=palette,
        sizes=(100, 800),
        alpha=0.7,
        ax=axes[1],
        edgecolor="black"
    )
    axes[1].set_title(f"{title}", fontsize=16)
    axes[1].set_xlabel("Gene Number")
    axes[1].legend([],[], frameon=False)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, filename), dpi=300)
    plt.close()
    print(f"Saved: {filename}")


plot_bar_and_bubble(up_df, "Top KEGG Pathways (UP genes)", "top_pathways_UP.png")
plot_bar_and_bubble(down_df, "Top KEGG Pathways (DOWN genes)", "top_pathways_DOWN.png")