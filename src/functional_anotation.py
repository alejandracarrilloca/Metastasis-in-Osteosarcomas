"""
24/11/2025
Functional Anotation
Authors: 
  - Dulce Alejandra Carrillo Carlos 
  - Andrea Villarruel Garcia 
  - Antonio Platas Renteral 

Functional annotation pipeline for DE genes.

This script performs:
1. Conversion from Ensembl IDs to Entrez IDs
2. Retrieval of KEGG pathways for DE genes
3. Generation of pathway summary tables and barplots

Usage:
    python annotation_pipeline.py \
        --deseq ../results/DESeq2_significant_genes.tsv \
        --outdir ../results/anot_func \
        --organism hsa
"""
import os 
import io 
import time
import argparse
import requests 
import pandas as pd
from Bio.KEGG import REST

# Function to filter DESeq2 results into UP and DOWN regulated genes
def filter_DEGs(input_file, output_dir):
    """
    Filter DESeq2 results to keep only DE genes labeled as "UP" or "DOWN".
    Creates separate files for UP and DOWN genes.

    Args:
        input_file (str): Path to DESeq2 results TSV file.
        output_dir (str): Directory to save filtered files.

    Returns:
        tuple: DataFrames of UP and DOWN genes.
    """
    os.makedirs(output_dir, exist_ok=True)  # Create output directory if it doesn't exist

    df = pd.read_csv(input_file, sep="\t")  # Read DESeq2 results

    # Clean the "Expresion" column and remove whitespace and carriage returns
    df["Expresion"] = df["Expresion"].astype(str).str.strip().str.replace("\r", "")

    # Filter UP and DOWN regulated genes
    df_up = df[df["Expresion"] == "UP"]
    df_down = df[df["Expresion"] == "DOWN"]

    # Prepare output file paths
    up_file = os.path.join(output_dir, "DESeq2_UP_genes.tsv")
    down_file = os.path.join(output_dir, "DESeq2_DOWN_genes.tsv")

    # Save filtered genes
    df_up.to_csv(up_file, sep="\t", index=False)
    df_down.to_csv(down_file, sep="\t", index=False)

    print(f"Filtered UP genes saved: {up_file} ({len(df_up)} genes)")
    print(f"Filtered DOWN genes saved: {down_file} ({len(df_down)} genes)")

    return df_up, df_down

# Function to map Ensembl IDs to Entrez IDs using Ensembl REST API
def ensembl_to_entrez(ensembl_ids, sleep_time=0.05):
    """
    Map Ensembl Gene IDs to Entrez IDs.

    Args:
        ensembl_ids (list): List of Ensembl IDs without version suffix.
        sleep_time (float): Delay between API calls.

    Returns:
        pd.DataFrame: Columns ['Ensembl', 'Entrez', 'Symbol'].
    """
    server = "https://rest.ensembl.org"
    headers = {"Content-Type": "application/json"}
    mapping = []

    for gene in ensembl_ids:
        url = f"{server}/xrefs/id/{gene}?external_db=EntrezGene"  # Ensembl endpoint
        try:
            res = requests.get(url, headers=headers, timeout=10).json()  # Call API
            if res:  # If mapping exists
                mapping.append({
                    "Ensembl": gene,
                    "Entrez": res[0]["primary_id"],
                    "Symbol": res[0].get("display_id", "")
                })
        except Exception as e:
            print(f"Error mapping {gene}: {e}")  # Handle errors
        time.sleep(sleep_time)  # To void overloading the server

    return pd.DataFrame(mapping)

def get_kegg_pathways(entrez_ids, organism="hsa", batch_size=200, sleep_time=0.05):
    """
    Retrieve KEGG pathways for a list of Entrez Gene IDs.

    Args:
        entrez_ids (list): List of Entrez IDs as strings.
        organism (str): KEGG organism code (default 'hsa').
        batch_size (int): Number of genes per batch request.
        sleep_time (float): Delay between batch requests.

    Returns:
        pd.DataFrame: Columns ['Gene', 'Pathway'].
    """
    pathways_all = []

    for i in range(0, len(entrez_ids), batch_size):
        batch = entrez_ids[i:i + batch_size]  # Slice batch
        kegg_ids = [f"{organism}:{g}" for g in batch]  # Format IDs for KEGG

        try:
            result = REST.kegg_link("pathway", "+".join(kegg_ids)).read()  # Get pathways
            if result.strip():  # If result is not empty
                df = pd.read_table(io.StringIO(result), header=None)  # Convert to DF
                df.columns = ["Gene", "Pathway"]
                pathways_all.append(df)
        except Exception as e:
            print(f"KEGG error for batch {i}-{i+batch_size}: {e}")  # For batch errors

        time.sleep(sleep_time)

    if pathways_all:
        return pd.concat(pathways_all, ignore_index=True)  # Merge batches

    return pd.DataFrame(columns=["Gene", "Pathway"])  # Return empty DF if no results

def kegg_id_to_name(pathway_ids, sleep_time=0.05):
    """
    Convert KEGG pathway IDs to names.

    Args:
        pathway_ids (list): KEGG pathway IDs.
        sleep_time (float): Delay between API calls.

    Returns:
        list: Pathway names corresponding to IDs.
    """
    names = []
    for pid in pathway_ids:
        try:
            res = REST.kegg_get(pid).read()  # Get KEGG pathway info
            first_line = res.split("\n")[1]  # Pathway "NAME" line
            names.append(first_line.replace("NAME", "").strip())  # Extract name
        except:
            names.append(pid)  # If not found keep the ID
        time.sleep(sleep_time)
    return names

def main():
    parser = argparse.ArgumentParser(description="Functional annotation pipeline for DE genes.")
    parser.add_argument("--deseq", required=True, help="DESeq2 results file (TSV)")
    parser.add_argument("--outdir", default="../results/anot_func", help="Output directory for annotations and plots")
    parser.add_argument("--organism", default="hsa", help="KEGG organism code (default: hsa)")

    args = parser.parse_args()

    outdir = args.outdir
    figs = os.path.join(outdir, "Figures")
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(figs, exist_ok=True)

    #Filter DEGs
    up_genes, down_genes = filter_DEGs(args.deseq, outdir)

    #Remove version numbers from Ensembl IDs
    up_list = [g.split(".")[0] for g in up_genes["ID"]]
    down_list = [g.split(".")[0] for g in down_genes["ID"]]

    #Map Ensembl IDs to Entrez IDs
    df_up = ensembl_to_entrez(up_list)
    df_down = ensembl_to_entrez(down_list)

    # Save mapping results
    df_up.to_csv(os.path.join(outdir, "UP_genes_Entrez.csv"), index=False)
    df_down.to_csv(os.path.join(outdir, "DOWN_genes_Entrez.csv"), index=False)

    # Retrieve KEGG pathways
    up_pw = get_kegg_pathways(df_up["Entrez"], organism=args.organism)
    down_pw = get_kegg_pathways(df_down["Entrez"], organism=args.organism)

    #Convert KEGG IDs to names
    if not up_pw.empty:
        up_pw["Pathway_Name"] = kegg_id_to_name(up_pw["Pathway"])
    if not down_pw.empty:
        down_pw["Pathway_Name"] = kegg_id_to_name(down_pw["Pathway"])

    up_pw.to_csv(os.path.join(outdir, "UP_genes_pathways.tsv"), sep="\t", index=False)
    down_pw.to_csv(os.path.join(outdir, "DOWN_genes_pathways.tsv"), sep="\t", index=False)
    
    print("\nFunctional annotation pipeline completed.\n")

if __name__ == "__main__":
    main()
