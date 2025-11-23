from Bio import Entrez
import pandas as pd
import os
from ftplib import FTP
from urllib.parse import urlparse
import urllib.request
import gzip

def verify_ref_genome(organism="Homo sapiens"):
    """
    Search and verify available reference or representative genomes
    for a given organism in NCBI Assembly.
    Returns a DataFrame with valid assemblies.
    """
    print(f"\nSearching reference genome for: {organism}")

    handle = Entrez.esearch(db="assembly", term=f"{organism}[Organism]", retmax=50)
    search_results = Entrez.read(handle)
    handle.close()

    if not search_results["IdList"]:
        print("No assemblies found for this organism.")
        return pd.DataFrame()

    assembly_data = []
    for uid in search_results["IdList"]:
        try:
            handle = Entrez.esummary(db="assembly", id=uid)
            summary = Entrez.read(handle)
            handle.close()
            docsum = summary['DocumentSummarySet']['DocumentSummary'][0]
            assembly_data.append(docsum)
        except Exception as e:
            print(f"Error retrieving assembly UID {uid}: {e}")
            continue

    if not assembly_data:
        print("No assembly summaries could be retrieved.")
        return pd.DataFrame()

    df = pd.DataFrame(assembly_data)

    if "RefSeq_category" not in df.columns:
        print("RefSeq_category not found in assembly summary data.")
        return pd.DataFrame()

    # Filter for reference or representative genomes
    df_filtered = df[df["RefSeq_category"].isin(["reference genome", "representative genome"])]

    if df_filtered.empty:
        print("No reference or representative genome found, returning first available assembly instead.")
        df_filtered = df.head(1)

    print("\nAvailable genomes:")
    print(df_filtered[["Organism", "AssemblyName", "RefSeq_category", "AssemblyAccession"]].head(5))

    return df_filtered

def download_ref_genome(df_filtered, output_dir):
    """
    Download the reference genome (FASTA) and annotation (GFF) files from NCBI FTP.
    Downloads the first RefSeq genome available in df_filtered.
    """

    if df_filtered.empty:
        print("No reference genomes available to download.")
        return

    ftp_url = df_filtered.iloc[0]["FtpPath_GenBank"]
    genome_name = df_filtered.iloc[0]["AssemblyName"]

    out_dir = os.path.join(output_dir, "genome", genome_name)
    os.makedirs(out_dir, exist_ok=True)

    parsed = urlparse(ftp_url)
    ftp_server = parsed.hostname
    ftp_path = parsed.path

    try:
        ftp = FTP(ftp_server)
        ftp.login()
        ftp.cwd(ftp_path)
        files = ftp.nlst()
        ftp.quit()
    except Exception as e:
        print(f"FTP connection failed: {e}")
        return

    fasta_files = [f for f in files if f.endswith("_genomic.fna.gz")]
    gff_files = [f for f in files if f.endswith("_genomic.gff.gz")]

    if not fasta_files and not gff_files:
        print("No FASTA or GFF files found on FTP.")
        return

    if fasta_files:
        fasta_url = ftp_url + "/" + fasta_files[0]
        fasta_out = os.path.join(out_dir, fasta_files[0])
        try:
            print(f"Downloading FASTA: {fasta_url}")
            urllib.request.urlretrieve(fasta_url, fasta_out)
            print(f"FASTA saved to: {fasta_out}")
        except Exception as e:
            print(f"Failed to download FASTA: {e}")

    if gff_files:
        gff_url = ftp_url + "/" + gff_files[0]
        gff_out = os.path.join(out_dir, gff_files[0])
        try:
            print(f"Downloading GFF: {gff_url}")
            urllib.request.urlretrieve(gff_url, gff_out)
            print(f"GFF saved to: {gff_out}")
        except Exception as e:
            print(f"Failed to download GFF: {e}")


if __name__ == "__main__":
    OUTPUT_DIR = "./test_ref_download"
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df_filtered = verify_ref_genome("Homo sapiens")
    download_ref_genome(df_filtered, OUTPUT_DIR)
