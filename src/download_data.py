"""
07/10/2025
Download data v.3
Dulce Alejandra Carrillo Carlos 
Andrea Villarruel Garcia 
Antonio Platas Renteral 

BioProject Downloader and Reference Genome Retriever

This script automates the retrieval of sequencing data (SRA) and reference genomes (NCBI Assembly)
for one or more BioProjects. It:
1. Queries BioProject metadata through NCBI Entrez.
2. Retrieves associated SRA and GEO records.
3. Downloads FASTQ data using SRA Toolkit (prefetch + fasterq-dump + gzip).
4. Runs a Bash script to concatenate FASTQ files by sample.
5. Downloads the reference genome (FASTA + GFF) for the specified organism.
"""

import argparse
from Bio import Entrez
import pandas as pd
from io import StringIO
import subprocess
import os
import sys
import urllib.request
from ftplib import FTP
from urllib.parse import urlparse


def parse_args():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments with the following fields:
            - project: BioProject accession(s)
            - output: Output directory for all results
            - concat_script: Path to Bash script for concatenating samples
            - organism: Organism for reference genome (default: "homo sapiens")
            - email: Email required for NCBI Entrez access
    """
    parser = argparse.ArgumentParser(description="Download BioProject SRA data and reference genomes")
    parser.add_argument("--project", "-p", nargs="+", required=True, help="BioProject accession(s)")
    parser.add_argument("--output", "-o", required=True, help="Output directory")
    parser.add_argument("--concat-script", "-c", required=True, help="Path to Bash script for concatenating samples")
    parser.add_argument("--organism", "-g", default="homo sapiens", help="Organism for reference genome")
    parser.add_argument("--email", "-e", required=True, help="Email for NCBI Entrez")
    return parser.parse_args()


def verify_bash(bash_path):
    """
    Verify that the provided Bash script exists and is accessible.

    Args:
        bash_path (str): Path to the Bash script.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    if not os.path.isfile(bash_path):
        raise FileNotFoundError(f"Bash script not found: {bash_path}")
    print(f"Bash script found: {bash_path}")


def check_sra_tools():
    """
    Check for the presence of essential SRA Toolkit tools in the system PATH.

    Tools verified:
        - prefetch
        - fasterq-dump
        - gzip

    Raises:
        Exception: If any tool is not found or not executable.
    """
    for tool in ["prefetch", "fasterq-dump", "gzip"]:
        found = False
        # Iterate through each directory in PATH
        for path_dir in os.environ["PATH"].split(os.pathsep):
            tool_path = os.path.join(path_dir, tool)
            # Check if the file exists and is executable
            if os.path.isfile(tool_path) and os.access(tool_path, os.X_OK):
                found = True
                break
        if not found:
            raise Exception(f"'{tool}' not found in PATH.")
    print("All SRA Toolkit tools found.")


def bioproject_uid(accession):
    """
    Retrieve the internal NCBI UID of a BioProject.

    Args:
        accession (str): BioProject accession.

    Returns:
        str or None: The BioProject UID, or None if not found.
    """
    handle = Entrez.esearch(db="bioproject", term=accession)
    res = Entrez.read(handle)
    handle.close()
    if not res["IdList"]:
        print(f"No BioProject found for {accession}")
        return None
    return res["IdList"][0]


def bioproject_summary(uid):
    """
    Retrieve metadata for a given BioProject UID.

    Args:
        uid (str): NCBI BioProject UID.

    Returns:
        dict: Summary fields including title, organism, and description.
    """
    handle = Entrez.esummary(db="bioproject", id=uid)
    rec = Entrez.read(handle)
    handle.close()
    return rec["DocumentSummarySet"]["DocumentSummary"][0]


def linked_ids(uid, dbtarget):
    """
    Retrieve all linked IDs between a BioProject and another NCBI database.

    Args:
        uid (str): BioProject UID.
        dbtarget (str): Target NCBI database (e.g., 'sra', 'biosample', 'assembly').

    Returns:
        list: List of linked IDs.
    """
    handle = Entrez.elink(dbfrom="bioproject", db=dbtarget, id=uid)
    links = Entrez.read(handle)
    handle.close()
    if not links or not links[0].get("LinkSetDb"): #Verify if empty
        return []
    return [l["Id"] for l in links[0]["LinkSetDb"][0]["Link"]] #Extract from dicts into a list of Ids


def geo_summary(geo_uid):
    """
    Retrieve summary information from a GEO dataset linked to a BioProject.

    Args:
        geo_uid (str): GEO dataset UID.

    Returns:
        dict: GEO record with accession, title, number of samples, and sample list.
    """
    handle = Entrez.esummary(db="gds", id=geo_uid)
    summary = Entrez.read(handle)
    handle.close()

    if not summary:
        print("No GEO summary returned")
        return {}

    geo_record = summary[0]
    print("\nGEO Summary")
    print("Accession:", geo_record.get("Accession"))
    print("Title:", geo_record.get("title"))
    print("Num Samples:", geo_record.get("n_samples"))

    print("\nSamples:")
    for sample in geo_record.get("Samples", []):
        print(sample.get('Accession'), ":", sample.get('Title'))

    return geo_record


def sra_runinfo(sra_ids):
    """
    Retrieve the run information table for a list of SRA UIDs.

    Args:
        sra_ids (list): List of SRA UIDs.

    Returns:
        pandas.DataFrame: Table containing run information.
    """
    try:
        handle = Entrez.efetch(db="sra", id=",".join(sra_ids), rettype="runinfo", retmode="text")
        data = handle.read()
        handle.close()

        if isinstance(data, bytes):
            data = data.decode()

        if not data.strip():
            return pd.DataFrame()
        return pd.read_csv(StringIO(data))
    
    except Exception as e:
        print(f"Could not get info for {sra_ids}: {e}")
        return pd.DataFrame()


def sra_table(uid):
    """
    Retrieve and combine SRA run tables for all runs associated with a BioProject.

    Args:
        uid (str): BioProject UID.

    Returns:
        pandas.DataFrame: Combined SRA run table or None if no runs are found.
    """
    sra_ids = linked_ids(uid, "sra")
    if not sra_ids:
        return None

    df_list = []
    batch = 50
    # Retrieve data in batches of 50 IDs
    for i in range(0, len(sra_ids), batch):
        sub = sra_ids[i:i + batch]
        print(f"Getting SRA runinfo batch {i + 1}-{i + len(sub)}...")
        df_list.append(sra_runinfo(sub))

    df_list = [df_chunk for df_chunk in df_list if not df_chunk.empty]
    if not df_list:
        print("No valid SRA runinfo tables were found.")
        return None

    df = pd.concat(df_list, ignore_index=True)
    return df.drop_duplicates(subset="Run")



def download_srr(srr_list, output_dir):
    """
    Download SRA runs, convert them to FASTQ, and validate output files.

    Args:
        srr_list (list): List of SRR accessions.
        output_dir (str): Output directory to store downloaded data.
    """
    for srr in srr_list:
        print(f"\nDownloading {srr}")

        # Step 1: Download the .sra file
        subprocess.run(["prefetch", "--output-directory", output_dir, srr], check=True)

        # Step 2: Convert to FASTQ
        srr_dir = os.path.join(output_dir, srr)
        os.makedirs(srr_dir, exist_ok=True)
        subprocess.run(["fasterq-dump", srr, "-O", srr_dir, "--split-files", "--gzip"], check=True)

        print(f"Finished {srr}. Files saved in {srr_dir}")

        # Step 3: Validate FASTQ files
        fastq_files = [f for f in os.listdir(srr_dir) if f.endswith(('.fastq', '.fastq.gz'))]
        if not fastq_files:
            print(f"{srr}: No FASTQ files found.")
            continue

        # Check file sizes (no empty files)
        for f in fastq_files:
            path = os.path.join(srr_dir, f)
            if os.path.getsize(path) == 0:
                print(f"Warning: {f} is empty.")


def verify_ref_genome(organism="Homo sapiens"):
    """
    Retrieve available reference genomes for a given organism from NCBI Assembly.

    Args:
        organism (str): Organism name. Default is "Homo sapiens".

    Returns:
        pandas.DataFrame: Filtered table with RefSeq reference genomes.
    """
    handle = Entrez.esearch(db="assembly", term=organism, retmax=150)
    search_results = Entrez.read(handle)
    handle.close()

    assembly_uids = search_results["IdList"]
    assembly_data = []

    for uid in assembly_uids:
        handle = Entrez.esummary(db="assembly", id=uid)
        summary = Entrez.read(handle)
        handle.close()
        docsum = summary['DocumentSummarySet']['DocumentSummary'][0]
        assembly_data.append(docsum)

    df_assembly = pd.DataFrame(assembly_data)
    df_selected = df_assembly[[
        "Organism", "AssemblyName", "AssemblyStatus", "AssemblyAccession", "RefSeq_category", "FtpPath_GenBank"
    ]]

    # Keep only RefSeq reference genomes
    df_filtered = df_selected[df_selected["RefSeq_category"] != "na"]
    return df_filtered


def download_ref_genome(df_filtered, output_dir):
    """
    Download the reference genome (FASTA and GFF files) from NCBI FTP.

    Args:
        df_filtered (pandas.DataFrame): Output from verify_ref_genome().
        output_dir (str): Output directory for saving genome files.
    """
    if df_filtered.empty:
        print("No reference genomes found.")
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

    # Download files
    for file_list, label in [(fasta_files, "FASTA"), (gff_files, "GFF")]:
        if not file_list:
            print(f"No {label} files found.")
            continue
        file_url = ftp_url + "/" + file_list[0]
        out_file = os.path.join(out_dir, file_list[0])
        print(f"Downloading {label}: {file_url}")
        try:
            urllib.request.urlretrieve(file_url, out_file)
            print(f"{label} saved to: {out_file}")
        except Exception as e:
            print(f"Failed to download {label}: {e}")


def process_bioproject(project_acc, output_dir, organism, concat_script):
    """
    Main function that orchestrates the download and processing for one BioProject.

    Args:
        project_acc (str): BioProject accession.
        output_dir (str): Output directory.
        organism (str): Organism for genome retrieval.
        concat_script (str): Path to Bash script for concatenation.
    """
    print(f"\nStarting project {project_acc}")

    uid = bioproject_uid(project_acc)
    if not uid:
        return

    # Get and print BioProject summary
    summary = bioproject_summary(uid)
    print("\nBioProject Summary")
    print("Accession:", summary.get("Project_Acc"))
    print("Title:", summary.get("Project_Title"))
    print("Organism:", summary.get("Organism_Name"))
    print("Description:", summary.get("Project_Description", "N/A"))

    # Retrieve GEO dataset if available
    geo_ids = linked_ids(uid, "gds")
    if geo_ids:
        geo_data = geo_summary(geo_ids[0])
        samples = geo_data.get("Samples", [])
        print(f"Retrieved {len(samples)} GEO samples.")
    else:
        print("No GEO datasets linked to this BioProject.")

    # Show linked database entries
    for db in ["sra", "biosample", "assembly"]:
        linked = linked_ids(uid, db)
        print(f"Entries in {db}: {len(linked)}")

    # Get SRA run table
    df = sra_table(uid)
    if df is None or df.empty:
        print("No SRR data found for this project.")
        return

    print(f"\nFound {len(df)} runs for {df['SampleName'].nunique()} samples.")
    print(df[["Run", "SampleName", "LibraryLayout"]].head(10))

    # Save run table
    outdir = os.path.join(output_dir, project_acc)
    os.makedirs(outdir, exist_ok=True)
    run_table_path = os.path.join(outdir, "sra_run_table.csv")
    df.to_csv(run_table_path, index=False)
    print(f"SRA run table saved at {run_table_path}")

    # Download sequencing data
    srr_list = df['Run'].tolist()
    download_srr(srr_list, output_dir=outdir)

    # Run concatenation script
    verify_bash(concat_script)
    subprocess.run([concat_script, run_table_path, outdir], check=True)

    # Download reference genome
    df_filtered = verify_ref_genome(organism=organism)
    download_ref_genome(df_filtered, output_dir=outdir)

    print(f"\nFinished BioProject {project_acc}")


if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_args()

    # Set NCBI Entrez email
    Entrez.email = args.email

    # Create output directory if not existing
    OUTPUT_DIR = args.output
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Process each BioProject accession provided
    for project in args.project:
        process_bioproject(
            project_acc=project,
            output_dir=OUTPUT_DIR,
            organism=args.organism,
            concat_script=args.concat_script
        )
