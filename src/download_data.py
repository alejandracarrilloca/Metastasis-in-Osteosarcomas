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
    parser = argparse.ArgumentParser(description="Download BioProject SRA data and reference genomes")
    parser.add_argument("--project", "-p", nargs="+", required=True, help="BioProject accession")
    parser.add_argument("--output", "-o", required=True, help="Output directory")
    parser.add_argument("--concat-script", "-c", required=True, help="Path to Bash script for concatenating samples")
    parser.add_argument("--organism", "-g", default="homo sapiens", help="Organism for reference genome")
    parser.add_argument("--email", "-e", required=True, help="Email for NCBI Entrez")
    return parser.parse_args()

def verify_bash(bash_path):
    """
    Verify that the Bash script exists
    """
    if not os.path.isfile(bash_path):
        raise FileNotFoundError(f"Bash script not found: {bash_path}")
    print(f"Bash script found: {bash_path}")

def check_sra_tools():
    """Check for required SRA Toolkit tools"""
    for tool in ["prefetch", "fasterq-dump", "gzip"]:
        found = False
        for path_dir in os.environ["PATH"].split(os.pathsep):
            tool_path = os.path.join(path_dir, tool)
            if os.path.isfile(tool_path) and os.access(tool_path, os.X_OK):
                found = True
                break
        if not found:
            raise Exception(f"'{tool}' not found in PATH.")

def bioproject_uid(accession):
    """Get internal UID of a BioProject"""
    handle = Entrez.esearch(database="bioproject", term=accession)
    res = Entrez.read(handle)
    handle.close()
    if not res["IdList"]:
        print(f"No BioProject found for {accession}")
        return None
    return res["IdList"][0]

def bioproject_summary(uid):
    """Retrieve BioProject metadata"""
    handle = Entrez.esummary(database="bioproject", id=uid)
    rec = Entrez.read(handle)
    handle.close()
    return rec["DocumentSummarySet"]["DocumentSummary"][0]

def geo_summary(geo_uid):
    """Retrieve GEO dataset summary information"""
    handle = Entrez.esummary(database="gds", id=geo_uid)
    summary = Entrez.read(handle)
    handle.close()
    
    if not summary:
        print("No GEO summary returned")
        return {}

    geo_record = summary[0]

    print("\n GEO Summary ")
    print("Accession:", geo_record.get("Accession"))
    print("Title:", geo_record.get("title"))
    print("Num Samples:", geo_record.get("n_samples"))

    print("\n Samples ")
    for sample in geo_record.get("Samples", []):
        print(sample.get('Accession'), ":", sample.get('Title'))
    
    return geo_record

def linked_ids(uid, dbtarget):
    """Return IDs linked to a BioProject in another NCBI database"""
    handle = Entrez.elink(dbfrom="bioproject", database=dbtarget, id=uid)
    links = Entrez.read(handle)
    handle.close()
    if not links or not links[0].get("LinkSetDb"):
        return []
    return [l["Id"] for l in links[0]["LinkSetDb"][0]["Link"]]

def sra_runinfo(sra_ids):
    """Get the runinfo table for a list of SRA UIDs"""
    try:
        handle = Entrez.efetch(database="sra", id=",".join(sra_ids), rettype="runinfo", retmode="text")
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
    """Retrieve SRR table for the BioProject."""
    sra_ids = linked_ids(uid, "sra")
    if not sra_ids:
        return None
    df_list = []
    batch = 50
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
    
    for srr in srr_list:
        print(f"\nDownloading {srr}")
        # Prefetch
        subprocess.run(["prefetch", "--output-directory", output_dir, srr], check=True)
        # Fastq-dump
        srr_dir = os.path.join(output_dir, srr)
        os.makedirs(srr_dir, exist_ok=True)
        subprocess.run([
            "fasterq-dump", srr, "-O", srr_dir, "--split-files", "--gzip"
        ], check=True)
        print(f"Finished {srr}. Files saved in {srr_dir}")

        # Validate FASTQ files
        fastq_files = [f for f in os.listdir(srr_dir) if f.endswith(('.fastq', '.fastq.gz'))]
        if not fastq_files:
            print(f"{srr}: No FASTQ files found.")
            continue
        
        #Validate files in directory
        for f in os.listdir(srr_dir):
            path = os.path.join(srr_dir, f)
            if os.path.getsize(path) == 0:
                print(f"{f} is empty")
                
def verify_ref_genome(organism = "Homo sapiens"):
    """
    Verify available reference genomes for Homo sapiens in NCBI Assembly.
    Returns a filtered DataFrame with genomes that have RefSeq references.
    """
    handle = Entrez.esearch(database="assembly", term=organism, retmax=150)
    search_results = Entrez.read(handle)
    handle.close()

    # Print the assembly UIDs found
    assembly_uids = search_results["IdList"]
    print("Assembly UIDs found:", assembly_uids)

    # Retrieve detailed information for each UID
    assembly_data = []
    for uid in assembly_uids:
        handle = Entrez.esummary(database="assembly", id=uid)
        summary = Entrez.read(handle)
        handle.close()
        docsum = summary['DocumentSummarySet']['DocumentSummary'][0]
        assembly_data.append(docsum)

    # Create DataFrame with all assembly info
    df_assembly = pd.DataFrame(assembly_data)
    print("Available columns:", df_assembly.columns.tolist())

    # Select only relevant columns
    df_selected = df_assembly[[
        "Organism",
        "AssemblyName",
        "AssemblyStatus",
        "AssemblyAccession",
        "RefSeq_category",
        "FtpPath_GenBank"
    ]]

    print("First 10 rows of the dataframe:")
    print(df_selected.head(10))

    # Filter only reference genomes
    df_filtered = df_selected[df_selected["RefSeq_category"] != "na"]
    print("Reference genomes:")
    print(df_filtered.head(10))

    return df_filtered

def download_ref_genome(df_filtered, output_dir):
    """
    Download the reference genome (FASTA) and annotation (GFF) files from NCBI FTP.
    Downloads the first RefSeq genome available in df_filtered.
    """

    if df_filtered.empty:
        print("No reference genomes found.")
        return

    # Select the first filtered genome
    ftp_url = df_filtered.iloc[0]["FtpPath_GenBank"]
    genome_name = df_filtered.iloc[0]["AssemblyName"]

    # Create output folder
    out_dir = os.path.join(output_dir, "genome", genome_name)
    os.makedirs(out_dir, exist_ok=True)

    # List available files from FTP
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

    # Identify FASTA and GFF files
    fasta_files = [f for f in files if f.endswith("_genomic.fna.gz")]
    gff_files = [f for f in files if f.endswith("_genomic.gff.gz")]

    if not fasta_files and not gff_files:
        print("No FASTA or GFF files found on FTP.")
        return

    # Download FASTA
    if fasta_files:
        fasta_url = ftp_url + "/" + fasta_files[0]
        fasta_out = os.path.join(out_dir, fasta_files[0])
        try:
            print(f"Downloading FASTA: {fasta_url}")
            urllib.request.urlretrieve(fasta_url, fasta_out)
            print(f"FASTA saved to: {fasta_out}")
        except Exception as e:
            print(f"Failed to download FASTA: {e}")

    # Download GFF
    if gff_files:
        gff_url = ftp_url + "/" + gff_files[0]
        gff_out = os.path.join(out_dir, gff_files[0])
        try:
            print(f"Downloading GFF: {gff_url}")
            urllib.request.urlretrieve(gff_url, gff_out)
            print(f"GFF saved to: {gff_out}")
        except Exception as e:
            print(f"Failed to download GFF: {e}")

def process_bioproject(project_acc, output_dir, organism, concat_script):

    print(f"\nStarting project {project_acc}")
    uid = bioproject_uid(project_acc)
    if not uid:
        return

    summary = bioproject_summary(uid)
    print("\nBioProject Summary")
    print("Accession:", summary.get("Project_Acc"))
    print("Title:", summary.get("Project_Title"))
    print("Organism:", summary.get("Organism_Name"))
    print("Description:", summary.get("Project_Description", "N/A"))

    # Retrieve GEO samples
    geo_ids = linked_ids(uid, "gds")
    if geo_ids:
        geo_summarys = geo_summary(geo_ids[0])
        samples = geo_summarys.get('Samples', [])
        if samples:
            print(f"Retrieved {len(samples)} GEO samples.")
        else:
            print("No GEO samples found in this GEO dataset.")
    else:
        print("No GEO datasets linked to this BioProject.")

    # Associated databases
    for database in ["sra", "biosample", "assembly"]:
        linked = linked_ids(uid, database)
        print(f"Entries in {database}: {len(linked)}")

    # Retrieve SRA table
    df = sra_table(uid)
    if df is None or df.empty:
        print("No SRR data found for this project.")
        return

    print(f"\nFound {len(df)} runs for {df['SampleName'].nunique()} samples.")
    print(df[["Run", "SampleName", "LibraryLayout"]].head(10))

    outdir = os.path.join(output_dir, project_acc)
    os.makedirs(outdir, exist_ok=True)
    run_table_path = os.path.join(outdir, "sra_run_table.csv")
    df.to_csv(run_table_path, index=False)
    print(f"SRA run table saved at {run_table_path}")

    srr_list = df['Run'].tolist()
    print(f"SRR list for project {project_acc}:", srr_list)

    if srr_list:
        download_srr(srr_list, output_dir=outdir)
    
    verify_bash(concat_script)
    subprocess.run([concat_script, run_table_path, outdir], check=True)

    df_filtered = verify_ref_genome(organism=organism)
    download_ref_genome(df_filtered, output_dir=outdir)

    print(f"\nFinished BioProject {project_acc}")

if __name__ == "__main__":
    args = parse_args()

    Entrez.email = args.email

    OUTPUT_DIR = args.output
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for p in args.project:
        process_bioproject(
            p,
            output_dir=OUTPUT_DIR,
            organism=args.organism,
            concat_script=args.concat_script
        )