#!/bin/bash

# ==========================================================
# Script: Concatenate FASTQ files by sample
# For each sample in the SRA run table CSV, 
# concatenate all corresponding FASTQ files (_1 and _2)
# into a single file per read direction.
# Input:
#   $1 - Path to sra_run_table.csv
#   $2 - Output directory containing SRR subdirectories
# ==========================================================

CSV="$1"
OUTPUT_DIR="$2"

# Column number in CSV that contains the sample name
SAMPLE_COL=30

# -----------------------------
# Check for required arguments
# -----------------------------
if [[ -z "$CSV" || -z "$OUTPUT_DIR" ]]; then
    echo "Usage: $0 <sra_run_table.csv> <output_directory>"
    exit 1
fi

# -----------------------------
# Extract unique sample names
# -----------------------------
# Skip header with tail, extract sample column, sort and get unique names
SAMPLES=$(tail -n +2 "$CSV" | cut -d, -f$SAMPLE_COL)
UNIQUE_SAMPLES=$(echo "$SAMPLES" | sort | uniq)

# -----------------------------
# Loop through each sample
# -----------------------------
for SAMPLE in $UNIQUE_SAMPLES; do
    # Count how many times the sample appears
    COUNT=$(echo "$SAMPLES" | grep -c "^$SAMPLE$")
    
    # If sample appears only once, nothing to concatenate
    if [[ $COUNT -le 1 ]]; then
        echo "Sample ID $SAMPLE appears only once. Nothing to concatenate."
        continue
    fi

    # Get all SRR IDs associated with this sample from the CSV
    SRRS=$(grep ",$SAMPLE," "$CSV" | cut -d, -f1)

    # Arrays to hold FASTQ file paths
    FILES_R1=()
    FILES_R2=()
    
    # -----------------------------
    # Collect existing FASTQ files
    # -----------------------------
    for SRR in $SRRS; do
        [[ -f "$OUTPUT_DIR/$SRR/${SRR}_1.fastq.gz" ]] && FILES_R1+=("$OUTPUT_DIR/$SRR/${SRR}_1.fastq.gz")
        [[ -f "$OUTPUT_DIR/$SRR/${SRR}_2.fastq.gz" ]] && FILES_R2+=("$OUTPUT_DIR/$SRR/${SRR}_2.fastq.gz")
    done

    # Define output file paths for concatenated FASTQ
    OUT_R1="$OUTPUT_DIR/${SAMPLE}_1.fastq.gz"
    OUT_R2="$OUTPUT_DIR/${SAMPLE}_2.fastq.gz"

    # -----------------------------
    # Concatenate read 1 files (_1)
    # -----------------------------
    if [[ ${#FILES_R1[@]} -gt 0 ]]; then
        echo "Concatenating read 1 (_1) files for Sample ID $SAMPLE..."
        zcat "${FILES_R1[@]}" | gzip > "$OUT_R1"
    fi

    # -----------------------------
    # Concatenate read 2 files (_2)
    # -----------------------------
    if [[ ${#FILES_R2[@]} -gt 0 ]]; then
        echo "Concatenating read 2 (_2) files for Sample ID $SAMPLE..."
        zcat "${FILES_R2[@]}" | gzip > "$OUT_R2"
    fi
done

echo "All concatenations completed successfully."
