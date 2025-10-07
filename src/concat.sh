#!/bin/bash

CSV="$1"
OUTPUT_DIR="$2"

SAMPLE_COL=30
if [[ -z "$CSV" || -z "$OUTPUT_DIR" ]]; then
    echo "Usage: $0 <sra_run_table.csv> <output_directory>"
    exit 1
fi

SAMPLES=$(tail -n +2 "$CSV" | cut -d, -f$SAMPLE_COL)
UNIQUE_SAMPLES=$(echo "$SAMPLES" | sort | uniq)

for SAMPLE in $UNIQUE_SAMPLES; do
    COUNT=$(echo "$SAMPLES" | grep -c "^$SAMPLE$")
    
    if [[ $COUNT -le 1 ]]; then
        echo "Sample ID $SAMPLE appears only once. Nothing to concatenate."
        continue
    fi

    SRRS=$(grep ",$SAMPLE," "$CSV" | cut -d, -f1)

    FILES_R1=()
    FILES_R2=()
    
    for SRR in $SRRS; do
        [[ -f "$OUTPUT_DIR/$SRR/${SRR}_1.fastq.gz" ]] && FILES_R1+=("$OUTPUT_DIR/$SRR/${SRR}_1.fastq.gz")
        [[ -f "$OUTPUT_DIR/$SRR/${SRR}_2.fastq.gz" ]] && FILES_R2+=("$OUTPUT_DIR/$SRR/${SRR}_2.fastq.gz")
    done

    OUT_R1="$OUTPUT_DIR/${SAMPLE}_1.fastq.gz"
    OUT_R2="$OUTPUT_DIR/${SAMPLE}_2.fastq.gz"

    if [[ ${#FILES_R1[@]} -gt 0 ]]; then
        echo "Concatenating read 1 (_1) files for Sample ID $SAMPLE..."
        zcat "${FILES_R1[@]}" | gzip > "$OUT_R1"
    fi

    if [[ ${#FILES_R2[@]} -gt 0 ]]; then
        echo "Concatenating read 2 (_2) files for Sample ID $SAMPLE..."
        zcat "${FILES_R2[@]}" | gzip > "$OUT_R2"
    fi
done

echo "All concatenations completed successfully."
