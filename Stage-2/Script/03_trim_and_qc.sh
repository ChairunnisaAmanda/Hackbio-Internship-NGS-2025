#!/usr/bin/bash
# Script: 03_trim_and_qc.sh
# Description: Perform trimming with fastp followed by quality control (FastQC + MultiQC)
#              for single-end (transcriptome) FASTQ files.

# Directories
RAW_DATA_DIR="raw_fastq"
TRIMMED_DATA_DIR="trimmed"
FASTP_REPORT_DIR="results/fastp_reports"
TRIMMED_QC_DIR="results/trimmed_fastqc_reports"

# Create output directories
echo "Creating output directories..."
mkdir -p "$TRIMMED_DATA_DIR" "$FASTP_REPORT_DIR" "$TRIMMED_QC_DIR"

# Check for input files
if [ -z "$(ls -A $RAW_DATA_DIR/*.fastq.gz 2>/dev/null)" ]; then
    echo "Error: No FASTQ files found in $RAW_DATA_DIR!"
    exit 1
fi

echo "STEP 1: TRIMMING WITH FASTP"
echo "Input directory : $RAW_DATA_DIR"
echo "Output directory: $TRIMMED_DATA_DIR"

sample_count=0
for fq in "$RAW_DATA_DIR"/*.fastq.gz; do
    base_name=$(basename "$fq" .fastq.gz)
    sample_count=$((sample_count + 1))
    echo "Processing sample $sample_count: $base_name"

    # Run fastp for single-end data
    fastp \
        -i "$fq" \
        -o "${TRIMMED_DATA_DIR}/${base_name}_trimmed.fastq.gz" \
        --html "${FASTP_REPORT_DIR}/${base_name}_fastp.html" \
        --json "${FASTP_REPORT_DIR}/${base_name}_fastp.json" \
        --thread 4

    echo "Completed: $base_name"
done

echo ""
echo "Trimming completed! Processed $sample_count samples."
echo "Trimmed files saved to: $TRIMMED_DATA_DIR"

# ── Quality control on trimmed reads ─────────────────────────────────────────
if [ -z "$(ls -A $TRIMMED_DATA_DIR/*.fastq.gz 2>/dev/null)" ]; then
    echo "Error: No trimmed files were produced!"
    exit 1
fi

echo ""
echo "STEP 2: QUALITY CONTROL ON TRIMMED DATA"
fastqc "$TRIMMED_DATA_DIR"/*.fastq.gz \
  --outdir "$TRIMMED_QC_DIR" \
  --threads 4

# ── Consolidate reports with MultiQC ─────────────────────────────────────────
echo "Generating MultiQC report..."
multiqc "$TRIMMED_QC_DIR" \
  --outdir "$TRIMMED_QC_DIR" \
  --filename "multiqc_report_trimmed.html"

echo ""
echo "=== SUMMARY ==="
echo "✓ Trimming completed: $sample_count samples processed"
echo "✓ Quality assessment completed on trimmed data"
echo ""
echo "Output directories:"
echo "  - Trimmed FASTQ files : $TRIMMED_DATA_DIR"
echo "  - fastp reports       : $FASTP_REPORT_DIR"
echo "  - FastQC/MultiQC QC   : $TRIMMED_QC_DIR"
echo ""
echo "Script completed successfully!"
