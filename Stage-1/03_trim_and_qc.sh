#!/usr/bin/bash
# Script: 03_trim_and_qc.sh
# Description: Perform trimming with fastp followed by quality control on trimmed data.

# Set directories
RAW_DATA_DIR="../data/raw_data_50"
TRIMMED_DATA_DIR="../data/trimmed_data"
FASTP_REPORT_DIR="../data/results/fastp_reports"
TRIMMED_QC_DIR="../data/results/trimmed_fastqc_reports"

# Create output directories
echo "Creating output directories..."
mkdir -p "$TRIMMED_DATA_DIR"
mkdir -p "$FASTP_REPORT_DIR"
mkdir -p "$TRIMMED_QC_DIR"

# Check if raw data exists
if [ -z "$(ls -A $RAW_DATA_DIR/*.fastq.gz 2>/dev/null)" ]; then
    echo "Error: No FASTQ files found in $RAW_DATA_DIR!"
    exit 1
fi

echo "TRIMMING WITH FASTP"
echo "Input: $RAW_DATA_DIR"
echo "Output: $TRIMMED_DATA_DIR"

# Process each sample with fastp
sample_count=0
for r1 in "$RAW_DATA_DIR"/*_1.fastq.gz; do
    base_name=$(basename "$r1" _1.fastq.gz)
    r2="${RAW_DATA_DIR}/${base_name}_2.fastq.gz"
    
    if [ ! -f "$r2" ]; then
        echo "Warning: R2 file not found for $base_name, skipping..."
        continue
    fi
    
    sample_count=$((sample_count + 1))
    echo "Processing sample $sample_count: $base_name"
    
    # Run fastp with basic settings
    fastp \
        -i "$r1" \
        -I "$r2" \
        -o "${TRIMMED_DATA_DIR}/${base_name}_1_trimmed.fastq.gz" \
        -O "${TRIMMED_DATA_DIR}/${base_name}_2_trimmed.fastq.gz" \
        --html "${FASTP_REPORT_DIR}/${base_name}_fastp.html" \
        --json "${FASTP_REPORT_DIR}/${base_name}_fastp.json" \
        --thread 2
    
    echo "Completed: $base_name"
done

echo ""
echo "✓ Trimming completed! Processed $sample_count samples."
echo "Trimmed files saved to: $TRIMMED_DATA_DIR"

# Check if trimming produced any files
if [ -z "$(ls -A $TRIMMED_DATA_DIR/*.fastq.gz 2>/dev/null)" ]; then
    echo "Error: No trimmed files were produced!"
    exit 1
fi

echo ""
echo "=== STEP 2: QUALITY CONTROL ON TRIMMED DATA ==="
echo "Starting FastQC analysis on trimmed data..."
echo "Input: $TRIMMED_DATA_DIR"
echo "Output: $TRIMMED_QC_DIR"

# Run FastQC on all trimmed files
fastqc "$TRIMMED_DATA_DIR"/*.fastq.gz \
  --outdir "$TRIMMED_QC_DIR" \
  --threads 4

echo "FastQC analysis complete."

# Generate consolidated MultiQC report
echo "Generating MultiQC report for trimmed data..."
multiqc "$TRIMMED_QC_DIR" \
  --outdir "$TRIMMED_QC_DIR" \
  --filename "multiqc_report_trimmed.html"

echo ""
echo "=== SUMMARY ==="
echo "✓ Trimming completed: $sample_count samples processed"
echo "✓ Quality assessment completed on trimmed data"
echo ""
echo "Output directories:"
echo "  - Trimmed data: $TRIMMED_DATA_DIR"
echo "  - FASTP reports: $FASTP_REPORT_DIR"
echo "  - Trimmed QC reports: $TRIMMED_QC_DIR"
echo ""
echo "Important files to review:"
echo "  - Individual FASTP reports: $FASTP_REPORT_DIR/*_fastp.html"
echo "  - MultiQC report: $TRIMMED_QC_DIR/multiqc_report_trimmed.html"
echo ""
echo "Compare with raw data QC: ../results/raw_fastqc_reports/multiqc_report_raw.html"
echo ""
echo "Script completed successfully!"
