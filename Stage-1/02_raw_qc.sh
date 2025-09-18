#!/usr/bin/bash
# Script: 02_raw_qc.sh
# Description: Run FastQC on raw sequencing data for initial quality assessment.

# Set directories for clarity
RAW_DATA_DIR="../data/raw_data_50"
QC_OUTPUT_DIR="../results/raw_fastqc_reports"

# Create the output directory if it doesn't exist
echo "Creating output directory: $QC_OUTPUT_DIR"
mkdir -p "$QC_OUTPUT_DIR"

# Check if the raw data directory exists and has files
if [ ! -d "$RAW_DATA_DIR" ]; then
    echo "Error: Raw data directory $RAW_DATA_DIR not found!"
    exit 1
fi

if [ -z "$(ls -A $RAW_DATA_DIR/*.fastq.gz 2>/dev/null)" ]; then
    echo "Error: No FASTQ files found in $RAW_DATA_DIR!"
    exit 1
fi

echo "Starting FastQC analysis on raw data in $RAW_DATA_DIR..."
echo "Output will be saved to: $QC_OUTPUT_DIR"

# Run FastQC on all gzipped FASTQ files in the raw data directory
# Using 4 threads for faster processing
fastqc "$RAW_DATA_DIR"/*.fastq.gz \
  --outdir "$QC_OUTPUT_DIR" \
  --threads 4 \
  --quiet

echo "FastQC analysis complete."
echo "Individual reports are in: $QC_OUTPUT_DIR"

# Generate a consolidated MultiQC report for easy viewing
echo "Generating MultiQC report..."
multiqc "$QC_OUTPUT_DIR" \
  --outdir "$QC_OUTPUT_DIR" \
  --filename "multiqc_report_raw.html" \
  --quiet

echo "MultiQC report generated: $QC_OUTPUT_DIR/multiqc_report_raw.html"
echo "Raw data quality assessment step finished successfully."
