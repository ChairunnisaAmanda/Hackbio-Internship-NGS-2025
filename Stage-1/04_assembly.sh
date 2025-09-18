#!/usr/bin/bash
# Script: 04_assembly.sh
# Description: Assemble trimmed reads into genomes using SPAdes.

# Set directories
TRIMMED_DATA_DIR="../data/trimmed_data"
ASSEMBLY_DIR="../results/assembly"

# Create output directory
echo "Creating output directory: $ASSEMBLY_DIR"
mkdir -p "$ASSEMBLY_DIR"

# Check if trimmed data exists
if [ -z "$(ls -A $TRIMMED_DATA_DIR/*.fastq.gz 2>/dev/null)" ]; then
    echo "Error: No trimmed FASTQ files found in $TRIMMED_DATA_DIR!"
    echo "Please run 03_trim_and_qc.sh first."
    exit 1
fi

echo "=== GENOME ASSEMBLY WITH SPADES ==="
echo "Input: $TRIMMED_DATA_DIR"
echo "Output: $ASSEMBLY_DIR"

# Process each sample
sample_count=0
success_count=0
for r1 in "$TRIMMED_DATA_DIR"/*_1_trimmed.fastq.gz; do
    base_name=$(basename "$r1" _1_trimmed.fastq.gz)
    r2="${TRIMMED_DATA_DIR}/${base_name}_2_trimmed.fastq.gz"
    
    sample_count=$((sample_count + 1))
    echo ""
    echo "Processing sample $sample_count: $base_name"
    
    # Create output directory for this sample
    sample_outdir="${ASSEMBLY_DIR}/${base_name}"
    mkdir -p "$sample_outdir"
    
    # Run SPAdes assembly
    echo "Running SPAdes assembly..."
    spades.py \
        -1 "$r1" \
        -2 "$r2" \
        -o "$sample_outdir" \
        --careful \  tries to reduce number of mismatches and short indels
        -t 4 \
        --memory 32 \
        --isolate --phred-offset 33
                # optimized for single isolate assembly
    
    # Check if assembly was successful
    if [ -f "${sample_outdir}/contigs.fasta" ] && [ -s "${sample_outdir}/contigs.fasta" ]; then
        echo "✓ Assembly successful: ${sample_outdir}/contigs.fasta"
        success_count=$((success_count + 1))
        
        # Check assembly statistics
        echo "Assembly statistics for $base_name:"
        echo "Number of contigs: $(grep -c '^>' ${sample_outdir}/contigs.fasta)"
        echo "Total length: $(awk '/^>/ {if (seqlen) print seqlen; seqlen=0; next} {seqlen += length($0)} END {print seqlen}' ${sample_outdir}/contigs.fasta) bp"
    else
        echo "✗ Assembly failed for: $base_name"
    fi
done

echo ""
echo "=== ASSEMBLY SUMMARY ==="
echo "Total samples processed: $sample_count"
echo "Successful assemblies: $success_count"
echo "Assembly results saved to: $ASSEMBLY_DIR"
echo ""
echo "Next step: Run 05_blast_id.sh for species identification"
