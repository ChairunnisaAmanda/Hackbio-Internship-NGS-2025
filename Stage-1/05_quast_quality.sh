#!/bin/bash
# Script: 05_quast_quality.sh
# Description: Run QUAST quality assessment on assembled genomes.

# Set directories
ASSEMBLY_DIR="../results/assembly"
QUAST_DIR="../results/quast_reports"

# Create output directory
echo "Creating output directory: $QUAST_DIR"
mkdir -p "$QUAST_DIR"

# Check if assemblies exist
if [ -z "$(ls -A $ASSEMBLY_DIR/*/contigs.fasta 2>/dev/null)" ]; then
    echo "Error: No assembly files found in $ASSEMBLY_DIR!"
    echo "Please run 04_assembly.sh first."
    exit 1
fi

echo "=== GENOME QUALITY ASSESSMENT WITH QUAST ==="
echo "Input: $ASSEMBLY_DIR"
echo "Output: $QUAST_DIR"

# Process each successful assembly
success_count=0
total_count=0

for assembly_dir in "$ASSEMBLY_DIR"/*; do
    sample_name=$(basename "$assembly_dir")
    contigs_file="$assembly_dir/contigs.fasta"
    
    total_count=$((total_count + 1))
    
    if [ -f "$contigs_file" ] && [ -s "$contigs_file" ]; then
        echo "Running QUAST for sample: $sample_name"
        
        # Run QUAST quality assessment
        quast.py \
            -o "$QUAST_DIR/$sample_name" \
            "$contigs_file" \
            --threads 2 \
            --silent
        
        success_count=$((success_count + 1))
        echo "✓ QUAST completed for $sample_name"
    else
        echo "✗ No contigs file found for $sample_name, skipping QUAST"
    fi
done

echo ""
echo "=== QUAST SUMMARY ==="
echo "Total assemblies checked: $total_count"
echo "Successful QUAST reports: $success_count"
echo "QUAST reports saved to: $QUAST_DIR"
echo ""
echo "Next step: Run 06_amr_toxin_analysis.sh"
