#!/usr/bin/bash
# Script: 06_amr_toxin_analysis.sh
# Description: Run ABRicate for AMR and toxin gene detection.

# Set directories
ASSEMBLY_DIR="../results/assembly"
ABRICATE_DIR="../results/abricate_results"

# Create output directories
echo "Creating output directories..."
mkdir -p "$ABRICATE_DIR/amr"
mkdir -p "$ABRICATE_DIR/toxin"
mkdir -p "$ABRICATE_DIR/summary"

# Check if assemblies exist
if [ -z "$(ls -A $ASSEMBLY_DIR/*/contigs.fasta 2>/dev/null)" ]; then
    echo "Error: No assembly files found in $ASSEMBLY_DIR!"
    echo "Please run 04_assembly.sh first."
    exit 1
fi

echo "=== AMR AND TOXIN GENE DETECTION WITH ABRICATE ==="
echo "Input: $ASSEMBLY_DIR"
echo "Output: $ABRICATE_DIR"

# Process each successful assembly
success_count=0
total_count=0

for assembly_dir in "$ASSEMBLY_DIR"/*; do
    sample_name=$(basename "$assembly_dir")
    contigs_file="$assembly_dir/contigs.fasta"
    
    total_count=$((total_count + 1))
    
    if [ -f "$contigs_file" ] && [ -s "$contigs_file" ]; then
        echo "Processing sample: $sample_name"
        
        # Run ABRicate for AMR genes (CARD database)
        echo "  Detecting AMR genes..."
        abricate --db card \
            --quiet \
            "$contigs_file" > "$ABRICATE_DIR/amr/${sample_name}_amr.tsv"
        
        # Run ABRicate for toxin/virulence genes (VFDB database)
        echo "  Detecting toxin genes..."
        abricate --db vfdb \
            --quiet \
            "$contigs_file" > "$ABRICATE_DIR/toxin/${sample_name}_toxin.tsv"
        
        success_count=$((success_count + 1))
        echo "✓ ABRicate completed for $sample_name"
    else
        echo "✗ No contigs file found for $sample_name, skipping ABRicate"
    fi
done

# Generate summary reports
echo ""
echo "Generating summary reports..."

# Summarize AMR results
abricate --summary "$ABRICATE_DIR/amr"/*.tsv > "$ABRICATE_DIR/summary/amr_summary.csv"
# Summarize toxin results
abricate --summary "$ABRICATE_DIR/toxin"/*.tsv > "$ABRICATE_DIR/summary/toxin_summary.csv"

# Combine all results into single files (for easier analysis)
cat "$ABRICATE_DIR/amr"/*.tsv > "$ABRICATE_DIR/summary/all_amr_results.tsv"
cat "$ABRICATE_DIR/toxin"/*.tsv > "$ABRICATE_DIR/summary/all_toxin_results.tsv"

echo ""
echo "=== ABRICATE SUMMARY ==="
echo "Total assemblies checked: $total_count"
echo "Successful ABRicate analyses: $success_count"
echo "Results saved to: $ABRICATE_DIR"
echo ""
echo "Summary files created:"
echo "  - AMR summary: $ABRICATE_DIR/summary/amr_summary.csv"
echo "  - Toxin summary: $ABRICATE_DIR/summary/toxin_summary.csv"
echo "  - Combined AMR results: $ABRICATE_DIR/summary/all_amr_results.tsv"
echo "  - Combined toxin results: $ABRICATE_DIR/summary/all_toxin_results.tsv"
echo ""
echo "Next step: Analyze the results and generate final report"
