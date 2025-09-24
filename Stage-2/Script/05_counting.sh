#!/usr/bin/bash
# Script: 05_counting.sh
# Simple featureCounts script for Arabidopsis (single-end)
# Make sure subread/featureCounts is installed

# Create output directory
mkdir -p counts

# Run featureCounts
featureCounts -O -t gene -g ID -a a_thaliana.gff3 -o counts/counts.txt mapped/*_Aligned.sortedByCoord.out.bam
# Check the summary
cat counts/counts.txt.summary
