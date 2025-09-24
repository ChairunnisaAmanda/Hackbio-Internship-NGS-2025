#!/bin/bash
# Script: 04_mapping_star.sh
# Script for mapping Arabidopsis single-end reads with STAR
# and indexing BAM files for IGV

# Directories
GENOME_DIR="genome/genomeIndex"
TRIM_DIR="trimmed"
OUT_DIR="mapped"

mkdir -p $OUT_DIR

# 1. Run STAR alignment for each trimmed FASTQ
for infile in $TRIM_DIR/*.fq.gz ; do
    [ -e "$infile" ] || continue

    sample=$(basename "$infile" .fq.gz)
    echo ">>> Aligning sample: $sample"

    STAR --genomeDir $GENOME_DIR \
         --readFilesIn $infile \
         --readFilesCommand zcat \
         --outFileNamePrefix $OUT_DIR/${sample}_ \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes All
done

# 2. Index all BAMs for IGV
echo ">>> Indexing BAM files..."
for bam in $OUT_DIR/*Aligned.sortedByCoord.out.bam; do
    samtools index $bam
done

echo ">>> All alignments complete and indexed!"
