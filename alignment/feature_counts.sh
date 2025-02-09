#!/bin/bash

# Directories
GTF_FILE=~/rna-seqLDB/genome_gtf/Homo_sapiens.GRCh38.113.gtf
BAM_DIR=~/rna-seqLDB/datasets/alignment
OUTPUT_FILE=~/rna-seqLDB/datasets/alignment/gene_counts_matrix.txt

# Ensure proper permissions
chmod +rw $BAM_DIR/*.bam
chmod +rw $BAM_DIR

# Run featureCounts with paired-end option
featureCounts -T 4 \
  -a $GTF_FILE \
  -o $OUTPUT_FILE \
  -p --countReadPairs \
  $BAM_DIR/*.bam

# Check if featureCounts ran successfully
if [ $? -eq 0 ]; then
    echo "featureCounts completed successfully. Gene counts saved to $OUTPUT_FILE."
else
    echo "Error: featureCounts failed. Please check the input files and paths."
fi
