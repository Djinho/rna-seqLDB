#!/bin/bash

# Directories
GENOME_DIR=~/rna-seqLDB/genome_gtf/STAR_index
OUTPUT_DIR=~/rna-seqLDB/datasets/alignment
INPUT_DIR=~/rna-seqLDB/datasets/fastq_files/trimmed

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Iterate over paired files
for SAMPLE in $(ls $INPUT_DIR/*_1_paired.fastq.gz | xargs -n 1 basename | sed 's/_1_paired.fastq.gz//'); do
    BAM_FILE="$OUTPUT_DIR/${SAMPLE}_Aligned.sortedByCoord.out.bam"

    # Skip if the BAM file already exists and is valid
    if [ -f "$BAM_FILE" ] && samtools quickcheck "$BAM_FILE"; then
        echo "Skipping $SAMPLE (already completed)"
    else
        echo "Processing sample: $SAMPLE"
        STAR --runThreadN 4 \
             --genomeDir $GENOME_DIR \
             --readFilesIn $INPUT_DIR/${SAMPLE}_1_paired.fastq.gz $INPUT_DIR/${SAMPLE}_2_paired.fastq.gz \
             --readFilesCommand zcat \
             --outFileNamePrefix $OUTPUT_DIR/${SAMPLE}_ \
             --outSAMtype BAM SortedByCoordinate \
             --sjdbOverhang 75 \
             --genomeLoad NoSharedMemory  # Changed from LoadAndKeep to NoSharedMemory
    fi
done
