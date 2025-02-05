#!/bin/bash

# Define paths
GENOME_DIR="/home/djinho/rna-seqLDB/genome_gtf/STAR_index"
GENOME_FASTA="/home/djinho/rna-seqLDB/genome_gtf/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_FILE="/home/djinho/rna-seqLDB/genome_gtf/Homo_sapiens.GRCh38.113.gtf"
THREADS=8  # Number of threads to use
SJDB_OVERHANG=75  # Read length minus 1 (adjust if needed)

# Create the STAR index directory if it doesn't exist
mkdir -p "$GENOME_DIR"

# Generate the STAR genome index
echo "Generating STAR genome index..."
STAR --runThreadN $THREADS \
     --runMode genomeGenerate \
     --genomeDir "$GENOME_DIR" \
     --genomeFastaFiles "$GENOME_FASTA" \
     --sjdbGTFfile "$GTF_FILE" \
     --sjdbOverhang $SJDB_OVERHANG

# Check if the index was generated successfully
if [ $? -eq 0 ]; then
    echo "STAR genome index successfully generated in $GENOME_DIR."
else
    echo "Error: Failed to generate STAR genome index. Please check the inputs and paths."
fi
