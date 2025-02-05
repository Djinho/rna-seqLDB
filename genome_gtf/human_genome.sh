#!/bin/bash

# Define URLs for downloading GTF and FASTA files
GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz"
FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

# Download GTF file
echo "Downloading GTF file..."
wget -c $GTF_URL -O Homo_sapiens.GRCh38.113.gtf.gz

# Download FASTA file
echo "Downloading FASTA file..."
wget -c $FASTA_URL -O Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Decompress files
echo "Decompressing files..."
gunzip -f Homo_sapiens.GRCh38.113.gtf.gz
gunzip -f Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

echo "Download and extraction complete!"
