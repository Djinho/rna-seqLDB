#!/bin/bash

# Directory containing FASTQ files
fastq_dir="/home/djinho/rna-seqLDB/datasets/fastq_files"

# Output directory for FastQC results
output_dir="${fastq_dir}/fastqc_results"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Run FastQC on all FASTQ files in the directory
echo "Running FastQC on FASTQ files in $fastq_dir..."
for file in "$fastq_dir"/*.fastq.gz; do
    echo "Processing: $file"
    fastqc -o "$output_dir" "$file"
done

echo "FastQC analysis complete. Results saved to $output_dir."
