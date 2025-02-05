#!/bin/bash

# Define the input CSV file
input_file="SraAccList.csv"

# Directory to store downloaded data
output_dir="fastq_files"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Read the CSV file and skip the header
tail -n +2 "$input_file" | while IFS= read -r accession; do
    echo "Processing $accession..."

    # Download the SRA file using prefetch
    echo "Downloading $accession..."
    prefetch "$accession" -o "$output_dir/$accession.sra"

    # Check if the download succeeded
    if [ $? -ne 0 ]; then
        echo "Failed to download $accession. Skipping."
        continue
    fi

    # Convert SRA file to FASTQ format using fastq-dump
    echo "Converting $accession to FASTQ format..."
    fastq-dump --split-files --gzip --outdir "$output_dir" "$output_dir/$accession.sra"

    # Check if the conversion succeeded
    if [ $? -ne 0 ]; then
        echo "Failed to convert $accession. Skipping."
        continue
    fi

    # Optionally remove the .sra file to save space
    echo "Cleaning up $accession.sra..."
    rm "$output_dir/$accession.sra"
done

echo "All tasks completed."
