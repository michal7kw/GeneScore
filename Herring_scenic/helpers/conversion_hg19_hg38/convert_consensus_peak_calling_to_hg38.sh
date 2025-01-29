#!/bin/bash

# Set paths to liftOver binary and chain file
LIFTOVER="../liftOver"
CHAIN_FILE="../hg19ToHg38.over.chain"

cd /scratch/michal.kubacki/herring/output

# Create the output directory
mkdir -p consensus_peak_calling_hg38

# Iterate over files and directories in consensus_peak_calling
for item in consensus_peak_calling/*; do
    # Check if the item is a file
    if [ -f "$item" ]; then
        # Get the filename without the path
        filename=$(basename "$item")

        # Set output file paths
        output_file="consensus_peak_calling_hg38/$filename"
        unmapped_file="consensus_peak_calling_hg38/${filename%.*}_unmapped.bed"

        # Run liftOver
        "$LIFTOVER" "$item" "$CHAIN_FILE" "$output_file" "$unmapped_file"
    fi

    # Check if the item is a directory
    if [ -d "$item" ]; then
        # Create corresponding directory in consensus_peak_calling_hg38
        output_dir="consensus_peak_calling_hg38/$(basename "$item")"
        mkdir -p "$output_dir"

        # Iterate over files in the current directory
        for file in "$item"/*; do
            # Get the filename without the path
            filename=$(basename "$file")

            # Set output file paths
            output_file="$output_dir/$filename"
            unmapped_file="$output_dir/${filename%.*}_unmapped.bed"

            # Run liftOver
            "$LIFTOVER" "$file" "$CHAIN_FILE" "$output_file" "$unmapped_file"
        done
    fi
done

## clean up
# find consensus_peak_calling_hg38 -name "*_unmapped.bed" -exec rm {} \;