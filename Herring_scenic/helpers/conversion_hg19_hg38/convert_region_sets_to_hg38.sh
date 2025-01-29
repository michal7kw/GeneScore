#!/bin/bash

# Set paths to liftOver binary and chain file
LIFTOVER="/scratch/michal.kubacki/herring/liftOver"
CHAIN_FILE="/scratch/michal.kubacki/herring/hg19ToHg38.over.chain"

cd /scratch/michal.kubacki/herring/output

# Create the output directory
mkdir -p region_sets_hg38

# Iterate over subdirectories in region_sets
for subdir in region_sets/*; do
    # Create corresponding subdirectory in region_sets_hg38
    output_subdir="region_sets_hg38/$(basename "$subdir")"
    mkdir -p "$output_subdir"

    # Iterate over BED files in the current subdirectory
    for bed_file in "$subdir"/*.bed; do
        # Get the filename without the path
        filename=$(basename "$bed_file")

        # Set output file paths
        output_file="$output_subdir/$filename"
        unmapped_file="$output_subdir/${filename%.bed}_unmapped.bed"

        # Run liftOver
        "$LIFTOVER" "$bed_file" "$CHAIN_FILE" "$output_file" "$unmapped_file"
    done
done

## Rcelan up
# find region_sets_hg38 -name "*_unmapped.bed" -exec rm {} \;