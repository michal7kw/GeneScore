#!/bin/bash

# Set paths to liftOver binary and chain file
LIFTOVER="/scratch/michal.kubacki/herring/liftOver"
CHAIN_FILE="/scratch/michal.kubacki/herring/hg19ToHg38.over.chain"

cd /scratch/michal.kubacki/herring/Herring_data

# Create the output directory
mkdir -p GSE168408_RAW_hg38

# Function to process each snATAC file
process_file() {
    file="$1"
    filename=$(basename "$file")
    output_file="GSE168408_RAW_hg38/${filename}"
    unmapped_file="GSE168408_RAW_hg38/${filename%.tsv.gz}.unmapped.bed"
    temp_file="GSE168408_RAW_hg38/${filename%.tsv.gz}.temp"

    echo "Processing file: $file"

    # Uncompress the input file and extract the first 6 columns
    echo "Uncompressing and extracting columns..."
    zcat "$file" | awk -F '\t' 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6}' > "${temp_file}"

    # Perform the liftOver
    echo "Performing liftOver..."
    "${LIFTOVER}" "${temp_file}" "${CHAIN_FILE}" "${output_file}.tmp" "${unmapped_file}"

    # Compress the unmapped file if it exists
    if [ -s "${unmapped_file}" ]; then
        echo "Compressing unmapped file..."
        gzip -c "${unmapped_file}" > "${unmapped_file}.gz"
        rm "${unmapped_file}"
    fi

    # Sort and compress the output file if it exists
    if [ -s "${output_file}.tmp" ]; then
        echo "Sorting and compressing output file..."
        sort -k1,1 -k2,2n "${output_file}.tmp" | gzip -c > "${output_file}"
        rm "${output_file}.tmp"
    fi

    # Remove the temporary file
    echo "Removing temporary file..."
    rm "${temp_file}"

    echo "Finished processing file: $file"
}

export -f process_file
export LIFTOVER
export CHAIN_FILE

files=(
    "GSE168408_RAW/GSM5138529_RL1784_2y_snATAC_fragments.tsv.gz"
    "GSE168408_RAW/GSM5138532_RL2210_4y_snATAC_fragments.tsv.gz"
    "GSE168408_RAW/GSM5138510_RL2366_ga22_snATAC_fragments.tsv.gz"  
    "GSE168408_RAW/GSM5138512_RL2207_ga24_snATAC_fragments.tsv.gz"  
    "GSE168408_RAW/GSM5138515_RL2367_1m_snATAC_fragments.tsv.gz"    
    "GSE168408_RAW/GSM5138518_RL1914_3m_snATAC_fragments.tsv.gz"    
    "GSE168408_RAW/GSM5138521_RL2208_6m_snATAC_fragments.tsv.gz"    
    "GSE168408_RAW/GSM5138523_RL2371_10m_snATAC_fragments.tsv.gz"   
    "GSE168408_RAW/GSM5138526_RL2209_1y_snATAC_fragments.tsv.gz"    
    "GSE168408_RAW/GSM5138529_RL1784_2y_snATAC_fragments.tsv.gz"    
    "GSE168408_RAW/GSM5138532_RL2210_4y_snATAC_fragments.tsv.gz"
    "GSE168408_RAW/GSM5138534_RL2364_6y_snATAC_fragments.tsv.gz"
    "GSE168408_RAW/GSM5138536_RL1994_8y_snATAC_fragments.tsv.gz"
    "GSE168408_RAW/GSM5138539_RL2368_10y_snATAC_fragments.tsv.gz"
    "GSE168408_RAW/GSM5138542_RL2372_14y_snATAC_fragments.tsv.gz"
    "GSE168408_RAW/GSM5138544_RL1785_16y_snATAC_fragments.tsv.gz"
    "GSE168408_RAW/GSM5138548_RL2085_20y_snATAC_fragments.tsv.gz"
    "GSE168408_RAW/GSM5138550_RL2369_25y_snATAC_fragments.tsv.gz"
    "GSE168408_RAW/GSM5138552_RL2373_40y_snATAC_fragments.tsv.gz"
)

# Process files in parallel using GNU Parallel
printf "%s\n" "${files[@]}" | parallel --line-buffer -j 32 process_file