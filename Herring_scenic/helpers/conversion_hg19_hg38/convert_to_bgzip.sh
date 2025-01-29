#!/bin/bash

cd /scratch/michal.kubacki/herring/Herring_data

# Create the output directory
mkdir -p GSE168408_RAW_hg38_bgzip

# Function to process each snATAC file
process_file() {
    file="$1"
    filename=$(basename "$file")
    output_file="GSE168408_RAW_hg38_bgzip/${filename%.gz}.bgz"

    echo "Processing file: $file"

    # Uncompress the input file and compress it with bgzip
    echo "Uncompressing and compressing with bgzip..."
    zcat "$file" | bgzip -c > "${output_file}"

    # Index the output file with tabix
    echo "Indexing output file..."
    tabix -p bed "${output_file}"

    echo "Finished processing file: $file"
}

export -f process_file

files=(
    "GSE168408_RAW_hg38/GSM5138529_RL1784_2y_snATAC_fragments.tsv.gz.gz"
    "GSE168408_RAW_hg38/GSM5138532_RL2210_4y_snATAC_fragments.tsv.gz.gz"
    "GSE168408_RAW_hg38/GSM5138510_RL2366_ga22_snATAC_fragments.tsv.gz.gz"  
    "GSE168408_RAW_hg38/GSM5138512_RL2207_ga24_snATAC_fragments.tsv.gz.gz"  
    "GSE168408_RAW_hg38/GSM5138515_RL2367_1m_snATAC_fragments.tsv.gz.gz"    
    "GSE168408_RAW_hg38/GSM5138518_RL1914_3m_snATAC_fragments.tsv.gz.gz"    
    "GSE168408_RAW_hg38/GSM5138521_RL2208_6m_snATAC_fragments.tsv.gz.gz"    
    "GSE168408_RAW_hg38/GSM5138523_RL2371_10m_snATAC_fragments.tsv.gz.gz"   
    "GSE168408_RAW_hg38/GSM5138526_RL2209_1y_snATAC_fragments.tsv.gz.gz"    
    "GSE168408_RAW_hg38/GSM5138529_RL1784_2y_snATAC_fragments.tsv.gz.gz"    
    "GSE168408_RAW_hg38/GSM5138532_RL2210_4y_snATAC_fragments.tsv.gz.gz"
    "GSE168408_RAW_hg38/GSM5138534_RL2364_6y_snATAC_fragments.tsv.gz.gz"
    "GSE168408_RAW_hg38/GSM5138536_RL1994_8y_snATAC_fragments.tsv.gz.gz"
    "GSE168408_RAW_hg38/GSM5138539_RL2368_10y_snATAC_fragments.tsv.gz.gz"
    "GSE168408_RAW_hg38/GSM5138542_RL2372_14y_snATAC_fragments.tsv.gz.gz"
    "GSE168408_RAW_hg38/GSM5138544_RL1785_16y_snATAC_fragments.tsv.gz.gz"
    "GSE168408_RAW_hg38/GSM5138548_RL2085_20y_snATAC_fragments.tsv.gz.gz"
    "GSE168408_RAW_hg38/GSM5138550_RL2369_25y_snATAC_fragments.tsv.gz.gz"
    "GSE168408_RAW_hg38/GSM5138552_RL2373_40y_snATAC_fragments.tsv.gz.gz"
)

# Process each file sequentially
for file in "${files[@]}"; do
    process_file "$file"
done

# Process files in parallel using GNU Parallel
# parallel process_file ::: "${files[@]}"