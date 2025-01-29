#!/bin/bash
#SBATCH --job-name=pycistopic_run_full
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michal.kubacki@external.fht.org
#SBATCH --partition=cpuq
#SBATCH --time=12:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --mem=64G
#SBATCH --output=pycistopic_run_full_%j.log

# Source the conda profile script
source /home/michal.kubacki/new_miniconda/etc/profile.d/conda.sh
conda activate scenicplus

zcat /group/testa/michal.kubacki/herring/data/GSE168408_RAW/merged_output.tsv.gz | awk '{if(NF==1) {age=$1; next} else {print $0"\t"age}}' | gzip > /group/testa/michal.kubacki/herring/data/GSE168408_RAW/corrected_merged_output.tsv.gz
