#!/bin/bash
#SBATCH --job-name=integrate
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michal.kubacki@external.fht.org
#SBATCH --partition=cpuq
#SBATCH --time=12:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --mem=64G
#SBATCH --output=./logs/integrate_%j.log

set -e

# Source the conda profile script
source /home/michal.kubacki/new_miniconda/etc/profile.d/conda.sh
conda activate scenicplus

PYTHON_PATH=$(which python)

"$PYTHON_PATH" /home/michal.kubacki/Githubs/Re-MEND/code/External_Datasets/GeneSet_Derivation/Herring_celloracle/integrate_tss_with_coaccess_regions.py || { echo "Error running prepare_data_for_celloracle_peaks.py"; exit 1; }