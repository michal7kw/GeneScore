#!/bin/bash
#SBATCH --job-name=coaccess
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michal.kubacki@external.fht.org
#SBATCH --partition=cpuq
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=128G
#SBATCH --output=coaccess_%j.log

set -e

# Source the conda profile script
source /home/michal.kubacki/new_miniconda/etc/profile.d/conda.sh
conda activate scenicplus

PYTHON_PATH=$(which python)
BASE_PATH="/home/michal.kubacki/Githubs/Re-MEND/code/External_Datasets/GeneSet_Derivation/Herring_scenic/Main_analysis/Get_coaccess_values"

"$PYTHON_PATH" "$BASE_PATH/1_save_celltype_consensus_regions.py" || { echo "Error running save_celltype_consensus_regions.py"; exit 1; }
"$PYTHON_PATH" "$BASE_PATH/2_cpu_coaccess.py" || { echo "Error running cpu_coaccess.py"; exit 1; }
"$PYTHON_PATH" "$BASE_PATH/3_prepare_data_for_celloracle_peaks.py" || { echo "Error running prepare_data_for_celloracle_peaks.py"; exit 1; }
