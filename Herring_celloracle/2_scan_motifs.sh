#!/bin/bash
#SBATCH --job-name=motifs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michal.kubacki@external.fht.org
#SBATCH --partition=cpuq
#SBATCH --time=12:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64G
#SBATCH --output=./logs/motifs_%j.log

set -e

# Source the conda profile script
source /home/michal.kubacki/new_miniconda/etc/profile.d/conda.sh
conda activate scenicplus

PYTHON_PATH=$(which python)

"$PYTHON_PATH" /home/michal.kubacki/Githubs/Re-MEND/code/External_Datasets/GeneSet_Derivation/Herring_celloracle/scan_motifs_for_different_cell_types.py || { echo "Error running scan_motifs_for_different_cell_types.py"; exit 1; }