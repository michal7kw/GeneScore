#!/bin/bash
#SBATCH --job-name=2_herring_atac_processing
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michal.kubacki@external.fht.org
#SBATCH --partition=cpuq
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64G
#SBATCH --output=./logs/2_herring_atac_preprocessing.log

# Set base directories
WORK_DIR="/home/michal.kubacki/Githubs/GeneScore/trimmed_GRN_derivation"
cd $WORK_DIR

# Source the conda profile script
source /home/michal.kubacki/new_miniconda/etc/profile.d/conda.sh
conda activate scenicplus

PYTHON_PATH=$(which python)

"$PYTHON_PATH" "2_herring_atac_processing.py"