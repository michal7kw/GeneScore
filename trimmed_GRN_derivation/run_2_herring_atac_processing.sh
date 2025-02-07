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

# Load environment variables from .env file
set -a
source .env
set +a

# Change to work directory from .env
cd $WORK_DIR

# Source the conda profile script (use path from .env)
source $CONDA_PROFILE_PATH
conda activate scenicplus

PYTHON_PATH=$(which python)

"$PYTHON_PATH" "2_herring_atac_processing.py"