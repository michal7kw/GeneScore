#!/bin/bash
#SBATCH --job-name=1_Herring_celloracle_process_data
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michal.kubacki@external.fht.org
#SBATCH --partition=cpuq
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=64G
#SBATCH --output=./logs/1_herring_scRNA_processing.log

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

"$PYTHON_PATH" "1_herring_scRNA_processing.py"