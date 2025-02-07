#!/bin/bash
#SBATCH --job-name=8_celloracle_scan_motifs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michal.kubacki@external.fht.org
#SBATCH --partition=cpuq
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64G
#SBATCH --output=./logs/8_celloracle_scan_motifs.log

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

"$PYTHON_PATH" "8_celloracle_scan_motifs.py"