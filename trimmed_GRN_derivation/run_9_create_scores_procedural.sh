#!/bin/bash
#SBATCH --job-name=9_create_scores_procedural
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michal.kubacki@external.fht.org
#SBATCH --partition=cpuq
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=128G
#SBATCH --output=./logs/9_create_scores_procedural.log

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

"$PYTHON_PATH" "9_create_scores_procedural.py"