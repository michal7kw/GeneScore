#!/bin/bash
#SBATCH --job-name=lable_transfer
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michal.kubacki@external.fht.org
#SBATCH --partition=cpuq
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=128G
#SBATCH --output=./logs/lable_transfer_%j.log

# Source the conda profile script
source /home/michal.kubacki/new_miniconda/etc/profile.d/conda.sh
conda activate scenicplus

PYTHON_PATH=$(which python)

BASE_DIR="/home/michal.kubacki/Githubs/GeneScore/Herring_scenic/PycisTopic_Scenic"

"$PYTHON_PATH" "$BASE_DIR/label_transfer.py"
