#!/bin/bash
#SBATCH --job-name=7_scan_motifs_for_different_cell_types
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michal.kubacki@external.fht.org
#SBATCH --partition=cpuq
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=128G
#SBATCH --output=./logs/7_scan_motifs_for_different_cell_types.log

# Set base directories
WORK_DIR="/home/michal.kubacki/Githubs/GeneScore/trimmed_GRN_derivation"
cd $WORK_DIR

# Source the conda profile script
source /home/michal.kubacki/new_miniconda/etc/profile.d/conda.sh
conda activate scenicplus

PYTHON_PATH=$(which python)

"$PYTHON_PATH" "7_scan_motifs_for_different_cell_types.py"