#!/bin/bash
#SBATCH --job-name=2_Herring_celloracle_scan_motifs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michal.kubacki@external.fht.org
#SBATCH --partition=cpuq
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=256G
#SBATCH --output=./logs/2_Herring_celloracle_scan_motifs.log

# Set base directories
WORK_DIR="/home/michal.kubacki/Githubs/GeneScore/Herring_celloracle/Main_analysis/Excitatory_neurons/scripts"
cd $WORK_DIR

# Source the conda profile script
source /home/michal.kubacki/new_miniconda/etc/profile.d/conda.sh
conda activate scenicplus

PYTHON_PATH=$(which python)

"$PYTHON_PATH" "2_Herring_celloracle_scan_motifs.py"