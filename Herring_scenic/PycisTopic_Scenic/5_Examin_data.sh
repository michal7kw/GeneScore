#!/bin/bash
#SBATCH --job-name=examin_data
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michal.kubacki@external.fht.org
#SBATCH --partition=cpuq
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=256G
#SBATCH --output=./logs/examin_data_%j.log

# Source the conda profile script
source /home/michal.kubacki/new_miniconda/etc/profile.d/conda.sh
conda activate scenicplus

PYTHON_PATH=$(which python)

"$PYTHON_PATH" /home/michal.kubacki/Githubs/Re-MEND/code/External_Datasets/GeneSet_Derivation/Herring_scenic/PycisTopic_Scenic/Examin_data.py