#!/bin/bash
#SBATCH --job-name=scenicplus
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michal.kubacki@external.fht.org
#SBATCH --partition=cpuq
#SBATCH --time=32:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=512G
#SBATCH --output=./logs/scenicplus_%j.log

cd /home/michal.kubacki/Githubs/Re-MEND/code/External_Datasets/GeneSet_Derivation/Herring_scenic/scplus_pipeline/Snakemake

# Source the conda.sh script directly
source /home/michal.kubacki/new_miniconda/etc/profile.d/conda.sh
source /home/michal.kubacki/new_miniconda/bin/activate scenicplusnew

# Set the full path to the snakemake executable within the scenicplus environment
SNAKEMAKE_PATH=$(which snakemake)
# Set the base path for config files
CONFIG_BASE_PATH="/home/michal.kubacki/Githubs/Re-MEND/code/External_Datasets/GeneSet_Derivation/Herring_scenic/scplus_pipeline/Snakemake/config"

# Check if a config file name is provided as an argument
if [ $# -eq 0 ]; then
    echo "Error: No config file name provided."
    echo "Usage: sbatch $0 <config_file_name>"
    exit 1
fi

# Construct the full path to the config file
CONFIG_FILE="${CONFIG_BASE_PATH}/${1}"

# Check if the config file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Config file not found: $CONFIG_FILE"
    exit 1
fi

# --forceall  --unlock
"$SNAKEMAKE_PATH" --cores 96 --forceall --configfile "$CONFIG_FILE"


######## Use Example #######
# sbatch 6_scenic_plus.sh config_all_excitatory.yaml
# sbatch 6_scenic_plus.sh config_all_excitatory_all_ages.yaml

# sbatch 6_scenic_plus.sh config_all_inhibitory.yaml
# sbatch 6_scenic_plus.sh config_all_inhibitory_all_ages.yaml

# sbatch 6_scenic_plus.sh config_all_cells.yaml
# sbatch 6_scenic_plus.sh config_all_excitatory_all_ages_relaxed.yaml

# sbatch 6_scenic_plus.sh /home/michal.kubacki/Githubs/Re-MEND/code/External_Datasets/GeneSet_Derivation/Herring_scenic/scplus_pipeline/Snakemake/config/config_all_excitatory_all_ages.yaml