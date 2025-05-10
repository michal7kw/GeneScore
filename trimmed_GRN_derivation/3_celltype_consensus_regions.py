# %%
import os
import sys
import pandas as pd
from pycisTopic import *

# Set working directory
# work_dir = '/home/michal.kubacki/Githubs/GeneScore/trimmed_GRN_derivation'
# work_dir = 'D:/Github/GeneScore/trimmed_GRN_derivation'
work_dir = '/mnt/d/Github/GeneScore/trimmed_GRN_derivation'

os.chdir(work_dir)

# Load environment variables from .env file
from dotenv import load_dotenv

# Explicitly specify the path to the .env file
env_path = os.path.join(work_dir, '.env')
load_dotenv(env_path)

# Get environment variables with error handling
project_functions_path = os.getenv('PROJECT_FUNCTIONS_PATH')
if not project_functions_path:
    raise ValueError("PROJECT_FUNCTIONS_PATH environment variable not found in .env file")

print(f"Using PROJECT_FUNCTIONS_PATH: {project_functions_path}")
sys.path.insert(0, project_functions_path)

from grn_helpers import set_output_folders

# %%
root_dir = os.getenv('BASE_PATH')

# %%
n_cpus = 8
# neurons_set = "L2-3_CUX2"
neurons_set = "all_ex"
# neurons_set = "all_ex_all_ages"

# %%
# all_excitatory ex_neurons ex_neurons_combined ex_progenitors all_inhibitory
out_dir, in_dir, root_dir, tmp_dir, data_folder = set_output_folders(root_dir, neurons_set)

# Read the consensus_regions.bed file using pandas
consensus_regions = pd.read_csv(os.path.join(out_dir, 'consensus_peak_calling/consensus_regions.bed'), sep='\t', header=None)
consensus_regions.columns = ['chrom', 'start', 'end', 'peak_id', 'score', 'strand']

# Define the cell types of interest
cells_dict = {
    "all_ex"            :   ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "all_ex_all_ages"   :   ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "L2-3_CUX2"         :   ['L2-3_CUX2'],
    "all_ex_comb"       :   ['ex_neurons']
}
cell_types = cells_dict[neurons_set]

# %%
# Create the peak_name column
consensus_regions['peak_name'] = consensus_regions['chrom'].astype(str) + '_' + consensus_regions['start'].astype(str) + '_' + consensus_regions['end'].astype(str)

# Define a score threshold for selecting valid peaks
score_threshold = 1.0

# Create a directory to store the cell type-specific consensus regions
cell_type_dir = os.path.join(out_dir, "cell_type_consensus_regions")
os.makedirs(cell_type_dir, exist_ok=True)

# Save the valid consensus regions for each cell type
for cell_type in cell_types:
    cell_type_regions = consensus_regions[consensus_regions['peak_id'].str.contains(cell_type) & (consensus_regions['score'] >= score_threshold)]
    cell_type_regions.loc[:, 'cell_type'] = cell_type
    cell_type_regions_file = os.path.join(cell_type_dir, f"{cell_type}_consensus_regions.bed")
    cell_type_regions.to_csv(cell_type_regions_file, sep='\t', header=False, index=False)
    print(f"Saved valid consensus regions for cell type {cell_type} to {cell_type_regions_file}")