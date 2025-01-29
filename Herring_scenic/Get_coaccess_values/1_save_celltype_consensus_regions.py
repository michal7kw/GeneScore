# %%
import os
import sys
import importlib
import pandas as pd
from pycisTopic import *

sys.path.insert(0, "/home/michal.kubacki/Githubs/Re-MEND/code/External_Datasets/GeneSet_Derivation/Herring_scenic/helpers")
import config
importlib.reload(config)
from config import *
n_cpus = 32

reference = "hg19"
get_set = "all_inhibitory_all_ages"

# all_excitatory ex_neurons ex_neurons_combined ex_progenitors all_inhibitory
out_dir, in_dir, root_dir, tmp_dir, data_folder = set_output_folders(reference, get_set)

# Read the consensus_regions.bed file using pandasex_neurons_combined
consensus_regions = pd.read_csv(os.path.join(out_dir, 'consensus_peak_calling/consensus_regions.bed'), sep='\t', header=None)
consensus_regions.columns = ['chrom', 'start', 'end', 'peak_id', 'score', 'strand']

# Define the cell types of interest
cell_types_dict = {
    "all_excitatory" : ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "all_inhibitory" : ['SST', 'VIP', 'MGE_dev'],
    "all_cell_types" : ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev', 'SST', 'VIP'],
    "all_inhibitory_all_ages" : ['VIP', 'SST', 'PV', 'MGE_dev']
}
cell_types = cell_types_dict[get_set]

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