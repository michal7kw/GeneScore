# %%
import os
import sys
import importlib
import pandas as pd
from pycisTopic import *

# %%
n_cpus = 8
neurons_set = "all_ex"
root_dir = "/group/testa/michal.kubacki/herring_minimal"

# %%
def set_output_folders(suffix):
    tmp_dir = os.path.join(root_dir, "tmp")
    in_dir = os.path.join(root_dir, "data")
    out_dir = os.path.join(root_dir, suffix)

    print(f"root_dir: {root_dir}")
    print(f"out_dir: {out_dir}")
    print(f"in_dir: {in_dir}")
    print(f"tmp_dir: {tmp_dir}")

    data_folder = "GSE168408_RAW"

    os.makedirs(out_dir, exist_ok = True)
    return out_dir, in_dir, root_dir, tmp_dir, data_folder

# %%
# all_excitatory ex_neurons ex_neurons_combined ex_progenitors all_inhibitory
out_dir, in_dir, root_dir, tmp_dir, data_folder = set_output_folders(neurons_set)

# Read the consensus_regions.bed file using pandas
consensus_regions = pd.read_csv(os.path.join(out_dir, 'consensus_peak_calling/consensus_regions.bed'), sep='\t', header=None)
consensus_regions.columns = ['chrom', 'start', 'end', 'peak_id', 'score', 'strand']

# Define the cell types of interest
cell_types_dict = {
    "all_ex" : ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
}
cell_types = cell_types_dict[neurons_set]

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