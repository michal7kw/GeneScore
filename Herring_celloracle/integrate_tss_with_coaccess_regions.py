#  %%
print("Loading packages")
import os
import re
import gc
import sys
import importlib
import pandas as pd

from celloracle import motif_analysis as ma
import celloracle as co

sys.path.insert(0, "/home/michal.kubacki/Githubs/Re-MEND/code/External_Datasets/GeneSet_Derivation/Herring_celloracle/helpers")
import config
importlib.reload(config)
from config import *
n_cpus = 32

# %%
reference = "hg19"

# neurons_set = "all_excitatory"
# neurons_set = "all_inhibitory"
neurons_set = "all_excitatory_all_ages"
# neurons_set = "all_inhibitory_all_ages"

output_dir, input_dir, root_dir, tmp_dir, in_dir_from_scenic = set_custom_folders("hg19", neurons_set)
    
celltypes_dict = {
    "all_excitatory"                : ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "all_inhibitory"                : ['SST', 'VIP', 'MGE_dev'],
    "all_excitatory_all_ages"       : ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "all_inhibitory_all_ages"       : ['VIP', 'SST', 'PV', 'MGE_dev']
}

cell_types = celltypes_dict[neurons_set]

coaccess_files = [f"{cell_type}_coaccess.csv" for cell_type in cell_types]

peaks_files = [f"{cell_type}_peaks.csv" for cell_type in cell_types]

def extract_cell_type(file_name):
    pattern = r'^(.+)_coaccess\.csv$'
    match = re.match(pattern, file_name)
    if match:
        return match.group(1)
    else:
        return None

# # %%
# print("Loading tss annotations")
# tss_path = os.path.join(output_dir, 'tss_annotated.csv')
# tss_annotated = pd.read_csv(tss_path)

# %%
print("Processing loop")
for coaccess_file, peaks_file in zip(coaccess_files, peaks_files):
    print(f"Processing {coaccess_file}, {peaks_file}")
    
    print("Loading connections")
    coaccess_path = os.path.join(in_dir_from_scenic, coaccess_file)
    cicero_connections = pd.read_csv(coaccess_path)
    print("Debugging: Checking the contents of cicero_connections")
    print(cicero_connections.head())
    
    print("Loading peaks")
    peak_path = os.path.join(in_dir_from_scenic, peaks_file)
    with open(peak_path, 'r') as file:
        peaks = file.read().split()
    print("Debugging: Checking the contents of peaks")
    print(peaks[:5])
    
    print("Formating peaks")
    # peaks = [peak.strip('"') for peak in peaks]

    tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome="hg19")
    
    print("TSS integration")
    integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss_annotated,
                                                   cicero_connections=cicero_connections)
    
    print("Filtering peaks")
    cell_type = extract_cell_type(coaccess_file)
    peak_filtered = integrated[integrated.coaccess >= 0.8]
    peak_filtered = peak_filtered[["peak_id", "gene_short_name"]].reset_index(drop=True)
    
    print("Debugging: Checking the contents of peak_filtered")
    print(peak_filtered.head())
    print(f"Number of rows in peak_filtered: {len(peak_filtered)}")
    
    print(f"Saving results to {output_dir}")
    peak_filtered_path = os.path.join(output_dir, f'processed_peak_file_{cell_type}.csv')
    peak_filtered.to_csv(peak_filtered_path, index=False)
    
    print(f"Processed peak file saved for cell type {cell_type}")
    gc.collect()