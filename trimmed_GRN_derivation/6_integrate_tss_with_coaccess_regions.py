#  %%
print("Loading packages")
import os
import re
import gc
import sys
import pandas as pd
from celloracle import motif_analysis as ma

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

from grn_helpers import *

# %%
neurons_set = "L2-3_CUX2"
# neurons_set = "all_ex"
# neurons_set = "all_ex_all_ages"
root_dir = os.getenv('BASE_PATH')

# %%
output_dir, input_dir, root_dir, tmp_dir, in_dir_from_scenic = set_custom_folders(root_dir, neurons_set)
    
celltypes_dict = {
    "all_ex"                : ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "all_ex_all_ages"       : ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "L2-3_CUX2"             : ['L2-3_CUX2']
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
    peak_filtered = integrated[integrated.coaccess >= 0.6]
    peak_filtered = peak_filtered[["peak_id", "gene_short_name"]].reset_index(drop=True)
    
    print("Debugging: Checking the contents of peak_filtered")
    print(peak_filtered.head())
    print(f"Number of rows in peak_filtered: {len(peak_filtered)}")
    
    print(f"Saving results to {output_dir}")
    peak_filtered_path = os.path.join(output_dir, f'processed_peak_file_{cell_type}.csv')
    peak_filtered.to_csv(peak_filtered_path, index=False)
    
    print(f"Processed peak file saved for cell type {cell_type}")
    gc.collect()
# %%
