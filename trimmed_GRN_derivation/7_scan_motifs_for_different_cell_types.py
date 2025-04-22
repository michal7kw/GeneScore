# %%
import gc
import os
import sys
import importlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import requests
import io
import gzip
import shutil
from pathlib import Path

import celloracle as co
from celloracle import motif_analysis as ma
co.__version__

# Set working directory
work_dir = '/home/michal.kubacki/Githubs/GeneScore/trimmed_GRN_derivation'
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

# Try to import from project_functions
try:
    from grn_helpers import *
except ImportError:
    print("Warning: Could not import from project_functions path, trying absolute path")
    # Try absolute import path as fallback
    sys.path.insert(0, '/home/michal.kubacki/Githubs/GeneScore/project_functions')
    from grn_helpers import *

# %%
n_cpus = 8
neurons_set = "L2-3_CUX2"
# neurons_set = "all_ex"
# neurons_set = "all_ex_all_ages"
reference = "hg19"
root_dir = os.getenv('BASE_PATH')

# %%
output_dir, input_dir, root_dir, tmp_dir, in_dir_from_scenic = set_custom_folders(root_dir, neurons_set)
    
celltypes_dict = {
    "all_ex"                : ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "all_ex_all_ages"       : ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "L2-3_CUX2"             : ['L2-3_CUX2']
}

sel_celltypes = celltypes_dict[neurons_set]

# %% Download and prepare additional motif databases
motif_dir = os.path.join(output_dir, "motif_databases")
os.makedirs(motif_dir, exist_ok=True)

def download_file(url, output_path):
    response = requests.get(url)
    if response.status_code == 200:
        with open(output_path, 'wb') as f:
            f.write(response.content)
        print(f"Downloaded {url} to {output_path}")
        return True
    else:
        print(f"Failed to download {url}, status code: {response.status_code}")
        return False

# Download JASPAR motifs
jaspar_path = os.path.join(motif_dir, "JASPAR2022_CORE_vertebrates_non-redundant_pfms.txt")
if not os.path.exists(jaspar_path):
    download_file(
        "https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_vertebrates_non-redundant_pfms.txt", 
        jaspar_path
    )

# Download HOCOMOCO motifs
hocomoco_path = os.path.join(motif_dir, "HOCOMOCOv11_core_HUMAN_mono_meme_format.meme")
if not os.path.exists(hocomoco_path):
    download_file(
        "https://hocomoco11.autosome.org/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme",
        hocomoco_path
    )

# HOMER motifs - Assuming you have HOMER installed
homer_path = os.path.join(motif_dir, "homer_motifs")
os.makedirs(homer_path, exist_ok=True)

# %%
# Load motifs from multiple databases and combine them
# Load original CisBP motifs
cisbp_motifs = ma.load_motifs("CisBP_ver2_Homo_sapiens.pfm")
print(f"Loaded {len(cisbp_motifs)} CisBP motifs")

# Load JASPAR motifs
try:
    jaspar_motifs = ma.load_motifs(jaspar_path)
    print(f"Loaded {len(jaspar_motifs)} JASPAR motifs")
except Exception as e:
    print(f"Error loading JASPAR motifs: {e}")
    jaspar_motifs = []

# Load HOCOMOCO motifs (may need conversion depending on format support)
try:
    hocomoco_motifs = ma.load_motifs(hocomoco_path)
    print(f"Loaded {len(hocomoco_motifs)} HOCOMOCO motifs")
except Exception as e:
    print(f"Error loading HOCOMOCO motifs: {e}")
    hocomoco_motifs = []

# Combine all motifs
all_motifs = cisbp_motifs + jaspar_motifs + hocomoco_motifs
print(f"Total combined motifs: {len(all_motifs)}")

# Load the base GRN for reference
base_GRN = pd.read_parquet(os.path.join(input_dir, "2023_11_tfi.celloracle.parquet"), engine='pyarrow')
base_GRN_non_zero = base_GRN.iloc[:, 2:].astype(bool).sum().sum()

# %%
# Modified parameters for more comprehensive GRNs
fpr_threshold = 0.05  # Increased from 0.01
motif_score_threshold = 6  # Lowered from 8

# %%
for cell_type in sel_celltypes:
    print(f"\nProcessing cell type: {cell_type}")
    peaks_path = os.path.join(output_dir, f'processed_peak_file_{cell_type}.csv')

    peaks = pd.read_csv(os.path.join(output_dir, peaks_path), index_col=0)
    print(f"Loaded {len(peaks)} peaks for {cell_type}")
    
    # Create TFinfo object
    tfi = ma.TFinfo(peak_data_frame=peaks, 
                    ref_genome=reference,
                    genomes_dir=None) 

    gc.collect()

    # Scan with the combined motifs and more permissive FPR
    print(f"Scanning with {len(all_motifs)} motifs (FPR: {fpr_threshold})")
    tfi.scan(fpr=fpr_threshold, 
            motifs=all_motifs,
            verbose=True, n_cpus=n_cpus)

    # Save TF info
    file_name = os.path.join(output_dir, f"{cell_type}.celloracle.tfinfo")
    tfi.to_hdf5(file_path=file_name)

    # Display sample of scanned motifs
    print("Sample of scanned motifs:")
    print(tfi.scanned_df.head())

    # Reset and apply more permissive filtering
    tfi.reset_filtering()
    tfi.filter_motifs_by_score(threshold=motif_score_threshold)
    
    # Create TF info dataframe
    tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)

    # Save results
    file_path = os.path.join(output_dir, f"{cell_type}.celloracle.parquet")
    df = tfi.to_dataframe()
    df.to_parquet(file_path)
    
    # Calculate and print statistics
    print(f"Shape of final dataframe: {df.shape}")
    GRN_non_zero = df.iloc[:, 2:].astype(bool).sum().sum()
    print(f"Non-zero elements: {GRN_non_zero}")
    print(f"Ratio to base GRN: {GRN_non_zero/base_GRN_non_zero:.2f}")
    
    # Save summary statistics
    stats = {
        "cell_type": cell_type,
        "peaks_count": len(peaks),
        "motifs_count": len(all_motifs),
        "fpr_threshold": fpr_threshold,
        "motif_score_threshold": motif_score_threshold,
        "final_shape": df.shape,
        "non_zero_elements": GRN_non_zero,
        "ratio_to_base": GRN_non_zero/base_GRN_non_zero
    }
    
    stats_df = pd.DataFrame([stats])
    stats_df.to_csv(os.path.join(output_dir, f"{cell_type}_motif_scan_stats.csv"), index=False)
