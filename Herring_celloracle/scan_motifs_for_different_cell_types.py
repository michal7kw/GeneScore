# %%
import gc
import os
import sys
import importlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import celloracle as co
from celloracle import motif_analysis as ma
co.__version__

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

sel_celltypes = celltypes_dict[neurons_set]


# Load motifs from celloracle dataset or any other custom TFs set
motifs = ma.load_motifs("CisBP_ver2_Homo_sapiens.pfm")
base_GRN = pd.read_parquet(os.path.join(input_dir, "2023_11_tfi.celloracle.parquet"), engine='pyarrow')

base_GRN_non_zero = base_GRN.iloc[:, 2:].astype(bool).sum().sum()

# %%
for cell_type in sel_celltypes:
    peaks_path = os.path.join(output_dir, f'processed_peak_file_{cell_type}.csv')

    peaks = pd.read_csv(os.path.join(output_dir, peaks_path), index_col=0)
    peaks.head()

    tfi = ma.TFinfo(peak_data_frame=peaks, 
                    ref_genome=reference,
                    genomes_dir=None) 

    gc.collect()

    tfi.scan(fpr=0.01, 
            motifs=motifs, #None
            verbose=True, n_cpus=n_cpus)

    file_name = os.path.join(output_dir, f"{cell_type}.celloracle.tfinfo")
    tfi.to_hdf5(file_path = file_name)

    tfi.scanned_df.head()

    tfi.reset_filtering()

    tfi.filter_motifs_by_score(threshold=8)

    tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)

    tfi.scanned_df.head()

    file_path = os.path.join(output_dir, f"{cell_type}.celloracle.parquet")

    df = tfi.to_dataframe()
    df.to_parquet(file_path)

    df.shape

    GRN_non_zero = df.iloc[:, 2:].astype(bool).sum().sum()

    print(GRN_non_zero/base_GRN_non_zero)
