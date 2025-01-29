# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.6
#   kernelspec:
#     display_name: Python (scenicplus)
#     language: python
#     name: scenicplus
# ---

# 1. Data Exploration
# 2. Importing and Formating snATAC data
# 3. Getting pseudobulk profiles from cell annotations
# 4. Infering consensus peaks

# # Environment

# +
# Standard library imports
import os
import gc
import sys
import pickle
import importlib

# Data manipulation imports
import pandas as pd

# Visualization imports
import matplotlib.pyplot as plt
import seaborn as sns

# pycisTopic imports
import pycisTopic
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, peak_calling
from pycisTopic.iterative_peak_calling import get_consensus_peaks

importlib.reload(pycisTopic)
from pycisTopic import *
pycisTopic.__version__

sys.path.insert(0, "/home/michal.kubacki/Githubs/Re-MEND/code/External_Datasets/GeneSet_Derivation/Herring_scenic/helpers")
import config
importlib.reload(config)
from config import *
n_cpu = 32


# +
#################################################################
reference = "hg19"


# neurons_set = "all_excitatory"
# neurons_set = "all_inhibitory"
neurons_set = "all_excitatory_all_ages"
# neurons_set = "all_inhibitory_all_ages"

cells_dict = {
    "all_inhibitory"            :   ['SST', 'VIP', 'MGE_dev'],
    "all_inhibitory_all_ages"   :   ['VIP', 'SST', 'PV', 'MGE_dev'],
    "all_excitatory"            :   ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "all_excitatory_all_ages"   :   ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev']
}

ages_dict = {
    "all_inhibitory"            :   ['1m','3m','6m','10m','1y','2y','4y','ga22','ga24'],
    "all_inhibitory_all_ages"   :   ['1m','3m','6m','10m','1y','2y','4y','6y','10y','16y','20y','40y','ga22','ga24'],
    "all_excitatory"            :   ['1m','3m','6m','10m','1y','2y','4y','ga22','ga24'],
    "all_excitatory_all_ages"   :   ['1m','3m','6m','10m','1y','2y','4y','6y','10y','16y','20y','40y','ga22','ga24']
}

out_dir, in_dir, root_dir, tmp_dir, data_folder = set_output_folders(reference, neurons_set)

sel_celltypes  = cells_dict[neurons_set]
sel_ages = ages_dict[neurons_set]

#################################################################
# -

ATAC_metadata_path = os.path.join(in_dir, "Processed_data_ATAC_BCs-meta-data.csv")

fragments_dict = select_files(reference, selected_fragments = sel_ages)

fragments_dict

# # Load ATAC metadata

pd.set_option('display.max_columns', None)
cells_data = pd.read_csv(ATAC_metadata_path, sep=",", index_col = 0)
cells_data.head()


# ## Add new columns

# +
def process_age(name):
    parts = name.split('-')[1]
    processed_name = parts.split('_')[1]
    return processed_name

def process_name(name):
    parts = name.split('_')[:2]
    processed_name = '_'.join(parts)
    processed_name = processed_name.replace('/', '-')
    return processed_name

def process_chem(name):
    processed_name = name.split('_')[-1]
    return processed_name

def map_major_clust(name, mapped_names):
    for prefix, name_list in mapped_names.items():
        if name in name_list:
            return prefix
    return process_name(name)

cells_data['age'] = cells_data.predictedCell.apply(process_age)
cells_data['chem'] = cells_data.predictedCell.apply(process_chem)

unique_names = cells_data.predictedGroup.unique()

prefixes = sel_celltypes
mapped_names = {}

for prefix in prefixes:
    mapped_names[prefix] = []
    for name in unique_names:
        if name.startswith(prefix):
            mapped_names[prefix].append(name)

cells_data['major_clust'] = cells_data.predictedGroup.apply(lambda x: map_major_clust(x, mapped_names))

# +
# def process_age(name):
#     parts = name.split('-')[1]
#     processed_name = parts.split('_')[1]
#     return processed_name

# def process_name(name):
#     parts = name.split('_')[:2]
#     processed_name = '_'.join(parts)
#     processed_name = processed_name.replace('/', '-')
#     return processed_name

# def process_chem(name):
#     processed_name = name.split('_')[-1]
#     return processed_name

# cells_data['age'] = cells_data.predictedCell.apply(process_age)
# cells_data['chem'] = cells_data.predictedCell.apply(process_chem)

# unique_names = cells_data.predictedGroup.unique()
# name_mapping = {name: process_name(name) for name in unique_names}
# cells_data['major_clust'] = cells_data.predictedGroup.map(name_mapping)
# -

# ## Match `ages` from the fragments files

# +
mapping = {
    '2d': '1m',
    '34d': '1m',
    '86d': '3m',
    '118d': '3m',
    '179d': '6m',
    '301d': '10m',
    '422d': '1y',
    '2yr': '2y',
    '627d': '2y',
    '3yr': '4y',
    '4yr': '4y',
    '6yr': '6y',
    '8yr': '8y',
    '10yr': '10y',
    '12yr': '14y',
    '14yr': '14y',
    '16yr': '16y',
    '17yr': '16y',
    '20yr': '20y',
    '25yr': '25y',
    '40yr': '40y',
    'ga22': 'ga22',
    'ga24': 'ga24',
    'ga34': 'ga24'
}

cells_data["age_mapped"] = [mapping.get(age, age) for age in cells_data.age]
cells_data["age_mapped"].unique()
# -

# # Filter ATAC metadata

print(cells_data.shape)

cells_data = cells_data[cells_data['major_clust'].isin(sel_celltypes)]
print(cells_data.shape)

print(cells_data['chem'].value_counts())
print(cells_data['PassQC'].value_counts())

cells_data = cells_data[cells_data['chem']=="v3"]
print(cells_data.shape)

cells_data = cells_data[cells_data['age_mapped'].isin(sel_ages)]
print(cells_data.shape)


# # Format indexes

# +
def format_index(index):
    parts = index.split("#")
    formatted_index = parts[1]
    return formatted_index

cells_data = cells_data.rename(index=format_index)
# -

cells_data["old_index"] = cells_data.index
cells_data.index = cells_data.index + "-" + cells_data["age_mapped"]

cells_data.head()

cells_data = cells_data[cells_data.age_mapped.isin(fragments_dict)]
cells_data.head()

print(cells_data.shape)

for cell in cells_data['major_clust'].unique():
    print(f"{cell}: {(cells_data['major_clust']==cell).sum()}")

cells_data.to_csv(os.path.join(out_dir, 'cells_data.csv'), index=True)

pd.reset_option('display.max_columns')

# # Getting pseudobulk profiles from cell annotations

cell_data = pd.read_csv(os.path.join(out_dir, 'cells_data.csv'), index_col = 0)
cell_data.head()

chromsizes = pd.read_table(
    os.path.join(in_dir, f"{reference}.chrom.sizes"),
    header = None,
    names = ["Chromosome", "End"]
)
chromsizes.insert(1, "Start", 0)
chromsizes.head()

# potentialy beforehand you  might want to run `for file in *fragments.tsv.gz; do tabix -p bed "$file"; done`

gc.collect()

fragments_dict

# +
#########  Requied modification to the export_pseudobulk function ###

# bed_paths = {}
# for cell_type in cell_data[variable].unique():
#     _bed_fname = os.path.join(
#         bed_path,
#         f"{_santize_string_for_filename(cell_type)}.fragments.tsv.gz")
#     if os.path.exists(_bed_fname):
#         bed_paths[cell_type] = _bed_fname
#     else:
#         log.warning(f"Missing fragments for {cell_type}!")

# # log.info("generating bigwig files")
# # joblib.Parallel(n_jobs=n_cpu)(
# #     joblib.delayed(_generate_bigwig)
# #     (
# #         path_to_fragments = bed_paths[cell_type],
# #         chromsizes = chromsizes_dict,
# #         normalize_bigwig = normalize_bigwig,
# #         bw_filename = os.path.join(bigwig_path, f"{_santize_string_for_filename(cell_type)}.bw"),
# #         log = log
# #     )
# #     for cell_type in bed_paths.keys()
# # )
# # bw_paths = {}
# # for cell_type in cell_data[variable].unique():
# #     _bw_fname = os.path.join(
# #         bigwig_path,
# #         f"{_santize_string_for_filename(cell_type)}.bw")
# #     if os.path.exists(_bw_fname):
# #         bw_paths[cell_type] = _bw_fname
# #     else:
# #         log.warning(f"Missing bigwig for {cell_type}!")

# # return bw_paths, bed_paths
# return  bed_paths

# +
# # %%script false --no-raise-error

os.makedirs(os.path.join(out_dir, "consensus_peak_calling"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bw_files"), exist_ok = True)

paths = export_pseudobulk(
    input_data = cells_data,
    variable = "major_clust",
    sample_id_col = "age_mapped",
    chromsizes = chromsizes,
    bigwig_path = os.path.join(root_dir, "", "consensus_peak_calling/pseudobulk_bw_files"),
    bed_path = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"),
    path_to_fragments = fragments_dict,
    n_cpu = 1,
    temp_dir = tmp_dir,
    split_pattern = "-"
)
# -

bw_path, bed_paths = paths

# # %%script false --no-raise-error
with open(os.path.join(out_dir, "consensus_peak_calling/bed_paths.tsv"), "wt") as f:
    for v in bed_paths:
        _ = f.write(f"{v}\t{bed_paths[v]}\n")

# +
# # %%script false --no-raise-error

# directory = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files")

# # Get the list of file names in the directory
# file_names = os.listdir(directory)
# # Create a dictionary to store the file paths
# file_paths = {}

# # Iterate over the file names and store their paths in the dictionary
# for file_name in file_names:
#     file_path = os.path.join(directory, file_name)
#     file_paths[file_name] = file_path

# # Specify the output directory and file name
# output_file = "bed_paths.tsv"

# # Write the file paths to the output file
# with open(os.path.join(out_dir, "consensus_peak_calling/", output_file), "wt") as f:
#     for file_name, file_path in file_paths.items():
#         f.write(f"{file_name}\t{file_path}\n")
# -

# # Inferring consensus peaks

bed_paths = {}
with open(os.path.join(out_dir, "consensus_peak_calling/bed_paths.tsv")) as f:
    for line in f:
        v, p = line.strip().split("\t")
        bed_paths.update({v: p})

bed_paths

# +
# # %%script false --no-raise-error 
import logging

macs_path = "/home/michal.kubacki/.conda/envs/scenicplus/bin/macs2" 
# macs_path = "macs2"

os.makedirs(os.path.join(out_dir, "consensus_peak_calling/MACS"), exist_ok = True)

narrow_peak_dict = peak_calling(
    macs_path = macs_path,
    bed_paths = bed_paths,
    outdir = os.path.join(os.path.join(out_dir, "consensus_peak_calling/MACS")),
    genome_size = 'hs',
    n_cpu = 1, # n_cpu,
    input_format = 'BEDPE',
    shift = 73,
    ext_size = 146,
    keep_dup = 'all',
    q_value = 0.05,
    _temp_dir = '/tmp',
    skip_empty_peaks=True,
    logging_level=logging.DEBUG
)

# +
# # %%script false --no-raise-error

# Other param
peak_half_width=250
path_to_blacklist=os.path.join(in_dir, f"{reference}-blacklist.v2.bed")
# Get consensus peaks
consensus_peaks = get_consensus_peaks(
    narrow_peaks_dict = narrow_peak_dict,
    peak_half_width = peak_half_width,
    chromsizes = chromsizes,
    path_to_blacklist = path_to_blacklist)

# +
# # %%script false --no-raise-error

consensus_peaks.to_bed(
    path = os.path.join(out_dir, "consensus_peak_calling/consensus_regions.bed"),
    keep = True,
    compression = 'infer',
    chain = False)
# -








