# %% [markdown]
# 1. Data Exploration
# 2. Importing and Formating snATAC data
# 3. Getting pseudobulk profiles from cell annotations
# 4. Infering consensus peaks

# 1. input:
#     1. `Processed_data_ATAC_BCs-meta-data.csv`
# 2. output:
#     1. `consensus_regions.bed`
#     2. `cistopic_obj.pkl`

# %% # Load libraries
# Standard library imports
import os
import gc
import sys
import pickle
import logging

# Data manipulation imports
import pandas as pd
import polars as pl

# pycisTopic imports
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, peak_calling
from pycisTopic.iterative_peak_calling import get_consensus_peaks
from pycisTopic.plotting.qc_plot import plot_sample_stats, plot_barcode_stats
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments, merge
from pycisTopic.qc import get_barcodes_passing_qc_for_sample

import matplotlib.pyplot as plt
import scrublet as scr

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

from grn_helpers import set_output_folders, select_files

# %% # Load environment variables
root_dir = os.getenv('BASE_PATH')

# %% # Define parameters
n_cpu = 20
# neurons_set = "L2-3_CUX2"
# neurons_set = "all_ex"
neurons_set = "all_ex_comb"
# neurons_set = "all_ex_all_ages"

# %% # Define parameters
cells_dict = {
    "all_ex"            :   ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "all_ex_all_ages"   :   ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "L2-3_CUX2"         :   ['L2-3_CUX2'],
    "all_ex_comb"       :   ['ex_neurons']
}

ages_dict = {
    "all_ex"           :   ['1m','3m','6m','10m','1y','2y','4y','ga22','ga24'],
    "all_ex_all_ages"  :   ['1m','3m','6m','10m','1y','2y','4y','6y','10y','16y','20y','40y','ga22','ga24'],
    "L2-3_CUX2"        :   ['1m','3m','6m','10m','1y','2y','4y','ga22','ga24'],
    "all_ex_comb"      :   ['1m','3m','6m','10m','1y','2y','4y','ga22','ga24']
}

out_dir, in_dir, root_dir, tmp_dir, data_folder = set_output_folders(root_dir, neurons_set)
sel_celltypes = cells_dict[neurons_set]
sel_ages = ages_dict[neurons_set]

# %% # Load ATAC metadata
ATAC_metadata_path = os.path.join(in_dir, "Processed_data_ATAC_BCs-meta-data.csv")

# %% # load fragments data
fragments_dict = select_files(root_dir, selected_fragments = sel_ages)

# %% # load cells data
cells_data = pd.read_csv(ATAC_metadata_path, sep=",", index_col=0)

# %% # process cells data
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
print(f"DEBUG: Unique values in cells_data.predictedGroup after load: {sorted(list(unique_names))}")
cells_data['original_clust'] = cells_data.predictedGroup.apply(process_name)

if neurons_set == "all_ex_comb":
    # Specific mapping for "all_ex_comb"
    ex_neuron_original_types = ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS']
    
    # Create a new column to store the original cell type
    # The line `cells_data['original_clust'] = cells_data.predictedGroup.apply(process_name)` was moved before this `if` block.
    
    
    # Map all ex_neuron_original_types to 'ex_neurons' for filtering
    def map_func_all_ex_comb(name_val):
        processed_name_val = process_name(name_val) # First process the name
        if processed_name_val in ex_neuron_original_types: # Then check if the processed name is in the list
            return "ex_neurons"
        return processed_name_val # Return the processed name if not in the list (e.g. for non-excitatory types)
    
    cells_data['major_clust'] = cells_data.predictedGroup.apply(map_func_all_ex_comb)
else:
    # Original mapping logic
    prefixes = sel_celltypes
    mapped_names = {}
    for prefix in prefixes:
        mapped_names[prefix] = []
        for name_iter in unique_names:
            if name_iter.startswith(prefix):
                mapped_names[prefix].append(name_iter)
    cells_data['major_clust'] = cells_data.predictedGroup.apply(lambda x: map_major_clust(x, mapped_names))

# %% # mapping ages
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

# %% # filter cells data
print(cells_data.shape)
print(f"DEBUG: Unique values in cells_data.major_clust before filtering by sel_celltypes: {sorted(list(cells_data['major_clust'].unique()))}")
print(f"DEBUG: sel_celltypes for this run: {sel_celltypes}")

cells_data = cells_data[cells_data['major_clust'].isin(sel_celltypes)]
print(cells_data.shape)

print(cells_data['chem'].value_counts())
print(cells_data['PassQC'].value_counts())

cells_data = cells_data[cells_data['chem']=="v3"]
print(cells_data.shape)

cells_data = cells_data[cells_data['age_mapped'].isin(sel_ages)]
print(cells_data.shape)

# %% # format indexes
def format_index(index):
    parts = index.split("#")
    formatted_index = parts[1]
    return formatted_index

cells_data = cells_data.rename(index=format_index)

cells_data["old_index"] = cells_data.index
cells_data.index = cells_data.index + "-" + cells_data["age_mapped"]
cells_data = cells_data[cells_data.age_mapped.isin(fragments_dict)]

print(cells_data.shape)

for cell in cells_data['major_clust'].unique():
    print(f"{cell}: {(cells_data['major_clust']==cell).sum()}")

# %% # Create a dictionary to store the count of each cell type
cell_counts = {}

# Count the occurrences of each cell type
for cell in cells_data['major_clust'].unique():
    count = (cells_data['major_clust'] == cell).sum()
    cell_counts[cell] = count
    print(f"{cell}: {count}")

# Create a mask to filter cell types with count >= 1000
mask = cells_data['major_clust'].isin([cell for cell, count in cell_counts.items() if count >= 800])

# Subset cells_data based on the mask
cells_data = cells_data[mask]

# %% # print cell counts
for cell in cells_data['major_clust'].unique():
    print(f"{cell}: {(cells_data['major_clust']==cell).sum()}")

# cells_data = pd.read_csv(os.path.join(out_dir, 'cells_data.csv'), index_col = 0)
cells_data.head()

# %% # load chromsizes
chromsizes = pd.read_table(
    os.path.join(in_dir, f"hg19.chrom.sizes"),
    header = None,
    names = ["Chromosome", "End"]
)
chromsizes.insert(1, "Start", 0)
 
# %% # collect garbage
gc.collect()

# %% # create output folders
os.makedirs(os.path.join(out_dir, "consensus_peak_calling"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bw_files"), exist_ok = True)

# %% # export pseudobulk profiles
# beforehand one might need to run `for file in *fragments.tsv.gz; do tabix -p bed "$file"; done`
# For all_ex_comb, use original_clust instead of major_clust for export_pseudobulk
# This preserves the original cell type information for the function
variable_to_use = "original_clust" if neurons_set == "all_ex_comb" else "major_clust"

# Ensure fragments_dict keys match the samples (ages) actually present in cells_data
# after all filtering. The sample_id_col used in export_pseudobulk is 'age_mapped'.
final_present_ages_in_cells_data = set(cells_data['age_mapped'].unique())
original_fragments_dict_keys = sorted(list(fragments_dict.keys())) # For logging
fragments_dict = {
    age: path
    for age, path in fragments_dict.items()
    if age in final_present_ages_in_cells_data
}
# Logging the filtering action
print(f"Original fragments_dict had keys: {original_fragments_dict_keys}")
print(f"After filtering based on cells_data['age_mapped'], fragments_dict now has keys: {sorted(list(fragments_dict.keys()))}")
print(f"Unique ages in final cells_data['age_mapped']: {sorted(list(final_present_ages_in_cells_data))}")

if not fragments_dict:
    raise ValueError("fragments_dict is empty after filtering against cells_data. "
                     "This means no common samples (ages) are left after processing cells_data. "
                     "Check filtering steps for cells_data and initial sample lists.")

paths = export_pseudobulk(
    input_data = cells_data,
    variable = variable_to_use,  # Use original_clust for all_ex_comb
    sample_id_col = "age_mapped",
    chromsizes = chromsizes,
    bigwig_path = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bw_files"),
    bed_path = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"),
    path_to_fragments = fragments_dict,
    n_cpu = 1,
    temp_dir = '/tmp',
    split_pattern = "-"
)

# %% # create bed paths
bw_paths, bed_paths = paths

with open(os.path.join(out_dir, "consensus_peak_calling/bed_paths.tsv"), "wt") as f:
    for v in bed_paths:
        f.write(f"{v}\t{bed_paths[v]}\n")

# %% # infer consensus peaks
bed_paths = {}
with open(os.path.join(out_dir, "consensus_peak_calling/bed_paths.tsv")) as f:
    for line in f:
        v, p = line.strip().split("\t")
        bed_paths.update({v: p})

#%% # peak_calling
# macs_path = "/home/michal.kubacki/.conda/envs/scenicplus/bin/macs2" 
macs_path = "/home/michal/miniforge3/envs/scenicplus-new/bin/macs2"
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/MACS"), exist_ok=True)

narrow_peak_dict = peak_calling(
    macs_path = macs_path,
    bed_paths = bed_paths,
    outdir = os.path.join(out_dir, "consensus_peak_calling/MACS"),
    genome_size = 'hs',
    n_cpu = 1,
    input_format = 'BEDPE',
    shift = 73,
    ext_size = 146,
    keep_dup = 'all',
    q_value = 0.05,
    _temp_dir = '/tmp',
    skip_empty_peaks=True,
    logging_level=logging.DEBUG
)

# %% # get consensus peaks
peak_half_width = 250
path_to_blacklist = os.path.join(in_dir, "hg19-blacklist.v2.bed")
consensus_peaks = get_consensus_peaks(
    narrow_peaks_dict=narrow_peak_dict,
    peak_half_width=peak_half_width,
    chromsizes=chromsizes,
    path_to_blacklist=path_to_blacklist)

# Pre-compute column types and optimize DataFrame operations
df = consensus_peaks.df

# Create a mask for object columns to avoid repeated dtype checks
object_cols = df.select_dtypes(include=['object']).columns

# Process all columns in one pass using vectorized operations where possible
for col in object_cols:
    # Check if first non-null value is a set
    first_valid = df[col].first_valid_index()
    if first_valid is not None and isinstance(df.at[first_valid, col], set):
        df[col] = df[col].apply(list)
    else:
        # Use vectorized operation for checking sets in column
        has_sets = df[col].apply(lambda x: isinstance(x, set) if pd.notnull(x) else False).any()
        if has_sets:
            df[col] = df[col].apply(lambda x: list(x) if isinstance(x, set) else x)

# Process list columns more efficiently
list_columns = ['Name', 'Score', 'Strand', 'ThickStart', 'ThickEnd', 'ItemRGB', 'BlockCount', 'BlockSizes', 'BlockStarts']
list_cols_present = [col for col in list_columns if col in df.columns]

if list_cols_present:
    for col in list_cols_present:
        df[col] = df[col].apply(lambda x: [] if pd.isna(x) else (
            list(x) if isinstance(x, (set, tuple)) else (
                [x] if isinstance(x, (int, float, str)) and not isinstance(x, (list, tuple)) else x
            )
        ))

# Handle index sets if present (rare case)
if isinstance(df.index, pd.Index) and any(isinstance(x, set) for x in df.index):
    df.index = pd.Index([list(x) if isinstance(x, set) else x for x in df.index])

# Convert any remaining sets in all columns
for col in df.columns:
    if df[col].dtype == 'object':
        df[col] = df[col].apply(lambda x: list(x) if isinstance(x, set) else x)

# Update columns in place instead of replacing the entire DataFrame
for col in df.columns:
    consensus_peaks.df[col] = df[col]

# Force conversion of any remaining sets in noncanonical columns
noncanonical = list(set(consensus_peaks.df.columns) - {'Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand', 'ThickStart', 'ThickEnd', 'ItemRGB', 'BlockCount', 'BlockSizes', 'BlockStarts'})
for col in noncanonical:
    if consensus_peaks.df[col].dtype == 'object':
        consensus_peaks.df[col] = consensus_peaks.df[col].apply(lambda x: list(x) if isinstance(x, set) else x)

del df  # Clean up reference

consensus_peaks.to_bed(
    path=os.path.join(out_dir, "consensus_peak_calling/consensus_regions.bed"),
    keep=noncanonical,  # Pass the list instead of set
    compression='infer',
    chain=False)

# %% # collect garbage
gc.collect()

# %% # create regions bed filename
regions_bed_filename = os.path.join(out_dir, "consensus_peak_calling", "consensus_regions.bed")
tss_bed_path = os.path.join(out_dir, "qc")
os.makedirs(tss_bed_path, exist_ok=True)
tss_bed_filename = os.path.join(tss_bed_path, "tss.bed")

pycistopic_qc_commands_filename = "pycistopic_qc_commands.txt"

# Create text file with all pycistopic qc command lines.
with open(pycistopic_qc_commands_filename, "w") as fh:
    for sample, fragment_filename in fragments_dict.items():
        print(
            "pycistopic qc",
            f"--fragments {fragment_filename}",
            f"--regions {regions_bed_filename}",
            f"--tss {tss_bed_filename}",
            f"--output {os.path.join(out_dir, 'qc', sample)}",
            sep=" ",
            file=fh,
        )

# %% # run pycistopic qc
import subprocess

with open("pycistopic_qc_commands.txt", "r") as file:
    commands = file.readlines()
    for command in commands:
        process = subprocess.Popen(command.strip().split(" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            print(f"Error executing command: {command.strip()}")
            print(f"Error message: {stderr.decode('utf-8')}")

# %% # collect garbage
gc.collect()

# %% # plot sample stats
for sample_id in fragments_dict:
    fig = plot_sample_stats(
        sample_id=sample_id,
        pycistopic_qc_output_dir=os.path.join(out_dir, "qc")
    )
    plt.savefig(os.path.join(out_dir, f"sample_stats_{sample_id}.png"))
    plt.close(fig)

# %% # get barcodes passing filters
sample_id_to_barcodes_passing_filters = {}
sample_id_to_thresholds = {}
for sample_id in fragments_dict:
    (
        sample_id_to_barcodes_passing_filters[sample_id],
        sample_id_to_thresholds[sample_id]
    ) = get_barcodes_passing_qc_for_sample(
            sample_id=sample_id,
            pycistopic_qc_output_dir=os.path.join(out_dir, "qc"),
            unique_fragments_threshold=None,
            tss_enrichment_threshold=None,
            frip_threshold=0,
            use_automatic_thresholds=True,
    )

# %% # save barcodes passing filters
with open(os.path.join(out_dir, "sample_id_to_barcodes_passing_filters.pkl"), "wb") as file:
    pickle.dump(sample_id_to_barcodes_passing_filters, file)

# %% # plot barcode stats
for sample_id in fragments_dict:
    fig = plot_barcode_stats(
        sample_id=sample_id,
        pycistopic_qc_output_dir=os.path.join(out_dir, "qc"),
        bc_passing_filters=sample_id_to_barcodes_passing_filters[sample_id],
        detailed_title=False,
        **sample_id_to_thresholds[sample_id]
    )
    plt.savefig(os.path.join(out_dir, f"barcode_stats_{sample_id}.png"))
    plt.close(fig)

# %% # collect garbage
gc.collect()

# %% # load barcodes passing filters
with open(os.path.join(out_dir,"sample_id_to_barcodes_passing_filters.pkl"), "rb") as file:
    sample_id_to_barcodes_passing_filters = pickle.load(file)

# %% # create cistopic object
path_to_regions = os.path.join(out_dir, "consensus_peak_calling/consensus_regions.bed")
path_to_blacklist = os.path.join(in_dir, "hg19-blacklist.v2.bed")
pycistopic_qc_output_dir = os.path.join(out_dir,"qc")

cistopic_obj_list = []
for sample_id in fragments_dict:
    gc.collect()
    sample_metrics = pl.read_parquet(
        os.path.join(pycistopic_qc_output_dir, f'{sample_id}.fragments_stats_per_cb.parquet')
    ).to_pandas().set_index("CB").loc[sample_id_to_barcodes_passing_filters[sample_id]]
    gc.collect()
    cistopic_obj = create_cistopic_object_from_fragments(
        path_to_fragments = fragments_dict[sample_id],
        path_to_regions = path_to_regions,
        path_to_blacklist = path_to_blacklist,
        metrics = sample_metrics,
        valid_bc = sample_id_to_barcodes_passing_filters[sample_id],
        n_cpu = 1,
        project = str(sample_id),
        split_pattern = '-'
    )
    gc.collect()
    cistopic_obj_list.append(cistopic_obj)

# %% # collect garbage
gc.collect()

# %% # create single cistopic object
cistopic_obj = cistopic_obj_list[0]
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj_single.pkl"), "wb")
)

# %% # merge cistopic objects
merged_cistopic_obj = merge(cistopic_obj_list, project="cisTopic_merge", split_pattern="-")

# %% # save merged cistopic object
pickle.dump(
    merged_cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj_merged.pkl"), "wb")
)

# %% # collect garbage
gc.collect()

# %% # adding metadata to a cisTopic object
file_path = os.path.join(out_dir, "cistopic_obj_merged.pkl")

with open(file_path, "rb") as file:
    cistopic_obj = pickle.load(file)

# %% # clean cell data index
cistopic_obj.cell_data.index = cistopic_obj.cell_data.index.str.replace(r'___.*', '', regex=True)

# %% # print cell data head
cistopic_obj.cell_data.head()

# %% # load cells data
# cell_data = pd.read_csv(os.path.join(out_dir, 'cells_data.csv'), index_col = 0)

# %% # print cells data shape
print(cells_data.shape)
cells_data = cells_data[~cells_data.index.duplicated()]
print(cells_data.shape)

# %% # print cistopic object cell data index
print(cistopic_obj.cell_data.index.is_unique)
print(cells_data.index.is_unique)
gc.collect()

# %% # add cell data to cistopic object
cistopic_obj.add_cell_data(cells_data, split_pattern='-')

cistopic_obj.cell_data.head()

# %% # Running scrublet
scrub = scr.Scrublet(cistopic_obj.fragment_matrix.T, expected_doublet_rate=0.1)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram();
scrub.call_doublets(threshold=0.22)
scrub.plot_histogram();
scrublet = pd.DataFrame([scrub.doublet_scores_obs_, scrub.predicted_doublets_], columns=cistopic_obj.cell_names, index=['Doublet_scores_fragments', 'Predicted_doublets_fragments']).T

# %% # add scrublet to cistopic object
cistopic_obj.add_cell_data(scrublet, split_pattern = '-')
sum(cistopic_obj.cell_data.Predicted_doublets_fragments == True)

# %% # remove doublets
singlets = cistopic_obj.cell_data[cistopic_obj.cell_data.Predicted_doublets_fragments == False].index.tolist()
# Subset cisTopic object
cistopic_obj_noDBL = cistopic_obj.subset(singlets, copy=True, split_pattern='-')
print(cistopic_obj_noDBL)
gc.collect()

# %% # save cistopic object
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj.pkl"), "wb")
)
# %%
