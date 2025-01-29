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

# %%
# Standard library imports
import os
import gc
import sys
import pickle
import logging

# Data manipulation imports
import pandas as pd

# pycisTopic imports
import pycisTopic
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, peak_calling
from pycisTopic.iterative_peak_calling import get_consensus_peaks
from pycisTopic.plotting.qc_plot import plot_sample_stats, plot_barcode_stats
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments, merge
from pycisTopic.qc import get_barcodes_passing_qc_for_sample

import matplotlib.pyplot as plt
import scrublet as scr

# %%
n_cpu = 8
neurons_set = "all_ex"
# neurons_set = "all_ex_all_ages"
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

def select_files(selected_fragments = ['2y', '4y']):
    input_dir = os.path.join(root_dir, "data")
    data_folder = "GSE168408_RAW"
    filenames = """
    GSM5138510_RL2366_ga22_snATAC_fragments.tsv GSM5138526_RL2209_1y_snATAC_fragments.tsv GSM5138542_RL2372_14y_snATAC_fragments.tsv GSM5138512_RL2207_ga24_snATAC_fragments.tsv GSM5138529_RL1784_2y_snATAC_fragments.tsv GSM5138544_RL1785_16y_snATAC_fragments.tsv GSM5138515_RL2367_1m_snATAC_fragments.tsv GSM5138532_RL2210_4y_snATAC_fragments.tsv GSM5138548_RL2085_20y_snATAC_fragments.tsv GSM5138518_RL1914_3m_snATAC_fragments.tsv GSM5138534_RL2364_6y_snATAC_fragments.tsv GSM5138550_RL2369_25y_snATAC_fragments.tsv GSM5138521_RL2208_6m_snATAC_fragments.tsv GSM5138536_RL1994_8y_snATAC_fragments.tsv GSM5138552_RL2373_40y_snATAC_fragments.tsv GSM5138523_RL2371_10m_snATAC_fragments.tsv GSM5138539_RL2368_10y_snATAC_fragments.tsv
    """.strip().split()

    fragments_dict = {}
    for filename in filenames:
        sample_id = filename.split('_')[2]
        fragments_dict[sample_id] = os.path.join(input_dir, data_folder, f"{filename}.gz")  

    print(f"All fragments: {fragments_dict}")
    print(f"Selected fragments: {selected_fragments}") 
    fragments_dict = {key: fragments_dict[key] for key in selected_fragments if key in fragments_dict}
    return fragments_dict

# %%
cells_dict = {
    "all_ex"            :   ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "all_ex_all_ages"   :   ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev']
}

ages_dict = {
    "all_ex"            :   ['1m','3m','6m','10m','1y','2y','4y','ga22','ga24'],
    "all_ex_all_ages"   :   ['1m','3m','6m','10m','1y','2y','4y','6y','10y','16y','20y','40y','ga22','ga24']
}

out_dir, in_dir, root_dir, tmp_dir, data_folder = set_output_folders(neurons_set)
sel_celltypes = cells_dict[neurons_set]
sel_ages = ages_dict[neurons_set]

# %%
ATAC_metadata_path = os.path.join(in_dir, "Processed_data_ATAC_BCs-meta-data.csv")

# %%
fragments_dict = select_files(selected_fragments = sel_ages)

# %%
cells_data = pd.read_csv(ATAC_metadata_path, sep=",", index_col=0)

# %%
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

# %%
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

# %% [markdown]
# # Filter ATAC metadata

# %%
print(cells_data.shape)

# %%
cells_data = cells_data[cells_data['major_clust'].isin(sel_celltypes)]
print(cells_data.shape)

# %%
print(cells_data['chem'].value_counts())
print(cells_data['PassQC'].value_counts())

# %%
cells_data = cells_data[cells_data['chem']=="v3"]
print(cells_data.shape)

# %%
cells_data = cells_data[cells_data['age_mapped'].isin(sel_ages)]
print(cells_data.shape)

# %% [markdown]
# # Format indexes

# %%
def format_index(index):
    parts = index.split("#")
    formatted_index = parts[1]
    return formatted_index

cells_data = cells_data.rename(index=format_index)

# %%
cells_data["old_index"] = cells_data.index
cells_data.index = cells_data.index + "-" + cells_data["age_mapped"]
cells_data = cells_data[cells_data.age_mapped.isin(fragments_dict)]

# %%
print(cells_data.shape)

# %%
for cell in cells_data['major_clust'].unique():
    print(f"{cell}: {(cells_data['major_clust']==cell).sum()}")

# %%
# Create a dictionary to store the count of each cell type
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

# %%
for cell in cells_data['major_clust'].unique():
    print(f"{cell}: {(cells_data['major_clust']==cell).sum()}")


# %% [markdown]
# # Getting pseudobulk profiles from cell annotations

# %%
# cells_data = pd.read_csv(os.path.join(out_dir, 'cells_data.csv'), index_col = 0)
cells_data.head()

# %%
chromsizes = pd.read_table(
    os.path.join(in_dir, f"hg19.chrom.sizes"),
    header = None,
    names = ["Chromosome", "End"]
)
chromsizes.insert(1, "Start", 0)

# %% [markdown]
# potentialy beforehand you  might want to run `for file in *fragments.tsv.gz; do tabix -p bed "$file"; done`

# %%
gc.collect()

# %%
os.makedirs(os.path.join(out_dir, "consensus_peak_calling"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bw_files"), exist_ok = True)

paths = export_pseudobulk(
    input_data = cells_data,
    variable = "major_clust",
    sample_id_col = "age_mapped",
    chromsizes = chromsizes,
    bigwig_path = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bw_files"),
    bed_path = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"),
    path_to_fragments = fragments_dict,
    n_cpu = 1,
    temp_dir = '/home/michal.kubacki/Githubs/Re-MEND/code/External_Datasets/GeneSet_Derivation/Herring_scenic/PycisTopic_Scenic/tmp',
    split_pattern = "-"
)


# %%
bw_paths, bed_paths = paths

# %%
with open(os.path.join(out_dir, "consensus_peak_calling/bed_paths.tsv"), "wt") as f:
    for v in bed_paths:
        f.write(f"{v}\t{bed_paths[v]}\n")

# %% [markdown]
# # Inferring consensus peaks

# %%
bed_paths = {}
with open(os.path.join(out_dir, "consensus_peak_calling/bed_paths.tsv")) as f:
    for line in f:
        v, p = line.strip().split("\t")
        bed_paths.update({v: p})

#%%
macs_path = "/home/michal.kubacki/.conda/envs/scenicplus/bin/macs2" 
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

# %%
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

# Update columns in place instead of replacing the entire DataFrame
for col in df.columns:
    consensus_peaks.df[col] = df[col]

del df  # Clean up reference

consensus_peaks.to_bed(
    path=os.path.join(out_dir, "consensus_peak_calling/consensus_regions.bed"),
    keep=True,
    compression='infer',
    chain=False)

# %%
gc.collect()

# %%
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

# %%
import subprocess

with open("pycistopic_qc_commands.txt", "r") as file:
    commands = file.readlines()
    for command in commands:
        process = subprocess.Popen(command.strip().split(" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            print(f"Error executing command: {command.strip()}")
            print(f"Error message: {stderr.decode('utf-8')}")

# %%
gc.collect()

# %%
for sample_id in fragments_dict:
    fig = plot_sample_stats(
        sample_id=sample_id,
        pycistopic_qc_output_dir=os.path.join(out_dir, "qc")
    )
    plt.savefig(os.path.join(out_dir, f"sample_stats_{sample_id}.png"))
    plt.close(fig)

# %%
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

# %%
with open(os.path.join(out_dir, "sample_id_to_barcodes_passing_filters.pkl"), "wb") as file:
    pickle.dump(sample_id_to_barcodes_passing_filters, file)

# %%
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

# %%
gc.collect()

# %%
with open(os.path.join(out_dir,"sample_id_to_barcodes_passing_filters.pkl"), "rb") as file:
    sample_id_to_barcodes_passing_filters = pickle.load(file)

# %%
gc.collect()
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

# %%
gc.collect()

cistopic_obj = cistopic_obj_list[0]
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj_single.pkl"), "wb")
)

merged_cistopic_obj = merge(cistopic_obj_list, project="cisTopic_merge", split_pattern="-")

pickle.dump(
    merged_cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj_merged.pkl"), "wb")
)

# %%
gc.collect()


# %% [markdown]
# # Adding metadata to a cisTopic object

# %%
file_path = os.path.join(out_dir, "cistopic_obj_merged.pkl")

with open(file_path, "rb") as file:
    cistopic_obj = pickle.load(file)

# %%
cistopic_obj.cell_data.index = cistopic_obj.cell_data.index.str.replace(r'___.*', '', regex=True)

# %%
cistopic_obj.cell_data.head()

# %%
# cell_data = pd.read_csv(os.path.join(out_dir, 'cells_data.csv'), index_col = 0)

# %%
print(cells_data.shape)
cells_data = cells_data[~cells_data.index.duplicated()]
print(cells_data.shape)

# %%
print(cistopic_obj.cell_data.index.is_unique)
print(cells_data.index.is_unique)
gc.collect()

# %%
cistopic_obj.add_cell_data(cells_data, split_pattern='-')

# %%
cistopic_obj.cell_data.head()

# %% [markdown]
# # Running scrublet

# %%
scrub = scr.Scrublet(cistopic_obj.fragment_matrix.T, expected_doublet_rate=0.1)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram();
scrub.call_doublets(threshold=0.22)
scrub.plot_histogram();
scrublet = pd.DataFrame([scrub.doublet_scores_obs_, scrub.predicted_doublets_], columns=cistopic_obj.cell_names, index=['Doublet_scores_fragments', 'Predicted_doublets_fragments']).T

# %%
cistopic_obj.add_cell_data(scrublet, split_pattern = '-')
sum(cistopic_obj.cell_data.Predicted_doublets_fragments == True)

# %%
# Remove doublets
singlets = cistopic_obj.cell_data[cistopic_obj.cell_data.Predicted_doublets_fragments == False].index.tolist()
# Subset cisTopic object
cistopic_obj_noDBL = cistopic_obj.subset(singlets, copy=True, split_pattern='-')
print(cistopic_obj_noDBL)
gc.collect()

# %%
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj.pkl"), "wb")
)
