# %%
import mudata
import anndata as ad
import pickle
import pandas as pd
import sys
import os
sys.path.insert(0, "/home/michal.kubacki/Githubs/Re-MEND/code/External_Datasets/GeneSet_Derivation/Herring_scenic/helpers")
import config
from config import *
n_cpu = 32

# %%
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

# %%
#################################################################
reference = "hg19"

neurons_set = "all_cells"
# neurons_set = "all_excitatory"
# neurons_set = "all_inhibitory"
# neurons_set = "all_excitatory_all_ages"
# neurons_set = "all_inhibitory_all_ages"

cells_dict = {
    "all_inhibitory"            :   ['SST', 'VIP', 'MGE_dev'],
    "all_inhibitory_all_ages"   :   ['VIP', 'SST', 'PV', 'MGE_dev'],
    "all_excitatory"            :   ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "all_excitatory_all_ages"   :   ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "all_cells"                 :   ['L4_RORB', 'L2-3_CUX2', 'SST', 'Astro', 'L5-6_TLE4', 'L5-6_THEMIS', 'VIP', 'OPC', 'PV', 'Oligo', 'Micro', 'MGE_dev', 'PN_dev']
}

ages_dict = {
    "all_inhibitory"            :   ['1m','3m','6m','10m','1y','2y','4y','ga22','ga24'],
    "all_inhibitory_all_ages"   :   ['1m','3m','6m','10m','1y','2y','4y','6y','10y','16y','20y','40y','ga22','ga24'],
    "all_excitatory"            :   ['1m','3m','6m','10m','1y','2y','4y','ga22','ga24'],
    "all_excitatory_all_ages"   :   ['1m','3m','6m','10m','1y','2y','4y','6y','10y','16y','20y','40y','ga22','ga24'],
    "all_cells"                 :   ['1m','3m','6m','10m','1y','2y','4y','6y','10y','16y','20y','40y','ga22','ga24']
}


out_dir, in_dir, root_dir, tmp_dir, data_folder = set_output_folders(reference, neurons_set)

sel_celltypes  = cells_dict[neurons_set]
sel_ages = ages_dict[neurons_set]

print(f"out_dir: {out_dir}")
#################################################################

# %%
############ Load data
ATAC_metadata_path = os.path.join(in_dir, "Processed_data_ATAC_BCs-meta-data.csv")

# %%
fragments_dict = select_files(reference, selected_fragments = sel_ages)

# %%
pd.set_option('display.max_columns', None)
cells_data = pd.read_csv(ATAC_metadata_path, sep=",", index_col = 0)
cells_data.head()


# %%
############ Add new columns

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
############ Map ages
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

# %%
print(cells_data.shape)

# %%
############ Filter data
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

# %%
############ Format data
def format_index(index):
    parts = index.split("#")
    formatted_index = parts[1]
    return formatted_index

cells_data = cells_data.rename(index=format_index)

# %%
cells_data["old_index"] = cells_data.index
cells_data.index = cells_data.index + "-" + cells_data["age_mapped"]

# %%
cells_data.head()

# %%
cells_data = cells_data[cells_data.age_mapped.isin(fragments_dict)]
cells_data.head()

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

# Create a mask to filter cell types with too small counts number
mask = cells_data['major_clust'].isin([cell for cell, count in cell_counts.items() if count >= 800])

# Subset cells_data based on the mask
cells_data = cells_data[mask]

# %%
for cell in cells_data['major_clust'].unique():
    print(f"{cell}: {(cells_data['major_clust']==cell).sum()}")

# %%
cells_data['id'] = cells_data.index 

# %%

# Function to convert predictedCell to match index format
def convert_predicted_cell(predicted_cell, index):
    cell_part = predicted_cell.split('-')[0]
    index_suffix = index.split('-')[2]
    return f"{cell_part}-1-{index_suffix}"

# Convert columns to string
cells_data['id'] = cells_data['id'].astype(str)
cells_data['predictedCell'] = cells_data['predictedCell'].astype(str)

# Apply function to the column using pandas apply
cells_data['predictedCellFormatted'] = cells_data.apply(
    lambda row: convert_predicted_cell(row['predictedCell'], row.name), 
    axis=1
)

# %%
# formatted_to_original = dict(zip(cells_data['predictedCellFormatted'], cells_data.index))
# original_index = cells_data.index

# cells_data = cells_data.set_index('predictedCellFormatted')
# cells_data = cells_data.rename_axis(None)

# %%
############ Save cells data

cells_data.to_csv(os.path.join(out_dir, 'cells_data.csv'), index=True)


# %%
############ Inspect PycisTopuc model

# mdata = mudata.read("/group/testa/michal.kubacki/herring/output_hg19_all_excitatory/ACC_GEX.h5mu")
# print(mdata["scRNA"].shape)

# %%
# mdata

# %%
# mdata.obs.head()

# %% [markdown]
# ## adata

# %%
############ Inspect adata
adata = ad.read_h5ad(os.path.join(out_dir, "adata.h5ad"))

# %%
adata

# %%
adata.obs.head()

# %%
############ Inspect cistopic_obj
file_path = os.path.join(out_dir, "cistopic_obj.pkl")

with open(file_path, "rb") as f:
    cistopic_obj = pickle.load(f)

print(type(cistopic_obj))

# %%
cistopic_obj.cell_data.head()


# %%
############ Inspect cell_data
cell_data = pd.read_csv(os.path.join(out_dir, "cells_data.csv"), index_col = 0)

# %%
cell_data.head()

# %%
############ Integrate data
cell_data.columns

# %%
cell_data['predictedCell'].head()

# %%
print("Data type of 'predictedCell' column:", cell_data['predictedCell'].dtype)

# find and print non-string elements
non_string_elements = cell_data[cell_data['predictedCell'].apply(lambda x: not isinstance(x, str))]

if non_string_elements.empty:
    print("All elements in 'predictedCell' are strings.")
else:
    print("Elements in 'predictedCell' that are not strings:")
    print(non_string_elements['predictedCell'])
    
    # Print additional information about these elements
    print("\nDetailed information about non-string elements:")
    for index, value in non_string_elements['predictedCell'].items():
        print(f"Index: {index}, Value: {value}, Type: {type(value)}")

# Print the number of non-string elements
print(f"\nTotal number of non-string elements: {len(non_string_elements)}")

# If there are non-string elements, print a sample of string elements for comparison
if not non_string_elements.empty:
    print("\nSample of string elements for comparison:")
    string_elements = cell_data[cell_data['predictedCell'].apply(lambda x: isinstance(x, str))]['predictedCell'].head()
    for index, value in string_elements.items():
        print(f"Index: {index}, Value: {value}, Type: {type(value)}")

# %%
cell_data['id'] = cell_data.index 

# %%
# Function to convert predictedCell to match index format
def convert_predicted_cell(predicted_cell, index):
    cell_part = predicted_cell.split('-')[0]
    index_suffix = index.split('-')[2]
    return f"{cell_part}-1-{index_suffix}"

# Convert columns to string
cell_data['id'] = cell_data['id'].astype(str)
cell_data['predictedCell'] = cell_data['predictedCell'].astype(str)

# Apply function to the column
cell_data['predictedCellFormatted'] = cell_data.apply(
    lambda row: convert_predicted_cell(row['predictedCell'], row.name), 
    axis=1
)

# %%
cell_data.predictedCellFormatted.head()

# %%
cell_data.id.head()

# %%
cell_data.age_mapped.head()

# %%
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

# %%
type(cistopic_obj.cell_data.index)

# %%
cistopic_obj.cell_data.index.intersection(adata.obs.index)

# %%
cistopic_obj.cell_data.index.intersection(adata.obs.sample_id)

# %%
matching_cells = set(cell_data.predictedCellFormatted).intersection(set(adata.obs.sample_id))

# Calculate the proportion
total_cells = len(adata.obs.index)
matching_proportion = len(matching_cells) / total_cells

# Print the results
print(f"Number of matching cells: {len(matching_cells)}")
print(f"Total number of cells in adata: {total_cells}")
print(f"Proportion of matching cells: {matching_proportion:.2%}")

# %%
# Get the intersection of the two sets
matching_cells = set(cell_data.predictedCellFormatted).intersection(set(adata.obs.index))

# Calculate the proportion
total_cells = len(adata.obs.index)
matching_proportion = len(matching_cells) / total_cells

print(f"Number of matching cells: {len(matching_cells)}")
print(f"Total number of cells in adata: {total_cells}")
print(f"Proportion of matching cells: {matching_proportion:.2%}")

# %%
# Create a dictionary mapping predictedCellFormatted to index in cell_data
cell_data_dict = dict(zip(cell_data.predictedCellFormatted, cell_data.index))

# Function to replace index if it exists in cell_data_dict
def replace_index(idx):
    return cell_data_dict.get(idx, idx)

# Create a copy of the original AnnData object
adata_new = adata.copy()

# Apply the replacement function to the index
adata_new.obs.index = adata_new.obs.index.map(replace_index)

adata_new.var_names_make_unique()
adata_new.obs_names_make_unique()

# Verify the changes
total_cells = len(adata_new.obs.index)
changed_cells = sum(adata_new.obs.index != adata.obs.index)

print(f"Total number of cells: {total_cells}")
print(f"Number of cells with changed index: {changed_cells}")
print(f"Proportion of cells changed: {changed_cells / total_cells:.2%}")

############ Save integrated data
adata_new.write_h5ad(os.path.join(out_dir, "modified_adata.h5ad"))