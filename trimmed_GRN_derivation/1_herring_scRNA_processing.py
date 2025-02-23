# %% [markdown]

# 1. **Data Loading and Initial Setup**:
#     - Loads single-cell RNA data using Scanpy
#     - Sets up configurations for specific neuron types and age groups
#     - Handles both full counts and downsampled CPM (Counts Per Million) data
# 2. **Data Filtering and Quality Control**:
#     - Filters cells based on chemistry version (v3)
#     - Selects specific cell types and age groups
#     - Applies quality control metrics:
#         - Number of genes per cell
#         - Number of counts per cell
#         - Mitochondrial gene percentage
#     - Removes outliers based on percentile thresholds
# 3. **Data Processing**:
#     - Performs downsampling to 30,000 cells if the dataset is larger
#     - Normalizes the data
#     - Applies log transformation
#     - Preserves raw data and creates different data layers

# 1. input:
#     1. `Processed_data_RNA-all_full-counts-and-downsampled-CPM.h5ad`
#     2. `cells_data.csv`
# 2. output
#     1. `subseted_rna_andata.h5ad`

# # Environment

# %%
import gc
import os
import sys
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt 

from dotenv import load_dotenv
load_dotenv()
sys.path.insert(0, os.getenv('PROJECT_FUNCTIONS_PATH'))

from grn_helpers import set_custom_folders

# %%
root_dir = os.getenv('BASE_PATH')

# %%
neurons_set = "L2-3_CUX2"
# neurons_set = "all_ex"
# neurons_set = "all_ex_all_ages"

# gois = ["AR", "THRB", "ESR2", "NR1H3", "NR1H2", "RARA", "RARG", "AHR", "NR3C1"]
# gois = ['AHR', 'AR', 'NR1I2', 'NR1I3', 'NR3C1', 'NR3C2', 'ESR1', 'RARA', 'ESR2', 'THRB', 'THRA']
gois = ['FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FGFRL1'] # FGF pathway
gois = gois + ['PTCH1', 'SMO', 'GLI1', 'GLI2', 'GLI3', 'GLI4'] # SAG pathway
gois = gois + ['BMPR1A', 'BMPR1B'] # BMP4 pathway
gois = gois + ['ACVR1'] # BMP7 pathway
gois = gois + ['CTNNB1', 'WNT5A', 'WNT3A', 'WNT3', 'APC', 'WNT10B'] # WNT pathway ('WNT1' is missing)
gois = gois + ['RARA', 'RARB', 'RARG', 'RXRA', 'RXRB', 'RXRG'] # Retinoic Acid pathway
print(f"gois: {gois}")
# Available genes: ['FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FGFRL1', 'PTCH1', 'SMO', 'GLI1', 'GLI2', 'GLI3', 'GLI4', 'BMPR1A', 'BMPR1B', 'ACVR1', 'CTNNB1', 'WNT5A', 'WNT3A', 'WNT3', 'APC', 'WNT10B', 'RARA', 'RARB', 'RARG', 'RXRA', 'RXRB', 'RXRG']
# Missing genes: ['WNT1']

min_genes_percentile = 2
max_genes_percentile = 98
min_counts_percentile = 2
max_counts_percentile = 98
max_mito_percent = 30

# %%
cells_dict = {
    "all_ex"            :   ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "all_ex_all_ages"   :   ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "L2-3_CUX2"         :   ['L2-3_CUX2']
}

ages_dict = {
    "all_ex"           :   ['1m','3m','6m','10m','1y','2y','4y','ga22','ga24'],
    "all_ex_all_ages"  :   ['1m','3m','6m','10m','1y','2y','4y','6y','10y','16y','20y','40y','ga22','ga24'],
    "L2-3_CUX2"        :   ['1m','3m','6m','10m','1y','2y','4y','ga22','ga24']
}

output_dir, input_dir, root_dir, _, in_dir_from_scenic = set_custom_folders(root_dir, neurons_set)

sel_celltypes  = cells_dict[neurons_set]
sel_ages = ages_dict[neurons_set]

# %%

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300

# %% [markdown]
# # Load data

# %%
adata = sc.read_h5ad(os.path.join(input_dir, 'Processed_data_RNA-all_full-counts-and-downsampled-CPM.h5ad'))
adata_metadata = pd.read_csv(os.path.join(in_dir_from_scenic, 'cells_data.csv'), sep=",", index_col = 0)

# %%
adata_metadata.head()

# %%
print(f"adata_metadata.age: {adata_metadata.age[:5]}")
print(f"adata_metadata.age_mapped: {adata_metadata.age_mapped[:5]}")

# %%
print(adata.obs.age.head())

# %%
age_mapping = dict(zip(adata_metadata.age, adata_metadata.age_mapped))

# Apply the mapping to adata.obs.age
adata.obs['age_mapped'] = adata.obs.age.map(age_mapping)

# If there are any values in adata.obs.age that don't have a mapping,
# they will become NaN. To keep the original values for these cases, use:
adata.obs['age_mapped'] = adata.obs.age.map(age_mapping).fillna(adata.obs.age)

# Display the first few rows to verify the mapping
print(adata.obs[['age', 'age_mapped']].head())

# how many values were mapped vs. unmapped:
mapped_count = adata.obs.age_mapped.notna().sum()
total_count = len(adata.obs.age)
print(f"Mapped {mapped_count} out of {total_count} values")
print(f"Mapping rate: {mapped_count/total_count:.2%}")

# %%
adata.obs.age_mapped.unique()

# %%
print(adata.shape)
adata = adata[adata.obs['chem']=='v3']
print(adata.shape)

# %%
# sel_celltypes = adata[adata.obs.cell_type == "PN", :].obs.major_clust.unique()
# sel_celltypes

# %%
adata = adata[adata.obs.age_mapped.isin(sel_ages)]
print(adata.shape)

# %%
adata = adata[adata.obs['major_clust'].isin(sel_celltypes)]
print(adata.shape)

# %%
adata.obs.age.unique()

# %%
for cell in adata.obs['major_clust'].unique():
    print(f"{cell}: {(adata.obs['major_clust']==cell).sum()}")

# %% [markdown]
# # Filter data

# %%
fig, axs = plt.subplots(1, 3, figsize=(18, 4.5))
sc.pl.violin(adata, ["n_genes_by_counts", "n_counts", "percent_mito"], 
             jitter=0.4, multi_panel=True, ax=axs, show=False)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'violin_plots_before_filtering.png'))
plt.close()

# %%
print(f"Number of cells before filtering: {adata.n_obs}")

cell_type_to_exclude = 'none'

# percentiles for filtering criteria
min_genes = np.percentile(adata.obs['n_genes_by_counts'], min_genes_percentile)
max_genes = np.percentile(adata.obs['n_genes_by_counts'], max_genes_percentile)
min_counts = np.percentile(adata.obs['n_counts'], min_counts_percentile)
max_counts = np.percentile(adata.obs['n_counts'], max_counts_percentile)
max_mito = max_mito_percent

mask = (adata.obs['n_genes_by_counts'] >= min_genes) & (adata.obs['n_genes_by_counts'] <= max_genes) & \
        (adata.obs['n_counts'] >= min_counts) & (adata.obs['n_counts'] <= max_counts) & \
        (adata.obs['percent_mito'] <= max_mito) | (adata.obs.major_clust==cell_type_to_exclude)

adata = adata[mask, :]

print(f"Number of cells after filtering: {adata.n_obs}")

# %%
fig, axs = plt.subplots(1, 3, figsize=(18, 4.5))
sc.pl.violin(adata, ["n_genes_by_counts", "n_counts", "percent_mito"], 
             jitter=0.4, multi_panel=True, ax=axs, show=False)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'violin_plots_after_filtering.png'))
plt.close()

# %%
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
sc.pl.scatter(adata, x="n_counts", y="percent_mito", ax=ax1, show=False)
sc.pl.scatter(adata, x="n_counts", y="n_genes_by_counts", ax=ax2, show=False)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'scatter_plots.png'))
plt.close()

# %%
gc.collect()

# %% [markdown]
# ## Downsample data

# %%
adata.shape

# %%
print([f"{celltype}: {len(adata.obs.major_clust[adata.obs.major_clust==celltype])}" for celltype in adata.obs.major_clust.unique()])

# %%
n_cells_downsample = 30000
if adata.shape[0] > n_cells_downsample:
    sc.pp.subsample(adata, n_obs=n_cells_downsample, random_state=123)

# %%
print([f"{celltype}: {len(adata.obs.major_clust[adata.obs.major_clust==celltype])}" for celltype in adata.obs.major_clust.unique()])

# %%
types = list(adata.obs.major_clust.unique())
types

# %% [markdown]
# # Normalize data

# %%
adata.raw = adata

# %%
adata.layers['counts'] = adata.X.copy()
sc.pp.filter_genes(adata, min_cells = 1)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

print(adata.layers['counts'][:20, :20])
print(adata.X.data[:10])
print(adata.layers['ds_norm_cts'].data[:10])

# %%
gc.collect()

# %%
sc.pp.highly_variable_genes(adata, inplace = True, n_top_genes=3000) #subset = True
adata

# %%
np.sum(adata.var['highly_variable'])

# %% [markdown]
# # Pseudobulk

# %%
# Group cells by major_clust and calculate mean expression values
counts_data = adata.layers["ds_norm_cts"].toarray()
pseudobulk_df = pd.DataFrame(counts_data, index=adata.obs.index, columns=adata.var_names)
pseudobulk_df["major_clust"] = adata.obs["major_clust"]
pseudobulk_df = pseudobulk_df.groupby("major_clust").mean()

# Create a new AnnData object for pseudobulk data
pseudobulk_adata = sc.AnnData(X=pseudobulk_df.values, obs=pd.DataFrame(index=pseudobulk_df.index), var=pd.DataFrame(index=pseudobulk_df.columns))

# %%
pseudobulk_df.head()

# %%
pseudobulk_adata

# %% [markdown]
# ### scRNAseq - receptors expression

# %%
# Find available genes and their indices
gene_indexs = []
available_gois = []
missing_gois = []

for goi in gois:
    gene_matches = np.where(adata.var_names == goi)[0]
    if len(gene_matches) > 0:
        gene_indexs.append(gene_matches[0])
        available_gois.append(goi)
    else:
        missing_gois.append(goi)

print("\nAvailable genes:", available_gois)
print("Missing genes:", missing_gois)
print(f"Found {len(available_gois)} out of {len(gois)} genes")

if len(available_gois) == 0:
    print("\nNo genes of interest found in the dataset. Please check gene names.")
    sys.exit(1)

# Update subsequent code to use available_gois instead of gois
for i, goi in enumerate(available_gois):
    expression_bulk = pseudobulk_adata.X[:,gene_indexs[i]]
    print([f"Expression of {goi} in {cell}: {expression}" for cell, expression in zip(pseudobulk_adata.obs_names, expression_bulk)])

# %%
expression_data = [pseudobulk_adata.X[:,gene_indexs[i]] for i in range(0,len(available_gois))]
cell_types = pseudobulk_adata.obs_names

bar_width = 0.1
spacing = 0.01
fig, ax = plt.subplots(figsize=(16, 6))

x = np.arange(len(cell_types))
ax.set_xticks(x)
ax.set_xticklabels(cell_types, rotation=45, ha='right', fontsize=15)

for i, goi in enumerate(available_gois):
    offset = (i - (len(available_gois) - 1) / 2) * (bar_width + spacing)
    ax.bar(x + offset, expression_data[i], width=bar_width, label=goi)

ax.set_xlabel('Cell Types', fontsize=16)
ax.set_ylabel('Expression', fontsize=16)
ax.set_title('Gene Expression across Cell Types', fontsize=16)
ax.legend()

plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'gene_expression_across_cell_types.png'))
plt.close()

# %%
summed_counts = np.sum(counts_data, axis=0)
print(summed_counts.shape)
print(summed_counts)

# %%
for i, gene_index in enumerate(gene_indexs[:3]):
    count = summed_counts[gene_index]
                               
    plt.figure(figsize=(8, 6))
    n, bins, patches = plt.hist(summed_counts, bins=50, log=True)
    
    # Find the bin index where the gene count falls into
    bin_index = np.where(bins <= count)[0][-1]
    patches[bin_index].set_facecolor('red')
    
    plt.xlabel("Counts")
    plt.ylabel("Frequency (log scale)")
    plt.title(f"Histogram of Counts: {available_gois[i]}")
    plt.savefig(os.path.join(output_dir, f'histogram_of_counts_{available_gois[i]}.png'))
    plt.close()
    print(f"{available_gois[i]} count: {count}")

# %%
adata.shape

# %%
file_name = os.path.join(output_dir, f"subseted_rna_andata.h5ad")
adata.write(file_name)

# %%
output_dir
