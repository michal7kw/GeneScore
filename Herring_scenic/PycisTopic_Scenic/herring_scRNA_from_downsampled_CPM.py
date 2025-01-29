# %% [markdown]
# in: Herring_data/RNA-all_full-counts-and-downsampled-CPM.h5ad
# 
# out: Herring_data/adata.h5ad

# %% [markdown]
# # 1. Environment

# %%
import os
import sys
import scanpy as sc

sys.path.insert(0, "/home/michal.kubacki/Githubs/Re-MEND/code/External_Datasets/GeneSet_Derivation/Herring_scenic/helpers")
import config
from config import *

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

#################################################################

# %%
file = "Processed_data_RNA-all_full-counts-and-downsampled-CPM.h5ad"
file_path = os.path.join(in_dir, file)

# %%
adata = sc.read_h5ad(file_path)
adata.var_names_make_unique()
adata.obs_names_make_unique()

# %%
adata

# %%
print(adata.obs.cell_type.unique())
print(adata.obs.major_clust.unique())


# %%
adata[adata.obs.cell_type == "PN", :].obs.major_clust.unique()

# %%
adata.obs_names

# %%
adata.X.sum(axis = 1)[:5]

# %%
adata.layers['ds_norm_cts'].sum(axis = 1)[:5]

# %% [markdown]
# # 3. Map ages

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

adata.obs["age_mapped"] = [mapping.get(age, age) for age in adata.obs.age]
adata.obs["age_mapped"].unique()

# %%
adata.obs['old_id'] = adata.obs_names

# %%
def process_obs_name(name):
    parts = name.split('-')[0]
    return parts

adata.obs['sample_id'] = adata.obs.old_id.apply(process_obs_name)

# %%
adata.obs['sample_id'] = adata.obs['sample_id'] + "-1-" + adata.obs["age_mapped"]
adata.obs_names = adata.obs['sample_id']

# %%
adata.obs_names = adata.obs_names.rename(None).astype(str)
adata.obs_names

# %% [markdown]
# # 4. QC

# %%
adata.raw = adata

# %%
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)

# %%
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8, 6))

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=False,
    ax=ax
)

plt.show()

# %%
adata = adata[adata.obs.n_genes_by_counts < 8000, :]
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()

# %%
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# %%
sc.pl.highly_variable_genes(adata)

# %%
adata = adata[:, adata.var.highly_variable]

# %%
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])

# %%
sc.pp.scale(adata, max_value=10)

# %%
adata.var_names_make_unique()
adata.obs_names_make_unique()

# %%
adata.write_h5ad(os.path.join(in_dir,  "adata.h5ad"))

# %%
adata = adata[adata.obs.age_mapped.isin(sel_ages)]

# %%
adata.shape

# %%
adata = adata[adata.obs['major_clust'].isin(sel_celltypes)]

# %%
adata.write_h5ad(os.path.join(out_dir,  "adata.h5ad"))

# %%
adata.shape

# %%
out_dir

# %%
adata_new = sc.read_h5ad(os.path.join(out_dir,  "adata.h5ad"))
adata_new.obs['major_clust']

# %%
list(adata_new.obs.columns)


