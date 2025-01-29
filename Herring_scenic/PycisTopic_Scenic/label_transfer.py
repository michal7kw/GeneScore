# %% [markdown]
# # Label transfer
# %% 
from pycisTopic.label_transfer import label_transfer
from pycisTopic.clust_vis import plot_metadata
import seaborn as sns
import pickle
import os
import sys
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, "/home/michal.kubacki/Githubs/Re-MEND/code/External_Datasets/GeneSet_Derivation/Herring_scenic/helpers")
import config
from config import *

n_cpu = 32

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
cistopic_obj = pickle.load(open(os.path.join(out_dir, "cistopic_obj.pkl"), "rb"))

# %%
# Load gene_act from the file
with open(os.path.join(out_dir, 'gene_act.pkl'), 'rb') as file:
    gene_act = pickle.load(file)

# %%
rna_anndata = sc.read_h5ad(
    os.path.join(out_dir,  "modified_adata.h5ad")
).raw.to_adata()
atac_anndata = sc.AnnData(gene_act.mtx.T, obs = pd.DataFrame(index = gene_act.cell_names), var = pd.DataFrame(index = gene_act.feature_names))

# %%
atac_anndata.obs["sample_id"] = "multiome_brain"
rna_anndata.obs["sample_id"] = "multiome_brain"

# %%
label_dict = label_transfer(
    rna_anndata,
    atac_anndata,
    labels_to_transfer = ['major_clust'],
    variable_genes = True,
    methods = ['scanorama'],
    return_label_weights = False,
    _temp_dir= tmp_dir,
    n_cpu = 1
)

# %%
label_dict_x=[label_dict[key] for key in label_dict.keys()]
label_pd = pd.concat(label_dict_x, axis=1, sort=False)
label_pd.index = cistopic_obj.cell_names
label_pd.columns = ['pycisTopic_' + x for x in label_pd.columns]
cistopic_obj.add_cell_data(label_pd, split_pattern = '-')

# %%
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj_labeled.pkl"), "wb")
)

# %%
plot_metadata(
    cistopic_obj,
    reduction_name='UMAP',
    variables=['major_clust'] + label_pd.columns.to_list(),
    remove_nan=True,
    seed=555,
    num_columns=3)

# %%
fig, axs = plt.subplots(ncols = 3, nrows = 2, figsize = (3 * 5, 2 * 5))
for method, ax in zip(label_pd.columns.to_list(), axs.ravel()):
    conf_mat = pd.crosstab(cistopic_obj.cell_data["major_clust"], cistopic_obj.cell_data[method])
    conf_mat = conf_mat / conf_mat.sum()
    sns.heatmap(conf_mat.loc[conf_mat.columns], ax = ax)
fig.tight_layout()
fig.show()
# %%
