# %% [markdown]
# # Environment

# %%
import gc 
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from datetime import datetime
import importlib

import celloracle as co
co.__version__

sys.path.insert(0, "/home/michal.kubacki/Githubs/Re-MEND/code/External_Datasets/GeneSet_Derivation/Herring_celloracle/helpers")


import config
importlib.reload(config)
from config import *
n_cpus = 32

# %%
reference = "hg19"

neurons_set = "all_excitatory_all_ages"
# neurons_set = "all_excitatory"

cells_dict = {
    "all_excitatory"            :   ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "all_excitatory_all_ages"   :   ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev']
}

ages_dict = {
    "all_excitatory"            :   ['1m','3m','6m','10m','1y','2y','4y','ga22','ga24'],
    "all_excitatory_all_ages"   :   ['1m','3m','6m','10m','1y','2y','4y','6y','10y','16y','20y','40y','ga22','ga24']
}

output_dir, input_dir, root_dir, tmp_dir, in_dir_from_scenic = set_custom_folders(reference, neurons_set)

sel_celltypes  = cells_dict[neurons_set]
sel_ages = ages_dict[neurons_set]

# %%
%config InlineBackend.figure_format = 'retina'
%matplotlib inline

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300

# %% [markdown]
# # Load data

# %% [markdown]
# ## Load scRNA-seq data

# %%
adata = sc.read_h5ad(os.path.join(output_dir, 'subseted_rna_andata.h5ad'))
print(adata)
print(adata.var.shape)
print(adata.obs.shape)

# %%
print([f"{celltype}: {len(adata.obs.major_clust[adata.obs.major_clust==celltype])}" for celltype in adata.obs.major_clust.unique()])

# %% [markdown]
# ### Add genes of interests

# %%
toadd = ['RARA', 'ESR2', 'THRB']
print(len(toadd))

hvgs = list(adata.var_names[adata.var['highly_variable']])
print(len(hvgs))

# %%
diff = list(set(toadd).difference(set(hvgs)))
intersec = list(set(toadd).intersection(set(hvgs)))
print(len(diff))
print(len(intersec))

# %%
values_to_remove = list(set(toadd).difference(set(adata.var_names)))
print(len(values_to_remove))

toadd = [item for item in toadd if item not in values_to_remove]

# %%
hvgs.extend(toadd)
hvgs = pd.Series(hvgs).unique()
print(len(hvgs))


# %%
# adata = adata[:, adata.var['highly_variable']]
adata = adata[:, hvgs]
adata

# %%
gc.collect()

# %% [markdown]
# # CellOracle

# %% [markdown]
# ## Init CellOracle object

# %%
oracle = co.Oracle()

# %%
oracle.import_anndata_as_raw_count(adata,
                                   cluster_column_name="major_clust",
                                   embedding_name="X_umap")

# %%
oracle

# %%
# base_GRN = co.data.load_human_promoter_base_GRN() 
# base_GRN = pd.read_parquet("./data/2023_11_tfi.celloracle.parquet", engine='pyarrow')
# base_GRN = pd.read_parquet(os.path.join(output_dir, "Herring_motif_scan.celloracle.parquet"), engine='pyarrow')

# %%
cell_type = sel_celltypes[0]

file_path = os.path.join(output_dir, f"{cell_type}.celloracle.parquet")
base_GRN = pd.read_parquet(file_path, engine='pyarrow')


base_GRN.head()

# %%
oracle.import_TF_data(TF_info_matrix=base_GRN)

# %%
External_data = pd.read_csv(os.path.join(input_dir, "2023_11_CellOracleProof.tsv"),delimiter="\t")
External_data

# %%
TF_to_TG_dictionary = {}

for TF, TGs in zip(External_data.TF, External_data.Target_genes):
    # convert target gene to list
    TG_list = TGs.replace(" ", "").split(",")
    # store target gene list in a dictionary
    TF_to_TG_dictionary[TF] = TG_list

# We invert the dictionary above using a utility function in celloracle.
TG_to_TF_dictionary = co.utility.inverse_dictionary(TF_to_TG_dictionary)

oracle.addTFinfo_dictionary(TG_to_TF_dictionary)

# %%
oracle.perform_PCA()

plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
plt.axvline(n_comps, c="k")
plt.show()
print(n_comps)
n_comps = min(n_comps, 50)

# %%
gc.collect()

# %%
n_cell = oracle.adata.shape[0]
k = int(0.025*n_cell)
print(f"Auto-selected k is :{k}")
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=n_cpus)

# %%
sc.pp.neighbors(oracle.adata)

# %%
sc.tl.umap(oracle.adata) #, min_dist=0.3

# %%
sc.pl.umap(oracle.adata, color = "major_clust", layer="imputed_count")

# %%
# current_datetime = datetime.now().strftime("%Y%m%d_%H%M%S")
# file_name = os.path.join(output_dir, f"oracle_1.celloracle.oracle")
# oracle.to_hdf5(file_name)

# %%
gc.collect()

# %% [markdown]
# ## Inferring GRN

# %%
links = oracle.get_links(cluster_name_for_GRN_unit="major_clust", alpha=10,
                         verbose_level=10, n_jobs=n_cpus)

# %%
links.links_dict.keys()

# %%
# links.filter_links()
links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)


# %%
links.get_network_score()

# %%
links.merged_score.head()

# %%
# current_datetime = datetime.now().strftime("%Y%m%d_%H%M%S")
# file_name = os.path.join(output_dir, f"links_1.celloracle.links")
# links.to_hdf5(file_path=file_name)

# %%
# file_name = os.path.join(output_dir, "Herring_links_20240404_160414.celloracle.links")
# links = co.load_hdf5(file_name)

# %%
plt.rcParams["figure.figsize"] = [9, 4.5]

# %%
links.plot_degree_distributions(plot_model=True)

# %%
plt.rcParams["figure.figsize"] = [6, 4.5]

# %%
links.cluster

# %%
# Compare GRN score between two clusters
links.plot_score_comparison_2D(value="eigenvector_centrality",
                               cluster1="L2-3_CUX2", cluster2="L4_RORB",
                               percentile=98)

# %%
# Compare GRN score between two clusters
links.plot_score_comparison_2D(value="betweenness_centrality",
                               cluster1="L2-3_CUX2", cluster2="L4_RORB",
                               percentile=98)



