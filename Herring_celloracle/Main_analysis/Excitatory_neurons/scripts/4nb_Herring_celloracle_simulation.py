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
importlib.reload(co)
from celloracle import *

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
# # Load the data

# %%
oracle_file_name = os.path.join(output_dir, "oracle_1.celloracle.oracle")
links_file_name = os.path.join(output_dir, "links_1.celloracle.links")

# %%
oracle = co.load_hdf5(oracle_file_name)

# %%
links = co.load_hdf5(links_file_name)

# %%
links

# %%
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)

# %%
gc.collect()

# %% [markdown]
# # Fit GRNs

# %%
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)

# %%
# current_datetime = datetime.now().strftime("%Y%m%d_%H%M%S")
# current_datetime = ""
# file_name = os.path.join(output_dir, f"Herring_simulation.celloracle.oracle")

# %%
# oracle.to_hdf5(file_name)

# %%
# current_datetime = ""
# file_name = f'{folder}Herring_simulation.celloracle.oracle'
# oracle = co.load_hdf5(file_name)

# %%
oracle.adata.obs.major_clust.unique()

# %%
goi = "THRB"

sc.pl.umap(oracle.adata, color = [goi], layer="imputed_count")

# %%
sc.pl.umap(oracle.adata, color=[goi, oracle.cluster_column_name],
                 layer="imputed_count", use_raw=False, cmap="viridis")

# %%
gc.collect()

# %% [markdown]
# # Run simulation

# %%
goi = "THRB"

print(oracle.all_regulatory_genes_in_TFdict[:10])
print(len(oracle.all_regulatory_genes_in_TFdict))
print(goi in oracle.all_regulatory_genes_in_TFdict)

# %%
# Enter perturbation conditions to simulate signal propagation after the perturbation.
oracle.simulate_shift(perturb_condition={goi: 0.0},
                      n_propagation=3)

# %%
gc.collect()

# %%
# Get transition probability
oracle.estimate_transition_prob(n_neighbors=200,
                                knn_random=True,
                                sampled_fraction=1)

# %%
gc.collect()

# %%
# Calculate embedding
oracle.calculate_embedding_shift(sigma_corr=0.05)

# %%
# file_name = os.path.join(output_dir, f"Herring_simulation_{goi}.celloracle.oracle")
# oracle.to_hdf5(file_name)

# %%
fig, ax = plt.subplots(1, 2,  figsize=[13, 6])

scale = 50
# Show quiver plot
oracle.plot_quiver(scale=scale, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: {goi} KO")

# Show quiver plot that was calculated with randomized graph.
oracle.plot_quiver_random(scale=scale, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")

plt.show()

# %%
n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)

# %%
oracle.suggest_mass_thresholds(n_suggestion=12)

# %%
min_mass = 60
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

# %%
fig, ax = plt.subplots(1, 2,  figsize=[13, 6])

scale_simulation = 10
# Show quiver plot
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: {goi} KO")

# Show quiver plot that was calculated with randomized graph.
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")

plt.show()

# %%
# Plot vector field with cell cluster
fig, ax = plt.subplots(figsize=[8, 8])

oracle.plot_cluster_whole(ax=ax, s=5)
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax, show_background=False)

# %%
def plot_arrows_legend(oracle, labels=None, colormap='viridis', scale=1, data_random=False, points_size=5, filename=None):
    fig, ax = plt.subplots(figsize=[8, 8])

    embedding = oracle.adata.obsm['X_umap']
    cluster_labels = oracle.adata.obs[labels]
    cluster_categories = pd.Categorical(cluster_labels)
    cmap = plt.cm.get_cmap(colormap, len(cluster_categories.categories))

    scatter = ax.scatter(embedding[:, 0], embedding[:, 1], c=cluster_categories.codes, cmap=cmap, s=points_size)

    # Arrow selection
    if data_random:
        flow = oracle.flow_rndm
    else:
        flow = oracle.flow

    if hasattr(oracle, "mass_filter"):
        mass_filter = oracle.mass_filter
        gridpoints_coordinates = oracle.flow_grid
    else:
        mass_filter = np.zeros(flow.shape[0], dtype=bool)
        gridpoints_coordinates = embedding

    ax.quiver(gridpoints_coordinates[~mass_filter, 0],
              gridpoints_coordinates[~mass_filter, 1],
              flow[~mass_filter, 0],
              flow[~mass_filter, 1],
              scale=scale)

    ax.axis("off")

    if labels is not None:
        # Create legend elements based on the cluster categories
        legend_elements = [plt.Line2D([0], [0], marker='o', color='w', label=str(label),
                                      markerfacecolor=cmap(i), markersize=10)
                           for i, label in enumerate(cluster_categories.categories)]
        ax.legend(handles=legend_elements, loc='best')

    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
    else:
        plt.show()

# %%
plot_arrows_legend(oracle, labels='major_clust', scale=10, data_random=False, points_size = 10, filename=os.path.join(output_dir, "plot_simulation_flow_on_grid.png"))

# %%
#file_name = os.path.join(output_dir, f"simulation_{goi}.celloracle.oracle")

# %%
#Save checkpoint
#oracle.to_hdf5(file_name)

# %% [markdown]
# # Evaluate perturbation results

# %%
# file_name = os.path.join(output_dir, f"Herring_simulation_{goi}.celloracle.oracle")
# oracle = co.load_hdf5(file_name)

# %%
simulated_count = oracle.adata.layers["simulated_count"]

# %%
original_count = oracle.adata.X.toarray()
log_fold_change = np.log2(simulated_count + 1) - np.log2(original_count + 1)

# %%
print(log_fold_change.shape)
log_fold_change

# %%
top_n = 1000

# %%
cell_types = oracle.adata.obs['major_clust']
gene_names = oracle.adata.var_names

unique_cell_types = cell_types.unique()

table_data_list = []

for cell_type in unique_cell_types:
    # Get the indices of cells belonging to the current cell type
    cell_type_indices = np.where(cell_types == cell_type)[0]
    
    # Calculate the mean log fold change for each gene in the current cell type
    cell_type_log_fold_change = log_fold_change[cell_type_indices, :].mean(axis=0)
    
    # Calculate the absolute values of the mean log fold changes
    abs_cell_type_log_fold_change = np.abs(cell_type_log_fold_change)
    
    # Sort the absolute mean log fold changes in descending order
    sorted_indices = np.argsort(abs_cell_type_log_fold_change)[::-1]
    
    # Select the top_n genes based on the sorted indices
    top_gene_indices = sorted_indices[:top_n]
    top_genes = gene_names[top_gene_indices]
    
    top_log_fold_changes = cell_type_log_fold_change[top_gene_indices]
    
    cell_type_table_data = pd.DataFrame({
        'source': [cell_type] * top_n,
        'target': top_genes,
        'log_fold_change': top_log_fold_changes,
    })
    
    table_data_list.append(cell_type_table_data)

final_table_data_sim = pd.concat(table_data_list, ignore_index=True)

# %%
final_table_data_sim['fold_change'] = np.exp2(final_table_data_sim['log_fold_change'])

final_table_data_sim.head()

# %%
# final_table_data_sim.to_csv(os.path.join(output_dir, f'final_table_data_sim_{goi}.csv'), index=False)

# %%
# final_table_data_sim = pd.read_csv(os.path.join(output_dir, f'final_table_data_sim_{goi}.csv'))
# final_table_data_sim

# %%
print(list(final_table_data_sim.target))

# %%
cell_types = oracle.adata.obs['major_clust']
unique_cell_types = cell_types.unique()

tops = []

for cell_type in unique_cell_types:
    tops_df = final_table_data_sim[final_table_data_sim['source'] == cell_type][:5]
    tops.append(tops_df)

tops_df = pd.concat(tops, ignore_index=True)
tops_df

# %% [markdown]
# # Get scores from `links`

# %%
links_file_name = os.path.join(output_dir, "PN_dev.celloracle.links")

# %%
links = co.load_hdf5(links_file_name)

# %%
links

# %%
links.filtered_links

# %%
goi = "THRB"

# %%
links.plot_score_per_cluster(goi=goi)

# %%
all_table_data = []

# Iterate over all cell types in the links.filtered_links dictionary
for celltype in links.filtered_links:
    # Get the GRN data for the current cell type
    grn_data = links.filtered_links[celltype]
    
    grn_data = grn_data[grn_data["source"] == goi]
    
    # Calculate the score for each row
    grn_data["score"] = -np.log10(grn_data["p"])
    
    grn_data["celltype"] = celltype
    
    grn_data = grn_data.rename(columns={"-logp": "X.logp"})
    
    table_data = grn_data[["source", "target", "coef_mean", "coef_abs", "p", "X.logp", "score", "celltype"]]
    
    all_table_data.append(table_data)

# Concatenate the table data from all cell types into a single DataFrame
final_table_data = pd.concat(all_table_data, ignore_index=True)

final_table_data.head()

# %%
# final_table_data.to_csv(os.path.join(output_dir, f'final_table_data_GRN_{goi}.csv'), index=False)

# %%
# final_table_data = pd.read_csv(f'{folder}final_table_data_GRN.csv')

# %%
# links.filtered_links['L2-3_CUX2'].source


