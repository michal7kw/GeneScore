# %% [markdown]
# 1. Data Exploration
# 2. Importing and Formating snATAC data
# 3. Getting pseudobulk profiles from cell annotations
# 4. Infering consensus peaks

# %%
# Standard library imports
import os
import gc
import sys
import pickle
import logging
import importlib

# Data manipulation imports
import pandas as pd
import scanpy as sc
import scrublet as scr
import polars as pl

# Visualization imports
import matplotlib.pyplot as plt
import seaborn as sns

# pycisTopic imports
import pycisTopic
from pycisTopic.lda_models import run_cgs_models_mallet, evaluate_models
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, peak_calling
from pycisTopic.iterative_peak_calling import get_consensus_peaks
from pycisTopic.plotting.qc_plot import plot_sample_stats, plot_barcode_stats
from pycisTopic.topic_qc import compute_topic_metrics, plot_topic_qc, topic_annotation
from pycisTopic.utils import fig2img
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments, merge
from pycisTopic.qc import get_barcodes_passing_qc_for_sample
from pycisTopic.clust_vis import (
    find_clusters,
    run_umap,
    run_tsne,
    plot_metadata,
    plot_topic,
    cell_topic_heatmap
)

importlib.reload(pycisTopic)
from pycisTopic import *
pycisTopic.__version__

sys.path.insert(0, "/home/michal.kubacki/Githubs/Re-MEND/code/External_Datasets/GeneSet_Derivation/Herring_scenic/helpers")
import config
importlib.reload(config)
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
# %% [markdown]
# # Model selection

# %% [markdown]
# # Clustering and visualization

# %%
cistopic_obj = pickle.load(open(os.path.join(out_dir, "cistopic_obj.pkl"), "rb"))

# %%
find_clusters(
    cistopic_obj,
    target  = 'cell',
    k = 10,
    res = [0.6, 1.2, 3],
    prefix = 'pycisTopic_',
    scale = True,
    split_pattern = '-'
)

# %%
run_umap(
    cistopic_obj,
    target  = 'cell', scale=True)

# %%
run_tsne(
    cistopic_obj,
    target  = 'cell', scale=True)

# %%
cistopic_obj.cell_data.columns

# %%
# plot_metadata(
#     cistopic_obj,
#     reduction_name='UMAP',
#     variables=['pycisTopic_leiden_10_0.6', 'pycisTopic_leiden_10_1.2', 'pycisTopic_leiden_10_3'],
#     target='cell', num_columns=4,
#     text_size=10,
#     dot_size=20)

# %%
annot_dict = {}
for resolution in [0.6, 1.2, 3]:
    annot_dict[f"pycisTopic_leiden_10_{resolution}"] = {}
    for cluster in set(cistopic_obj.cell_data[f"pycisTopic_leiden_10_{resolution}"]):
        counts = cistopic_obj.cell_data.loc[
            cistopic_obj.cell_data.loc[cistopic_obj.cell_data[f"pycisTopic_leiden_10_{resolution}"] == cluster].index, "major_clust"].value_counts()
        if not counts.empty:
            annot_dict[f"pycisTopic_leiden_10_{resolution}"][cluster] = f"{counts.index[counts.argmax()]}({cluster})"
        else:
            annot_dict[f"pycisTopic_leiden_10_{resolution}"][cluster] = f"N/A({cluster})"

# %%
annot_dict

# %%
for resolution in [0.6, 1.2, 3]:
    cistopic_obj.cell_data[f'pycisTopic_leiden_10_{resolution}'] = [
        annot_dict[f'pycisTopic_leiden_10_{resolution}'][x] for x in cistopic_obj.cell_data[f'pycisTopic_leiden_10_{resolution}'].tolist()
    ]

# %%
# plot_metadata(
#     cistopic_obj,
#     reduction_name='UMAP',
#     variables=['major_clust', 'pycisTopic_leiden_10_0.6', 'pycisTopic_leiden_10_1.2', 'pycisTopic_leiden_10_3'],
#     target='cell', num_columns=4,
#     text_size=10,
#     dot_size=5)

# %%
# plot_metadata(
#     cistopic_obj,
#     reduction_name='UMAP',
#     variables=['log10_unique_fragments_count', 'tss_enrichment', 'Doublet_scores_fragments', 'fraction_of_fragments_in_peaks'],
#     target='cell', num_columns=4,
#     text_size=10,
#     dot_size=5)

# %%
# plot_topic(
#     cistopic_obj,
#     reduction_name = 'UMAP',
#     target = 'cell',
#     num_columns=5
# )

# %%

color_dict = {'major_clust': {}}
unique_values = cistopic_obj.cell_data['major_clust'].unique()
colors = sns.color_palette('husl', len(unique_values)).as_hex()
for value, color in zip(unique_values, colors):
    color_dict['major_clust'][value] = color

cell_topic_heatmap(
    cistopic_obj,
    variables=['major_clust'],
    scale=False,
    legend_loc_x=1.0,
    legend_loc_y=-1.2,
    legend_dist_y=-1,
    figsize=(10, 10),
    color_dictionary=color_dict
)

# %% [markdown]
# # Topic binarization & QC

# %%
region_bin_topics_top_3k = binarize_topics(
    cistopic_obj, method='ntop', ntop = 3_000,
    plot=False, num_columns=5
)

# %%
region_bin_topics_otsu = binarize_topics(
    cistopic_obj, method='otsu',
    plot=False, num_columns=5
)

# %%
binarized_cell_topic = binarize_topics(
    cistopic_obj,
    target='cell',
    method='li',
    plot=False,
    num_columns=5, nbins=100)

# %%
topic_qc_metrics = compute_topic_metrics(cistopic_obj)

# %%
fig_dict={}
fig_dict['CoherenceVSAssignments']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Log10_Assignments', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['AssignmentsVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Log10_Assignments', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSRegions_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Regions_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSMarginal_dist']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Marginal_topic_dist', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSGini_index']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Gini_index', var_color='Gini_index', plot=False, return_fig=True)
gc.collect()

# %%
# Plot topic stats in one figure
fig=plt.figure(figsize=(40, 43))
i = 1
for fig_ in fig_dict.keys():
    plt.subplot(2, 3, i)
    img = fig2img(fig_dict[fig_]) #To convert figures to png to plot together, see .utils.py. This converts the figure to png.
    plt.imshow(img)
    plt.axis('off')
    i += 1
plt.subplots_adjust(wspace=0, hspace=-0.70)
plt.show()

# %%
topic_annot = topic_annotation(
    cistopic_obj,
    annot_var='major_clust',
    binarized_cell_topic=binarized_cell_topic,
    general_topic_thr = 0.2
)

# %%
topic_annot

# %%
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj.pkl"), "wb")
)

# %% [markdown]
# # Differentially Accessible Regions (DARs)

# %%
cistopic_obj = pickle.load(open(os.path.join(out_dir, "cistopic_obj.pkl"), "rb"))

# %%
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)
import numpy as np

# %%
imputed_acc_obj = impute_accessibility(
    cistopic_obj,
    selected_cells=None,
    selected_regions=None,
    scale_factor=10**6
)

# %%
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)

# %%
variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj,
    min_disp = 0.05,
    min_mean = 0.0125,
    max_mean = 3,
    max_disp = np.inf,
    n_bins=20,
    n_top_features=None,
    plot=True
)
gc.collect()

# %%
len(variable_regions)

# %%
file_path = os.path.join(out_dir, "find_diff_features_cistopic_obj.pkl")
with open(file_path, "wb") as file:
    pickle.dump(cistopic_obj, file)

file_path = os.path.join(out_dir, "find_diff_features_imputed_acc_obj.pkl")
with open(file_path, "wb") as file:
    pickle.dump(imputed_acc_obj, file)

file_path = os.path.join(out_dir, "find_diff_features_variable_regions.pkl")
with open(file_path, "wb") as file:
    pickle.dump(variable_regions, file)

# %%
# file_path = os.path.join(out_dir, "find_diff_features_cistopic_obj.pkl")
# with open(file_path, "rb") as file:
#     cistopic_obj = pickle.load(file)

# file_path = os.path.join(out_dir, "find_diff_features_imputed_acc_obj.pkl")
# with open(file_path, "rb") as file:
#     imputed_acc_obj = pickle.load(file)

# file_path = os.path.join(out_dir, "find_diff_features_variable_regions.pkl")
# with open(file_path, "rb") as file:
#     variable_regions = pickle.load(file)

# %%
markers_dict= find_diff_features(
    cistopic_obj,
    imputed_acc_obj,
    variable='major_clust',
    var_features=variable_regions,
    contrasts=None,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=1,
    _temp_dir='/tmp',
    split_pattern = '-'
)

# %%
print(markers_dict.keys())

# %%
# for x in ['L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev']:
#     print(f"Cell type: {x}")
#     print(markers_dict[x])
#     print("---")

# %%
# markers_dict['L2-3_CUX2'].head()

# %%
from pycisTopic.clust_vis import plot_imputed_features

# %%
# plot_imputed_features(
#     cistopic_obj,
#     reduction_name='UMAP',
#     imputed_data=imputed_acc_obj,
#     features=[markers_dict[x].index.tolist()[0] for x in  ['L2-3_CUX2']],
#     scale=False,
#     num_columns=3
# )

# %%
print("Number of DARs found:")
print("---------------------")
for x in markers_dict:
    print(f"  {x}: {len(markers_dict[x])}")

# %%
os.makedirs(os.path.join(out_dir, "region_sets"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_otsu"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_top_3k"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "DARs_cell_type"), exist_ok = True)

# %%
from pycisTopic.utils import region_names_to_coordinates

# %%
for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "Topics_otsu", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )
gc.collect()

# %%
for topic in region_bin_topics_top_3k:
    region_names_to_coordinates(
        region_bin_topics_top_3k[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "Topics_top_3k", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

# %%
for cell_type in markers_dict:
    region_names_to_coordinates(
        markers_dict[cell_type].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "DARs_cell_type", f"{cell_type}.bed"),
        sep = "\t",
        header = False, index = False
    )

# %%
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj.pkl"), "wb")
)

# %% [markdown]
# # Gene activity

# %%
import pyranges as pr
from pycisTopic.gene_activity import get_gene_activity

# %%
chromsizes = pd.read_table(os.path.join(out_dir, "qc", f"{reference}.chrom_sizes_and_alias.tsv"))
chromsizes

# %%
chromsizes.rename({"# ucsc": "Chromosome", "length": "End"}, axis = 1, inplace = True)
chromsizes["Start"] = 0
chromsizes = pr.PyRanges(chromsizes[["Chromosome", "Start", "End"]])

# %%
chromsizes

# %%
pr_annotation = pd.read_table(
        os.path.join(out_dir, "qc", "tss.bed")
    ).rename(
        {"Name": "Gene", "# Chromosome": "Chromosome"}, axis = 1)
pr_annotation["Transcription_Start_Site"] = pr_annotation["Start"]
pr_annotation = pr.PyRanges(pr_annotation)
pr_annotation
gc.collect()

# %%
gene_act, weigths = get_gene_activity(
    imputed_acc_obj,
    pr_annotation,
    chromsizes,
    use_gene_boundaries=True, # Whether to use the whole search space or stop when encountering another gene
    upstream=[1000, 100000], # Search space upstream. The minimum means that even if there is a gene right next to it
                             # these bp will be taken (1kbp here)
    downstream=[1000,100000], # Search space downstream
    distance_weight=True, # Whether to add a distance weight (an exponential function, the weight will decrease with distance)
    decay_rate=1, # Exponent for the distance exponential funciton (the higher the faster will be the decrease)
    extend_gene_body_upstream=10000, # Number of bp upstream immune to the distance weight (their value will be maximum for
                          #this weight)
    extend_gene_body_downstream=500, # Number of bp downstream immune to the distance weight
    gene_size_weight=False, # Whether to add a weights based on the length of the gene
    gene_size_scale_factor='median', # Dividend to calculate the gene size weigth. Default is the median value of all genes
                          #in the genome
    remove_promoters=False, # Whether to remove promoters when computing gene activity scores
    average_scores=True, # Whether to divide by the total number of region assigned to a gene when calculating the gene
                          #activity score
    scale_factor=1, # Value to multiply for the final gene activity matrix
    extend_tss=[10,10], # Space to consider a promoter
    gini_weight = True, # Whether to add a gini index weigth. The more unique the region is, the higher this weight will be
    return_weights= True, # Whether to return the final weights
    project='Gene_activity') # Project name for the gene activity object
# %%
print(type(gene_act))
print(type(weigths))

# %%
# Save gene_act to a file
with open(os.path.join(out_dir, 'gene_act.pkl'), 'wb') as file:
    pickle.dump(gene_act, file)

# Save weights to a file
weigths.to_pickle(os.path.join(out_dir, 'weights.pkl'))
# %%
DAG_markers_dict= find_diff_features(
    cistopic_obj,
    gene_act,
    variable='major_clust',
    var_features=None,
    contrasts=None,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=1,
    _temp_dir='/tmp',
    split_pattern = '-')

# %%
# plot_imputed_features(
#     cistopic_obj,
#     reduction_name='UMAP',
#     imputed_data=gene_act,
#     features=['ZBED1', 'OLIG2', 'SOX10', 
#                'ENPP6', 'OLIG1',
#                'VIP', 'SST',
#               'NFIB', 'SOX9'], 
#     scale=True,
#     num_columns=4
# )

# %%
print("Number of DAGs found:")
print("---------------------")
for x in markers_dict:
    print(f"  {x}: {len(DAG_markers_dict[x])}")

# %%
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj.pkl"), "wb")
)
gc.collect()
