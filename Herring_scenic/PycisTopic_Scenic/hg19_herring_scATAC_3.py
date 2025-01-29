# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.6
#   kernelspec:
#     display_name: pycistopic
#     language: python
#     name: pycistopic
# ---

# 1. Creating a cisTopic object 

# # Environment

# +
# Standard library imports
import os
import gc
import pickle
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
n_cpus = 32

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

fragments_dict = select_files(reference, selected_fragments = sel_ages)

# # Creating a cisTopic object

# +
with open(os.path.join(out_dir,"sample_id_to_barcodes_passing_filters.pkl"), "rb") as file:
    sample_id_to_barcodes_passing_filters = pickle.load(file)

# print(sample_id_to_barcodes_passing_filters)

# +
# # %%script false --no-raise-error
gc.collect()
path_to_regions = os.path.join(out_dir, "consensus_peak_calling/consensus_regions.bed")
path_to_blacklist = os.path.join(in_dir, f"{reference}-blacklist.v2.bed")
pycistopic_qc_output_dir = os.path.join(out_dir,"qc")


cistopic_obj_list = []
for sample_id in fragments_dict:
    gc.collect()
    sample_metrics = pl.read_parquet(
        os.path.join(pycistopic_qc_output_dir, f'{sample_id}.fragments_stats_per_cb.parquet')
    ).to_pandas().set_index("CB").loc[ sample_id_to_barcodes_passing_filters[sample_id] ]
    gc.collect()
    cistopic_obj = create_cistopic_object_from_fragments(
        path_to_fragments = fragments_dict[sample_id],
        path_to_regions = path_to_regions,
        path_to_blacklist = path_to_blacklist,
        metrics = sample_metrics,
        valid_bc = sample_id_to_barcodes_passing_filters[sample_id],
        n_cpu = 1, #n_cpu
        project = str(sample_id),
        split_pattern = '-'
    )
    gc.collect()
    cistopic_obj_list.append(cistopic_obj)

# +
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
# -

gc.collect()
