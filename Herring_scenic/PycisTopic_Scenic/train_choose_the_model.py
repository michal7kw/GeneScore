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

# %%
cistopic_obj = pickle.load(open(os.path.join(out_dir, "cistopic_obj.pkl"), "rb"))

# %%
with open(os.path.join(out_dir, 'models.pkl'), 'rb') as file:
    models = pickle.load(file)

# %%
model = evaluate_models(
    models,
    select_model = 40,
    return_model = True
)

# %%
cistopic_obj.add_LDA_model(model)

# %%
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj.pkl"), "wb")
)
# %%
