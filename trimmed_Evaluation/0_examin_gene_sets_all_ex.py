# %% [markdown]
# # Environment

# %%
import os
import sys
import numpy as np
import pandas as pd

# Set working directory
# work_dir = '/home/michal.kubacki/Githubs/GeneScore/trimmed_Evaluation'
work_dir = 'D:/Github/GeneScore/trimmed_Evaluation'
# work_dir = '/mnt/d/Github/GeneScore/trimmed_Evaluation'
os.chdir(work_dir)

# Load environment variables from .env file
from dotenv import load_dotenv

# Explicitly specify the path to the .env file
env_path = os.path.join(work_dir, '.env')
load_dotenv(env_path)

# Get environment variables with error handling
project_functions_path = os.getenv('PROJECT_FUNCTIONS_PATH')
if not project_functions_path:
    raise ValueError("PROJECT_FUNCTIONS_PATH environment variable not found in .env file")

print(f"Using PROJECT_FUNCTIONS_PATH: {project_functions_path}")
sys.path.insert(0, project_functions_path)

import evaluated_helpers

import importlib
importlib.reload(evaluated_helpers)

from evaluated_helpers import *

# %%
base_path = os.getenv('BASE_PATH')
gene_set = "L2-3_CUX2"

# %% [markdown]
# # Load Gene Sets Data 

# %%
gere_sets_dict_celloracle, gene_sets_dict_cell_type_first_celloracle = load_GRNs_gene_sets(root_dir=base_path, gene_set_list = [gene_set], weights_list="scores_grn_all_from_comb_run_new.csv")

# %% [markdown]
# # Analyse separately

# %%
goi = "RARA"
score = "coef_mean"

# %% [markdown]
# ### Remove duplicates

# %%
print(gene_sets_dict_cell_type_first_celloracle[gene_set][gene_set].keys())
gois = gene_sets_dict_cell_type_first_celloracle[gene_set][gene_set].keys()

# %%
data={key: gene_sets_dict_cell_type_first_celloracle[gene_set][gene_set][key]["targets"] for key in gois}
duplicates = print_number_of_duplicate_genes(data)

# %%
if duplicates:
    gene_sets_dict_cell_type_first_celloracle = remove_duplicates_preserve_order_GRNs(gene_sets_dict_cell_type_first_celloracle)
    data={key: gene_sets_dict_cell_type_first_celloracle[gene_set][gene_set][key]["targets"] for key in gois}
    print_number_of_duplicate_genes(data)

# %% [markdown]
# ### Statistics

# %%
sets = list(gene_sets_dict_cell_type_first_celloracle.keys())
cell_types = list(gene_sets_dict_cell_type_first_celloracle[gene_set].keys())
scored_genes = list(gene_sets_dict_cell_type_first_celloracle[gene_set][gene_set].keys())
print(scored_genes)
print(len(gene_sets_dict_cell_type_first_celloracle[gene_set][gene_set][goi]['targets']))

# %%
for cell_type in cell_types:
    strings = [f"{scored_gene}: {len(gene_sets_dict_cell_type_first_celloracle[gene_set][cell_type][scored_gene]['targets'])}" for scored_gene in scored_genes]
    print(f'{cell_type}: {strings}') 


# %% [markdown]
# ### Heatmap

# %%
create_heatmap(gene_sets_dict_cell_type_first_celloracle, gene_set, scored_genes, cell_types)

# %% [markdown]
# ### Venn diagrams

# %%
scored_genes=["RARA", "THRB", "AR"]
cell_type = "L2-3_CUX2"
analyze_gene_sets_gene_set(gene_sets_dict_cell_type_first_celloracle, gene_set, cell_type=cell_type, scored_genes=scored_genes, mode="positive", printouts=True)

# %%
genes_to_mark = ['PAM', 'EPHA6', 'OXR1', 'RBFOX1', 'RASGRF2', 'CPVL']
better_hist_GRNs(gene_sets_dict_cell_type_first_celloracle, gene_set, gene_set, goi, score="coef_mean", bins=30, genes_to_mark=genes_to_mark)

# %%
genes_to_mark = ['BDNF', 'PRKG2', 'KIRREL3', 'CUX1', 'UBB']
better_hist_GRNs(gene_sets_dict_cell_type_first_celloracle, gene_set, gene_set, goi, score="coef_mean", bins=30, genes_to_mark=genes_to_mark)

# %%
scored_genes=["RARA", "THRB", "AR"]
cell_type = gene_set
analyze_gene_sets_gene_set(gene_sets_dict_cell_type_first_celloracle, gene_set, cell_type=cell_type, scored_genes=scored_genes, mode="negative")

# %%
cell_types=["L2-3_CUX2", "L4_RORB", "L5-6_THEMIS"]
scored_gene="NR3C2"

analyze_gene_sets_cell_types(gene_sets_dict_cell_type_first_celloracle, set_selected, scored_gene, cell_types, mode = "positive", printouts=True)

# %%
cell_types=["L2-3_CUX2", "L4_RORB", "L5-6_THEMIS"]
scored_gene="NR3C2"

analyze_gene_sets_cell_types(gene_sets_dict_cell_type_first_celloracle, set_selected, scored_gene, cell_types, mode = "negative")

# %% [markdown]
# ### Histograms

# %%
better_hist_GRNs(gene_sets_dict_cell_type_first_celloracle, set_selected, cell_type_selected, "AR", score="coef_mean", bins=30)

# %%
better_hist_GRNs(gene_sets_dict_cell_type_first_celloracle, set_selected, cell_type_selected, scored_gene_selected, score="coef_mean", bins=30)

# %%
visualize_weight_distributions(gene_sets_dict_cell_type_first_celloracle, set_selected, cell_type_selected, scored_genes, score_type="coef_mean")

# %% [markdown]
# ### Intersections

# %%
data={key: gene_sets_dict_cell_type_first_celloracle[set_selected][cell_type_selected][key]["targets"] for key in gois}

plot_gene_set_intersections(data, title="Intersection of Gene Sets (Ag)")

# %% [markdown]
# ### Histograms of the unique genes

# %% [markdown]
# Unique genes for each set:
# 
# L2-3_CUX2: {'SASH1', 'HTR1B', 'MICAL2', 'BCL9', 'RASGRF2'}
# 
# L4_RORB: {'GFRA1', 'RGS6', 'ATF5', 'GDF5', 'NCAM2', 'SESN3', 'FN1', 'KIF26B', 'NR4A3', 'PTPRZ1', 'TUBB2B', 'SOX5', 'BHLHE22', 'LINGO2'}
# 
# L5-6_THEMIS: {'HIST1H1E', 'COMMD3-BMI1', 'ETV5', 'TXNIP', 'CH25H', 'USP39', 'RARB', 'HIST1H2BN', 'LPL'}

# %%
genes_to_mark = ['SASH1', 'HTR1B', 'MICAL2', 'BCL9', 'RASGRF2']
better_hist_GRNs(gene_sets_dict_cell_type_first_celloracle, set_selected, "L2-3_CUX2", "NR3C2", score="coef_mean", bins=30, genes_to_mark=genes_to_mark)

# %%
genes_to_mark = ['GFRA1', 'RGS6', 'ATF5', 'GDF5', 'NCAM2', 'SESN3', 'FN1', 'KIF26B', 'NR4A3', 'PTPRZ1', 'TUBB2B', 'SOX5', 'BHLHE22', 'LINGO2']
better_hist_GRNs(gene_sets_dict_cell_type_first_celloracle, set_selected, "L4_RORB", "NR3C2", score="coef_mean", bins=30, genes_to_mark=genes_to_mark)

# %%
genes_to_mark = ['HIST1H1E', 'COMMD3-BMI1', 'ETV5', 'TXNIP', 'CH25H', 'USP39', 'RARB', 'HIST1H2BN', 'LPL']
better_hist_GRNs(gene_sets_dict_cell_type_first_celloracle, set_selected, "L5-6_THEMIS", "NR3C2", score="coef_mean", bins=30, genes_to_mark=genes_to_mark)


