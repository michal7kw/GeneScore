# %% [markdown]
# # Environment

# %%
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import scanpy as sc
from datetime import datetime

# Set working directory
# work_dir = '/home/michal.kubacki/Githubs/GeneScore/trimmed_GRN_derivation'
# work_dir = 'D:/Github/GeneScore/trimmed_GRN_derivation'
work_dir = '/mnt/d/Github/GeneScore/trimmed_GRN_derivation'

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

from grn_helpers import *
from results_analysis import *

# %%
# Configuration parameters
# neurons_set = "L2-3_CUX2"
neurons_set = "all_ex"
# neurons_set = "all_ex_all_ages" 

# Get base path from environment variables with error handling
root_dir = os.getenv('BASE_PATH')
if not root_dir:
    raise ValueError("BASE_PATH environment variable not found in .env file")

# %%
# Set up directories
output_dir, input_dir, root_dir, tmp_dir, in_dir_from_scenic = set_custom_folders(root_dir, neurons_set)

# Plot settings
plt.rcParams['figure.figsize'] = [10, 6]
plt.rcParams["savefig.dpi"] = 300

# %%
# Load the results from the previous script
print("Loading simulation scores...")
sim_scores_path = os.path.join(output_dir, 'scores_sim_all_new.csv')
sim_scores = pd.read_csv(sim_scores_path)
print(f"Loaded simulation scores with shape: {sim_scores.shape}")

print("\nLoading GRN scores...")
grn_scores_path = os.path.join(output_dir, 'scores_grn_all_from_comb_run_new.csv')
grn_scores = pd.read_csv(grn_scores_path)
print(f"Loaded GRN scores with shape: {grn_scores.shape}")

# %% [markdown]
# # Simulation results

# %%
# Display the first few rows of each dataset
print("\nSimulation scores sample:")
sim_scores.rename(columns={'gene': 'target_gene'}, inplace=True)
sim_scores['score'] = sim_scores['log_fold_change']
# print(sim_scores.head())
sim_scores.head()

# %%
# Get basic statistics for the datasets
print("\nSimulation scores statistics:")
print(sim_scores.describe())

# %%
fig, axes = plt.subplots(1, 2, figsize=(15, 5))

axes[0].hist(sim_scores['score'], bins=50)
axes[0].set_title('score')
axes[0].set_xlabel('score')
axes[0].set_ylabel('Frequency')

axes[1].hist(sim_scores['fold_change'], bins=50)
axes[1].set_title('fold_change')
axes[1].set_xlabel('fold_change')
axes[1].set_ylabel('Frequency')

plt.tight_layout()
plt.show()

print("Score:", f"Min: {min(sim_scores['score']):.6f}", f"Max: {max(sim_scores['score']):.6f}")
print("fold_change:", f"Min: {min(sim_scores['fold_change']):.6f}", f"Max: {max(sim_scores['fold_change']):.6f}")

# %%
# Get unique perturbed genes and cell types
perturbed_genes = sim_scores['goi'].unique()
print(f"Unique perturbed genes: {len(perturbed_genes)}")
print(perturbed_genes)

cell_types = sim_scores['cell_type'].unique()
print(f"\nUnique cell types: {len(cell_types)}")
print(cell_types)

# %% [markdown]
# # GRN results

# %%
print("\nGRN scores sample:")
grn_scores.rename(columns={'target': 'target_gene'}, inplace=True)
grn_scores.drop(columns=['score'], inplace=True)
grn_scores.rename(columns={'coef_mean': 'score'}, inplace=True)
# print(grn_scores.head())
grn_scores.head()

# %%
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

axes[0].hist(grn_scores['p'], bins=50)
axes[0].set_title('P Values')
axes[0].set_xlabel('p')
axes[0].set_ylabel('Frequency')

axes[1].hist(grn_scores['X.logp'], bins=50)
axes[1].set_title('X.logp')
axes[1].set_xlabel('X.logp')
axes[1].set_ylabel('Frequency')

axes[2].hist(grn_scores['score'], bins=50)
axes[2].set_title('Scores')
axes[2].set_xlabel('score')
axes[2].set_ylabel('Frequency')

plt.tight_layout()
plt.show()

print("P-values range:", f"Min: {min(grn_scores['p']):.6f}", f"Max: {max(grn_scores['p']):.6f}")
print("Coefficient means range:", f"Min: {min(grn_scores['X.logp']):.6f}", f"Max: {max(grn_scores['X.logp']):.6f}")
print("Scores range:", f"Min: {min(grn_scores['score']):.6f}", f"Max: {max(grn_scores['score']):.6f}")

# %%
print("\nGRN scores statistics:")
print(grn_scores.describe())

# %%
# Get unique perturbed genes and cell types
perturbed_genes = grn_scores['goi'].unique()
print(f"Unique perturbed genes: {len(perturbed_genes)}")
print(perturbed_genes)

cell_types = grn_scores['cell_type'].unique()
print(f"\nUnique cell types: {len(cell_types)}")
print(cell_types)

# %% [markdown]
# # Common analysis

# %%
example_gene = 'RARA'
sim_data, grn_data = analyze_gene_perturbation(sim_scores, grn_scores, example_gene)

# %%
# Find common target genes
common_genes = find_common_targets(sim_scores, grn_scores, example_gene, top_n=50)

# %%
# Print execution time
start_time = datetime.now()
print(f"Script started at {start_time}")
end_time = datetime.now()
print(f"Script ended at {end_time}")
print(f"Total execution time: {end_time - start_time}")


