# %% [markdown]
# # Environment

# %%
import gc
import os
import importlib
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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

# %%
n_cpus = 8
neurons_set = "L2-3_CUX2"
# neurons_set = "all_ex"
# neurons_set = "all_ex_all_ages"
root_dir = os.getenv('BASE_PATH')

# %%
cells_dict = {
    "all_ex"            :   ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "all_ex_all_ages"   :   ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "L2-3_CUX2"         :   ['L2-3_CUX2'],
    "all_ex_comb"       :   ['ex_neurons']
}

ages_dict = {
    "all_ex"            :   ['1m','3m','6m','10m','1y','2y','4y','ga22','ga24'],
    "all_ex_all_ages"   :   ['1m','3m','6m','10m','1y','2y','4y','6y','10y','16y','20y','40y','ga22','ga24'],
    "L2-3_CUX2"         :   ['1m','3m','6m','10m','1y','2y','4y','ga22','ga24'],
    "all_ex_comb"       :   ['1m','3m','6m','10m','1y','2y','4y','ga22','ga24']
}

output_dir, input_dir, root_dir, tmp_dir, in_dir_from_scenic = set_custom_folders(root_dir, neurons_set)

sel_celltypes  = cells_dict[neurons_set]
sel_ages = ages_dict[neurons_set]

# %%

plt.rcParams['figure.figsize'] = (15,7)
plt.rcParams["savefig.dpi"] = 600

# %% [markdown]
# # Check results

# %%
cell_type = sel_celltypes[0]

# %%
file_path = os.path.join(output_dir, f"{cell_type}.celloracle.parquet")

# %%
df = pd.read_parquet(file_path)

# %%
df.shape

# %%
# total number of non zero elements
result = df.iloc[:, 2:].astype(bool).sum().sum()
print(result)

# %%
selected_columns = df.iloc[:, 2:]
row_counts = selected_columns.apply(lambda row: len(row[row != 0.0]), axis=1)
print(list(row_counts))

first_row = df.iloc[0]
non_zero_elements = first_row[first_row != 0]
print(non_zero_elements)

# %% [markdown]
# ## Compare with `2023_11_tfi.celloracle.parquet`

# %%
base_GRN = pd.read_parquet(os.path.join(input_dir, "2023_11_tfi.celloracle.parquet"), engine='pyarrow')

# %%
base_GRN.head()

# %%
base_GRN.shape

# %%
selected_columns = base_GRN.iloc[:, 2:]
row_counts = selected_columns.apply(lambda row: len(row[row != 0.0]), axis=1)
print(list(row_counts[:100]))

first_row = base_GRN.iloc[0]
non_zero_elements = first_row[first_row != 0]
print(non_zero_elements)

# %%
result2 = base_GRN.iloc[:, 2:].astype(bool).sum().sum()
print(result2)

# %%
print(result/result2)


