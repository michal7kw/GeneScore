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

import celloracle as co
co.__version__

sys.path.insert(0, "/home/michal.kubacki/Githubs/GeneScore/Herring_celloracle/helpers")

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


