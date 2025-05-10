# %% [markdown]
# # Environment

# %%
from IPython.display import display
import os
import sys
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ipywidgets as widgets

from dotenv import load_dotenv
load_dotenv()
sys.path.insert(0, os.getenv('PROJECT_FUNCTIONS_PATH'))

from evaluated_helpers import process_gene_sets, boxplot_EDCs_GRN_scores_parameters_local

from gene_scoring import score_genes

# %%
gpu_support = False
recompute = True
plotting = True

# %%
base_path = os.getenv('BASE_PATH')
root_dir = base_path
data_path = os.path.join(base_path, "data")
output_path = os.path.join(base_path, "all_ex", "results")

# %% [markdown]
# # Load Precomputed Scores

# %%
%%capture
if not recompute:
    file_name = os.path.join(output_path, f"EDCs_andata_scored_sim.loom")
    adata = ad.read_loom(file_name, sparse=False)

    adata.var['original_var_names'] = adata.var.index
    adata.var_names = adata.var['var_names']

    adata.var_names_make_unique()
    adata.var.set_index('var_names', inplace=True)

# %% [markdown]
# # Load Gene Sets Data

# %%
gene_set_list = [
        "all_ex"
    ]

# %%
gene_sets = {}

for gene_set in gene_set_list:
    path = os.path.join(root_dir, f"{gene_set}", "celloracle")
    print(path)

    # Load scores_sim_all.csv
    sim_1_path = os.path.join(path, 'scores_sim_all_new.csv')
    if os.path.exists(sim_1_path):
        gene_sets[gene_set] = pd.read_csv(sim_1_path)
    else:
        print(f"File not found: {sim_1_path}")

# %%
gene_sets.keys()

# %%
gene_sets['all_ex'].head()

# %%
from project_functions.load_gene_sets import process_gene_sets


gene_sets_dict_cell_type_first = process_gene_sets(gene_sets)

# %%
sets = list(gene_sets.keys())
print(sets)

set_selected = sets[0]
cell_types = list(gene_sets_dict_cell_type_first[set_selected].keys())
print(cell_types)

cell_type_selected = cell_types[0]
scored_genes = list(gene_sets_dict_cell_type_first[set_selected][cell_type_selected].keys())
print(scored_genes)

scored_gene_selected = scored_genes[0]
print(len(gene_sets_dict_cell_type_first[set_selected][cell_type_selected][scored_gene_selected]['targets']))

# %% [markdown]
# # Load Expression Data

# %%
if recompute:
    adata = ad.read_h5ad(os.path.join(data_path,'CTL04_EDCs.h5ad'))

# %%
adata

# %%
adata.obs.columns

# %%
adata.obs['condition_concentraion'] = [item1 + '_' + item2 for item1, item2 in zip(list(adata.obs.Condition), list(adata.obs.Concentration))]

# %%
adata.obs.condition_concentraion.unique()

# %%
samples_count = [f"{unique}: {np.sum(adata.obs.condition_concentraion == unique)}" for unique in list(adata.obs.condition_concentraion.unique())]
print(samples_count)

# %% [markdown]
# # Scoring

# %%
n_genes = 80
print(list(gene_sets_dict_cell_type_first["all_ex"].keys()))
print(list(gene_sets_dict_cell_type_first["all_ex"]['L4_RORB'].keys()))

# %%
gois = ['AHR', 'AR', 'NR1I2', 'NR1I3', 'NR3C1', 'NR3C2', 'ESR1', 'RARA', 'ESR2', 'THRB', 'THRA']
gene_sets = ['all_ex']
cell_types = ['L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'L5-6_TLE4', 'PN_dev']

# %%
if recompute:
    for control in [True, False]:
        for control_condition in ['DMSO', None]:
            for normalize_weights in [True, False]:
                for scaling_only_based_on_control in [True, False]: 
                    for scale_by_variance in [True, False]:
                        for gene_set in gene_sets:
                            for cell_type in list(gene_sets_dict_cell_type_first[gene_set].keys()):
                                for goi in list(gene_sets_dict_cell_type_first[gene_set][cell_type].keys()):
                                    score_genes(
                                        adata,
                                        gene_list=gene_sets_dict_cell_type_first[gene_set][cell_type][goi]['targets'][:n_genes], 
                                        gene_weights=gene_sets_dict_cell_type_first[gene_set][cell_type][goi]['coef_mean'][:n_genes],   
                                        score_name = (
                                            f'gene_score_{gene_set}_{cell_type}_{goi}_{control}_'
                                            f'normalized_{normalize_weights}_'
                                            f'scaled_{scale_by_variance}_'
                                            f'cc_{control_condition}_'
                                            f'sc_{scaling_only_based_on_control}'
                                        ),                                    
                                        ctrl_size=50,
                                        gene_pool=None,
                                        n_bins=25,
                                        random_state=0,
                                        copy=False,
                                        used_layer='cpm',
                                        return_scores=False,
                                        control=control,
                                        weighted=True,
                                        abs_diff=False,
                                        gpu=gpu_support,
                                        chunk_size=10000,
                                        disable_chunking=True,
                                        scale_by_variance=scale_by_variance,
                                        normalize_weights=normalize_weights,
                                        conditions_labels='Condition',
                                        control_condition=control_condition,
                                        debug=False,
                                        scaling_only_based_on_control=scaling_only_based_on_control
                                )

# %% [markdown]
# # Plotting

# %% [markdown]
# ## Configure

# %%
%%capture
if plotting:
    condition_dropdown = widgets.Dropdown(
        options=list(adata.obs.Condition.unique()),
        value=list(adata.obs.Condition.unique())[0],
        description='Condition:',
        disabled=False,
    )

    gene_set_dropdown = widgets.Dropdown(
        options=list(gene_sets_dict_cell_type_first.keys()),
        value=list(gene_sets_dict_cell_type_first.keys())[0],
        description='Gene Set:',
        disabled=False,
    )

    control_dropdown = widgets.Dropdown(
        options=list(['True', 'False']),
        value=list(['True', 'False'])[0],
        description='Control:',
        disabled=False,
    )

    control_condition_dropdown = widgets.Dropdown(
        options=list(['DMSO', "None"]),
        value=list(['DMSO', "None"])[0],
        description='Condition Control:',
        disabled=False,
    )

    normalized_dropdown = widgets.Dropdown(
        options=list(['True', 'False']),
        value=list(['True', 'False'])[0],
        description='Normalized weights:',
        disabled=False,
    )

    scaled_dropdown = widgets.Dropdown(
        options=list(['True', 'False']),
        value=list(['True', 'False'])[0],
        description='Scale by variance:',
        disabled=False,
    )

    scaling_only_based_on_control_dropdown = widgets.Dropdown(
        options=list(['True', 'False']),
        value=list(['True', 'False'])[0],
        description='Scale only with Control:',
        disabled=False,
    )

    cell_type_dropdown = widgets.Dropdown(
        options=cell_types,
        value=cell_types[0],
        description='Cell Type:',
        disabled=False,
    )

    scored_gene_dropdown = widgets.Dropdown(
        options=gois,
        value=gois[0],
        description='Scored Gene:',
        disabled=False,
    )

# %% [markdown]
# ## Display

# %%
if plotting:
    interactive_plot = widgets.interactive(boxplot_EDCs_GRN_scores_parameters_local,
                                    adata=widgets.fixed(adata), # type: ignore
                                    conditions=condition_dropdown, # type: ignore
                                    gene_set=gene_set_dropdown, # type: ignore
                                    cell_type=cell_type_dropdown, # type: ignore
                                    goi=scored_gene_dropdown, # type: ignore
                                    control=control_dropdown, # type: ignore
                                    normalize_weights=normalized_dropdown, # type: ignore
                                    scale_by_variance=scaled_dropdown, # type: ignore
                                    control_condition=control_condition_dropdown, # type: ignore
                                    scaling_only_based_on_control=scaling_only_based_on_control_dropdown # type: ignore
                                    )

    display(interactive_plot)

# %% [markdown]
# # Save results

# %%
if recompute:
    file_name = os.path.join(output_path, f"EDCs_andata_scored_sim.loom")
    adata.write_loom(file_name)


