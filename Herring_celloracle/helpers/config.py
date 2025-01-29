import os
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


root_dir = "/group/testa/michal.kubacki/herring"
# root_dir = "/scratch/michal.kubacki/herring"

def setup_config_folders(time_select_by, cell_select_by):
    ### time_select_by = ['stage_id', 'age', 'day']     ###
    ### cell_select_by = ['cell_type', 'major_clust']   ###

    if cell_select_by == "cell_type":
        selected_celltypes = ['PN'] # ['PN', 'IN', 'Non-Neu', 'Poor-Quality']
    elif cell_select_by == "major_clust":
        selected_celltypes = ['L2-3_CUX2', 'L4_RORB', 'L5-6_TLE4', 'L5-6_THEMIS']   # ['L4_RORB', 'L2-3_CUX2', 'SST', 'Astro', 'L5-6_TLE4', 'L5-6_THEMIS', 'VIP', 'OPC', 'PV', 'PV_SCUBE3', 'ID2', 'Oligo', 'LAMP5_NOS1', 'Micro', 'Vas', 'MGE_dev', 'Poor-Quality', 'PN_dev', 'CGE_dev']

    else:
        print("Not valid cells selector")
        return

    if time_select_by == 'stage_id':
        selected_ages = ['Neonatal']     # ['Adolescence', 'Adult', 'Childhood', 'Fetal', 'Infancy', 'Neonatal']
    elif time_select_by == 'age':
        selected_ages = ['ga22', 'ga24'] # ['34d', '2yr', '8yr', '2d', '86d', '16yr', 'ga22', '118d', '627d', '6yr', 'ga24', '179d', '4yr', '10yr', 'ga34', '301d', '20yr', '40yr', '422d', '12yr', '3yr', '14yr', '17yr', '25yr']
    elif time_select_by == 'day':
        selected_ages = 1001             # 0 - to select all days
    else:
        print("Not valid time selector")
        return
    
    if time_select_by != 'day': 
        ages_string = '_'.join(selected_ages) 
    else: 
        ages_string = str(selected_ages)

    celltypes_string = '_'.join(selected_celltypes)

    output_dir = f"{root_dir}/output/celloracle/age_{str(ages_string)}_celltypes_{celltypes_string}"
    input_dir = f"{root_dir}/data"

    if not os.path.exists(output_dir): 
        print("Creating:", end=" ")
        os.mkdir(output_dir)

    print("output_dir: ",output_dir)
    print("input_dir: ", input_dir)

    return output_dir, input_dir, selected_ages, selected_celltypes


def set_custom_folders(genome, suffix):
    tmp_dir = os.path.join(root_dir, "celloracle", "tmp")
    in_dir = os.path.join(root_dir, "data")
    if(genome == "hg38"):
        out_dir = os.path.join(root_dir, f"output_hg38_{suffix} ", "celloracle")
        in_dir_from_scenic = os.path.join(root_dir, f"output_hg38_{suffix}")
    else:
        out_dir = os.path.join(root_dir, f"output_hg19_{suffix}", "celloracle")
        in_dir_from_scenic = os.path.join(root_dir, f"output_hg19_{suffix}")

    print(f"root_dir: {root_dir}")
    print(f"out_dir: {out_dir}")
    print(f"in_dir: {in_dir}")
    print(f"tmp_dir: {tmp_dir}")

    os.makedirs(out_dir, exist_ok = True)
    return out_dir, in_dir, root_dir, tmp_dir, in_dir_from_scenic

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