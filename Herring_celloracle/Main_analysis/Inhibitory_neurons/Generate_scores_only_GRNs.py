import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from datetime import datetime
import importlib
import celloracle as co

sys.path.insert(0, "/home/michal.kubacki/Githubs/Re-MEND/code/External_Datasets/GeneSet_Derivation/Herring_celloracle/helpers")

import config
importlib.reload(config)
from config import *

reference = "hg19"

neurons_set = "all_inhibitory_all_ages"
# neurons_set = "all_inhibitory"

cells_dict = {
    "all_inhibitory"            :   ['SST', 'VIP', 'MGE_dev'],
    "all_inhibitory_all_ages"   :   ['VIP', 'SST', 'PV', 'MGE_dev']
}

ages_dict = {
    "all_inhibitory"            :   ['1m','3m','6m','10m','1y','2y','4y','ga22','ga24'],
    "all_inhibitory_all_ages"   :   ['1m','3m','6m','10m','1y','2y','4y','6y','10y','16y','20y','40y','ga22','ga24']
}

output_dir, input_dir, root_dir, tmp_dir, in_dir_from_scenic = set_custom_folders(reference, neurons_set)

sel_celltypes  = cells_dict[neurons_set]
sel_ages = ages_dict[neurons_set]

single_file = True
plotting = False


# Load scRNA-seq data
print("Load scRNA-seq data")
adata = sc.read_h5ad(os.path.join(output_dir, 'subseted_rna_andata.h5ad'))

# Add genes of interest
gois = ["RET", "NTRK1", "NTRK2", "NTRK3", "GFRA1", "GFRA2", "GFRA3", "GFRA4",
        "AHR", "ARNT", "ARNT2", "CLOCK",
        "AR",
        "NR1I2", "NR1I3",
        "NR3C1", "NR3C2",
        "ESR1", "GPER1",
        "DIO3", "DIO2",
        'RARA', 'ESR2', 'THRB',
        "THRA", "THRSP", "THRAP3"]

# Filter gois to only include genes present in the dataset
gois_present = [gene for gene in gois if gene in adata.var_names]

print(f"Genes of interest found: {len(gois_present)}")
print(f"Genes of interest not found: {set(gois) - set(gois_present)}")

all_grn = []

# Infer and fit GRN for each cell type
for cell_type in sel_celltypes:
    print(f"Processing cell type: {cell_type}")
    
    # Save links for the current cell type
    file_name = os.path.join(output_dir, f"{cell_type}.celloracle.links")
    
    try:
        links = co.load_hdf5(file_name)
    except FileNotFoundError:
        print(f"File not found: {file_name}")
        continue
    except Exception as e:
        print(f"Error reading file {file_name}: {str(e)}")
        continue
    
    for goi in gois_present:
        print(f"Processing {goi} for cell type {cell_type}")
        # Get scores from links
        for celltype in links.filtered_links:
            grn_data = links.filtered_links[celltype]
            grn_data = grn_data[grn_data["source"] == goi]
            if not grn_data.empty:
                grn_data["score"] = -np.log10(grn_data["p"])
                grn_data["celltype"] = celltype
                grn_data = grn_data.rename(columns={"-logp": "X.logp"})
                table_data = grn_data[["source", "target", "coef_mean", "coef_abs", "p", "X.logp", "score", "celltype"]]
                table_data['goi'] = goi
                all_grn.append(table_data)
                
                if not single_file:
                    try:
                        table_data.to_csv(os.path.join(output_dir, f'scores_grn_{goi}_{cell_type}.csv'), index=False)
                    except Exception as e:
                        print(f"Error writing file for {goi} in {cell_type}: {str(e)}")
            else:
                print(f"No data found for gene {goi} in cell type {celltype}")

# Concatenate all final_table_data DataFrames
if all_grn:
    try:
        all_grn_df = pd.concat(all_grn, ignore_index=True)
        all_grn_df.to_csv(os.path.join(output_dir, f'scores_grn_all.csv'), index=False)
        print(f"Combined GRN data saved to {os.path.join(output_dir, 'scores_grn_all.csv')}")
    except Exception as e:
        print(f"Error saving combined GRN data: {str(e)}")
else:
    print("No GRN data to save.")