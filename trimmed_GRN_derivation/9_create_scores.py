import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import celloracle as co
import importlib
from datetime import datetime

from dotenv import load_dotenv
load_dotenv()
sys.path.insert(0, os.getenv('PROJECT_FUNCTIONS_PATH'))

from grn_helpers import set_custom_folders

# %%
# Configuration
n_cpus = 16
single_file = True
plotting = True
neurons_set = "all_ex"
# neurons_set = "all_ex_all_ages"
root_dir = os.getenv('BASE_PATH')

# %%
cells_dict = {
    "all_ex"            :   ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "all_ex_all_ages"   :   ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev']
}

ages_dict = {
    "all_ex"            :   ['1m','3m','6m','10m','1y','2y','4y','ga22','ga24'],
    "all_ex_all_ages"   :   ['1m','3m','6m','10m','1y','2y','4y','6y','10y','16y','20y','40y','ga22','ga24']
}

motif_scan_files = {
    "all_ex"            : {'L2-3_CUX2': 'L2-3_CUX2.celloracle.parquet',
                                   'L4_RORB': 'L4_RORB.celloracle.parquet',
                                   'L5-6_THEMIS': 'L5-6_THEMIS.celloracle.parquet',
                                   'L5-6_TLE4': '5-6_TLE4.celloracle.parquet',
                                   'PN_dev': 'PN_dev.celloracle.parquet'},
    "all_ex_all_ages"   : {'L2-3_CUX2': 'L2-3_CUX2.celloracle.parquet',
                                   'L4_RORB': 'L4_RORB.celloracle.parquet',
                                   'L5-6_THEMIS': 'L5-6_THEMIS.celloracle.parquet',
                                   'L5-6_TLE4': '5-6_TLE4.celloracle.parquet',
                                   'PN_dev': 'PN_dev.celloracle.parquet'}
}

# Setup directories
output_dir, input_dir, root_dir, tmp_dir, in_dir_from_scenic = set_custom_folders(root_dir, neurons_set)

sel_celltypes = cells_dict[neurons_set]
sel_ages = ages_dict[neurons_set]
motif_scan_files = motif_scan_files[neurons_set]

# Plot settings
plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300

# Genes of interest
gois = ["RET", "NTRK1", "NTRK2", "NTRK3", "GFRA1", "GFRA2", "GFRA3", "GFRA4",
        "AHR", "ARNT", "ARNT2", "CLOCK", "AR", "NR1I2", "NR1I3", "NR3C1", "NR3C2",
        "ESR1", "GPER1", "DIO3", "DIO2", 'RARA', 'ESR2', 'THRB', "THRA", "THRSP", "THRAP3"]

def load_and_process_data():
    print("Loading scRNA-seq data")
    try:
        adata = sc.read_h5ad(os.path.join(output_dir, 'subseted_rna_andata.h5ad'))
    except FileNotFoundError:
        print(f"Error: File 'subseted_rna_andata.h5ad' not found in {output_dir}")
        sys.exit(1)

    hvgs = list(adata.var_names[adata.var['highly_variable']])
    gois_present = [gene for gene in gois if gene in adata.var_names]
    combined_genes = pd.Series(hvgs + gois_present).unique()
    adata = adata[:, combined_genes]

    print(f"Number of cells: {adata.n_obs}")
    print(f"Number of genes: {adata.n_vars}")
    print(f"Number of genes of interest found: {len(gois_present)}")
    print(f"Genes of interest not found: {set(gois) - set(gois_present)}")
    print(f"Unique cell types: {adata.obs['major_clust'].unique()}")

    return adata, gois_present

def initialize_oracle(adata):
    print("Initializing CellOracle object")
    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(adata, cluster_column_name="major_clust", embedding_name="X_umap")

    print("Adding external data")
    try:
        external_data = pd.read_csv(os.path.join(input_dir, "2023_11_CellOracleProof.tsv"), delimiter="\t")
        TF_to_TG_dictionary = {TF: TGs.replace(" ", "").split(",") for TF, TGs in zip(external_data.TF, external_data.Target_genes)}
        TG_to_TF_dictionary = co.utility.inverse_dictionary(TF_to_TG_dictionary)
        oracle.addTFinfo_dictionary(TG_to_TF_dictionary)
    except FileNotFoundError:
        print(f"Error: File '2023_11_CellOracleProof.tsv' not found in {input_dir}")
        sys.exit(1)

    print("Performing PCA and KNN imputation")
    oracle.perform_PCA()
    n_comps = min(np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0], 50)
    n_cell = oracle.adata.shape[0]
    k = int(0.025*n_cell)
    oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8, b_maxl=k*4, n_jobs=n_cpus)

    print("Calculating UMAP")
    sc.pp.neighbors(oracle.adata)
    sc.tl.umap(oracle.adata)

    return oracle

def process_cell_type(oracle, cell_type, motif_scan_file, gois_present):
    print(f"Processing cell type: {cell_type}")
    
    file_path = os.path.join(output_dir, motif_scan_file)
    try:
        base_GRN = pd.read_parquet(file_path, engine='pyarrow')
        oracle.import_TF_data(TF_info_matrix=base_GRN)
    except FileNotFoundError:
        print(f"Error: File '{motif_scan_file}' not found in {output_dir}")
        return None, None

    links = oracle.get_links(cluster_name_for_GRN_unit="major_clust", alpha=10, verbose_level=10, n_jobs=n_cpus)
    # links.filter_links()
    links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)

    links.get_network_score()
    
    # Save links for the current cell type
    file_name = os.path.join(output_dir, f"{cell_type}.celloracle.links")
    links.to_hdf5(file_path=file_name)

    if plotting:
        links.plot_degree_distributions(plot_model=True)
        plt.savefig(os.path.join(output_dir, f"degree_distributions_{cell_type}.png"), bbox_inches='tight')
        plt.close()

    oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
    oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)

    return process_gois(oracle, cell_type, gois_present, links)

def process_gois(oracle, cell_type, gois_present, links):
    all_sim = []
    all_grn = []

    for goi in gois_present:
        try:
            if goi in oracle.adata.var_names:
                print(f"Processing {goi} for cell type {cell_type}")
                
                if plotting:
                    plot_gene_expression(oracle, goi, cell_type)

                oracle.simulate_shift(perturb_condition={goi: 0.0}, n_propagation=3)
                oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
                oracle.calculate_embedding_shift(sigma_corr=0.05)

                if plotting:
                    plot_quiver(oracle, goi, cell_type)

                oracle.calculate_p_mass(smooth=0.8, n_grid=40, n_neighbors=200)
                oracle.calculate_mass_filter(min_mass=60, plot=plotting)

                if plotting:
                    plot_simulation_flow(oracle, goi, cell_type)

                sim_data = process_simulation_results(oracle, goi, cell_type)
                all_sim.extend(sim_data)

                grn_data = process_grn_data(links, goi, cell_type)
                all_grn.extend(grn_data)

        except Exception as e:
            print(f"Error processing gene {goi} for cell type {cell_type}: {str(e)}")

    return all_sim, all_grn

def process_simulation_results(oracle, goi, cell_type):
    simulated_count = oracle.adata.layers["simulated_count"]
    original_count = oracle.adata.X.toarray()
    log_fold_change = np.log2(simulated_count + 1) - np.log2(original_count + 1)
    
    top_n = 1000
    local_cell_types = oracle.adata.obs['major_clust']
    gene_names = oracle.adata.var_names
    unique_cell_types = local_cell_types.unique()

    all_data = []
    for local_cell_type in unique_cell_types:
        cell_type_indices = np.where(local_cell_types == local_cell_type)[0]
        cell_type_log_fold_change = log_fold_change[cell_type_indices, :].mean(axis=0)
        abs_cell_type_log_fold_change = np.abs(cell_type_log_fold_change)
        
        sorted_indices = np.argsort(abs_cell_type_log_fold_change)[::-1]
        top_gene_indices = sorted_indices[:top_n]
        top_genes = gene_names[top_gene_indices]
        top_log_fold_changes = cell_type_log_fold_change[top_gene_indices]
        
        cell_type_data = pd.DataFrame({
            'local_cell_type': [local_cell_type] * top_n,
            'gene': top_genes,
            'log_fold_change': top_log_fold_changes,
            'goi': [goi] * top_n,
            'fold_change': np.exp2(top_log_fold_changes),
            'cell_type': [cell_type] * top_n
        })
        all_data.append(cell_type_data)

    return all_data

def process_grn_data(links, goi, cell_type):
    all_grn_data = []
    for celltype in links.filtered_links:
        grn_data = links.filtered_links[celltype]
        grn_data = grn_data[grn_data["source"] == goi]
        if not grn_data.empty:
            grn_data["score"] = -np.log10(grn_data["p"])
            grn_data["celltype"] = celltype
            grn_data = grn_data.rename(columns={"-logp": "X.logp"})
            table_data = grn_data[["source", "target", "coef_mean", "coef_abs", "p", "X.logp", "score", "celltype"]]
            table_data['goi'] = goi
            table_data['cell_type'] = cell_type
            all_grn_data.append(table_data)
    return all_grn_data

def plot_gene_expression(oracle, goi, cell_type):
    sc.pl.umap(oracle.adata, color=[goi, oracle.cluster_column_name], layer="imputed_count", use_raw=False, cmap="viridis")
    plt.savefig(os.path.join(output_dir, f"gene_expression_{goi}_{cell_type}.png"), bbox_inches='tight')
    plt.close()

def plot_quiver(oracle, goi, cell_type):
    fig, ax = plt.subplots(1, 2, figsize=[13, 6])
    scale = 50
    oracle.plot_quiver(scale=scale, ax=ax[0])
    ax[0].set_title(f"Simulated cell identity shift vector: {goi} KO")
    oracle.plot_quiver_random(scale=scale, ax=ax[1])
    ax[1].set_title(f"Randomized simulation vector")
    plt.savefig(os.path.join(output_dir, f"quiver_plot_{goi}_{cell_type}.png"), bbox_inches='tight')
    plt.close()

def plot_simulation_flow(oracle, goi, cell_type):
    fig, ax = plt.subplots(1, 2, figsize=[13, 6])
    scale_simulation = 10
    oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
    ax[0].set_title(f"Simulated cell identity shift vector: {goi} KO")
    oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
    ax[1].set_title(f"Randomized simulation vector")
    plt.savefig(os.path.join(output_dir, f"simulation_flow_{goi}_{cell_type}.png"), bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=[8, 8])
    oracle.plot_cluster_whole(ax=ax, s=5)
    oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax, show_background=False)
    plt.savefig(os.path.join(output_dir, f"cluster_specific_simulation_flow_{goi}_{cell_type}.png"), bbox_inches='tight')
    plt.close()

def main():
    adata, gois_present = load_and_process_data()
    oracle = initialize_oracle(adata)

    all_sim_top = []
    all_grn_combined = []

    for cell_type, motif_scan_file in motif_scan_files.items():
        all_sim, all_grn = process_cell_type(oracle, cell_type, motif_scan_file, gois_present)
        if all_sim:
            all_sim_cell = pd.concat(all_sim, ignore_index=True)
            all_sim_top.append(all_sim_cell)
        if all_grn:
            all_grn_combined.extend(all_grn)

    if all_sim_top:
        all_sim_save = pd.concat(all_sim_top, ignore_index=True)
        all_sim_save.to_csv(os.path.join(output_dir, 'scores_sim_all_new.csv'), index=False)
        print(f"Simulation scores saved to {os.path.join(output_dir, 'scores_sim_all_new.csv')}")

    if all_grn_combined:
        all_grn_save = pd.concat(all_grn_combined, ignore_index=True)
        all_grn_save.to_csv(os.path.join(output_dir, 'scores_grn_all_from_comb_run_new.csv'), index=False)
        print(f"GRN scores saved to {os.path.join(output_dir, 'scores_grn_all_from_comb_run_new.csv')}")

if __name__ == "__main__":
    start_time = datetime.now()
    print(f"Script started at {start_time}")

    try:
        main()
    except Exception as e:
        print(f"An error occurred during script execution: {str(e)}")
    finally:
        end_time = datetime.now()
        print(f"Script ended at {end_time}")
        print(f"Total execution time: {end_time - start_time}")