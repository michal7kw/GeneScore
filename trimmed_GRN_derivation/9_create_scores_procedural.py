# %%
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import celloracle as co
import importlib
import pickle
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

# Try to import from project_functions
try:
    from grn_helpers import *
except ImportError:
    print("Warning: Could not import from project_functions path, trying absolute path")
    # Try absolute import path as fallback
    # sys.path.insert(0, '/home/michal.kubacki/Githubs/GeneScore/project_functions')
    # sys.path.insert(0, 'D:/Github/GeneScore/project_functions')
    sys.path.insert(0,'/mnt/d/Github/GeneScore/project_functions')
    from grn_helpers import *

# %%
n_cpus = 20
single_file = True
plotting = True
neurons_set = "L2-3_CUX2"
# neurons_set = "all_ex"
# neurons_set = "all_ex_all_ages"
root_dir = os.getenv('BASE_PATH')

# %%
cells_dict = {
    "all_ex"            :   ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "all_ex_all_ages"   :   ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "L2-3_CUX2"         :   ['L2-3_CUX2']
}

ages_dict = {
    "all_ex"            :   ['1m','3m','6m','10m','1y','2y','4y','ga22','ga24'],
    "all_ex_all_ages"   :   ['1m','3m','6m','10m','1y','2y','4y','6y','10y','16y','20y','40y','ga22','ga24'],
    "L2-3_CUX2"         :   ['1m','3m','6m','10m','1y','2y','4y','ga22','ga24']
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
                                   'PN_dev': 'PN_dev.celloracle.parquet'},
    "L2-3_CUX2"         : {'L2-3_CUX2': 'L2-3_CUX2.celloracle.parquet'}
}

# %%
output_dir, input_dir, root_dir, tmp_dir, in_dir_from_scenic = set_custom_folders(root_dir, neurons_set)

sel_celltypes = cells_dict[neurons_set]
sel_ages = ages_dict[neurons_set]
motif_scan_files = motif_scan_files[neurons_set]

# Plot settings
plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300

# %%
gois = ['FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FGFRL1'] # FGF pathway
gois = gois + ['PTCH1', 'SMO', 'GLI1', 'GLI2', 'GLI3', 'GLI4'] # SAG pathway
gois = gois + ['BMPR1A', 'BMPR1B'] # BMP4 pathway
gois = gois + ['ACVR1'] # BMP7 pathway
gois = gois + ['CTNNB1', 'WNT5A', 'WNT3A', 'WNT3', 'APC', 'WNT10B'] # WNT pathway
gois = gois + ['RARA', 'RARB', 'RARG', 'RXRA', 'RXRB', 'RXRG'] # Retinoic Acid pathway
print(f"gois: {gois}")

# %%
print("Loading scRNA-seq data")
adata = sc.read_h5ad(os.path.join(output_dir, 'subseted_rna_andata.h5ad'))

hvgs = list(adata.var_names[adata.var['highly_variable']])
gois_present = [gene for gene in gois if gene in adata.var_names]
combined_genes = pd.Series(hvgs + gois_present).unique()
adata = adata[:, combined_genes]

print(f"Number of cells: {adata.n_obs}")
print(f"Number of genes: {adata.n_vars}")
print(f"Number of genes of interest found: {len(gois_present)}")
print(f"Genes of interest not found: {set(gois) - set(gois_present)}")
print(f"Unique cell types: {adata.obs['major_clust'].unique()}")

# %%
gois_present

# %%
oracle = co.Oracle()
oracle.import_anndata_as_raw_count(adata, cluster_column_name="major_clust", embedding_name="X_umap")

# %% [markdown]
# # Enhance TF-TG dictionary

# %% [markdown]
# ## 2023_11_CellOracleProof.tsv

# %%
df_grouped = pd.read_csv(os.path.join(input_dir, "2023_11_CellOracleProof.tsv"), delimiter="\t")

# %%
print(df_grouped.shape)
df_grouped.head()

# %%
for _, row in df_grouped[:5].iterrows():
    print(f"{row['TF']} n_targets: {len(row['Target_genes'].split(','))}")

# %%
# Find intersection between TFs and genes of interest
tf_array = df_grouped.TF.unique()
gois_present_in_tfs = np.intersect1d(tf_array, gois_present)
gois_present_in_tfs

# %%
# TF_to_TG_dictionary = {TF: TGs.replace(" ", "").split(",") for TF, TGs in zip(df_grouped.TF, df_grouped.Target_genes)}
# TG_to_TF_dictionary = co.utility.inverse_dictionary(TF_to_TG_dictionary)
# oracle.addTFinfo_dictionary(TG_to_TF_dictionary)

# %% [markdown]
# ## trrust_rawdata.human.tsv

# %%
df = pd.read_csv("./TF_TG/trrust_rawdata.human.tsv", sep="\t", header=None, 
                 names=["TF", "Target", "Mode", "PMID"])

# Group by TF and aggregate target genes into a comma-separated string.
df_grouped = df.groupby("TF")["Target"].apply(lambda genes: ",".join(genes)).reset_index()

# Rename the aggregated column to match the desired format.
df_grouped.rename(columns={"Target": "Target_genes"}, inplace=True)

# %%
print(df_grouped.shape)
df_grouped.head()

# %%
for _, row in df_grouped[:5].iterrows():
    print(f"{row['TF']} n_targets: {len(row['Target_genes'].split(','))}")

# %%
# Find intersection between TFs and genes of interest
tf_array = df_grouped.TF.unique()
gois_present_in_tfs = np.intersect1d(tf_array, gois_present)
gois_present_in_tfs

# %%
TF_to_TG_dictionary = {TF: TGs.replace(" ", "").split(",") for TF, TGs in zip(df_grouped.TF, df_grouped.Target_genes)}
TG_to_TF_dictionary = co.utility.inverse_dictionary(TF_to_TG_dictionary)
oracle.addTFinfo_dictionary(TG_to_TF_dictionary)

# %% [markdown]
# ## Brain_GTEx-regulons.txt

# %%
df = pd.read_csv("./TF_TG/Brain_GTEx-regulons.txt", sep="\t")

# Group by TF and aggregate the 'gene' column into a comma-separated string.
df_grouped = df.groupby("TF")["gene"].apply(lambda genes: ",".join(genes)).reset_index()

# Rename the aggregated column to 'Target_genes'.
df_grouped.rename(columns={"gene": "Target_genes"}, inplace=True)


# %%
print(df_grouped.shape)
df_grouped.head()

# %%
for _, row in df_grouped[:5].iterrows():
    print(f"{row['TF']} n_targets: {len(row['Target_genes'].split(','))}")

# %%
# Find intersection between TFs and genes of interest
tf_array = df_grouped.TF.unique()
gois_present_in_tfs = np.intersect1d(tf_array, gois_present)
gois_present_in_tfs

# %%
TF_to_TG_dictionary = {TF: TGs.replace(" ", "").split(",") for TF, TGs in zip(df_grouped.TF, df_grouped.Target_genes)}
TG_to_TF_dictionary = co.utility.inverse_dictionary(TF_to_TG_dictionary)
oracle.addTFinfo_dictionary(TG_to_TF_dictionary)

# %% [markdown]
# ## Fetal-Brain-regulons.txt

# %%
df = pd.read_csv("./TF_TG/Fetal-Brain-regulons.txt", sep="\t")

# Group by TF and aggregate the 'gene' column into a comma-separated string.
df_grouped = df.groupby("TF")["gene"].apply(lambda genes: ",".join(genes)).reset_index()

# Rename the aggregated column to 'Target_genes'.
df_grouped.rename(columns={"gene": "Target_genes"}, inplace=True)

# %%
print(df_grouped.shape)
df_grouped.head()

# %%
for _, row in df_grouped[:5].iterrows():
    print(f"{row['TF']} n_targets: {len(row['Target_genes'].split(','))}")

# %%
# Find intersection between TFs and genes of interest
tf_array = df_grouped.TF.unique()
gois_present_in_tfs = np.intersect1d(tf_array, gois_present)
gois_present_in_tfs

# %%
# TF_to_TG_dictionary = {TF: TGs.replace(" ", "").split(",") for TF, TGs in zip(df_grouped.TF, df_grouped.Target_genes)}
# TG_to_TF_dictionary = co.utility.inverse_dictionary(TF_to_TG_dictionary)
# oracle.addTFinfo_dictionary(TG_to_TF_dictionary)

# %% [markdown]
# ## TFLink_Homo_sapiens_interactions_SS_simpleFormat_v1.0.tsv

# %%
df = pd.read_csv("./TF_TG/TFLink_Homo_sapiens_interactions_SS_simpleFormat_v1.0.tsv", sep="\t")

# Group by the transcription factor column ("Name.TF") and aggregate the "Name.Target" column.
df_grouped = df.groupby("Name.TF")["Name.Target"].apply(lambda targets: ",".join(targets)).reset_index()

# Rename the columns to match the desired output.
df_grouped.rename(columns={"Name.TF": "TF", "Name.Target": "Target_genes"}, inplace=True)

# %%
print(df_grouped.shape)
df_grouped.head()

# %%
for _, row in df_grouped[:5].iterrows():
    print(f"{row['TF']} n_targets: {len(row['Target_genes'].split(','))}")

# %%
# Find intersection between TFs and genes of interest
tf_array = df_grouped.TF.unique()
gois_present_in_tfs = np.intersect1d(tf_array, gois_present)
gois_present_in_tfs

# %%
TF_to_TG_dictionary = {TF: TGs.replace(" ", "").split(",") for TF, TGs in zip(df_grouped.TF, df_grouped.Target_genes)}
TG_to_TF_dictionary = co.utility.inverse_dictionary(TF_to_TG_dictionary)
oracle.addTFinfo_dictionary(TG_to_TF_dictionary)

# %% [markdown]
# ## TFLink_Homo_sapiens_interactions_LS_simpleFormat_v1.0.tsv

# %%
df = pd.read_csv("./TF_TG/TFLink_Homo_sapiens_interactions_LS_simpleFormat_v1.0.tsv", sep="\t")

# Group by the transcription factor column ("Name.TF") and aggregate the "Name.Target" column.
df_grouped = df.groupby("Name.TF")["Name.Target"].apply(lambda targets: ",".join(targets)).reset_index()

# Rename the columns to match the desired output.
df_grouped.rename(columns={"Name.TF": "TF", "Name.Target": "Target_genes"}, inplace=True)

# %%
print(df_grouped.shape)
df_grouped.head()

# %%
for _, row in df_grouped[:5].iterrows():
    print(f"{row['TF']} n_targets: {len(row['Target_genes'].split(','))}")

# %%
# Find intersection between TFs and genes of interest
tf_array = df_grouped.TF.unique()
gois_present_in_tfs = np.intersect1d(tf_array, gois_present)
gois_present_in_tfs

# %%
TF_to_TG_dictionary = {TF: TGs.replace(" ", "").split(",") for TF, TGs in zip(df_grouped.TF, df_grouped.Target_genes)}
TG_to_TF_dictionary = co.utility.inverse_dictionary(TF_to_TG_dictionary)
oracle.addTFinfo_dictionary(TG_to_TF_dictionary)

# %% [markdown]
# ## TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv

# %%
df = pd.read_csv("./TF_TG/TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv", sep="\t")

# Group by the transcription factor column ("Name.TF") and aggregate the "Name.Target" column.
df_grouped = df.groupby("Name.TF")["Name.Target"].apply(lambda targets: ",".join(targets)).reset_index()

# Rename the columns to match the desired output.
df_grouped.rename(columns={"Name.TF": "TF", "Name.Target": "Target_genes"}, inplace=True)

# %%
print(df_grouped.shape)
df_grouped.head()

# %%
for _, row in df_grouped[:5].iterrows():
    print(f"{row['TF']} n_targets: {len(row['Target_genes'].split(','))}")

# %%
# Find intersection between TFs and genes of interest
tf_array = df_grouped.TF.unique()
gois_present_in_tfs = np.intersect1d(tf_array, gois_present)
gois_present_in_tfs

# %%
TF_to_TG_dictionary = {TF: TGs.replace(" ", "").split(",") for TF, TGs in zip(df_grouped.TF, df_grouped.Target_genes)}
TG_to_TF_dictionary = co.utility.inverse_dictionary(TF_to_TG_dictionary)
oracle.addTFinfo_dictionary(TG_to_TF_dictionary)

# %% [markdown]
# # DIM reduction

# %%
import pickle
# Check if saved data exists and load it, otherwise perform PCA
oracle_save_path = os.path.join(output_dir, 'oracle_after_pca.pkl')
pca_results_path = os.path.join(output_dir, 'pca_results.npz')

if os.path.exists(oracle_save_path) and os.path.exists(pca_results_path):
    print("Loading saved PCA results and oracle object...")
    with open(oracle_save_path, 'rb') as f:
        oracle = pickle.load(f)
    pca_data = np.load(pca_results_path)
    print("Successfully loaded saved data")
else:
    print("No saved data found. Performing PCA...")
    oracle.perform_PCA()
    
    # Save PCA results and oracle object
    print("Saving PCA results and oracle object...")
    with open(oracle_save_path, 'wb') as f:
        pickle.dump(oracle, f)
    
    # Save the PCA transformed data separately for quick access
    pca_data = {
        'pca_transformed': oracle.pca.transform(oracle.adata.X),
        'explained_variance_ratio': oracle.pca.explained_variance_ratio_,
        'components': oracle.pca.components_
    }
    np.savez(pca_results_path, **pca_data)
    
    print(f"Saved oracle object to: {oracle_save_path}")
    print("Saved PCA results to: pca_results.npz")

n_comps = min(np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0], 50)
# %%
plt.axvline(n_comps, c="k")
plt.savefig(os.path.join(output_dir, "pca_elbow.png"), bbox_inches='tight')
plt.close()
print(n_comps)
n_comps = min(n_comps, 50)

# %%
n_cell = oracle.adata.shape[0]
k = int(0.025*n_cell)
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8, b_maxl=k*4, n_jobs=n_cpus)

# %%
sc.pp.neighbors(oracle.adata)
sc.tl.umap(oracle.adata)

# %%
all_sim_top = []
all_grn_combined = []

# Iterate over the cell types and their corresponding motif scan files
for cell_type, motif_scan_file in motif_scan_files.items():
    print(f"Processing cell type: {cell_type}")

    # Load base GRN
    base_GRN_path = os.path.join(output_dir, motif_scan_file)
    print(f"Loading base GRN from: {base_GRN_path}")
    if not os.path.exists(base_GRN_path):
        print(f"Warning: Base GRN file not found at {base_GRN_path}. Skipping cell type {cell_type}.")
        continue  # Skip to the next cell type if file not found
    
    try:
        base_GRN = pd.read_parquet(base_GRN_path, engine='pyarrow')
        oracle.import_TF_data(TF_info_matrix=base_GRN)
    except Exception as e:
        print(f"Error loading or importing base GRN for {cell_type} from {base_GRN_path}: {e}")
        continue # Skip to the next cell type on error

    # Get links
    try:
        links = oracle.get_links(cluster_name_for_GRN_unit="major_clust", alpha=10, verbose_level=10, n_jobs=n_cpus)
        links.filter_links(p=0.05, weight="coef_abs", threshold_number=2000)
        links.get_network_score()
    except Exception as e:
        print(f"Error processing links for {cell_type}: {e}")
        continue # Skip to the next cell type on error

    # Save links
    try:
        file_name = os.path.join(output_dir, f"{cell_type}.celloracle.links")
        links.to_hdf5(file_path=file_name)
        print(f"Saved links for {cell_type} to {file_name}")
    except Exception as e:
        print(f"Error saving links for {cell_type}: {e}")
        # Continue processing even if saving fails

    if plotting:
        try:
            links.plot_degree_distributions(plot_model=True)
            plt.savefig(os.path.join(output_dir, f"degree_distributions_{cell_type}.png"), bbox_inches='tight')
            plt.close()
        except Exception as e:
            print(f"Error plotting degree distributions for {cell_type}: {e}")

    try:
        oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
        oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)
    except Exception as e:
        print(f"Error fitting GRN for simulation for {cell_type}: {e}")
        continue # Skip to the next cell type on error

    # Process each gene of interest
    for goi in gois_present:
        if goi in oracle.adata.var_names:
            print(f"Processing {goi} for cell type {cell_type}")
            
            if plotting:
                try:
                    sc.pl.umap(oracle.adata, color=[goi, oracle.cluster_column_name], layer="imputed_count", use_raw=False, cmap="viridis", show=False)
                    plt.savefig(os.path.join(output_dir, f"gene_expression_{goi}_{cell_type}.png"), bbox_inches='tight')
                    plt.close()
                except Exception as e:
                    print(f"Error plotting gene expression for {goi} in {cell_type}: {e}")

            try:
                # Simulate perturbation
                oracle.simulate_shift(perturb_condition={goi: 0.0}, n_propagation=3)
                oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
                oracle.calculate_embedding_shift(sigma_corr=0.05)

                # Get simulation scores
                sim_scores = oracle.get_simulation_score()
                sim_scores['cell_type'] = cell_type
                sim_scores['perturbed_gene'] = goi
                all_sim_top.append(sim_scores)

                # Get GRN scores
                grn_scores = links.get_network_score_for_each_target_gene()
                grn_scores['cell_type'] = cell_type
                grn_scores['perturbed_gene'] = goi
                all_grn_combined.append(grn_scores)

                if plotting:
                    oracle.plot_simulation_results(show=False)
                    plt.savefig(os.path.join(output_dir, f"simulation_results_{goi}_{cell_type}.png"), bbox_inches='tight')
                    plt.close()
            except Exception as e:
                print(f"Error during simulation or scoring for {goi} in {cell_type}: {e}")
                # Decide whether to continue with the next gene or cell type
                # continue # Uncomment to skip to the next gene on error

# %%
if all_sim_top:
    all_sim_save = pd.concat(all_sim_top, ignore_index=True)
    all_sim_save.to_csv(os.path.join(output_dir, 'scores_sim_all_new.csv'), index=False)

if all_grn_combined:
    all_grn_save = pd.concat(all_grn_combined, ignore_index=True)
    all_grn_save.to_csv(os.path.join(output_dir, 'scores_grn_all_from_comb_run_new.csv'), index=False)

# %%
start_time = datetime.now()
print(f"Script started at {start_time}")
end_time = datetime.now()
print(f"Script ended at {end_time}")
print(f"Total execution time: {end_time - start_time}")


