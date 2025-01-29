import os
import sys
import pickle
from itertools import islice
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, vstack
from pycisTopic import *
import importlib
from sklearn.metrics.pairwise import cosine_similarity

sys.path.insert(0, "/home/michal.kubacki/Githubs/Re-MEND/code/External_Datasets/GeneSet_Derivation/Herring_scenic/helpers")
import config
importlib.reload(config)
from config import *
n_cpus = 32

reference = "hg19"
cell_types_set = "all_inhibitory_all_ages"
# # all_excitatory all_inhibitory ex_neurons ex_neurons_combined ex_progenitors
out_dir, in_dir, root_dir, tmp_dir, data_folder = set_output_folders(reference, cell_types_set)

print("Load cistopic_obj")
file_path = os.path.join(out_dir, "cistopic_obj.pkl")
with open(file_path, "rb") as file:
    cistopic_obj = pickle.load(file)

print("Get the fragment matrix and peak names from the cisTopic object")
fragment_matrix_ori = cistopic_obj.fragment_matrix
peak_names = cistopic_obj.region_names
peak_names = [f"{name.split(':')[0]}_{name.split(':')[1].replace('-', '_')}" for name in peak_names]
print(f"fragment_matrix_ori.shape: {fragment_matrix_ori.shape}")
print(F"peak_names: {peak_names[:10]}")

    
celltypes_dict = {
    "all_excitatory" : ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev'],
    "all_inhibitory" : ['SST', 'VIP', 'MGE_dev'],
    "all_cell_types" : ['L5-6_TLE4', 'L2-3_CUX2', 'L4_RORB', 'L5-6_THEMIS', 'PN_dev', 'SST', 'VIP'],
    "all_inhibitory_all_ages" : ['VIP', 'SST', 'PV', 'MGE_dev']
}

cell_types = celltypes_dict[cell_types_set]

cell_type_dir = os.path.join(out_dir, "cell_type_consensus_regions")
cell_type_fragments_dir = os.path.join(out_dir, "cell_type_fragments_data")

cell_types = list(reversed(cell_types))

for cell_type in cell_types:
    print(f"Processing cell type: {cell_type}")

    print(f"Loading consensus regions for {cell_type}")
    cell_type_regions_file = os.path.join(cell_type_dir, f"{cell_type}_consensus_regions.bed")
    cell_type_regions = pd.read_csv(cell_type_regions_file, sep='\t', header=None)
    cell_type_regions.columns = ['chrom', 'start', 'end', 'peak_id', 'score', 'strand', 'peak_name', 'cell_type']
    cell_type_regions.head()

    print(f"Creating peak_id to peak_name mapping for {cell_type}")
    peak_id_to_name = {}
    for peak_id, peak_name in zip(cell_type_regions['peak_id'], cell_type_regions['peak_name']):
        if ',' in peak_id:
            peak_ids = peak_id.split(',')
            for p_id in peak_ids:
                peak_id_to_name[p_id] = peak_name
        else:
            peak_id_to_name[peak_id] = peak_name
    first_10_elements = dict(islice(peak_id_to_name.items(), 10))
    print(first_10_elements)

    print(f"Loading cell_type_peak_indices for {cell_type}")
    cell_type_peak_indices_file = os.path.join(out_dir, f"{cell_type}_peak_indices.pkl")
    if os.path.exists(cell_type_peak_indices_file):
        print("Retriving from a file")
        with open(cell_type_peak_indices_file, "rb") as file:
            cell_type_peak_indices = pickle.load(file)
    else:
        print("Generating")
        print(f"cell_type_regions['peak_id']: {cell_type_regions['peak_id'][:10]}")
        cell_type_peak_indices = []
        for peak_id in cell_type_regions['peak_id']:
            if peak_id in peak_id_to_name and peak_id_to_name[peak_id] in peak_names:
                cell_type_peak_indices.append(peak_names.index(peak_id_to_name[peak_id]))
            # if ',' in peak_id:
            #     peak_ids = peak_id.split(',')
            #     for p_id in peak_ids:
            #         if p_id in peak_id_to_name and peak_id_to_name[p_id] in peak_names:
            #             cell_type_peak_indices.append(peak_names.index(peak_id_to_name[p_id]))
            # else:
            #     if peak_id in peak_id_to_name and peak_id_to_name[peak_id] in peak_names:
            #         cell_type_peak_indices.append(peak_names.index(peak_id_to_name[peak_id]))
        with open(cell_type_peak_indices_file, "wb") as file:
            pickle.dump(cell_type_peak_indices, file)

    print(f"Processing {len(cell_type_peak_indices)} peak indices for {cell_type}")

    max_value = max(cell_type_peak_indices)
    min_value = min(cell_type_peak_indices)

    print("Maximum value:", max_value)
    print("Minimum value:", min_value)

    valid_indices = np.array(cell_type_peak_indices)

    if len(valid_indices) > 0:
        chunk_size = 10000
        num_chunks = len(valid_indices) // chunk_size + 1
        coaccess_results = []

        for i in range(num_chunks):
            start_idx = i * chunk_size
            end_idx = min((i + 1) * chunk_size, len(valid_indices))
            chunk_indices = valid_indices[start_idx:end_idx]

            max_value = np.max(chunk_indices)
            min_value = np.min(chunk_indices)

            print("chunk_indices Maximum value:", max_value)
            print("chunk_indices Minimum value:", min_value)

            fragment_matrix_chunk = fragment_matrix_ori[chunk_indices, :]
            fragment_matrix_chunk = (fragment_matrix_chunk > 0).astype(np.bool_)
            corr_matrix_chunk = cosine_similarity(fragment_matrix_chunk)

            upper_tri_indices = np.triu_indices(corr_matrix_chunk.shape[0], k=1)

            peak1_indices_chunk = upper_tri_indices[0]
            peak2_indices_chunk = upper_tri_indices[1]
            correlations_chunk = corr_matrix_chunk[peak1_indices_chunk, peak2_indices_chunk]

            mask_chunk = np.abs(correlations_chunk) > 0.2

            peak1_indices_chunk = peak1_indices_chunk[mask_chunk]
            peak2_indices_chunk = peak2_indices_chunk[mask_chunk]
            correlations_chunk = correlations_chunk[mask_chunk]

            coaccess_results_chunk = np.column_stack((chunk_indices[peak1_indices_chunk], chunk_indices[peak2_indices_chunk], correlations_chunk))
            coaccess_results.append(coaccess_results_chunk)

        coaccess_results = np.concatenate(coaccess_results, axis=0)

        coaccess_df = pd.DataFrame(coaccess_results, columns=['Peak1_idx', 'Peak2_idx', 'coaccess'])

        print(f"Mapping peak indices to peak names for {cell_type}")
        coaccess_df['Peak1'] = coaccess_df['Peak1_idx'].astype(int).map(lambda x: peak_names[x].replace(':', '_').replace('-', '_'))
        coaccess_df['Peak2'] = coaccess_df['Peak2_idx'].astype(int).map(lambda x: peak_names[x].replace(':', '_').replace('-', '_'))

        coaccess_df = coaccess_df.drop(columns=['Peak1_idx', 'Peak2_idx'])

        print(f"Saving co-accessibility results for {cell_type}")
        coaccess_df.to_csv(os.path.join(out_dir, f"{cell_type}_coaccess.csv"), index=False)

        print(f"Co-accessibility results for cell type {cell_type} saved to {cell_type}_coaccess.csv")
    else:
        print(f"Skipping cell type {cell_type} due to out-of-bounds indices")

















# # %%
# import os
# import pickle
# from itertools import islice
# import numpy as np
# import pandas as pd
# from scipy.sparse import csr_matrix
# from pycisTopic import *
# import importlib
# import sys
# current_dir = os.getcwd()
# parent_dir = os.path.dirname(current_dir)
# if parent_dir not in sys.path:
#     sys.path.append(parent_dir)

# # import config_pc
# # importlib.reload(config_pc)
# # from config_pc import select_files, set_output_folders
# # n_cpu = 8

# import config
# importlib.reload(config)
# from config import select_files, set_output_folders
# n_cpu = 32

# reference = "hg19"
# # all_excitatory ex_neurons ex_neurons_combined ex_progenitors
# out_dir, in_dir, root_dir, tmp_dir, data_folder = set_output_folders(reference, "ex_progenitors")
# # %%
# print("Load cistopic_obj")
# file_path = os.path.join(out_dir, "cistopic_obj.pkl")
# with open(file_path, "rb") as file:
#     cistopic_obj = pickle.load(file)
# # %%
# print("Get the fragment matrix and peak names from the cisTopic object")
# fragment_matrix_ori = cistopic_obj.fragment_matrix
# peak_names = cistopic_obj.region_names
# peak_names = [f"{name.split(':')[0]}_{name.split(':')[1].replace('-', '_')}" for name in peak_names]
# print(f"fragment_matrix_ori.shape: {fragment_matrix_ori.shape}")
# print(F"peak_names: {peak_names[:10]}")
# # %%
# # Read the consensus_regions.bed file
# consensus_regions = pd.read_csv(os.path.join(out_dir, 'consensus_peak_calling/consensus_regions.bed'), sep='\t', header=None)
# consensus_regions.columns = ['chrom', 'start', 'end', 'peak_id', 'score', 'strand']

# # %%
# # Extract the cell type information from the peak_id column
# consensus_regions['cell_type'] = consensus_regions['peak_id'].str.split('_', n=1, expand=True)[0]
# # Get the unique cell types
# cell_types = consensus_regions['cell_type'].unique().tolist()

# # %%
# cell_type_dir = os.path.join(out_dir, "cell_type_consensus_regions")
# cell_type_fragments_dir = os.path.join(out_dir, "cell_type_fragments_data")

# cell_types = list(reversed(cell_types))

# for cell_type in cell_types:
#     print(f"Processing cell type: {cell_type}")

#     print(f"Loading consensus regions for {cell_type}")
#     cell_type_regions_file = os.path.join(cell_type_dir, f"{cell_type}_consensus_regions.bed")
#     cell_type_regions = pd.read_csv(cell_type_regions_file, sep='\t', header=None)
#     cell_type_regions.columns = ['chrom', 'start', 'end', 'peak_id', 'score', 'strand', 'cell_type', 'peak_name']
#     cell_type_regions.head()

#     print(f"Creating peak_id to peak_name mapping for {cell_type}")
#     peak_id_to_name = {}
#     for peak_id, peak_name in zip(cell_type_regions['peak_id'], cell_type_regions['peak_name']):
#         if ',' in peak_id:
#             peak_ids = peak_id.split(',')
#             for p_id in peak_ids:
#                 peak_id_to_name[p_id] = peak_name
#         else:
#             peak_id_to_name[peak_id] = peak_name
#     first_10_elements = dict(islice(peak_id_to_name.items(), 10))
#     print(first_10_elements)

#     print(f"Loading cell_type_peak_indices for {cell_type}")
#     cell_type_peak_indices_file = os.path.join(out_dir, f"{cell_type}_peak_indices.pkl")
#     if os.path.exists(cell_type_peak_indices_file):
#         with open(cell_type_peak_indices_file, "rb") as file:
#             cell_type_peak_indices = pickle.load(file)
#     else:
#         cell_type_peak_indices = []
#         for peak_id in cell_type_regions['peak_id']:
#             if ',' in peak_id:
#                 peak_ids = peak_id.split(',')
#                 for p_id in peak_ids:
#                     if p_id in peak_id_to_name and peak_id_to_name[p_id] in peak_names:
#                         cell_type_peak_indices.append(peak_names.index(peak_id_to_name[p_id]))
#             else:
#                 if peak_id in peak_id_to_name and peak_id_to_name[peak_id] in peak_names:
#                     cell_type_peak_indices.append(peak_names.index(peak_id_to_name[peak_id]))
#         with open(cell_type_peak_indices_file, "wb") as file:
#             pickle.dump(cell_type_peak_indices, file)

#     print(f"Processing {len(cell_type_peak_indices)} peak indices for {cell_type}")

#     max_value = max(cell_type_peak_indices)
#     min_value = min(cell_type_peak_indices)

#     print("Maximum value:", max_value)
#     print("Minimum value:", min_value)

#     valid_indices = np.array(cell_type_peak_indices)

#     from sklearn.metrics.pairwise import cosine_similarity

#     if len(valid_indices) > 0:

#         fragment_matrix = fragment_matrix_ori[valid_indices, :]
#         fragment_matrix = (fragment_matrix > 0).astype(np.bool_)  # Convert to binary matrix
#         corr_matrix = cosine_similarity(fragment_matrix)

#         upper_tri_indices = np.triu_indices(corr_matrix.shape[0], k=1)

#         peak1_indices = upper_tri_indices[0]
#         peak2_indices = upper_tri_indices[1]
#         correlations = corr_matrix[peak1_indices, peak2_indices]

#         # Create a mask to select rows where correlation is higher than 0.1
#         mask = np.abs(correlations) > 0.1

#         # Apply the mask to select only the rows with correlation higher than threshold
#         peak1_indices = peak1_indices[mask]
#         peak2_indices = peak2_indices[mask]
#         correlations = correlations[mask]

#         # Create a NumPy array to store the co-accessibility results
#         coaccess_results = np.column_stack((valid_indices[peak1_indices], valid_indices[peak2_indices], correlations))

#         # Convert the NumPy array to a pandas DataFrame
#         coaccess_df = pd.DataFrame(coaccess_results, columns=['Peak1_idx', 'Peak2_idx', 'coaccess'])

#         print(f"Mapping peak indices to peak names for {cell_type}")
#         # Modify the mapping step to format peak names as chr1_24175168_24175668
#         coaccess_df['Peak1'] = coaccess_df['Peak1_idx'].astype(int).map(lambda x: peak_names[x].replace(':', '_').replace('-', '_'))
#         coaccess_df['Peak2'] = coaccess_df['Peak2_idx'].astype(int).map(lambda x: peak_names[x].replace(':', '_').replace('-', '_'))

#         coaccess_df = coaccess_df.drop(columns=['Peak1_idx', 'Peak2_idx'])

#         print(f"Saving co-accessibility results for {cell_type}")
#         coaccess_df.to_csv(os.path.join(out_dir, f"{cell_type}_coaccess.csv"), index=False)

#         print(f"Co-accessibility results for cell type {cell_type} saved to {cell_type}_coaccess.csv")
#     else:
#         print(f"Skipping cell type {cell_type} due to out-of-bounds indices")