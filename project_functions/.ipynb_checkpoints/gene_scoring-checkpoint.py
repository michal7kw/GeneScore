from scanpy._utils import _check_use_raw
from scipy.sparse import issparse
import numpy as np
import os
import pandas as pd
from joblib import Parallel, delayed
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad
import scipy.sparse as sp
import numpy as np
import pandas as pd
from scipy import sparse
import numba
import warnings
from typing import Optional, Union, List
try:
    import cupy as cp
except ImportError:
    cp = None

def _sparse_nanmean(X, axis):
    """
    np.nanmean equivalent for sparse matrices
    """
    if not issparse(X):
        raise TypeError("X must be a sparse matrix")

    # count the number of nan elements per row/column (dep. on axis)
    Z = X.copy()
    Z.data = np.isnan(Z.data)
    Z.eliminate_zeros()
    n_elements = Z.shape[axis] - Z.sum(axis)

    # set the nans to 0, so that a normal .sum() works
    Y = X.copy()
    Y.data[np.isnan(Y.data)] = 0
    Y.eliminate_zeros()

    # the average
    s = Y.sum(axis, dtype="float64")  # float64 for score_genes function compatibility)
    m = s / n_elements

    return m

def map_cell_types(adata, endocrine_rec_networks, slot1, slot2, cell_type_mapping):
    adata_cell_types = adata.obs[slot1].unique()
    
    endocrine_cell_types = endocrine_rec_networks[slot2].unique()
    
    endocrine_rec_networks['mapped_cell_types'] = endocrine_rec_networks[slot2].apply(
        lambda x: [key for key, value in cell_type_mapping.items() if value and x in value]
    )

    # If no match is found, the value will be an empty list, replace it with None
    endocrine_rec_networks['mapped_cell_types'] = endocrine_rec_networks['mapped_cell_types'].apply(
        lambda x: x if x else None
    )
    
    return endocrine_rec_networks

def process_network_celltype_mapping(df, source, cell_type_mapping):
    # Initialize the dictionaries to store the results
    genes = {}
    effects = {}
    
    for celltype in cell_type_mapping:
        print(f"key: {celltype}, value: {cell_type_mapping[celltype]}")
        for mapped_cell in cell_type_mapping[celltype]:
            # Ensure the filtering matches 'mapped_cell' with the 'celltype' column in the dataframe
            filtered_network = df[(df['source'] == source) & (df['celltype'] == celltype)]
            # Sort by 'p' value
            filtered_network = filtered_network.sort_values(by='p')
            # Store the results in the dictionaries
            genes[mapped_cell] = filtered_network['target'].tolist()
            effects[mapped_cell] = filtered_network['coef_mean'].tolist()
            
    return genes, effects

#######################################################################################
@numba.njit(parallel=True)
def fast_control_avgs(X, gene_indices, cut_indices, n_cells):
    """
    Calculate control averages for each gene using Numba for acceleration.
    
    This function computes the average expression of control genes for each target gene.
    Control genes are those in the same expression bin as the target gene, excluding the target gene itself.
    """
    control_avgs = np.zeros((len(gene_indices), n_cells))
    for i in numba.prange(len(gene_indices)):
        mask = cut_indices == cut_indices[i]
        mask[gene_indices[i]] = False
        masked_X = X[:, mask]
        for j in range(n_cells):
            row = masked_X[j]
            control_avgs[i, j] = np.nanmean(row)
    return control_avgs

def calculate_gene_variances(X, gene_indices):
    """
    Calculate the variance of each gene in the gene_list.
    """
    if sparse.issparse(X):
        variances = X.power(2).mean(axis=0).A1 - np.power(X.mean(axis=0).A1, 2)
    else:
        variances = np.nanvar(X, axis=0)
    return variances[gene_indices]

def score_genes_weighted_optimized(
    adata,
    gene_list,
    gene_weights=None,
    ctrl_size=50,
    gene_pool=None,
    n_bins=25,
    score_name="score",
    random_state=0,
    copy=False,
    use_raw=False,
    return_scores=False,
    control=True,
    weighted=True,
    abs_diff=False,
    gpu=True,
    chunk_size=10000,
    disable_chunking=True,
    scale_by_variance=False,
    normalize_weights=False    
):
    """
    Calculate gene set scores for cells in an AnnData object.
    
    This function computes scores based on the expression of a given list of genes,
    optionally comparing against a set of control genes.
    """
    
    # Check if CuPy is available for GPU acceleration
    if gpu and cp is None:
        raise ImportError("CuPy is required for GPU acceleration. Install it with 'pip install cupy'.")

    # Create a copy of the AnnData object if requested
    adata = adata.copy() if copy else adata

    # Set random seed for reproducibility
    if random_state is not None:
        np.random.seed(random_state)
    
    # Get variable names (gene names) from the appropriate layer
    var_names = adata.raw.var_names if use_raw else adata.var_names
    
    # Ensure gene_list is a pandas Index object
    gene_list = pd.Index([gene_list] if isinstance(gene_list, str) else gene_list)
    
    # Identify genes in gene_list that are not in var_names
    genes_to_ignore = gene_list.difference(var_names, sort=False)

    # Keep only genes that are in both gene_list and var_names
    gene_list = gene_list.intersection(var_names)
    if len(gene_list) == 0:
        raise ValueError("No valid genes were passed for scoring.")
    
    # Set up gene pool for control gene selection
    gene_pool = pd.Index(var_names, dtype="string") if gene_pool is None else pd.Index(gene_pool, dtype="string").intersection(var_names)
    if len(gene_pool) == 0:
        raise ValueError("No valid genes were passed for reference set.")
    
    # Get the appropriate data matrix
    _adata = adata.raw if use_raw else adata
    X = _adata.X
    
    # Ensure X is in CSR format if sparse
    if sparse.issparse(X):
        X = X.tocsr()
    
    # Calculate average gene expression across cells
    if sparse.issparse(X):
        obs_avg = pd.Series(np.array(X.mean(axis=0)).flatten(), index=gene_pool)
    else:
        obs_avg = pd.Series(np.nanmean(X, axis=0), index=gene_pool)
    
    # Remove genes with non-finite average expression
    obs_avg = obs_avg[np.isfinite(obs_avg)]
    
    # Bin genes based on their average expression
    n_items = int(np.round(len(obs_avg) / (n_bins - 1)))
    obs_cut = obs_avg.rank(method="min") // n_items
    
    # Select control genes
    control_genes = pd.Index([], dtype="string")
    for cut in np.unique(obs_cut.loc[gene_list]):
        r_genes = obs_cut[obs_cut == cut].index
        r_genes = r_genes.to_series().sample(ctrl_size).index if ctrl_size < len(r_genes) else r_genes
        control_genes = control_genes.union(r_genes.difference(gene_list))
    
    # Set up gene weights
    if gene_weights is None:
        weights = np.ones(len(gene_list))
    else:
        gene_list = gene_list.tolist()
        weights = np.array([gene_weights[gene_list.index(gene)] if gene in gene_list else 0 for gene in gene_list])
        if np.any(weights == 0):
            warnings.warn("Some genes in gene_list were not assigned a weight and will be ignored.")
    
    # Get indices of genes in gene_list
    gene_indices = var_names.get_indexer(gene_list)

    # Calculate gene variances and use them to adjust weights
    if scale_by_variance:
        gene_variances = calculate_gene_variances(X, gene_indices)
        variance_scaling = 1 / np.sqrt(gene_variances + 1e-8)  # Add small constant to avoid division by zero
        
        if gene_weights is None:
            weights = variance_scaling
        else:
            weights = np.array([gene_weights[gene_list.index(gene)] if gene in gene_list else 1 for gene in gene_list])
            weights *= variance_scaling
    else:
        if gene_weights is None:
            weights = np.ones(len(gene_list))
        else:
            weights = np.array([gene_weights[gene_list.index(gene)] if gene in gene_list else 1 for gene in gene_list])
    
    # Normalize weights
    if normalize_weights:
        weights /= np.sum(weights)

    # Calculate gene set scores
    if control:
        cut_indices = obs_cut.iloc[var_names.get_indexer(gene_pool)].values
        
        score = np.zeros(X.shape[0])

        # Set up chunking for large datasets
        if disable_chunking:
            chunk_starts = [0]
            chunk_ends = [X.shape[0]]
        else:
            chunk_starts = range(0, X.shape[0], chunk_size)
            chunk_ends = [min(start + chunk_size, X.shape[0]) for start in chunk_starts]

        for chunk_start, chunk_end in zip(chunk_starts, chunk_ends):
            
            X_chunk = X[chunk_start:chunk_end]
            
            if gpu:
                # GPU-accelerated calculation
                X_chunk = cp.asarray(X_chunk.toarray() if sparse.issparse(X_chunk) else X_chunk)
                gene_indices_gpu = cp.asarray(gene_indices)
                cut_indices_gpu = cp.asarray(cut_indices)
                weights_gpu = cp.asarray(weights)
                
                empty_mask_count = 0
                total_genes = len(gene_indices_gpu)

                control_avgs = cp.zeros((len(gene_indices), X_chunk.shape[0]))
                for i, gene_idx in enumerate(gene_indices_gpu):
                    mask = cut_indices_gpu == cut_indices_gpu[gene_idx]
                    mask[gene_idx] = False
                    if cp.any(mask):
                        control_avgs[i] = cp.nanmean(X_chunk[:, mask], axis=1)
                    else:
                        empty_mask_count += 1
                        # Use global mean as fallback
                        control_avgs[i] = cp.mean(X_chunk[:, gene_idx])

                if empty_mask_count > 0:
                    warnings.warn(f"{empty_mask_count} out of {total_genes} genes had no control genes in their expression bin. Results may be less reliable.")

                gene_avgs = X_chunk[:, gene_indices_gpu].T
                diff = gene_avgs - control_avgs
                if abs_diff:
                    diff = cp.abs(diff)
                
                if weighted:
                    chunk_score = cp.sum(diff * weights_gpu.reshape(-1, 1), axis=0)
                else:
                    chunk_score = cp.sum(diff, axis=0)
                
                score[chunk_start:chunk_end] = cp.asnumpy(chunk_score)
            else:
                # CPU calculation
                if sparse.issparse(X_chunk):
                    control_avgs = np.zeros((len(gene_indices), X_chunk.shape[0]))
                    for i, gene_idx in enumerate(gene_indices):
                        mask = cut_indices == cut_indices[gene_idx]
                        mask[gene_idx] = False
                        control_avgs[i] = X_chunk[:, mask].mean(axis=1).A1
                else:
                    control_avgs = fast_control_avgs(X_chunk, gene_indices, cut_indices, X_chunk.shape[0])
                
                gene_avgs = X_chunk[:, gene_indices].T.toarray() if sparse.issparse(X_chunk) else X_chunk[:, gene_indices].T
                diff = gene_avgs - control_avgs
                if abs_diff:
                    diff = np.abs(diff)
                
                if weighted:
                    chunk_score = np.sum(diff * weights.reshape(-1, 1), axis=0)
                else:
                    chunk_score = np.sum(diff, axis=0)
                
                score[chunk_start:chunk_end] = chunk_score
    else:
        # Calculate scores without control genes
        score = np.zeros(X.shape[0])
        for chunk_start in range(0, X.shape[0], chunk_size):
            chunk_end = min(chunk_start + chunk_size, X.shape[0])
            X_chunk = X[chunk_start:chunk_end]
            
            if gpu:
                # GPU-accelerated calculation
                X_chunk = cp.asarray(X_chunk.toarray() if sparse.issparse(X_chunk) else X_chunk)
                gene_indices_gpu = cp.asarray(gene_indices)
                weights_gpu = cp.asarray(weights)
                
                gene_avgs = X_chunk[:, gene_indices_gpu].T
                if weighted:
                    chunk_score = cp.sum(gene_avgs * weights_gpu.reshape(-1, 1), axis=0)
                else:
                    chunk_score = cp.sum(gene_avgs, axis=0)
                
                score[chunk_start:chunk_end] = cp.asnumpy(chunk_score)
            else:
                # CPU calculation
                gene_avgs = X_chunk[:, gene_indices].T.toarray() if sparse.issparse(X_chunk) else X_chunk[:, gene_indices].T
                if weighted:
                    chunk_score = np.sum(gene_avgs * weights.reshape(-1, 1), axis=0)
                else:
                    chunk_score = np.sum(gene_avgs, axis=0)
                
                score[chunk_start:chunk_end] = chunk_score
    
    # Create a new DataFrame with all columns at once
    result_df = pd.DataFrame({
        'obs_names': adata.obs_names,
        score_name: score
    })
    
    if return_scores:
        return result_df.set_index('obs_names')[score_name]
    else:
        # Update adata.obs with the new scores
        if score_name in adata.obs.columns:
            adata.obs = adata.obs.drop(columns=[score_name])
        adata.obs = adata.obs.join(result_df.set_index('obs_names'))
        return adata if copy else None

#######################################################################################
def score_genes_weighted_with_labels(
    adata,
    gene_list,
    gene_weights=None,
    ctrl_size=50,
    gene_pool=None,
    n_bins=25,
    score_name="score",
    random_state=0,
    copy=False,
    use_raw=False,
    return_scores=False,
    control=True,
    weighted=True,
    abs_diff=False,
    gpu=None,
    chunk_size=10000,
    disable_chunking=True,
    scale_by_variance=False,
    normalize_weights=False,
    labels=None,
    control_condition=None,
    debug=False,
    scaling_only_based_on_control=True
):

    if debug:
        print(f"Initial gene_list length: {len(gene_list)}")
        print(f"Initial gene_weights length: {len(gene_weights) if gene_weights is not None else 'None'}")
    
    # Check GPU availability
    if gpu!=False:
        try:
            import cupy as cp
            gpu = True
        except ImportError:
            gpu = False
        except cp.cuda.runtime.CUDARuntimeError:
            gpu = False

    # Create a copy of the AnnData object if requested
    adata = adata.copy() if copy else adata

    # Set random seed for reproducibility
    if random_state is not None:
        np.random.seed(random_state)
    
    # Get variable names (gene names) from the appropriate layer
    var_names = adata.raw.var_names if use_raw else adata.var_names
    
    # Ensure gene_list is a pandas Index object
    gene_list = pd.Index(gene_list)
    
    # Keep only genes that are in both gene_list and var_names
    valid_genes = gene_list.intersection(var_names)
    removed_genes = gene_list.difference(valid_genes)
    if debug: 
        print(f"Number of valid genes: {len(valid_genes)}")
        print(f"Removed genes: {list(removed_genes)}")
    if len(valid_genes) == 0:
        raise ValueError("No valid genes were passed for scoring.")
    
    # Align gene weights with valid genes
    if gene_weights is not None:
        gene_weights_dict = {gene: weight for gene, weight in zip(gene_list, gene_weights) if gene in valid_genes}
        aligned_weights = np.array([gene_weights_dict[gene] for gene in valid_genes])
        if debug: print(f"Aligned weights length: {len(aligned_weights)}")
    else:
        aligned_weights = np.ones(len(valid_genes))
    
    # Get indices of genes in valid_genes
    gene_indices = var_names.get_indexer(valid_genes)
    if debug: print(f"Number of gene indices: {len(gene_indices)}")
    
    # Set up gene pool for control gene selection
    gene_pool = pd.Index(var_names, dtype="string") if gene_pool is None else pd.Index(gene_pool, dtype="string").intersection(var_names)
    if len(gene_pool) == 0:
        raise ValueError("No valid genes were passed for reference set.")
    
    # Get the appropriate data matrix
    _adata = adata.raw if use_raw else adata
    X = _adata.X
    
    # Ensure X is in CSR format if sparse
    if sparse.issparse(X):
        X = X.tocsr()
        
    # Check if labels are provided and valid
    if labels is not None:
        if labels not in adata.obs.columns:
            raise ValueError(f"Label column '{labels}' not found in adata.obs")
        if control_condition is None:
            raise ValueError("Control condition must be specified when using labels")
        if control_condition not in adata.obs[labels].unique():
            raise ValueError(f"Control condition '{control_condition}' not found in the label column")
        
        # Create a mask for control cells
        control_mask = adata.obs[labels] == control_condition
    else:
        control_mask = np.ones(X.shape[0], dtype=bool)
    
    # Calculate average gene expression across control cells
    if sparse.issparse(X):
        obs_avg = pd.Series(np.array(X[control_mask].mean(axis=0)).flatten(), index=gene_pool)
    else:
        obs_avg = pd.Series(np.nanmean(X[control_mask], axis=0), index=gene_pool)
    
    # Remove genes with non-finite average expression
    obs_avg = obs_avg[np.isfinite(obs_avg)]
    
    # Bin genes based on their average expression in control cells
    n_items = int(np.round(len(obs_avg) / (n_bins - 1)))
    obs_cut = obs_avg.rank(method="min") // n_items
    
    # Select control genes based on control cells
    if control:
        control_genes = pd.Index([], dtype="string")
        for cut in np.unique(obs_cut.loc[valid_genes]):
            r_genes = obs_cut[obs_cut == cut].index
            r_genes = r_genes.to_series().sample(ctrl_size).index if ctrl_size < len(r_genes) else r_genes
            control_genes = control_genes.union(r_genes.difference(valid_genes))
        control_gene_indices = var_names.get_indexer(control_genes)
    
    # Initialize the score array
    score = np.zeros(X.shape[0])

    # Get unique conditions
    conditions = [None] if labels is None else adata.obs[labels].unique()

    # Calculate gene variances for the control condition if scaling_only_based_on_control is True
    if scale_by_variance and scaling_only_based_on_control:
        control_gene_variances = calculate_gene_variances(X[control_mask], gene_indices)
        control_variance_scaling = 1 / np.sqrt(control_gene_variances + 1e-8)

    for condition in conditions:
        if labels is not None:
            condition_mask = adata.obs[labels] == condition
        else:
            condition_mask = np.ones(X.shape[0], dtype=bool)
        
        X_condition = X[condition_mask]
        
        # Normalize weights for the current condition
        if normalize_weights:
            aligned_weights /= np.sum(aligned_weights)

        # Calculate gene variances for the current condition or use control condition variances
        if scale_by_variance:
            if scaling_only_based_on_control:
                variance_scaling = control_variance_scaling
            else:
                gene_variances = calculate_gene_variances(X_condition, gene_indices)
                variance_scaling = 1 / np.sqrt(gene_variances + 1e-8)
            
            if debug: 
                print(f"Variance scaling shape: {variance_scaling.shape}")
                print(f"Aligned weights shape: {aligned_weights.shape}")
            condition_weights = aligned_weights * variance_scaling
        else:
            condition_weights = aligned_weights

        if debug: print(f"Condition weights shape: {condition_weights.shape}")

        # Calculate gene set scores for the current condition
        condition_score = np.zeros(X_condition.shape[0])

        # Set up chunking for large datasets
        if disable_chunking:
            chunk_starts = [0]
            chunk_ends = [X_condition.shape[0]]
        else:
            chunk_starts = range(0, X_condition.shape[0], chunk_size)
            chunk_ends = [min(start + chunk_size, X_condition.shape[0]) for start in chunk_starts]

        for chunk_start, chunk_end in zip(chunk_starts, chunk_ends):
            X_chunk = X_condition[chunk_start:chunk_end]
            
            if gpu:
                # GPU-accelerated calculation
                X_chunk = cp.asarray(X_chunk.toarray() if sparse.issparse(X_chunk) else X_chunk)
                gene_indices_gpu = cp.asarray(gene_indices)
                weights_gpu = cp.asarray(condition_weights)
                
                gene_avgs = cp.mean(X_chunk[:, gene_indices_gpu], axis=1)
                
                if control:
                    control_gene_indices_gpu = cp.asarray(control_gene_indices)
                    control_avgs = cp.mean(X_chunk[:, control_gene_indices_gpu], axis=1)
                    diff = gene_avgs - control_avgs
                else:
                    diff = gene_avgs
                
                if abs_diff:
                    diff = cp.abs(diff)
                
                if weighted:
                    chunk_score = cp.sum(diff.reshape(-1, 1) * weights_gpu.reshape(1, -1), axis=1)
                else:
                    chunk_score = cp.sum(diff.reshape(-1, 1), axis=1)
                
                condition_score[chunk_start:chunk_end] = cp.asnumpy(chunk_score)
            else:
                # CPU calculation
                gene_avgs = np.mean(X_chunk[:, gene_indices].toarray() if sparse.issparse(X_chunk) else X_chunk[:, gene_indices], axis=1)
                
                if control:
                    control_avgs = np.mean(X_chunk[:, control_gene_indices].toarray() if sparse.issparse(X_chunk) else X_chunk[:, control_gene_indices], axis=1)
                    diff = gene_avgs - control_avgs
                else:
                    diff = gene_avgs
                
                if abs_diff:
                    diff = np.abs(diff)
                
                if weighted:
                    chunk_score = np.sum(diff.reshape(-1, 1) * condition_weights.reshape(1, -1), axis=1)
                else:
                    chunk_score = np.sum(diff.reshape(-1, 1), axis=1)
                
                condition_score[chunk_start:chunk_end] = chunk_score

        # Assign the condition_score to the corresponding cells in the main score array
        score[condition_mask] = condition_score
    
    # Create a new DataFrame with all columns at once
    result_df = pd.DataFrame({
        'obs_names': adata.obs_names,
        score_name: score
    })
    
    if return_scores:
        return result_df.set_index('obs_names')[score_name]
    else:
        # Update adata.obs with the new scores
        if score_name in adata.obs.columns:
            adata.obs = adata.obs.drop(columns=[score_name])
        adata.obs = adata.obs.join(result_df.set_index('obs_names'))
        return adata if copy else None
