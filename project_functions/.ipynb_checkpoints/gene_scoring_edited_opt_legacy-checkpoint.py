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

# def calculate_gene_variances(X, gene_indices):
#     """
#     Calculate the variance of each gene in the gene_list.
#     """
#     if sparse.issparse(X):
#         variances = X.power(2).mean(axis=0).A1 - np.power(X.mean(axis=0).A1, 2)
#     else:
#         variances = np.nanvar(X, axis=0)
#     return variances[gene_indices]

def calculate_gene_variances(X, gene_indices):
    """
    Calculate the variance of each gene in the gene_list, handling edge cases.
    """
    if sparse.issparse(X):
        mean = X.mean(axis=0).A1
        mean_sq = X.power(2).mean(axis=0).A1
        variances = mean_sq - np.power(mean, 2)
    else:
        # Count non-NaN values for each gene
        non_nan_count = np.sum(~np.isnan(X), axis=0)
        
        # For genes with more than one non-NaN value, use ddof=1
        # For genes with only one non-NaN value, use ddof=0
        # For genes with no non-NaN values, the result will be NaN
        variances = np.where(non_nan_count > 1,
                             np.nanvar(X, axis=0, ddof=1),
                             np.nanvar(X, axis=0, ddof=0))

    # Replace negative variances (due to numerical issues) with zero
    variances = np.maximum(variances, 0)
    
    # Handle cases where variance is undefined (e.g., all values are NaN)
    variances = np.nan_to_num(variances, nan=0.0)
    
    return variances[gene_indices]

def score_genes(
    adata,
    gene_list,
    gene_weights=None,
    score_name="score",
    gene_pool=None,
    n_bins=25,
    random_state=0,
    copy=False,
    used_layer='cpm',
    return_scores=False,
    control=True,
    weighted=True,
    abs_diff=False,
    gpu=True,
    chunk_size=10000,
    disable_chunking=True,
    normalize_weights=False,
    conditions_labels='Condition',
    control_condition='DMSO'
):
    
    # Check if CuPy is available for GPU acceleration
    if gpu and cp is None:
        raise ImportError("CuPy is required for GPU acceleration. Install it with 'pip install cupy'.")

    # Create a copy of the AnnData object if requested
    adata = adata.copy() if copy else adata

    # Set random seed for reproducibility
    if random_state is not None:
        np.random.seed(random_state)
    
    # Get variable names (gene names) from the appropriate layer
    # var_names = adata.raw.var_names if use_raw else adata.var_names
    var_names = adata.var_names
    
    # Ensure gene_list is a pandas Index object
    gene_list = pd.Index([gene_list] if isinstance(gene_list, str) else gene_list)
    
    # Keep only genes that are in both gene_list and var_names
    gene_list = gene_list.intersection(var_names)
    if len(gene_list) == 0:
        raise ValueError("No valid genes were passed for scoring.")
    
    # Set up gene pool for control gene selection
    gene_pool = pd.Index(var_names, dtype="string") if gene_pool is None else pd.Index(gene_pool, dtype="string").intersection(var_names)
    if len(gene_pool) == 0:
        raise ValueError("No valid genes were passed for reference set.")
    
    # Get the appropriate data matrix
    if used_layer == 'cpm':
        X = adata.layers['cpm']
    elif used_layer == 'raw':
        X = adata.layers['counts']
    else:
        X = adata.X
    
    # Ensure X is in CSR format if sparse
    if sparse.issparse(X):
        X = X.tocsr()
    
    # Calculate average gene expression across cells
    if conditions_labels is not None and control_condition is not None:
        control_mask = adata.obs[conditions_labels] == control_condition
        if sparse.issparse(X):
            obs_avg = pd.Series(np.array(X[control_mask].mean(axis=0)).flatten(), index=gene_pool)
        else:
            obs_avg = pd.Series(np.nanmean(X[control_mask], axis=0), index=gene_pool)
    else:
        if sparse.issparse(X):
            obs_avg = pd.Series(np.array(X.mean(axis=0)).flatten(), index=gene_pool)
        else:
            obs_avg = pd.Series(np.nanmean(X, axis=0), index=gene_pool)
    
    # Remove genes with non-finite average expression
    obs_avg = obs_avg[np.isfinite(obs_avg)]
    
    # Bin genes based on their average expression
    n_items = int(np.round(len(obs_avg) / (n_bins - 1)))
    obs_cut = obs_avg.rank(method="min") // n_items
    
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
