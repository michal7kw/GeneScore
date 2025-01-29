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
    ctrl_size=50,
    gene_pool=None,
    n_bins=25,
    random_state=0,
    copy=False,
    used_layer='cpm',
    return_scores=False,
    control=True,
    weighted=True,
    abs_diff=False,
    gpu=None,
    chunk_size=10000,
    disable_chunking=True,
    scale_by_variance=False,
    normalize_weights=False,
    conditions_labels='Condition',
    control_condition='DMSO',
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
    # var_names = adata.raw.var_names if use_raw else adata.var_names
    var_names = adata.var_names
    
    # Ensure gene_list is a pandas Index object
    gene_list = pd.Index([gene_list] if isinstance(gene_list, str) else gene_list)
    
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
    conditions = [None] if conditions_labels is None else adata.obs[conditions_labels].unique()

    # Calculate gene variances for the control condition if scaling_only_based_on_control is True
    if conditions_labels is not None and control_condition is not None and scale_by_variance and scaling_only_based_on_control:
        control_gene_variances = calculate_gene_variances(X[control_mask], gene_indices)
        control_variance_scaling = 1 / np.sqrt(control_gene_variances + 1e-8)

    for condition in conditions:
        if conditions_labels is not None:
            condition_mask = adata.obs[conditions_labels] == condition
        else:
            condition_mask = np.ones(X.shape[0], dtype=bool)
        
        X_condition = X[condition_mask]
        
        # Normalize weights for the current condition
        if normalize_weights:
            aligned_weights /= np.sum(aligned_weights)

        # Calculate gene variances for the current condition or use control condition variances
        if scale_by_variance:
            if conditions_labels is not None and control_condition is not None and scaling_only_based_on_control:
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

