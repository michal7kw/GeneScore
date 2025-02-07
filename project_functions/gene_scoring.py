import numpy as np
import pandas as pd
from scipy import sparse
import warnings
from typing import Optional, Union, List

try:
    import cupy as cp
except ImportError:
    cp = None

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
        
        # Initialize variances array with zeros
        variances = np.zeros(X.shape[1])
        
        # Calculate variance only for genes with sufficient non-NaN values
        mask_multiple = non_nan_count > 1
        mask_single = non_nan_count == 1
        
        if np.any(mask_multiple):
            variances[mask_multiple] = np.nanvar(X[:, mask_multiple], axis=0, ddof=1)
        if np.any(mask_single):
            variances[mask_single] = np.nanvar(X[:, mask_single], axis=0, ddof=0)
        # Genes with no valid values remain zero

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
    used_layer=None,
    return_scores=False,
    control=True,
    weighted=True,
    abs_diff=False,
    gpu=True,
    chunk_size=10000,
    disable_chunking=True,
    scale_by_variance=False,
    normalize_weights=False,
    conditions_labels=None,
    control_condition=None,
    debug=False,
    scaling_only_based_on_control=False
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
        # print("Warnig: Using 'X' attribute, make sure that you are using the right layer")
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
    
    # Get indices of genes in gene_list
    gene_indices = var_names.get_indexer(gene_list)

    # Sample control genes for each target gene
    control_indices = []
    for gene_idx in gene_indices:
        gene_cut = obs_cut.iloc[gene_idx]
        bin_genes = obs_cut[obs_cut == gene_cut].index
        bin_genes = bin_genes.drop(var_names[gene_idx])  # Remove the target gene
        if len(bin_genes) > ctrl_size:
            control_genes = np.random.choice(bin_genes, size=ctrl_size, replace=False)
        else:
            control_genes = bin_genes
        control_indices.append(var_names.get_indexer(control_genes))

    # Set up gene weights
    if gene_weights is None:
        weights = np.ones(len(gene_list))
    else:
        gene_list = gene_list.tolist()
        weights = np.array([gene_weights[gene_list.index(gene)] if gene in gene_list else 0 for gene in gene_list])
        if np.any(weights == 0):
            warnings.warn("Some genes in gene_list were not assigned a weight and will be ignored.")

    # Initialize the score array
    score = np.zeros(X.shape[0])

    # Get unique conditions
    conditions = [None] if conditions_labels is None else adata.obs[conditions_labels].unique()

    # Calculate gene variances for the control condition if scaling_only_based_on_control is True
    if scale_by_variance and scaling_only_based_on_control:
        control_mask = adata.obs[conditions_labels] == control_condition
        control_gene_variances = calculate_gene_variances(X[control_mask], gene_indices)
        control_variance_scaling = 1 / np.sqrt(control_gene_variances + 1e-8)
        # Replace infinite values with 1 (no scaling)
        control_variance_scaling[~np.isfinite(control_variance_scaling)] = 1

    scores = {}
    for condition in conditions:
        if conditions_labels is not None:
            condition_mask = adata.obs[conditions_labels] == condition
        else:
            condition_mask = np.ones(X.shape[0], dtype=bool)
    
        X_condition = X[condition_mask]

        # Calculate gene variances for the current condition or use control condition variances
        if scale_by_variance:
            if scaling_only_based_on_control:
                variance_scaling = control_variance_scaling
            else:
                gene_variances = calculate_gene_variances(X_condition, gene_indices)
                variance_scaling = 1 / np.sqrt(gene_variances + 1e-8)
                # Replace infinite values with 1 (no scaling)
                variance_scaling[~np.isfinite(variance_scaling)] = 1
            condition_weights = weights * variance_scaling
        else:
            condition_weights = weights
    
        # Normalize weights for the current condition
        if normalize_weights:
            abs_sum_weights = np.sum(np.abs(condition_weights))
            if abs_sum_weights > 0:
                condition_weights /= abs_sum_weights
            else:
                warnings.warn("Sum of absolute weights is zero. Skipping normalization.")

        # Calculate gene set scores
        if control:
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
                    condition_weights_gpu = cp.asarray(condition_weights)
                    
                    chunk_score = cp.zeros(X_chunk.shape[0])
                    for i, (gene_idx, ctrl_indices) in enumerate(zip(gene_indices, control_indices)):
                        gene_expr = X_chunk[:, gene_idx]
                        ctrl_expr = X_chunk[:, ctrl_indices]
                        ctrl_avg = cp.nanmean(ctrl_expr, axis=1)
                        diff = gene_expr - ctrl_avg
                        if abs_diff:
                            diff = cp.abs(diff)
                        if weighted:
                            chunk_score += diff * condition_weights_gpu[i]
                        else:
                            chunk_score += diff
                    
                    condition_score[chunk_start:chunk_end] = cp.asnumpy(chunk_score)
                else:
                    # CPU calculation
                    chunk_score = np.zeros(X_chunk.shape[0])
                    for i, (gene_idx, ctrl_indices) in enumerate(zip(gene_indices, control_indices)):
                        gene_expr = X_chunk[:, gene_idx].toarray().flatten() if sparse.issparse(X_chunk) else X_chunk[:, gene_idx]
                        ctrl_expr = X_chunk[:, ctrl_indices].toarray() if sparse.issparse(X_chunk) else X_chunk[:, ctrl_indices]
                        ctrl_avg = np.nanmean(ctrl_expr, axis=1)
                        diff = gene_expr - ctrl_avg
                        if abs_diff:
                            diff = np.abs(diff)
                        if weighted:
                            chunk_score += diff * condition_weights[i]
                        else:
                            chunk_score += diff
                    
                    condition_score[chunk_start:chunk_end] = chunk_score
        else:
            # Calculate scores without control genes
            condition_score = np.zeros(X_condition.shape[0])
            for chunk_start in range(0, X_condition.shape[0], chunk_size):
                chunk_end = min(chunk_start + chunk_size, X_condition.shape[0])
                X_chunk = X_condition[chunk_start:chunk_end]
                
                if gpu:
                    # GPU-accelerated calculation
                    X_chunk = cp.asarray(X_chunk.toarray() if sparse.issparse(X_chunk) else X_chunk)
                    gene_indices_gpu = cp.asarray(gene_indices)
                    condition_weights_gpu = cp.asarray(condition_weights)
                    
                    gene_avgs = X_chunk[:, gene_indices_gpu].T
                    if weighted:
                        chunk_score = cp.sum(gene_avgs * condition_weights_gpu.reshape(-1, 1), axis=0)
                    else:
                        chunk_score = cp.sum(gene_avgs, axis=0)
                    
                    condition_score[chunk_start:chunk_end] = cp.asnumpy(chunk_score)
                else:
                    # CPU calculation
                    gene_avgs = X_chunk[:, gene_indices].T.toarray() if sparse.issparse(X_chunk) else X_chunk[:, gene_indices].T
                    if weighted:
                        chunk_score = np.sum(gene_avgs * condition_weights.reshape(-1, 1), axis=0)
                    else:
                        chunk_score = np.sum(gene_avgs, axis=0)
                    
                    condition_score[chunk_start:chunk_end] = chunk_score
        
        score[condition_mask] = condition_score

    # Create a new Series with the scores
    result_series = pd.Series(score, index=adata.obs_names, name=score_name)
    
    if return_scores:
        return result_series
    else:
        # Update adata.obs with the new scores
        if score_name in adata.obs.columns:
            # Create a copy of the dataframe without the score column
            adata.obs = adata.obs.drop(columns=[score_name])
        # Use more efficient assignment that avoids fragmentation
        adata.obs = adata.obs.assign(**{score_name: result_series})
        return adata if copy else None