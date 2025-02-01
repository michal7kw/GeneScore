import os
import csv
import matplotlib
import pandas as pd
import anndata as ad
import seaborn as sns
import numpy as np
import warnings
from collections import defaultdict
import matplotlib.pyplot as plt
from upsetplot import from_contents, UpSet
from typing import Optional, Union


import itertools
from itertools import combinations
from matplotlib_venn import venn3, venn2

#### Simulation scores functions ####

def load_GRNs_gene_sets(root_dir: str, gene_set_list: list = ["all_ex"], weights_list: str = "scores_grn_all_from_comb_run_new.csv") -> tuple:
    """
    Load and format gene regulatory network (GRN) data from CSV files.
    
    Args:
        root_dir (str): Root directory containing the gene set data
        gene_set_list (list): List of gene set names to load. Defaults to ["all_ex"]
        weights_list (str): Name of the weights CSV file. Defaults to "scores_grn_all_from_comb_run_new.csv"
    
    Returns:
        tuple: Two dictionaries containing formatted GRN data:
            - First dict organized by gene of interest -> cell type -> targets
            - Second dict organized by cell type -> gene of interest -> targets
            
    Raises:
        FileNotFoundError: If the specified CSV file is not found
        ValueError: If the CSV file is missing required columns
    """
    if not os.path.isdir(root_dir):
        raise ValueError(f"Root directory does not exist: {root_dir}")
        
    if not isinstance(gene_set_list, list):
        raise TypeError("gene_set_list must be a list")
        
    if not isinstance(weights_list, str):
        raise TypeError("weights_list must be a string")
        
    gene_sets = {}
    required_columns = ['source', 'target', 'score', 'coef_mean', 'celltype']

    # Load data for each gene set
    for gene_set in gene_set_list:
        path = os.path.join(root_dir, f"{gene_set}", "celloracle")
        file_path = os.path.join(path, weights_list)
        
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
            
        df = pd.read_csv(file_path)
        
        # Verify required columns exist
        missing_cols = [col for col in required_columns if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns in {file_path}: {missing_cols}")
            
        gene_sets[gene_set] = df

    # Sets Formatting
    gene_sets_dict = {}
    gene_sets_dict_cell_type_first = {}

    for key, value in gene_sets.items():
        gene_sets_dict[key] = {}
        gene_sets_dict_cell_type_first[key] = {}

        for _, row in value.iterrows():
            try:
                goi = str(row['source'])
                target = str(row['target'])
                score1 = float(row['score']) * float(row['coef_mean'])
                score2 = float(row['coef_mean'])
                source = str(row['celltype'])
            except ValueError as e:
                raise ValueError(f"Error converting data types in file for gene set {key}: {str(e)}")

            # Format 1: Gene of interest first
            if goi not in gene_sets_dict[key]:
                gene_sets_dict[key][goi] = {}

            if source not in gene_sets_dict[key][goi]:
                gene_sets_dict[key][goi][source] = {'targets': [], 'scored_coef_mean': [], 'coef_mean': []}

            gene_sets_dict[key][goi][source]['targets'].append(target)
            gene_sets_dict[key][goi][source]['scored_coef_mean'].append(score1)
            gene_sets_dict[key][goi][source]['coef_mean'].append(score2)

            # Format 2: Cell type first
            if source not in gene_sets_dict_cell_type_first[key]:
                gene_sets_dict_cell_type_first[key][source] = {}

            if goi not in gene_sets_dict_cell_type_first[key][source]:
                gene_sets_dict_cell_type_first[key][source][goi] = {'targets': [], 'scored_coef_mean': [], 'coef_mean': []}

            gene_sets_dict_cell_type_first[key][source][goi]['targets'].append(target)
            gene_sets_dict_cell_type_first[key][source][goi]['scored_coef_mean'].append(score1)
            gene_sets_dict_cell_type_first[key][source][goi]['coef_mean'].append(score2)

    return gene_sets_dict, gene_sets_dict_cell_type_first

def remove_duplicates_preserve_order_GRNs(data_dict: dict) -> dict:
    """
    Remove duplicate targets while preserving order and keeping the highest coefficient mean.
    
    Args:
        data_dict (dict): Dictionary containing GRN data
        
    Returns:
        dict: Dictionary with duplicates removed while preserving order
        
    Raises:
        TypeError: If input is not a dictionary
        ValueError: If dictionary structure is invalid
    """
    if not isinstance(data_dict, dict):
        raise TypeError("Input must be a dictionary")
        
    result = {}
    for set_selected, set_data in data_dict.items():
        if not isinstance(set_data, dict):
            raise ValueError(f"Invalid structure for set {set_selected}")
            
        result[set_selected] = {}
        for cell_type_selected, cell_type_data in set_data.items():
            if not isinstance(cell_type_data, dict):
                raise ValueError(f"Invalid structure for cell type {cell_type_selected} in set {set_selected}")
                
            result[set_selected][cell_type_selected] = {}
            for scored_gene_selected, gene_data in cell_type_data.items():
                if not isinstance(gene_data, dict) or not all(k in gene_data for k in ['targets', 'coef_mean', 'scored_coef_mean']):
                    raise ValueError(f"Invalid gene data structure for {scored_gene_selected}")
                    
                targets = gene_data['targets']
                coef_mean = gene_data['coef_mean']
                scored_coef_mean = gene_data['scored_coef_mean']

                # Create a dictionary to store the highest coef_mean for each target
                target_dict = {}
                for i, target in enumerate(targets):
                    if target not in target_dict or coef_mean[i] > target_dict[target][1]:
                        target_dict[target] = (i, coef_mean[i])

                # Create new lists without duplicates
                new_targets = []
                new_coef_mean = []
                new_scored_coef_mean = []
                for target, (index, _) in sorted(target_dict.items(), key=lambda x: x[1][0]):
                    new_targets.append(target)
                    new_coef_mean.append(coef_mean[index])
                    new_scored_coef_mean.append(scored_coef_mean[index])

                # Update the result dictionary
                result[set_selected][cell_type_selected][scored_gene_selected] = {
                    'targets': new_targets,
                    'coef_mean': new_coef_mean,
                    'scored_coef_mean': new_scored_coef_mean
                }

    return result

def process_gene_sets(gene_sets: dict) -> dict:
    """
    Process gene sets data and organize it by cell type and gene of interest.
    
    Args:
        gene_sets (dict): Dictionary containing gene set data, where each value is a pandas DataFrame
                         with required columns: ['local_cell_type', 'cell_type', 'goi', 'gene', 'log_fold_change']
    
    Returns:
        dict: Nested dictionary organized as:
              {gene_set -> cell_type -> gene_of_interest -> {'targets': list, 'coef_mean': list}}
              
    Raises:
        TypeError: If gene_sets is not a dictionary
        ValueError: If any DataFrame is missing required columns
        ValueError: If any DataFrame is empty
    """
    if not isinstance(gene_sets, dict):
        raise TypeError("gene_sets must be a dictionary")
        
    required_columns = ['local_cell_type', 'cell_type', 'goi', 'gene', 'log_fold_change']
    gene_sets_dict_cell_type_first = defaultdict(lambda: defaultdict(dict))

    for gene_set, df in gene_sets.items():
        if not isinstance(df, pd.DataFrame):
            raise TypeError(f"Value for gene_set '{gene_set}' must be a pandas DataFrame")
            
        if df.empty:
            raise ValueError(f"DataFrame for gene_set '{gene_set}' is empty")
            
        # Check for required columns
        missing_cols = [col for col in required_columns if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns in gene_set '{gene_set}': {missing_cols}")

        df = df.copy()
        # Filter rows where local_cell_type matches cell_type
        df = df[df['local_cell_type'] == df['cell_type']]
        
        if df.empty:
            continue  # Skip if no matching cell types
        
        # Convert to categorical for faster grouping
        df['goi'] = pd.Categorical(df['goi'])
        df['cell_type'] = pd.Categorical(df['cell_type'])
        
        # Group by goi and cell_type
        grouped = df.groupby(['goi', 'cell_type'], observed=True)
        
        for (goi, cell_type), group in grouped:
            if group.empty:
                continue  # Skip empty groups
                
            # Sort by absolute value of log_fold_change in ascending order
            sorted_group = group.reindex(group['log_fold_change'].abs().sort_values(ascending=False).index)
            
            targets = sorted_group['gene'].tolist()
            weights = sorted_group['log_fold_change'].tolist()
            
            gene_sets_dict_cell_type_first[gene_set][cell_type][goi] = {
                'targets': targets,
                'coef_mean': weights
            }

    return dict(gene_sets_dict_cell_type_first)  # Convert defaultdict to regular dict

def boxplot_Reference_GRN_scores_local(
    adata: ad.AnnData,
    control: str,
    control_condition: str,
    condition: list,
    gene_set: str,
    cell_type: str,
    goi: str,
    normalize_weights: Union[bool, str],
    scale_by_variance: Union[bool, str],
    scaling_only_based_on_control: Union[bool, str],
    prefix: str = ""
) -> plt.Figure:
    """
    Create a boxplot comparing gene scores across different conditions.
    
    Args:
        adata (AnnData): Annotated data matrix containing gene expression data
        control (str): Name of the control condition
        control_condition (str): Control condition identifier
        condition (list): List of 3 condition names to compare
        gene_set (str): Name of the gene set
        cell_type (str): Cell type to analyze
        goi (str): Gene of interest
        normalize_weights (Union[bool, str]): Whether weights are normalized
        scale_by_variance (Union[bool, str]): Whether to scale by variance
        scaling_only_based_on_control (Union[bool, str]): Whether scaling is based only on control
        prefix (str, optional): Prefix for the selection string. Defaults to ""
        
    Returns:
        matplotlib.figure.Figure: Generated boxplot figure, or None if no data available
    """
    # Input validation
    if not isinstance(adata, ad.AnnData):
        raise TypeError("adata must be an AnnData object")
    if not isinstance(condition, list) or len(condition) != 3:
        raise ValueError("condition must be a list containing exactly 3 conditions")
    if not all(isinstance(x, str) for x in [control, control_condition, gene_set, cell_type, goi, prefix]):
        raise TypeError("control, control_condition, gene_set, cell_type, goi and prefix must be strings")

    # Convert string booleans to actual booleans
    def str_to_bool(value):
        if isinstance(value, bool):
            return value
        if isinstance(value, str):
            return value.lower() == 'true'
        raise TypeError(f"Expected bool or str, got {type(value)}")

    normalize_weights = str_to_bool(normalize_weights)
    scale_by_variance = str_to_bool(scale_by_variance)
    scaling_only_based_on_control = str_to_bool(scaling_only_based_on_control)

    if 'sample_type' not in adata.obs.columns:
        raise ValueError("adata.obs must contain 'sample_type' column")

    plt.close('all')

    mask = adata.obs['sample_type'].isin(condition)
    adata_filtered = adata[mask]

    if len(adata_filtered) > 0:
        fig, ax = plt.subplots(figsize=(5, 6))
        selection = (
                    f'{prefix}'
                    f'gene_score_{gene_set}_{cell_type}_{goi}_{control}_'
                    f'normalized_{normalize_weights}_'
                    f'scaled_{scale_by_variance}_'
                    f'cc_{control_condition}_'
                    f'sc_{scaling_only_based_on_control}'
                    )
                    
        if selection not in adata_filtered.obs.columns:
            raise ValueError(f"Selection column '{selection}' not found in adata.obs")
            
        gene_scores_dmso = adata_filtered[adata_filtered.obs['sample_type'] == condition[0]].obs[selection].values
        gene_scores_ag = adata_filtered[adata_filtered.obs['sample_type'] == condition[1]].obs[selection].values
        gene_scores_inh = adata_filtered[adata_filtered.obs['sample_type'] == condition[2]].obs[selection].values
                
        sns.boxplot(data=[gene_scores_dmso, gene_scores_ag, gene_scores_inh], notch=False,
                    boxprops=dict(alpha=0.5),
                    ax=ax)
        sns.stripplot(data=[gene_scores_dmso, gene_scores_ag, gene_scores_inh], 
                        jitter=True, color=".3", linewidth=1, ax=ax)
        
        ax.set_xticks(range(3))
        ax.set_xticklabels([condition[0], condition[1], condition[2]], fontsize=12)
        ax.set_title(f'Gene Scores - {goi}\n cell_type: {cell_type}, control: {control}', fontsize=16)
        ax.set_ylabel("Gene Score", fontsize=12)

        plt.tight_layout()
        return fig
    else:
        print(f'No data to plot for the selected condition: {condition}')
        return None

#### Gene Set Analysis ###

def create_heatmap(data: dict, set_selected: str, scored_genes: list, cell_types: list) -> None:
    """
    Create a heatmap showing the number of target genes for each scored gene across cell types.
    
    Args:
        data (dict): Nested dictionary containing gene regulatory network data
        set_selected (str): Name of the selected gene set
        scored_genes (list): List of scored genes to analyze
        cell_types (list): List of cell types to analyze
        
    Returns:
        None: Displays a heatmap plot
        
    Raises:
        ValueError: If required data is missing or invalid
        KeyError: If specified keys are not found in data structure
    """
    if not isinstance(data, dict) or not data:
        raise ValueError("data must be a non-empty dictionary")
    if not isinstance(set_selected, str):
        raise TypeError("set_selected must be a string")
    if not isinstance(scored_genes, list) or not scored_genes:
        raise ValueError("scored_genes must be a non-empty list")
    if not isinstance(cell_types, list) or not cell_types:
        raise ValueError("cell_types must be a non-empty list")
        
    if set_selected not in data:
        raise KeyError(f"Selected set '{set_selected}' not found in data")
        
    gene_counts = {}
    for scored_gene in scored_genes:
        if scored_gene not in data[set_selected][cell_types[0]]:
            raise KeyError(f"Scored gene '{scored_gene}' not found in data")
            
        gene_counts[scored_gene] = {}
        for cell_type in cell_types:
            if cell_type not in data[set_selected]:
                raise KeyError(f"Cell type '{cell_type}' not found in data")
                
            try:
                gene_counts[scored_gene][cell_type] = len(data[set_selected][cell_type][scored_gene]['targets'])
            except KeyError:
                raise KeyError(f"Missing target data for gene '{scored_gene}' in cell type '{cell_type}'")
    
    df = pd.DataFrame(gene_counts)
    
    plt.figure(figsize=(12, 8))
    sns.heatmap(df, annot=True, fmt='d', cmap='YlOrRd')
    plt.title("Number of Target Genes for Each Scored Gene Across Cell Types")
    plt.ylabel("Cell Types")
    plt.xlabel("Scored Genes")
    plt.tight_layout()
    plt.show()

def create_venn_diagram_gene_set(gene_sets: dict, cell_type: str, mode: str = "positive", printouts: bool = False) -> None:
    """
    Create a Venn diagram showing the intersection of gene sets for a given cell type.
    
    Args:
        gene_sets (dict): Dictionary mapping gene set names to sets of genes
        cell_type (str): Name of the cell type being analyzed
        mode (str, optional): Analysis mode - "positive" or "negative". Defaults to "positive"
        printouts (bool, optional): Whether to print detailed intersection info. Defaults to False
        
    Returns:
        None: Displays a Venn diagram plot
        
    Raises:
        ValueError: If gene_sets is empty or invalid
        TypeError: If arguments are of wrong type
    """
    if not isinstance(gene_sets, dict) or not gene_sets:
        raise ValueError("gene_sets must be a non-empty dictionary")
    if not isinstance(cell_type, str):
        raise TypeError("cell_type must be a string")
    if not isinstance(mode, str) or mode not in ["positive", "negative"]:
        raise ValueError("mode must be either 'positive' or 'negative'")
    if not isinstance(printouts, bool):
        raise TypeError("printouts must be a boolean")

    set_names = list(gene_sets.keys())
    venn_sets = [set(genes) for genes in gene_sets.values()]
    
    if printouts:
        print(f"\nShared genes between the sets:")
        all_sets = set.union(*venn_sets)
        for i in range(2, len(venn_sets) + 1):
            for combination in itertools.combinations(range(len(venn_sets)), i):
                shared_genes = set.intersection(*(venn_sets[j] for j in combination))
                if shared_genes:
                    set_names_str = " & ".join(set_names[j] for j in combination)
                    print(f"Shared between {set_names_str}: {shared_genes}")

    plt.figure(figsize=(10, 10))
    if len(set_names) == 2:
        venn2(venn_sets, set_labels=set_names)
    elif len(set_names) == 3:
        venn3(venn_sets, set_labels=set_names)
    else:
        plt.close()
        print(f"More than 3 sets detected for {cell_type}, coefficients: {mode}. Displaying text-based representation.")
        for i in range(2, min(len(set_names) + 1, 4)):
            for combo in combinations(range(len(set_names)), i):
                common_genes = set.intersection(*[venn_sets[j] for j in combo])
                if common_genes:
                    print(f"Genes common to {', '.join([set_names[j] for j in combo])}: {common_genes}")
        return

    plt.title(f"Genes Intersection for {cell_type}, coefficients: {mode}", fontsize=16)
    plt.show()

def print_number_of_duplicate_genes(gene_sets: dict) -> None:
    """
    Analyze and print information about duplicate genes in each gene set.

    Args:
        gene_sets (dict): Dictionary where keys are set names and values are collections of genes.
            Values can be lists, sets, or any iterable containing genes.

    Returns:
        None: Prints duplicate gene information for each set to stdout.

    Raises:
        TypeError: If gene_sets is not a dictionary
        ValueError: If gene_sets is empty
        TypeError: If any gene set values are not iterable

    Example:
        gene_sets = {
            "set1": ["A", "B", "A", "C"],
            "set2": ["D", "E", "E", "F"]
        }
        print_number_of_duplicate_genes(gene_sets)
        # Output:
        # for set1 total duplicates: 1 out of 4
        # for set2 total duplicates: 1 out of 4
    """
    # Input validation
    if not isinstance(gene_sets, dict):
        raise TypeError("gene_sets must be a dictionary")
    if not gene_sets:
        raise ValueError("gene_sets cannot be empty")

    for set_name, gene_set in gene_sets.items():
        try:
            # Convert to a list to preserve order and ensure iterability
            genes = list(gene_set)
        except TypeError:
            raise TypeError(f"Gene set '{set_name}' must be an iterable collection of genes")

        original_length = len(genes)
        
        # Find duplicates
        seen = set()
        duplicates = []
        for gene in genes:
            if gene in seen:
                duplicates.append(gene)
            else:
                seen.add(gene)
        
        # Print results
        if duplicates:
            print(f"for {set_name} total duplicates: {len(duplicates)} out of {original_length}")

def analyze_gene_sets_gene_set(data: dict, set_selected: str, cell_type: str, scored_genes: list, mode: str = "positive", printouts: bool = False) -> None:
    """
    Analyze gene sets for a given cell type and create Venn diagrams of gene intersections.

    Args:
        data (dict): Dictionary containing gene expression data organized by sets and cell types
        set_selected (str): Name of the gene set to analyze
        cell_type (str): Cell type to analyze
        scored_genes (list): List of genes to analyze
        mode (str, optional): Filter mode - "positive" for positive coefficients, "negative" for negative. Defaults to "positive"
        printouts (bool, optional): Whether to print additional analysis info. Defaults to False

    Raises:
        TypeError: If arguments are not of expected types
        ValueError: If mode is not "positive" or "negative", or if data structure is invalid
        KeyError: If keys are not found in data dictionary
    """
    # Input validation
    if not isinstance(data, dict):
        raise TypeError("data must be a dictionary")
    if not isinstance(set_selected, str):
        raise TypeError("set_selected must be a string")
    if not isinstance(cell_type, str):
        raise TypeError("cell_type must be a string")
    if not isinstance(scored_genes, list):
        raise TypeError("scored_genes must be a list")
    if mode not in ["positive", "negative"]:
        raise ValueError("mode must be either 'positive' or 'negative'")

    # Validate data structure
    if set_selected not in data:
        raise KeyError(f"set_selected '{set_selected}' not found in data")
    if cell_type not in data[set_selected]:
        raise KeyError(f"cell_type '{cell_type}' not found in data[{set_selected}]")

    gene_sets = {}
    for scored_gene in scored_genes:
        if scored_gene not in data[set_selected][cell_type]:
            print(f"Warning: scored_gene '{scored_gene}' not found in data")
            continue

        targets = data[set_selected][cell_type][scored_gene].get('targets', [])
        coefs = data[set_selected][cell_type][scored_gene].get('coef_mean', [])
        
        # Ensure targets and coefs are the same length
        if len(targets) != len(coefs):
            raise ValueError(f"Mismatch in length of targets ({len(targets)}) and coefficients ({len(coefs)}) for {scored_gene}")
        
        # Filter targets based on coefficient signs
        if mode == "positive":
            filtered_targets = [target for target, coef in zip(targets, coefs) if coef > 0]
        else:
            filtered_targets = [target for target, coef in zip(targets, coefs) if coef < 0]
        
        # Only add to gene_sets if there are filtered targets
        if filtered_targets:
            gene_sets[scored_gene] = set(filtered_targets)
    
    if not gene_sets:
        print("Warning: No gene sets were created after filtering")
        return
        
    create_venn_diagram_gene_set(gene_sets, cell_type, mode=mode, printouts=printouts)

def create_venn_diagram_cell_types(gene_sets: dict, scored_gene: str, mode: str = "positive", printouts: bool = False) -> None:
    """
    Create Venn diagrams showing gene intersections between cell types.

    Args:
        gene_sets (dict): Dictionary mapping cell types to sets of genes
        scored_gene (str): Name of the scored gene being analyzed
        mode (str, optional): Analysis mode ("positive" or "negative"). Defaults to "positive"
        printouts (bool, optional): Whether to print additional analysis info. Defaults to False

    Raises:
        TypeError: If arguments are not of expected types
        ValueError: If mode is not "positive" or "negative"
    """
    # Input validation
    if not isinstance(gene_sets, dict):
        raise TypeError("gene_sets must be a dictionary")
    if not isinstance(scored_gene, str):
        raise TypeError("scored_gene must be a string")
    if mode not in ["positive", "negative"]:
        raise ValueError("mode must be either 'positive' or 'negative'")
    
    if not gene_sets:
        print("Warning: Empty gene_sets dictionary provided")
        return

    set_names = list(gene_sets.keys())
    venn_sets = [set(genes) for genes in gene_sets.values()]
    
    if printouts:
        print(f"\nUnique genes for each set:")
        for name, gene_set in zip(set_names, venn_sets):
            unique_genes = gene_set - set.union(*(s for s in venn_sets if s != gene_set))
            print(f"{name}: {unique_genes}")
    
    if len(set_names) == 2:
        plt.figure(figsize=(10, 10))
        venn2(venn_sets, set_labels=set_names)
        plt.title(f"Gene Intersection for {scored_gene}, coefficients: {mode}", fontsize=16)
        plt.show()
    elif len(set_names) == 3:
        plt.figure(figsize=(10, 10))
        venn3(venn_sets, set_labels=set_names)
        plt.title(f"Gene Intersection for {scored_gene}, coefficients: {mode}", fontsize=16)
        plt.show()
    else:
        print(f"More than 3 cell types detected for {scored_gene}, coefficients: {mode}. Displaying text-based representation.")
        for i in range(2, min(len(set_names) + 1, 4)):
            for combo in combinations(range(len(set_names)), i):
                common_genes = set.intersection(*[venn_sets[j] for j in combo])
                if common_genes:
                    print(f"Genes common to {', '.join([set_names[j] for j in combo])}: {common_genes}")
        
def analyze_gene_sets_cell_types(data: dict, set_selected: str, scored_gene: str, cell_types: list, mode: str = "positive", printouts: bool = False) -> None:
    """
    Analyze gene sets across multiple cell types.

    Args:
        data (dict): Dictionary containing gene expression data
        set_selected (str): Name of the gene set to analyze
        scored_gene (str): Name of the scored gene to analyze
        cell_types (list): List of cell types to analyze
        mode (str, optional): Filter mode - "positive" or "negative". Defaults to "positive"
        printouts (bool, optional): Whether to print additional analysis info. Defaults to False

    Raises:
        TypeError: If arguments are not of expected types
        ValueError: If mode is not "positive" or "negative"
        KeyError: If keys are not found in data dictionary
    """
    # Input validation
    if not isinstance(data, dict):
        raise TypeError("data must be a dictionary")
    if not isinstance(set_selected, str):
        raise TypeError("set_selected must be a string")
    if not isinstance(scored_gene, str):
        raise TypeError("scored_gene must be a string")
    if not isinstance(cell_types, list):
        raise TypeError("cell_types must be a list")
    if mode not in ["positive", "negative"]:
        raise ValueError("mode must be either 'positive' or 'negative'")

    gene_sets = {}
    for cell_type in cell_types:
        try:
            targets = data[set_selected][cell_type][scored_gene]['targets']
            coefs = data[set_selected][cell_type][scored_gene]['coef_mean']
        except KeyError as e:
            print(f"Warning: Could not find data for cell_type '{cell_type}': {str(e)}")
            continue
        
        if mode == "positive":
            filtered_targets = [target for target, coef in zip(targets, coefs) if coef > 0]
        else:
            filtered_targets = [target for target, coef in zip(targets, coefs) if coef < 0]
        
        if filtered_targets:
            gene_sets[cell_type] = set(filtered_targets)
    
    if not gene_sets:
        print("Warning: No gene sets were created after filtering")
        return
        
    create_venn_diagram_cell_types(gene_sets, scored_gene, mode=mode, printouts=printouts)

def better_hist_DEGs(data: dict, gene_set: str) -> None:
    """
    Create an enhanced histogram for Differentially Expressed Genes.

    Args:
        data (dict): Dictionary containing gene expression data
        gene_set (str): Name of the gene set to plot

    Raises:
        TypeError: If arguments are not of expected types
        KeyError: If gene_set is not found in data
    """
    if not isinstance(data, dict):
        raise TypeError("data must be a dictionary")
    if not isinstance(gene_set, str):
        raise TypeError("gene_set must be a string")
    if gene_set not in data:
        raise KeyError(f"gene_set '{gene_set}' not found in data")

    plt.style.use('seaborn')
    sns.set_palette("deep")

    fig, ax = plt.subplots(figsize=(10, 6))

    hist = ax.hist(data[gene_set], bins=20, edgecolor='black', alpha=0.7)

    ax.set_xlabel('Weighting Coefficient', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title(f'Gene Set: {gene_set}', fontsize=14, fontweight='bold')

    ax.grid(True, linestyle='--', alpha=0.7)
    ax.tick_params(axis='both', which='major', labelsize=10)

    plt.tight_layout()
    plt.show()

def better_hist_GRNs(data: dict, set_selected: str, cell_type_selected: str, scored_gene_selected: str, 
                     score: str, bins: int = 20, genes_to_mark: list = None) -> None:
    """
    Create an enhanced histogram for Gene Regulatory Networks.

    Args:
        data (dict): Dictionary containing gene expression data
        set_selected (str): Name of the gene set to analyze
        cell_type_selected (str): Cell type to analyze
        scored_gene_selected (str): Name of the scored gene to analyze
        score (str): Score type to plot
        bins (int, optional): Number of histogram bins. Defaults to 20
        genes_to_mark (list, optional): List of genes to highlight. Defaults to None

    Raises:
        TypeError: If arguments are not of expected types
        ValueError: If bins is less than 1
        KeyError: If keys are not found in data dictionary
    """
    # Input validation
    if not isinstance(data, dict):
        raise TypeError("data must be a dictionary")
    if not isinstance(bins, int) or bins < 1:
        raise ValueError("bins must be a positive integer")
    if genes_to_mark is not None and not isinstance(genes_to_mark, list):
        raise TypeError("genes_to_mark must be a list or None")

    plt.style.use('seaborn')
    sns.set_palette("deep")

    fig, ax = plt.subplots(figsize=(12, 6))

    try:
        hist_data = data[set_selected][cell_type_selected][scored_gene_selected][score]
    except KeyError as e:
        raise KeyError(f"Could not find data: {str(e)}")
    
    hist, bin_edges = np.histogram(hist_data, bins=bins)
    ax.hist(hist_data, bins=bins, edgecolor='black', alpha=0.7)

    ax.set_xlabel('Weighting Coefficient', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title(f'Gene Set: {scored_gene_selected}', fontsize=14, fontweight='bold')

    ax.grid(True, linestyle='--', alpha=0.7)
    ax.tick_params(axis='both', which='major', labelsize=10)

    if genes_to_mark:
        colors = plt.cm.rainbow(np.linspace(0, 1, len(genes_to_mark)))
        for gene, color in zip(genes_to_mark, colors):
            if gene in data[set_selected][cell_type_selected][scored_gene_selected]["targets"]:
                gene_index = data[set_selected][cell_type_selected][scored_gene_selected]["targets"].index(gene)
                gene_value = data[set_selected][cell_type_selected][scored_gene_selected][score][gene_index]
                ax.axvline(x=gene_value, color=color, linestyle='--', label=gene)
            else:
                print(f"Warning: Gene {gene} not found in the dataset.")

        ax.legend(title="Marked Genes", title_fontsize=10, fontsize=8, loc='upper right', bbox_to_anchor=(1.25, 1))

    plt.tight_layout()
    plt.show()

def print_duplicate_genes(gene_sets: dict, max_examples: int = 5) -> None:
    """
    Print information about duplicate genes in gene sets.

    Args:
        gene_sets (dict): Dictionary mapping set names to collections of genes
        max_examples (int, optional): Maximum number of duplicate examples to show. Defaults to 5

    Raises:
        TypeError: If gene_sets is not a dictionary or max_examples is not an integer
        ValueError: If max_examples is less than 1
    """
    if not isinstance(gene_sets, dict):
        raise TypeError("gene_sets must be a dictionary")
    if not isinstance(max_examples, int) or max_examples < 1:
        raise ValueError("max_examples must be a positive integer")

    for set_name, gene_set in gene_sets.items():
        original_length = len(gene_set)
        # Convert to a list to preserve order
        genes = list(gene_set)
        
        # Find duplicates
        seen = set()
        duplicates = []
        for gene in genes:
            if gene in seen:
                duplicates.append(gene)
            else:
                seen.add(gene)
        
        # Print results
        if duplicates:
            print(f"\nDuplicates found in {set_name}:")
            print(f"Total duplicates: {len(duplicates)} out of {original_length}")
            print(f"Examples (up to {max_examples}):")
            for gene in duplicates[:max_examples]:
                indices = [i for i, x in enumerate(genes) if x == gene]
                print(f"  Gene '{gene}' appears at indices: {indices}")
        else:
            print(f"\nNo duplicates found in {set_name}")

def plot_gene_set_intersections(gene_sets: dict, figsize: tuple = (12, 8), title: str = "", 
                              save_path: str = None, dpi: int = 300, maxlength: int = 100) -> None:
    """
    Create an UpSet plot showing intersections between multiple gene sets.

    Args:
        gene_sets (dict): Dictionary mapping set names to collections of genes
        figsize (tuple, optional): Figure size as (width, height). Defaults to (12, 8)
        title (str, optional): Title for the plot. Defaults to ""
        save_path (str, optional): Path to save the figure. If None, displays plot. Defaults to None
        dpi (int, optional): DPI for saved figure. Defaults to 300
        maxlength (int, optional): Maximum number of genes to include from each set. Defaults to 100

    Raises:
        TypeError: If arguments are of incorrect type
        ValueError: If numeric arguments are invalid or gene_sets is empty
        RuntimeError: If plotting fails

    Returns:
        None: Displays or saves the plot
    """
    # Type checking
    if not isinstance(gene_sets, dict):
        raise TypeError("gene_sets must be a dictionary")
    if not isinstance(figsize, tuple) or len(figsize) != 2:
        raise TypeError("figsize must be a tuple of length 2")
    if not isinstance(title, str):
        raise TypeError("title must be a string")
    if save_path is not None and not isinstance(save_path, str):
        raise TypeError("save_path must be None or a string")
    if not isinstance(dpi, int):
        raise TypeError("dpi must be an integer")
    if not isinstance(maxlength, int):
        raise TypeError("maxlength must be an integer")

    # Value checking
    if not gene_sets:
        raise ValueError("gene_sets cannot be empty")
    if any(x <= 0 for x in figsize):
        raise ValueError("figsize values must be positive")
    if dpi <= 0:
        raise ValueError("dpi must be positive")
    if maxlength <= 0:
        raise ValueError("maxlength must be positive")

    # Suppress FutureWarnings
    warnings.simplefilter(action='ignore', category=FutureWarning)
    
    try:
        # Trim gene sets to maxlength
        trimmed_gene_sets = {key: set(list(value)[:maxlength]) for key, value in gene_sets.items()}
        
        # Convert trimmed gene sets to a format suitable for UpSet plot
        data = from_contents(trimmed_gene_sets)
        
        # Create the UpSet plot
        upset = UpSet(data, show_counts=True, sort_by='cardinality', element_size=40)
        
        # Set up the plot
        fig = plt.figure(figsize=figsize)
        
        # Plot with warning suppression
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            upset.plot(fig)
        
        # Remove the default title
        plt.title('')
        
        # Add the title vertically on the left side
        if title:
            fig.text(0.1, 0.5, title, 
                    rotation=90, 
                    verticalalignment='center', 
                    horizontalalignment='center',
                    fontsize=12,
                    fontweight='bold')
        
        # Disable tight_layout for this figure
        fig.set_tight_layout(False)
        
        # Handle saving/displaying with warning suppression
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            
            if save_path:
                plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
                print(f"Figure saved to {save_path}")
            else:
                plt.show()
        
        plt.close(fig)

    except Exception as e:
        plt.close()  # Clean up in case of error
        raise RuntimeError(f"Error creating UpSet plot: {str(e)}")


### GRNs Scores funcitons ####

def boxplot_EDCs_GRN_scores_parameters_local(
        adata: ad.AnnData,
        conditions: str,
        gene_set: str,
        cell_type: str,
        goi: str,
        control_condition: str,
        control: str = "True",
        normalize_weights: str = "False",
        scale_by_variance: str = "False",
        scaling_only_based_on_control: str = "True",
        prefix: str = ""
) -> None:
    """
    Create a boxplot of GRN scores for different EDC concentrations with overlaid scatter points.

    Args:
        adata (AnnData): Annotated data matrix containing gene expression and metadata
        conditions (str): Name of condition to filter concentrations by
        gene_set (str): Name of the gene set to use for scoring
        cell_type (str): Cell type to analyze
        goi (str): Gene of interest
        control_condition (str): Name of control condition
        control (str, optional): Whether to include control ("True"/"False"). Defaults to "True"
        normalize_weights (str, optional): Whether to normalize weights ("True"/"False"). Defaults to "False"
        scale_by_variance (str, optional): Whether to scale by variance ("True"/"False"). Defaults to "False"
        scaling_only_based_on_control (str, optional): Whether to scale based only on control. Defaults to "True"
        prefix (str, optional): Prefix for column names. Defaults to ""

    Raises:
        TypeError: If arguments are not of expected types
        ValueError: If required columns are missing or arguments have invalid values
        RuntimeError: If no data is available for plotting
    """
    # Type checking
    if not isinstance(adata, ad.AnnData):
        raise TypeError("adata must be an AnnData object")
    if not all(isinstance(x, str) for x in [conditions, gene_set, cell_type, goi, control_condition, 
                                          normalize_weights, scale_by_variance, control, 
                                          scaling_only_based_on_control, prefix]):
        raise TypeError("String arguments must be strings")

    # Validate string boolean values
    for param_name, param_value in [("normalize_weights", normalize_weights), 
                                  ("scale_by_variance", scale_by_variance),
                                  ("control", control),
                                  ("scaling_only_based_on_control", scaling_only_based_on_control)]:
        if param_value not in ["True", "False"]:
            raise ValueError(f"{param_name} must be 'True' or 'False'")

    # Validate required columns exist
    if 'condition_concentraion' not in adata.obs.columns:
        raise ValueError("adata.obs must contain 'condition_concentraion' column")

    plt.close('all')

    # Get all concentrations for the selected conditions
    condition_concentrations = adata.obs['condition_concentraion'].unique()
    selected_concentrations = [conc for conc in condition_concentrations if conditions in conc]
    
    if not selected_concentrations:
        raise ValueError(f"No concentrations found containing condition: {conditions}")
        
    all_concentrations = ['DMSO_0.1'] + selected_concentrations

    mask = adata.obs['condition_concentraion'].isin(all_concentrations)
    adata_filtered = adata[mask]

    if adata_filtered.shape[0] == 0:
        raise RuntimeError(f'No data available for the selected conditions: {all_concentrations}')

    fig, ax = plt.subplots(figsize=(12, 8))

    gene_scores = {}

    selection = (
        f'{prefix}'
        f'gene_score_{gene_set}_{cell_type}_{goi}_{control}_'
        f'normalized_{normalize_weights}_'
        f'scaled_{scale_by_variance}_'
        f'cc_{control_condition}_'
        f'sc_{scaling_only_based_on_control}'
    )

    # Validate selection column exists
    if selection not in adata_filtered.obs.columns:
        raise ValueError(f"Column '{selection}' not found in adata.obs")

    for conc in all_concentrations:
        gene_scores[conc] = adata_filtered[adata_filtered.obs['condition_concentraion'] == conc].obs[selection].values.tolist()

    data = [gene_scores[conc] for conc in all_concentrations]

    positions = np.arange(len(all_concentrations))
    width = 0.5

    bp = ax.boxplot(data, positions=positions, widths=width, patch_artist=True,
                    boxprops=dict(facecolor='C0', color='None', alpha=0.5),
                    whiskerprops=dict(color='C0'),
                    capprops=dict(color='C0'),
                    showfliers=False,
                    medianprops=dict(color='white'))

    for i, scores in enumerate(data):
        if scores:  # Check if scores is not empty
            ax.scatter([positions[i]] * len(scores), scores, color='C0', alpha=0.7, s=20)

    ax.set_xticks(positions)
    ax.set_xticklabels(all_concentrations, rotation=45, ha='right', fontsize=18)
    ax.set_title(f'GRNs Gene Scores - {goi}\nCell type - {cell_type}', fontsize=18)
    ax.set_ylabel("Gene Score", fontsize=18)

    plt.tight_layout()
    plt.show()

def boxplot_Reference_GRN_scores_parameters_all_conditions(
        adata: "AnnData", # type: ignore
        conditions: list[list[str]], 
        gene_set: str,
        cell_type: str,
        goi: str,
        control_condition: str,
        control: bool = True,
        normalize_weights: str = "False",
        scale_by_variance: str = "False", 
        scaling_only_based_on_control: bool = True,
        prefix: str = ""
) -> "matplotlib.figure.Figure":
    """
    Create boxplots comparing gene scores across multiple conditions.

    Args:
        adata: AnnData object containing the data
        conditions: List of lists, where each inner list contains sample types to compare
        gene_set: Name of the gene set
        cell_type: Cell type to analyze
        goi: Gene of interest
        control_condition: Control condition name
        control: Whether to use control in score calculation
        normalize_weights: Whether to normalize weights ("True"/"False")
        scale_by_variance: Whether to scale by variance ("True"/"False")
        scaling_only_based_on_control: Whether to scale using only control data
        prefix: Optional prefix for score column names

    Returns:
        matplotlib Figure object containing the plots

    Raises:
        ValueError: If input parameters are invalid
        KeyError: If required columns are missing from adata
    """
    if not isinstance(conditions, list) or not all(isinstance(c, list) for c in conditions):
        raise ValueError("conditions must be a list of lists of strings")
    if not isinstance(adata.obs, pd.DataFrame):
        raise ValueError("adata.obs must be a pandas DataFrame")
    if 'sample_type' not in adata.obs.columns:
        raise KeyError("adata.obs must contain 'sample_type' column")

    plt.close('all')

    fig, axes = plt.subplots(1, len(conditions), figsize=(3 * len(conditions), 5), sharey=False)
    
    for idx, condition in enumerate(conditions):
        mask = adata.obs['sample_type'].isin(condition)
        adata_filtered = adata[mask]

        if len(adata_filtered) > 0:
            ax = axes[idx] if len(conditions) > 1 else axes
            
            gene_scores = []
            for cond in condition:
                selection = (
                    f'{prefix}'
                    f'gene_score_{gene_set}_{cell_type}_{goi}_{control}_'
                    f'normalized_{normalize_weights}_'
                    f'scaled_{scale_by_variance}_'
                    f'cc_{control_condition}_'
                    f'sc_{scaling_only_based_on_control}'
                )
                
                if selection not in adata_filtered.obs.columns:
                    raise KeyError(f"Column '{selection}' not found in adata.obs")
                    
                scores = adata_filtered[adata_filtered.obs['sample_type'] == cond].obs[selection].values
                gene_scores.append(scores)
            
            sns.boxplot(data=gene_scores, notch=False, boxprops=dict(alpha=0.5), ax=ax)
            sns.stripplot(data=gene_scores, jitter=True, color=".3", linewidth=1, ax=ax)
            
            ax.set_xticks(range(len(condition)))
            ax.set_xticklabels(condition, fontsize=14, rotation=45, ha='right')
            ax.set_title(f'{condition[1].split("_")[0]}', fontsize=14)
            ax.set_ylabel("Gene Score", fontsize=12)
            
            # Automatically adjust y-axis limits for each subplot
            if gene_scores:  # Check if we have any data
                y_min = min([min(scores) for scores in gene_scores if len(scores) > 0])
                y_max = max([max(scores) for scores in gene_scores if len(scores) > 0])
                y_range = y_max - y_min
                ax.set_ylim(y_min - 0.1 * y_range, y_max + 0.1 * y_range)
        else:
            print(f'No data to plot for the condition: {condition}')

    fig.suptitle(f'Gene Scores - {goi} ({cell_type})', fontsize=16)
    plt.tight_layout()
    return fig

def boxplot_Reference_GRN_scores_parameters_all_gois(
        adata: "AnnData", # type: ignore
        gene_set: str,
        condition: Union[str, list[str]],
        cell_type: str,
        gois: list[str],
        control_condition: str,
        control: bool = True,
        normalize_weights: str = "False",
        scale_by_variance: str = "False",
        scaling_only_based_on_control: bool = True,
        prefix: str = ""
) -> Optional["matplotlib.figure.Figure"]:
    """
    Create boxplots comparing gene scores across multiple genes of interest.

    Args:
        adata: AnnData object containing the data
        gene_set: Name of the gene set
        condition: Single condition or list of conditions to analyze
        cell_type: Cell type to analyze
        gois: List of genes of interest
        control_condition: Control condition name
        control: Whether to use control in score calculation
        normalize_weights: Whether to normalize weights ("True"/"False")
        scale_by_variance: Whether to scale by variance ("True"/"False")
        scaling_only_based_on_control: Whether to scale using only control data
        prefix: Optional prefix for score column names

    Returns:
        matplotlib Figure object containing the plots, or None if no data available

    Raises:
        ValueError: If input parameters are invalid
        KeyError: If required columns are missing from adata
    """
    if not isinstance(gois, list) or not all(isinstance(g, str) for g in gois):
        raise ValueError("gois must be a list of strings")
    if not isinstance(adata.obs, pd.DataFrame):
        raise ValueError("adata.obs must be a pandas DataFrame")
    if 'sample_type' not in adata.obs.columns:
        raise KeyError("adata.obs must contain 'sample_type' column")

    plt.close('all')

    nrows = 2
    ncols = (len(gois) + 1) // 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 6 * nrows), sharey=False)
    axes = axes.flatten() if len(gois) > 1 else [axes]

    conditions = [condition] if isinstance(condition, str) else condition
    
    mask = adata.obs['sample_type'].isin(conditions)
    adata_filtered = adata[mask]

    if len(adata_filtered) > 0:
        for idx, goi in enumerate(gois):
            ax = axes[idx]
            
            selection = (
                f'{prefix}'
                f'gene_score_{gene_set}_{cell_type}_{goi}_{control}_'
                f'normalized_{normalize_weights}_'
                f'scaled_{scale_by_variance}_'
                f'cc_{control_condition}_'
                f'sc_{scaling_only_based_on_control}'
            )
            
            if selection not in adata_filtered.obs.columns:
                raise KeyError(f"Column '{selection}' not found in adata.obs")
            
            gene_scores = [adata_filtered[adata_filtered.obs['sample_type'] == cond].obs[selection].values 
                         for cond in conditions]

            if any(len(scores) == 0 for scores in gene_scores):
                print(f"Warning: No data found for some conditions for gene {goi}")
                continue

            sns.boxplot(data=gene_scores, notch=False, boxprops=dict(alpha=0.5), ax=ax)
            sns.stripplot(data=gene_scores, jitter=True, color=".3", linewidth=1, ax=ax)
            
            ax.set_xticks(range(len(conditions)))
            ax.set_xticklabels(conditions, fontsize=12, rotation=45, ha='right')
            ax.set_title(f'Gene: {goi}', fontsize=14)
            ax.set_ylabel("Gene Score", fontsize=12)
            
            # Automatically adjust y-axis limits
            all_scores = np.concatenate(gene_scores)
            if len(all_scores) > 0:
                y_min = np.min(all_scores)
                y_max = np.max(all_scores)
                y_range = y_max - y_min
                ax.set_ylim(y_min - 0.1 * y_range, y_max + 0.1 * y_range)

        # Hide any unused axes
        for j in range(len(gois), len(axes)):
            fig.delaxes(axes[j])

        fig.suptitle(f'Gene Scores - {gene_set} ({cell_type})', fontsize=16)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        return fig
    else:
        print(f'No data to plot for the conditions: {conditions}')
        return None










################################## UNUSED ###################################################3









def load_Sim_gene_sets(gene_set_list, root_dir):
    gene_sets = {}

    # Load data for each gene set and sim file type
    for gene_set in gene_set_list:
        path = os.path.join(root_dir, f"output_hg19_{gene_set}", "celloracle")
        print(path)
        
        # Load scores_sim_all.csv
        sim_1_path = os.path.join(path, 'scores_sim_all.csv')
        if os.path.exists(sim_1_path):
            gene_sets[gene_set] = pd.read_csv(sim_1_path)
        else:
            print(f"File not found: {sim_1_path}")
        
        # # Load v2_scores_sim_all.csv
        # sim_2_path = os.path.join(path, 'scores_sim_all_to_test_old.csv')
        # if os.path.exists(sim_2_path):
        #     gene_sets[f"{gene_set}_sim_2"] = pd.read_csv(sim_2_path)
        # else:
        #     print(f"File not found: {sim_2_path}")

    # Format Data
    for key, df in gene_sets.items():
        df['abs_fold_change'] = df['fold_change'].abs()
        df = df.groupby('goi').apply(lambda x: x.nlargest(50, 'abs_fold_change')).reset_index(drop=True)
        df['fold_change_enh'] = np.power(df['fold_change'], 10)
        gene_sets[key] = df

    # Sets Formatting
    gene_sets_dict = {}
    gene_sets_dict_cell_type_first = {}

    for key, df in gene_sets.items():
        gene_sets_dict[key] = {}
        gene_sets_dict_cell_type_first[key] = {}

        for _, row in df.iterrows():
            local_cell_type = row['local_cell_type']
            gene = row['gene']
            goi = row['goi']
            fold_change = float(row['fold_change'])
            cell_type = row['cell_type']

            # Format 1: Gene of interest first
            if goi not in gene_sets_dict[key]:
                gene_sets_dict[key][goi] = {}

            if cell_type not in gene_sets_dict[key][goi]:
                gene_sets_dict[key][goi][cell_type] = {'targets': [], 'weights': [], 'local_cell_types': []}

            gene_sets_dict[key][goi][cell_type]['targets'].append(gene)
            gene_sets_dict[key][goi][cell_type]['weights'].append(fold_change)
            gene_sets_dict[key][goi][cell_type]['local_cell_types'].append(local_cell_type)

            # Format 2: Cell type first
            if cell_type not in gene_sets_dict_cell_type_first[key]:
                gene_sets_dict_cell_type_first[key][cell_type] = {}

            if goi not in gene_sets_dict_cell_type_first[key][cell_type]:
                gene_sets_dict_cell_type_first[key][cell_type][goi] = {'targets': [], 'weights': [], 'local_cell_types': []}

            gene_sets_dict_cell_type_first[key][cell_type][goi]['targets'].append(gene)
            gene_sets_dict_cell_type_first[key][cell_type][goi]['weights'].append(fold_change)
            gene_sets_dict_cell_type_first[key][cell_type][goi]['local_cell_types'].append(local_cell_type)

    print("Processed datasets:", gene_sets_dict_cell_type_first.keys())

    return gene_sets_dict, gene_sets_dict_cell_type_first

def load_expression_data(data_dir = "/group/testa/michal.kubacki/gene_score"):
    # Read the expression matrix
    raw_matrix = pd.read_csv(os.path.join(data_dir, "Organoids/CTL04/raw_matrix.csv"), index_col=0)
    logcpm_matrix = pd.read_csv(os.path.join(data_dir, "Organoids/CTL04/logcpm_matrix.csv"), index_col=0)
    cpm_matrix = pd.read_csv(os.path.join(data_dir, "Organoids/CTL04/cpm_matrix.csv"), index_col=0)

    # Read the row data (genes/features)
    row_data = pd.read_csv(os.path.join(data_dir, "Organoids/CTL04/row_data.csv"), index_col=0)

    # Read the column data (sample metadata)
    col_data = pd.read_csv(os.path.join(data_dir, "Organoids/CTL04/col_data.csv"), index_col=0)

    expr_matrix = logcpm_matrix

    adata = ad.AnnData(X=expr_matrix.transpose(), obs=col_data, var=row_data)

    print([var_name for var_name in adata.var_names if var_name.startswith("AHR")])
    print([var_name for var_name in adata.var_names if var_name.startswith("AR")])
    print([var_name for var_name in adata.var_names if "NR3" in var_name])

    expr_matrix_df = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
    expr_matrix_df

    adata.obs['sample_type'] = 'Other'
    adata.obs.loc[adata.obs.index.str.endswith('_DMSO'), 'sample_type'] = 'DMSO'
    adata.obs.loc[adata.obs.index.str.endswith('_CTL'), 'sample_type'] = 'CTL'

    adata.obs.loc[adata.obs.index.str.endswith('_Ret_Ag'), 'sample_type'] = 'Ret_Ag'
    adata.obs.loc[adata.obs.index.str.endswith('_Ret_Inh'), 'sample_type'] = 'Ret_Inh'

    adata.obs.loc[adata.obs.index.str.endswith('_AhHyd_Ag'), 'sample_type'] = 'AhHyd_Ag'
    adata.obs.loc[adata.obs.index.str.endswith('_AhHyd_Inh'), 'sample_type'] = 'AhHyd_Inh'

    adata.obs.loc[adata.obs.index.str.endswith('_Andr_Ag'), 'sample_type'] = 'Andr_Ag'
    adata.obs.loc[adata.obs.index.str.endswith('_Andr_Inh'), 'sample_type'] = 'Andr_Inh'

    adata.obs.loc[adata.obs.index.str.endswith('_LivX_Ag'), 'sample_type'] = 'LivX_Ag'
    adata.obs.loc[adata.obs.index.str.endswith('_LivX_Inh'), 'sample_type'] = 'LivX_Inh'

    adata.obs.loc[adata.obs.index.str.endswith('_GC_Ag'), 'sample_type'] = 'GC_Ag'
    adata.obs.loc[adata.obs.index.str.endswith('_GC_Inh'), 'sample_type'] = 'GC_Inh'

    adata.obs.loc[adata.obs.index.str.endswith('_Estr_Ag'), 'sample_type'] = 'Estr_Ag'
    adata.obs.loc[adata.obs.index.str.endswith('_Estr_Inh'), 'sample_type'] = 'Estr_Inh'

    adata.obs.loc[adata.obs.index.str.endswith('_Thyr_Ag'), 'sample_type'] = 'Thyr_Ag'
    adata.obs.loc[adata.obs.index.str.endswith('_Thyr_Inh'), 'sample_type'] = 'Thyr_Inh'
    return(adata)

def compare_target_intersections(gene_sets_dict, set_names, cell_types, genes):
    result = {set_name: {} for set_name in set_names}
    
    for set_name in set_names:
        if set_name not in gene_sets_dict:
            print(f"Error: Set '{set_name}' not found in the dictionary.")
            continue
        
        for cell_type in cell_types:
            if cell_type not in gene_sets_dict[set_name]:
                print(f"Warning: Cell type '{cell_type}' not found in the '{set_name}' set.")
                continue
            
            cell_type_data = gene_sets_dict[set_name][cell_type]
            intersection = set()
            first_gene = True
            
            for gene in genes:
                if gene not in cell_type_data:
                    print(f"Warning: Gene '{gene}' not found for cell type '{cell_type}' in set '{set_name}'.")
                    continue
                
                gene_targets = set(cell_type_data[gene]['targets'])
                
                if first_gene:
                    intersection = gene_targets
                    first_gene = False
                else:
                    intersection &= gene_targets
            
            result[set_name][cell_type] = list(intersection)
    
    return result



########### Legacy Functions ###############
def process_network(df, source, celltype):
    filtered_network = df[(df['source'] == source) & (df['celltype'] == celltype)]
    # Sort by 'p' value
    filtered_network = filtered_network.sort_values(by='p')
    genes = filtered_network['target'].tolist()
    effects = filtered_network['score'].tolist()
    return genes, effects

def load_pathways(data_dir = "/group/testa/michal.kubacki/gene_score"):
    # gene_sets_files = ["REACTOME_CHOLESTEROL_BIOSYNTHESIS.v2023.2.Hs.grp", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE.v2023.2.Hs.grp",
    #                   "BIOCARTA_PROTEASOME_PATHWAY.v2023.2.Hs.grp", "BIOCARTA_TEL_PATHWAY.v2023.2.Hs.grp",
    #                   "BIOCARTA_P53HYPOXIA_PATHWAY.v2023.2.Hs.grp", "BIOCARTA_CARM_ER_PATHWAY.v2023.2.Hs.grp",
    #                   "REACTOME_PROCESSING_OF_CAPPED_INTRON_CONTAINING_PRE_MRNA.v2023.2.Hs.grp", "REACTOME_DEUBIQUITINATION.v2023.2.Hs.grp",
    #                   "REACTOME_SARS_COV_1_INFECTION.v2023.2.Hs.grp", "REACTOME_SIGNALING_BY_ROBO_RECEPTORS.v2023.2.Hs.grp"]
    gene_sets_files = ["REACTOME_CHOLESTEROL_BIOSYNTHESIS.v2023.2.Hs.grp", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE.v2023.2.Hs.grp"]

    base_path = os.path.join(data_dir, "Pathways")
    gene_sets = {}
    weight_sets = {}
    for gene_set in gene_sets_files:
        genes = []
        with open(os.path.join(base_path, gene_set), "r") as file:
            # Skip the first two lines (header and URL)
            next(file)
            next(file)

            for line in file:
                genes.append(line.strip())
        gene_sets[gene_set]= genes
        weight_sets[gene_set] = None
    return gene_sets

def load_endocrineReceptorsGRNs_L_old_format(data_dir = "/group/testa/michal.kubacki/gene_score"):
    # Load the data
    endocrine_rec_networks = pd.read_csv(os.path.join(data_dir, "endocrineReceptorsGRNs_L.tsv"), sep="\t")
    print(endocrine_rec_networks.head())
    print(endocrine_rec_networks.source.unique())
    celltype_counts = endocrine_rec_networks['celltype'].value_counts().reset_index()
    celltype_counts.columns = ['celltype', 'count']
    print(celltype_counts)
    gene_sets = {}
    weight_sets = {}
    celltype = "Exc_Mig"

    # ESR2: estrogen signalling
    genes_ESR2, effects_ESR2 = process_network(endocrine_rec_networks, "ESR2", celltype)
    if len(genes_ESR2) > 0:
        gene_sets["ESR2"] = genes_ESR2
        weight_sets["ESR2"] = effects_ESR2

    # THRB: thyroid hormone signalling
    genes_THRB, effects_THRB = process_network(endocrine_rec_networks, "THRB", celltype)
    if len(genes_THRB) > 0:
        gene_sets["THRB"] = genes_THRB
        weight_sets["THRB"] = effects_THRB

    # RARA: retonoic acid signalling
    genes_RARA, effects_RARA = process_network(endocrine_rec_networks, "RARA", celltype)
    if len(genes_RARA) > 0:
        gene_sets["RARA"] = genes_RARA
        weight_sets["RARA"] = effects_RARA

    # NR2F1: Potentially related to Ret (Retinoid Receptors) or LivX (Liver X Receptor)
    genes_NR2F1, effects_NR2F1 = process_network(endocrine_rec_networks, "NR2F1", celltype)
    if len(genes_NR2F1) > 0:
        gene_sets["NR2F1"] = genes_NR2F1
        weight_sets["NR2F1"] = effects_NR2F1

    # RORA: Potentially related to Ret
    genes_RORA, effects_RORA = process_network(endocrine_rec_networks, "RORA", celltype)
    if len(genes_RORA) > 0:
        gene_sets["RORA"] = genes_RORA
        weight_sets["RORA"] = effects_RORA
    
    # NR4A1: Potentially related to Ret (Retinoid Receptors) or LivX (Liver X Receptor)
    genes_NR4A1, effects_NR4A1 = process_network(endocrine_rec_networks, "NR4A1", celltype)
    if len(genes_NR4A1) > 0:
        gene_sets["NR4A1"] = genes_NR4A1
        weight_sets["NR4A1"] = effects_NR4A1

    # NR2F2: Potentially related to Ret (Retinoid Receptors) or LivX (Liver X Receptor)
    genes_NR2F2, effects_NR2F2 = process_network(endocrine_rec_networks, "NR2F2", celltype)
    if len(genes_NR2F2) > 0:
        gene_sets["NR2F2"] = genes_NR2F2
        weight_sets["NR2F2"] = effects_NR2F2 

    # NR4A2: Potentially related to Ret (Retinoid Receptors) or LivX (Liver X Receptor)
    genes_NR4A2, effects_NR4A2 = process_network(endocrine_rec_networks, "NR4A2", celltype)
    if len(genes_NR4A2) > 0:
        gene_sets["NR4A2"] = genes_NR4A2
        weight_sets["NR4A2"] = effects_NR4A2 
    return gene_sets

def load_topDEGs(file='topDEGs_all_conditions.csv'):
    gene_sets = {}
    weight_sets = {}
    
    with open(file, 'r') as file:
        reader = csv.DictReader(file)
        
        for row in reader:
            for column, value in row.items():
                condition, comparison, data_type = column.split('.')
                
                if condition not in gene_sets:
                    gene_sets[condition] = []
                    weight_sets[condition] = []
                
                if data_type == 'gene':
                    gene_sets[condition].append(value)
                elif data_type == 'score':
                    score = float(value) if value else None
                    weight_sets[condition].append(score)
    return gene_sets, weight_sets

def convert_weight_sets(weight_sets):
    if isinstance(weight_sets, dict) and all(isinstance(weights, dict) for weights in weight_sets.values()):
        weight_sets_list = {}
        for condition, weights in weight_sets.items():
            weight_sets_list[condition] = list(weights.values())
        return weight_sets_list
    else:
        print("weight_sets has already correct format")

def load_topDEGs_unique(file='topDEGs_all_conditions.csv'):
    gene_sets_unique = {}
    weight_sets_unique = {}
    gene_conditions = {}
    gene_scores = {}
    
    with open(file, 'r') as file:
        reader = csv.DictReader(file)
        
        for row in reader:
            for column, value in row.items():
                condition, comparison, data_type = column.split('.')
                
                if condition not in gene_sets_unique:
                    gene_sets_unique[condition] = []
                    weight_sets_unique[condition] = {}
                    gene_scores[condition] = {}
                
                if data_type == 'gene':
                    gene = value
                    if gene not in gene_conditions:
                        gene_conditions[gene] = []
                    gene_conditions[gene].append(condition)
                    gene_sets_unique[condition].append(gene)
                elif data_type == 'score':
                    score = float(value) if value else None
                    gene_scores[condition][gene] = score
        
        for condition in gene_sets_unique:
            for gene in gene_sets_unique[condition]:
                if gene in gene_scores[condition]:
                    weight_sets_unique[condition][gene] = gene_scores[condition][gene]
        
        for gene, conditions in gene_conditions.items():
            max_score_condition = None
            max_score = float('-inf')
            
            for condition in conditions:
                if gene in weight_sets_unique[condition]:
                    score = weight_sets_unique[condition][gene]
                    if score is not None and score > max_score:
                        max_score = score
                        max_score_condition = condition
            
            if max_score_condition is not None:
                for condition in conditions:
                    if condition == max_score_condition:
                        weight_sets_unique[condition][gene] = max_score
                    else:
                        weight_sets_unique[condition][gene] = 0.0
    
    weight_sets_unique = convert_weight_sets(weight_sets_unique)
    return gene_sets_unique, weight_sets_unique
    