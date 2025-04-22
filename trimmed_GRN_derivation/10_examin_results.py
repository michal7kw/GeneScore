# %%
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import scanpy as sc
from datetime import datetime

# Set working directory
work_dir = '/home/michal.kubacki/Githubs/GeneScore/trimmed_GRN_derivation'
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
    sys.path.insert(0, '/home/michal.kubacki/Githubs/GeneScore/project_functions')
    from grn_helpers import *

# %%
# Configuration parameters
neurons_set = "L2-3_CUX2"  # Options: "all_ex", "all_ex_all_ages", "L2-3_CUX2"

# Get base path from environment variables with error handling
root_dir = os.getenv('BASE_PATH')
if not root_dir:
    raise ValueError("BASE_PATH environment variable not found in .env file")

# %%
# Set up directories
output_dir, input_dir, root_dir, tmp_dir, in_dir_from_scenic = set_custom_folders(root_dir, neurons_set)

# Plot settings
plt.rcParams['figure.figsize'] = [10, 6]
plt.rcParams["savefig.dpi"] = 300

# %%
# Load the results from the previous script
print("Loading simulation scores...")
sim_scores_path = os.path.join(output_dir, 'scores_sim_all_new.csv')
sim_scores = pd.read_csv(sim_scores_path)
print(f"Loaded simulation scores with shape: {sim_scores.shape}")

print("\nLoading GRN scores...")
grn_scores_path = os.path.join(output_dir, 'scores_grn_all_from_comb_run_new.csv')
grn_scores = pd.read_csv(grn_scores_path)
print(f"Loaded GRN scores with shape: {grn_scores.shape}")

# %%
# Display the first few rows of each dataset
print("\nSimulation scores sample:")
# print(sim_scores.head())
sim_scores.head()

# %%
print("\nGRN scores sample:")
# print(grn_scores.head())
grn_scores.head()

# %%
# Get basic statistics for the datasets
print("\nSimulation scores statistics:")
print(sim_scores.describe())

print("\nGRN scores statistics:")
print(grn_scores.describe())

# %%
plt.hist(grn_scores['p'], bins=50)
print(min(grn_scores['p']))
print(max(grn_scores['p']))

# %%
# Get unique perturbed genes and cell types
perturbed_genes = sim_scores['perturbed_gene'].unique()
print(f"Unique perturbed genes: {len(perturbed_genes)}")
print(perturbed_genes)

cell_types = sim_scores['cell_type'].unique()
print(f"\nUnique cell types: {len(cell_types)}")
print(cell_types)

# %%
# Function to analyze a specific gene perturbation
def analyze_gene_perturbation(gene, cell_type=None):
    """Analyze the effects of perturbing a specific gene.
    
    Args:
        gene (str): The gene that was perturbed
        cell_type (str, optional): Filter by cell type. If None, all cell types are included.
    """
    # Filter the data
    sim_filter = sim_scores['perturbed_gene'] == gene
    grn_filter = grn_scores['perturbed_gene'] == gene
    
    if cell_type is not None:
        sim_filter &= (sim_scores['cell_type'] == cell_type)
        grn_filter &= (grn_scores['cell_type'] == cell_type)
    
    sim_data = sim_scores[sim_filter]
    grn_data = grn_scores[grn_filter]
    
    print(f"Analysis for perturbation of gene: {gene}")
    if cell_type:
        print(f"Cell type: {cell_type}")
    
    # Display summary statistics
    print(f"\nSimulation data shape: {sim_data.shape}")
    print(f"GRN data shape: {grn_data.shape}")
    
    # Top affected genes by simulation scores
    if not sim_data.empty:
        print("\nTop 10 affected genes by simulation scores:")
        top_sim = sim_data.sort_values('score', ascending=False).head(10)
        print(top_sim[['target_gene', 'score']])
        
        # Plot distribution of simulation scores
        plt.figure(figsize=(10, 6))
        sns.histplot(sim_data['score'], kde=True)
        plt.title(f'Distribution of Simulation Scores for {gene} Perturbation')
        plt.xlabel('Score')
        plt.ylabel('Frequency')
        plt.savefig(os.path.join(output_dir, f"sim_score_dist_{gene}{'_'+cell_type if cell_type else ''}.png"))
        plt.show()
    
    # Top affected genes by GRN scores
    if not grn_data.empty:
        print("\nTop 10 affected genes by GRN scores:")
        top_grn = grn_data.sort_values('score', ascending=False).head(10)
        print(top_grn[['target_gene', 'score']])
        
        # Plot distribution of GRN scores
        plt.figure(figsize=(10, 6))
        sns.histplot(grn_data['score'], kde=True)
        plt.title(f'Distribution of GRN Scores for {gene} Perturbation')
        plt.xlabel('Score')
        plt.ylabel('Frequency')
        plt.savefig(os.path.join(output_dir, f"grn_score_dist_{gene}{'_'+cell_type if cell_type else ''}.png"))
        plt.show()
    
    return sim_data, grn_data

# %%
# Function to compare scores across different cell types for a gene
def compare_across_cell_types(gene):
    """Compare the effects of perturbing a gene across different cell types.
    
    Args:
        gene (str): The gene that was perturbed
    """
    # Filter the data for the specified gene
    sim_data = sim_scores[sim_scores['perturbed_gene'] == gene]
    grn_data = grn_scores[grn_scores['perturbed_gene'] == gene]
    
    if sim_data.empty or grn_data.empty:
        print(f"No data available for gene {gene}")
        return
    
    # Compute average scores by cell type
    sim_avg_by_cell = sim_data.groupby('cell_type')['score'].mean().reset_index()
    grn_avg_by_cell = grn_data.groupby('cell_type')['score'].mean().reset_index()
    
    # Plot comparisons
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Simulation scores by cell type
    sns.barplot(x='cell_type', y='score', data=sim_avg_by_cell, ax=ax1)
    ax1.set_title(f'Average Simulation Scores by Cell Type for {gene}')
    ax1.set_xlabel('Cell Type')
    ax1.set_ylabel('Average Score')
    ax1.tick_params(axis='x', rotation=45)
    
    # GRN scores by cell type
    sns.barplot(x='cell_type', y='score', data=grn_avg_by_cell, ax=ax2)
    ax2.set_title(f'Average GRN Scores by Cell Type for {gene}')
    ax2.set_xlabel('Cell Type')
    ax2.set_ylabel('Average Score')
    ax2.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"cell_type_comparison_{gene}.png"))
    plt.show()
    
    return sim_avg_by_cell, grn_avg_by_cell

# %%
# Function to find common target genes between simulation and GRN scores
def find_common_targets(gene, cell_type=None, top_n=50):
    """Find common target genes between simulation and GRN scores.
    
    Args:
        gene (str): The perturbed gene
        cell_type (str, optional): Filter by cell type. If None, all cell types are included.
        top_n (int): Number of top genes to consider from each dataset
    """
    # Filter the data
    sim_filter = sim_scores['perturbed_gene'] == gene
    grn_filter = grn_scores['perturbed_gene'] == gene
    
    if cell_type is not None:
        sim_filter &= (sim_scores['cell_type'] == cell_type)
        grn_filter &= (grn_scores['cell_type'] == cell_type)
    
    sim_data = sim_scores[sim_filter]
    grn_data = grn_scores[grn_filter]
    
    if sim_data.empty or grn_data.empty:
        print(f"No data available for gene {gene}")
        return None
    
    # Get top N target genes from each dataset
    top_sim_genes = set(sim_data.sort_values('score', ascending=False).head(top_n)['target_gene'])
    top_grn_genes = set(grn_data.sort_values('score', ascending=False).head(top_n)['target_gene'])
    
    # Find common genes
    common_genes = top_sim_genes.intersection(top_grn_genes)
    
    print(f"Analysis for perturbation of gene: {gene}")
    if cell_type:
        print(f"Cell type: {cell_type}")
    print(f"Top {top_n} genes in simulation scores: {len(top_sim_genes)}")
    print(f"Top {top_n} genes in GRN scores: {len(top_grn_genes)}")
    print(f"Common genes: {len(common_genes)} ({len(common_genes)/top_n*100:.2f}%)")
    print(f"Common genes: {sorted(common_genes)}")
    
    # Create a scatter plot comparing scores for common genes
    if common_genes:
        common_sim_data = sim_data[sim_data['target_gene'].isin(common_genes)]
        common_grn_data = grn_data[grn_data['target_gene'].isin(common_genes)]
        
        # Merge the datasets
        merged_data = pd.merge(common_sim_data, common_grn_data, 
                              on=['target_gene', 'perturbed_gene', 'cell_type'],
                              suffixes=('_sim', '_grn'))
        
        plt.figure(figsize=(10, 8))
        sns.scatterplot(x='score_sim', y='score_grn', data=merged_data)
        
        # Add gene labels to points
        for _, row in merged_data.iterrows():
            plt.annotate(row['target_gene'], 
                        (row['score_sim'], row['score_grn']),
                        xytext=(5, 5), textcoords='offset points')
        
        plt.title(f'Comparison of Simulation vs GRN Scores for {gene}')
        plt.xlabel('Simulation Score')
        plt.ylabel('GRN Score')
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.savefig(os.path.join(output_dir, f"score_comparison_{gene}{'_'+cell_type if cell_type else ''}.png"))
        plt.show()
    
    return common_genes

# %%
# Example usage
# Analyze a specific gene perturbation
if len(perturbed_genes) > 0:
    example_gene = perturbed_genes[0]
    print(f"\nAnalyzing example gene: {example_gene}")
    sim_data, grn_data = analyze_gene_perturbation(example_gene)
    
    # If multiple cell types exist, compare across them
    if len(cell_types) > 1:
        print(f"\nComparing {example_gene} across cell types")
        sim_avg, grn_avg = compare_across_cell_types(example_gene)
    
    # Find common target genes
    print(f"\nFinding common target genes for {example_gene}")
    common_genes = find_common_targets(example_gene, top_n=50)

# %%
# Print execution time
start_time = datetime.now()
print(f"Script started at {start_time}")
end_time = datetime.now()
print(f"Script ended at {end_time}")
print(f"Total execution time: {end_time - start_time}")
