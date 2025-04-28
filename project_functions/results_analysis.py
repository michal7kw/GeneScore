import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Function to analyze a specific gene perturbation
def analyze_gene_perturbation(sim_scores, grn_scores, gene, output_dir=None, cell_type=None):
    """Analyze the effects of perturbing a specific gene.
    
    Args:
        gene (str): The gene that was perturbed
        cell_type (str, optional): Filter by cell type. If None, all cell types are included.
    """
    # Filter the data
    sim_filter = sim_scores['goi'] == gene
    grn_filter = grn_scores['goi'] == gene
    
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
        plt.xlabel('score')
        plt.ylabel('Frequency')
        if output_dir:
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
        plt.xlabel('score')
        plt.ylabel('Frequency')
        if output_dir:
            plt.savefig(os.path.join(output_dir, f"grn_score_dist_{gene}{'_'+cell_type if cell_type else ''}.png"))
        plt.show()
    
    return sim_data, grn_data

# Function to compare scores across different cell types for a gene
def compare_across_cell_types(sim_scores, grn_scores, output_dir, gene):
    """Compare the effects of perturbing a gene across different cell types.
    
    Args:
        gene (str): The gene that was perturbed
    """
    # Filter the data for the specified gene
    sim_data = sim_scores[sim_scores['goi'] == gene]
    grn_data = grn_scores[grn_scores['goi'] == gene]
    
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

# Function to find common target genes between simulation and GRN scores
def find_common_targets(sim_scores, grn_scores, gene, cell_type=None, top_n=50, output_dir=None):
    """Find common target genes between simulation and GRN scores.
    
    Args:
        gene (str): The perturbed gene
        cell_type (str, optional): Filter by cell type. If None, all cell types are included.
        top_n (int): Number of top genes to consider from each dataset
    """
    # Filter the data
    sim_filter = sim_scores['goi'] == gene
    grn_filter = grn_scores['goi'] == gene
    
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
                              on=['target_gene', 'goi', 'cell_type'],
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
        if output_dir:
            plt.savefig(os.path.join(output_dir, f"score_comparison_{gene}{'_'+cell_type if cell_type else ''}.png"))
        plt.show()
    
    return common_genes