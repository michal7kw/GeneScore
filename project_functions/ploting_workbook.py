import fnmatch
import numpy as np
import pandas as pd
import scipy
import scipy.stats as stats

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns

import scanpy as sc

def plot_scores_for_conditions_separate(adata, conditions_set, gene_set):
    colors = {'Ag': 'red', 'Inh': 'blue', 'DMSO': 'green'}
    
    # Group conditions
    condition_groups = {}
    for condition in conditions_set:
        base_condition = condition.split('_')[0]
        if base_condition not in condition_groups:
            condition_groups[base_condition] = {'Ag': None, 'Inh': None}
        if condition.endswith('_Ag'):
            condition_groups[base_condition]['Ag'] = condition
        elif condition.endswith('_Inh'):
            condition_groups[base_condition]['Inh'] = condition

    num_groups = len(condition_groups)
    positions = np.arange(num_groups) * 2  # Adjust spacing between groups

    for plot_type in ['Ag', 'Inh']:
        plt.figure(figsize=(10, 5))
        
        for i, (base_condition, conditions) in enumerate(condition_groups.items()):
            scores = adata[adata.obs['sample_type'] == conditions[plot_type]].obs[f'gene_score_{gene_set}_{plot_type.lower()}'].values
            dmso_scores = adata[adata.obs['sample_type'] == 'DMSO'].obs[f'gene_score_{gene_set}_{plot_type.lower()}'].values

            bp = plt.boxplot([scores, dmso_scores], 
                             positions=[positions[i], positions[i]+0.5], 
                             widths=0.4, 
                             patch_artist=True)
            
            for patch, condition in zip(bp['boxes'], [plot_type, 'DMSO']):
                patch.set_facecolor(colors[condition])
                patch.set_alpha(0.5)

            # Add scatter plots
            plt.scatter(np.random.normal(positions[i], 0.05, size=len(scores)), scores, color=colors[plot_type], alpha=0.5, s=20)
            plt.scatter(np.random.normal(positions[i]+0.5, 0.05, size=len(dmso_scores)), dmso_scores, color=colors['DMSO'], alpha=0.5, s=20)

        plt.xticks(positions + 0.25, list(condition_groups.keys()), rotation=45, ha='right')
        plt.ylabel(f'Gene Score ({gene_set}_{plot_type.lower()})', fontsize=16)
        plt.title(f'{plot_type}onist Scores with corresponding DMSO', fontsize=16)
        
        # Create custom legend
        legend_elements = [Patch(facecolor=colors[plot_type], edgecolor='black', alpha=0.5, label=f'{plot_type}onist'),
                           Patch(facecolor=colors['DMSO'], edgecolor='black', alpha=0.5, label='DMSO')]
        plt.legend(handles=legend_elements, loc='best', fontsize=16)
        
        plt.tight_layout()
        plt.show()

def plot_scores_for_conditions(adata, conditions_set, gene_set):
    plt.figure(figsize=(12, 10))
    
    # Group conditions
    condition_groups = {}
    for condition in conditions_set:
        base_condition = condition.split('_')[0]
        if base_condition not in condition_groups:
            condition_groups[base_condition] = {'Ag': None, 'Inh': None}
        if condition.endswith('_Ag'):
            condition_groups[base_condition]['Ag'] = condition
        elif condition.endswith('_Inh'):
            condition_groups[base_condition]['Inh'] = condition

    num_groups = len(condition_groups)
    positions = np.arange(num_groups) * 3  # Adjust spacing between groups

    colors = {'Ag': 'red', 'Inh': 'blue', 'DMSO': 'green'}
    
    for i, (base_condition, conditions) in enumerate(condition_groups.items()):
        ag_scores = adata[adata.obs['sample_type'] == conditions['Ag']].obs[f'gene_score_{gene_set}_ag'].values
        dmso_ag_scores = adata[adata.obs['sample_type'] == 'DMSO'].obs[f'gene_score_{gene_set}_ag'].values
        inh_scores = adata[adata.obs['sample_type'] == conditions['Inh']].obs[f'gene_score_{gene_set}_inh'].values
        dmso_inh_scores = adata[adata.obs['sample_type'] == 'DMSO'].obs[f'gene_score_{gene_set}_inh'].values

        bp = plt.boxplot([ag_scores, dmso_ag_scores, inh_scores, dmso_inh_scores], 
                         positions=[positions[i], positions[i], positions[i]+1.0, positions[i]+1.0], 
                         widths=0.4, 
                         patch_artist=True)
        
        for patch, condition in zip(bp['boxes'], ['Ag', 'DMSO', 'Inh', 'DMSO']):
            patch.set_facecolor(colors[condition])
            patch.set_alpha(0.5)

        # Add scatter plots
        plt.scatter(np.random.normal(positions[i], 0.05, size=len(ag_scores)), ag_scores, color=colors['Ag'], alpha=0.5, s=20)
        plt.scatter(np.random.normal(positions[i], 0.05, size=len(dmso_ag_scores)), dmso_ag_scores, color=colors['DMSO'], alpha=0.5, s=20)
        plt.scatter(np.random.normal(positions[i]+1.0, 0.05, size=len(inh_scores)), inh_scores, color=colors['Inh'], alpha=0.5, s=20)
        plt.scatter(np.random.normal(positions[i]+1.0, 0.05, size=len(dmso_inh_scores)), dmso_inh_scores, color=colors['DMSO'], alpha=0.5, s=20)

    plt.xticks(positions + 1, list(condition_groups.keys()), rotation=45, ha='right')
    plt.ylabel(f'Gene Score ({gene_set})')
    
    # Create custom legend
    legend_elements = [Patch(facecolor=color, edgecolor='black', alpha=0.5, label=label)
                       for label, color in colors.items()]
    plt.legend(handles=legend_elements, loc='best')
    
    plt.tight_layout()
    plt.show()

def plot_scores_for_conditions_GRNs(adata, conditions_set, gene_set):
    plt.figure(figsize=(12, 10))
    
    # Group conditions
    condition_groups = {}
    for condition in conditions_set:
        base_condition = condition.split('_')[0]
        if base_condition not in condition_groups:
            condition_groups[base_condition] = {'Ag': None, 'Inh': None}
        if condition.endswith('_Ag'):
            condition_groups[base_condition]['Ag'] = condition
        elif condition.endswith('_Inh'):
            condition_groups[base_condition]['Inh'] = condition

    num_groups = len(condition_groups)
    positions = np.arange(num_groups) * 5  # Adjust spacing between groups

    colors = {'Ag': 'red', 'Inh': 'blue', 'DMSO': 'green'}
    
    for i, (base_condition, conditions) in enumerate(condition_groups.items()):
        ag_scores = adata[adata.obs['sample_type'] == conditions['Ag']].obs[f'gene_score_{gene_set}'].values
        dmso_scores = adata[adata.obs['sample_type'] == 'DMSO'].obs[f'gene_score_{gene_set}'].values
        inh_scores = adata[adata.obs['sample_type'] == conditions['Inh']].obs[f'gene_score_{gene_set}'].values

        bp = plt.boxplot([ag_scores, dmso_scores, inh_scores], 
                         positions=[positions[i], positions[i]+0.5, positions[i]+1.5], 
                         widths=0.4, 
                         patch_artist=True)
        
        for patch, condition in zip(bp['boxes'], ['Ag', 'DMSO', 'Inh', 'DMSO']):
            patch.set_facecolor(colors[condition])
            patch.set_alpha(0.5)

        # Add scatter plots
        plt.scatter(np.random.normal(positions[i], 0.05, size=len(ag_scores)), ag_scores, color=colors['Ag'], alpha=0.5, s=20)
        plt.scatter(np.random.normal(positions[i]+0.5, 0.05, size=len(dmso_scores)), dmso_scores, color=colors['DMSO'], alpha=0.5, s=20)
        plt.scatter(np.random.normal(positions[i]+1.5, 0.05, size=len(inh_scores)), inh_scores, color=colors['Inh'], alpha=0.5, s=20)

    plt.xticks(positions + 1, list(condition_groups.keys()), rotation=45, ha='right')
    plt.ylabel(f'Gene Score ({gene_set})')
    
    # Create custom legend
    legend_elements = [Patch(facecolor=color, edgecolor='black', alpha=0.5, label=label)
                       for label, color in colors.items()]
    plt.legend(handles=legend_elements, loc='best')
    
    plt.tight_layout()
    plt.show()

############################################################################# 06.2024 ####################################################################################

def boxplot_Reference_GRN_scores(adata, use_raw, control, control_condition, condition, gene_set, cell_type, goi, normalize_weights, scale_by_variance, sufix=""):
    plt.close('all')

    mask = adata.obs['sample_type'].isin(condition)
    adata_filtered = adata[mask]

    if len(adata_filtered) > 0:
        fig, ax = plt.subplots(figsize=(5, 6))
        gene_scores_dmso = adata_filtered[adata_filtered.obs['sample_type'] == condition[0]].obs[f'gene_score_{gene_set}_{cell_type}_{goi}_{control}_normalized_{normalize_weights}_scaled_{scale_by_variance}_raw_{use_raw}_cc_{control_condition}{sufix}'].values
        gene_scores_ag = adata_filtered[adata_filtered.obs['sample_type'] == condition[1]].obs[f'gene_score_{gene_set}_{cell_type}_{goi}_{control}_normalized_{normalize_weights}_scaled_{scale_by_variance}_raw_{use_raw}_cc_{control_condition}{sufix}'].values
        gene_scores_inh = adata_filtered[adata_filtered.obs['sample_type'] == condition[2]].obs[f'gene_score_{gene_set}_{cell_type}_{goi}_{control}_normalized_{normalize_weights}_scaled_{scale_by_variance}_raw_{use_raw}_cc_{control_condition}{sufix}'].values
                
        sns.boxplot(data=[gene_scores_dmso, gene_scores_ag, gene_scores_inh], notch=False,
                    boxprops=dict(alpha=0.5),
                    ax=ax)
        sns.stripplot(data=[gene_scores_dmso, gene_scores_ag, gene_scores_inh], 
                        jitter=True, color=".3", linewidth=1, ax=ax)
        
        ax.set_xticks(range(3))
        ax.set_xticklabels([condition[0], condition[1], condition[2]], fontsize=12)
        ax.set_title(f'Gene Scores - {goi} ({cell_type})', fontsize=16)
        ax.set_ylabel("Gene Score", fontsize=12)

        plt.tight_layout()
        return fig
    else:
        print(f'No data to plot for the selected condition: {condition}')

def boxplot_Reference_GRN_scores_parameters_Velmeshev(adata, condition, gene_set, cell_type, scored_gene, normalize_weights="False", scale_by_variance="False", labels=None):
    plt.close('all')

    mask = adata.obs['Og_diagnosis'].isin(condition)
    adata_filtered = adata[mask]

    if len(adata_filtered) > 0:
        fig, ax = plt.subplots(figsize=(5, 6))
        
        gene_scores_1 = adata_filtered[
            (adata_filtered.obs['Og_diagnosis'] == condition[0]) & 
            (adata_filtered.obs['cell_label_mapped'] == cell_type)
        ].obs[f'gene_score_{gene_set}_{cell_type}_{scored_gene}_labeled_{labels}_normalized_{normalize_weights}_scaled_{scale_by_variance}'].values

        gene_scores_2 = adata_filtered[
            (adata_filtered.obs['Og_diagnosis'] == condition[1]) & 
            (adata_filtered.obs['cell_label_mapped'] == cell_type)
        ].obs[f'gene_score_{gene_set}_{cell_type}_{scored_gene}_labeled_{labels}_normalized_{normalize_weights}_scaled_{scale_by_variance}'].values

        sns.boxplot(data=[gene_scores_1, gene_scores_2], notch=False,
                    boxprops=dict(alpha=0.5),
                    ax=ax)
        sns.stripplot(data=[gene_scores_1, gene_scores_2], 
                        jitter=True, color=".3", linewidth=1, ax=ax)
        
        ax.set_xticks(range(2))
        ax.set_xticklabels([condition[0], condition[1]], fontsize=12)
        ax.set_title(f'Gene Scores - {scored_gene} ({cell_type})', fontsize=16)
        ax.set_ylabel("Gene Score", fontsize=12)

        plt.tight_layout()
        return fig
    else:
        print(f'No data to plot for the selected condition: {condition}')

def boxplot_Reference_DEGs_scores_parameters_Velmeshev(adata, conditions, cell_type, key, measure, normalize_weights="False", scale_by_variance="False", labels=None):
    plt.close('all')

    mask = adata.obs['Og_diagnosis'].isin(conditions)
    adata_filtered = adata[mask]

    if len(adata_filtered) > 0:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
        
        scores = {condition: {} for condition in conditions}
        
        for condition in conditions:
            for score_type in ['ag', 'inh']:
                scores[condition][score_type] = adata_filtered[
                    (adata_filtered.obs['Og_diagnosis'] == condition) & 
                    (adata_filtered.obs['cell_label'] == cell_type)
                ].obs[f'gene_score_{measure}_{key}_{score_type}_{cell_type}_labeled_{labels}_normalized_{normalize_weights}_scaled_{scale_by_variance}'].values

        # Plot AG scores
        ag_data = [scores[condition]['ag'] for condition in conditions]
        sns.boxplot(data=ag_data, notch=False, boxprops=dict(alpha=0.5), ax=ax1)
        sns.stripplot(data=ag_data, jitter=True, color=".3", linewidth=1, ax=ax1)
        
        ax1.set_xticks(range(len(conditions)))
        ax1.set_xticklabels(conditions, fontsize=12)
        ax1.set_title(f'AG Scores - {cell_type}', fontsize=16)
        ax1.set_ylabel("Gene Score", fontsize=12)

        # Plot INH scores
        inh_data = [scores[condition]['inh'] for condition in conditions]
        sns.boxplot(data=inh_data, notch=False, boxprops=dict(alpha=0.5), ax=ax2)
        sns.stripplot(data=inh_data, jitter=True, color=".3", linewidth=1, ax=ax2)
        
        ax2.set_xticks(range(len(conditions)))
        ax2.set_xticklabels(conditions, fontsize=12)
        ax2.set_title(f'INH Scores - {cell_type}', fontsize=16)
        ax2.set_ylabel("Gene Score", fontsize=12)

        plt.tight_layout()
        return fig
    else:
        print(f'No data to plot for the selected conditions: {conditions}')

def boxplot_Reference_GRN_scores_parameters(adata, conditions, gene_set, cell_type, goi, control=True, normalize_weights="False", scale_by_variance="False", scaling_only_based_on_control=True, labels=None, prefix=""):
    plt.close('all')

    mask = adata.obs['sample_type'].isin(conditions)
    adata_filtered = adata[mask]

    if len(adata_filtered) > 0:
        fig, ax = plt.subplots(figsize=(5, 6))
        
        gene_scores_dmso = adata_filtered[adata_filtered.obs['sample_type'] == conditions[0]].obs[f'{prefix}gene_score_{gene_set}_{cell_type}_{goi}_control_{control}_labeled_{labels}_normalized_{normalize_weights}_scaled_{scale_by_variance}_ctrlscale_{scaling_only_based_on_control}'].values
        gene_scores_ag = adata_filtered[adata_filtered.obs['sample_type'] == conditions[1]].obs[f'{prefix}gene_score_{gene_set}_{cell_type}_{goi}_control_{control}_labeled_{labels}_normalized_{normalize_weights}_scaled_{scale_by_variance}_ctrlscale_{scaling_only_based_on_control}'].values
        gene_scores_inh = adata_filtered[adata_filtered.obs['sample_type'] == conditions[2]].obs[f'{prefix}gene_score_{gene_set}_{cell_type}_{goi}_control_{control}_labeled_{labels}_normalized_{normalize_weights}_scaled_{scale_by_variance}_ctrlscale_{scaling_only_based_on_control}'].values
                
        sns.boxplot(data=[gene_scores_dmso, gene_scores_ag, gene_scores_inh], notch=False,
                    boxprops=dict(alpha=0.5),
                    ax=ax)
        sns.stripplot(data=[gene_scores_dmso, gene_scores_ag, gene_scores_inh], 
                        jitter=True, color=".3", linewidth=1, ax=ax)
        
        ax.set_xticks(range(3))
        ax.set_xticklabels([conditions[0], conditions[1], conditions[2]], fontsize=12)
        ax.set_title(f'Gene Scores - {goi} ({cell_type})', fontsize=16)
        ax.set_ylabel("Gene Score", fontsize=12)

        plt.tight_layout()
        return fig
    else:
        print(f'No data to plot for the selected condition: {conditions}')


def boxplot_Reference_GRN_scores_parameters_all_gois(adata, gene_set, condition, cell_type, gois, control_condition, control=True, normalize_weights="False", scale_by_variance="False", scaling_only_based_on_control=True, prefix=""):
    plt.close('all')

    nrows = 2
    ncols = (len(gois) + 1) // 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 6 * nrows), sharey=False)
    axes = axes.flatten() if len(gois) > 1 else [axes]

    # Ensure condition is a list
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
            
            gene_scores = [adata_filtered[adata_filtered.obs['sample_type'] == cond].obs[selection].values for cond in conditions]

            sns.boxplot(data=gene_scores, notch=False, boxprops=dict(alpha=0.5), ax=ax)
            sns.stripplot(data=gene_scores, jitter=True, color=".3", linewidth=1, ax=ax)
            
            ax.set_xticks(range(len(conditions)))
            ax.set_xticklabels(conditions, fontsize=12, rotation=45, ha='right')
            ax.set_title(f'Gene: {goi}', fontsize=14)
            
            # Set y-axis label
            ax.set_ylabel("Gene Score", fontsize=12)
            
            # Automatically adjust y-axis limits
            all_scores = np.concatenate(gene_scores)
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


def boxplot_Reference_GRN_scores_parameters_all_conditions(adata, conditions, gene_set, cell_type, goi, control_condition, control=True, normalize_weights="False", scale_by_variance="False", scaling_only_based_on_control=True, prefix=""):
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
                scores = adata_filtered[adata_filtered.obs['sample_type'] == cond].obs[
                    selection
                ].values
                gene_scores.append(scores)
            
            sns.boxplot(data=gene_scores, notch=False, boxprops=dict(alpha=0.5), ax=ax)
            sns.stripplot(data=gene_scores, jitter=True, color=".3", linewidth=1, ax=ax)
            
            ax.set_xticks(range(len(condition)))
            ax.set_xticklabels(condition, fontsize=14, rotation=45, ha='right')
            ax.set_title(f'{condition[1].split("_")[0]}', fontsize=14)
            
            # Set y-axis label for all subplots
            ax.set_ylabel("Gene Score", fontsize=12)
            
            # Automatically adjust y-axis limits for each subplot
            y_min = min([min(scores) for scores in gene_scores])
            y_max = max([max(scores) for scores in gene_scores])
            y_range = y_max - y_min
            ax.set_ylim(y_min - 0.1 * y_range, y_max + 0.1 * y_range)

        else:
            print(f'No data to plot for the condition: {condition}')

    fig.suptitle(f'Gene Scores - {goi} ({cell_type})', fontsize=16)
    plt.tight_layout()
    return fig

def boxplot_EDCs_GRN_scores(adata, condition, gene_set, cell_type, scored_gene):
    plt.close('all')

    # Get all concentrations for the selected condition
    condition_concentrations = adata.obs['condition_concentraion'].unique()
    selected_concentrations = [conc for conc in condition_concentrations if condition in conc]
    all_concentrations = ['DMSO_0.1'] + selected_concentrations

    mask = adata.obs['condition_concentraion'].isin(all_concentrations)
    adata_filtered = adata[mask]

    if adata_filtered.shape[0] > 0:
        fig, ax = plt.subplots(figsize=(12, 8))

        gene_scores = {}

        for conc in all_concentrations:
            gene_scores[conc] = adata_filtered[adata_filtered.obs['condition_concentraion'] == conc].obs[f'gene_score_{gene_set}_{cell_type}_{scored_gene}'].values.tolist()

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
        ax.set_title(f'Gene Scores - {scored_gene}\n{gene_set} - {cell_type}', fontsize=18)
        ax.set_ylabel("Gene Score", fontsize=18)

        plt.tight_layout()
        plt.show()
    else:
        print(f'No data to plot for the selected conditions: {all_concentrations}')

def boxplot_EDCs_GRN_scores_parameters(
        adata, 
        conditions, 
        gene_set, 
        cell_type, 
        goi, 
        control=True,
        normalize_weights="False", 
        scale_by_variance="False", 
        labels=None,
        scaling_only_based_on_control=True
        ):
    plt.close('all')

    # Get all concentrations for the selected conditions
    condition_concentrations = adata.obs['condition_concentraion'].unique()
    selected_concentrations = [conc for conc in condition_concentrations if conditions in conc]
    all_concentrations = ['DMSO_0.1'] + selected_concentrations

    mask = adata.obs['condition_concentraion'].isin(all_concentrations)
    adata_filtered = adata[mask]

    if adata_filtered.shape[0] > 0:
        fig, ax = plt.subplots(figsize=(12, 8))

        gene_scores = {}

        for conc in all_concentrations:
            gene_scores[conc] = adata_filtered[adata_filtered.obs['condition_concentraion'] == conc].obs[f'gene_score_{gene_set}_{cell_type}_{goi}_control_{control}_labeled_{labels}_normalized_{normalize_weights}_scaled_{scale_by_variance}_ctrlscale_{scaling_only_based_on_control}'].values.tolist()

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
        ax.set_title(f'Gene Scores - {goi}\n{gene_set} - {cell_type}', fontsize=18)
        ax.set_ylabel("Gene Score", fontsize=18)

        plt.tight_layout()
        plt.show()
    else:
        print(f'No data to plot for the selected conditions: {all_concentrations}')


def boxplot_Reference_GRN_scores_parameters_sim(adata, conditions, gene_set, control, cell_type, goi, scaling_only_based_on_control, normalize_weights="False", scale_by_variance="False", labels=None, prefix=""):
    plt.close('all')

    mask = adata.obs['sample_type'].isin(conditions)
    adata_filtered = adata[mask]

    if len(adata_filtered) > 0:
        fig, ax = plt.subplots(figsize=(5, 6))
        
        gene_scores_dmso = adata_filtered[adata_filtered.obs['sample_type'] == conditions[0]].obs[f'{prefix}gene_score_{gene_set}_{cell_type}_{goi}_control_{control}_labeled_{labels}_normalized_{normalize_weights}_scaled_{scale_by_variance}_ctrlscale_{scaling_only_based_on_control}'].values

        gene_scores_ag = adata_filtered[adata_filtered.obs['sample_type'] == conditions[1]].obs[f'{prefix}gene_score_{gene_set}_{cell_type}_{goi}_control_{control}_labeled_{labels}_normalized_{normalize_weights}_scaled_{scale_by_variance}_ctrlscale_{scaling_only_based_on_control}'].values
        gene_scores_inh = adata_filtered[adata_filtered.obs['sample_type'] == conditions[2]].obs[f'{prefix}gene_score_{gene_set}_{cell_type}_{goi}_control_{control}_labeled_{labels}_normalized_{normalize_weights}_scaled_{scale_by_variance}_ctrlscale_{scaling_only_based_on_control}'].values
                
        sns.boxplot(data=[gene_scores_dmso, gene_scores_ag, gene_scores_inh], notch=False,
                    boxprops=dict(alpha=0.5),
                    ax=ax)
        sns.stripplot(data=[gene_scores_dmso, gene_scores_ag, gene_scores_inh], 
                        jitter=True, color=".3", linewidth=1, ax=ax)
        
        ax.set_xticks(range(3))
        ax.set_xticklabels([conditions[0], conditions[1], conditions[2]], fontsize=12)
        ax.set_title(f'Gene Scores - {goi} ({cell_type})', fontsize=16)
        ax.set_ylabel("Gene Score", fontsize=12)

        plt.tight_layout()
        return fig
    else:
        print(f'No data to plot for the selected conditions: {conditions}')

def visualize_weight_distributions(gene_sets_dict_cell_type_first, set_selected, cell_type_selected, scored_genes, score_type="scored_coef_mean"):
    data = []
    labels = []
    
    for scored_gene in scored_genes:
        weights = gene_sets_dict_cell_type_first[set_selected][cell_type_selected][scored_gene][score_type]
        data.append(weights)
        labels.append(scored_gene)
    
    # Create the box plot
    fig, ax = plt.subplots(figsize=(12, 8))
    bp = ax.boxplot(data, labels=labels, patch_artist=True)
    
    # Customize the plot
    ax.set_title(f'Distribution of Weights for Different Genes\nSet: {set_selected}, Cell Type: {cell_type_selected}')
    ax.set_xlabel('Genes')
    ax.set_ylabel('Weights')
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right')
    
    # Add grid for easier comparison
    ax.yaxis.grid(True)
    
    # Color each box
    colors = plt.cm.viridis(np.linspace(0, 1, len(scored_genes)))
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    
    plt.tight_layout()
    plt.show()

def visualize_weight_distributions_DEGs(gene_sets_dict_cell_type_first, scored_genes, score_type="coef_mean"):
    data = []
    labels = []
    
    for scored_gene in scored_genes:
        weights = gene_sets_dict_cell_type_first[scored_gene][score_type]
        data.append(weights)
        labels.append(scored_gene)
    
    # Create the box plot
    fig, ax = plt.subplots(figsize=(12, 8))
    bp = ax.boxplot(data, labels=labels, patch_artist=True)
    
    # Customize the plot
    ax.set_xlabel('Genes')
    ax.set_ylabel('Weights')
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right')
    
    # Add grid for easier comparison
    ax.yaxis.grid(True)
    
    # Color each box
    colors = plt.cm.viridis(np.linspace(0, 1, len(scored_genes)))
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    
    plt.tight_layout()
    plt.show()

################################################################################### Legacy ###########################################################################

def update_plots_boxplot_EDCs_with_scatter_separated_selection(adata, celltype, gene_set, gene, conditions):
    conditions = ['DMSO', conditions]
    plt.close('all')

    if isinstance(conditions, str):
        conditions = [conditions]

    mask = adata.obs['Condition'].isin(conditions)
    adata_filtered = adata[mask]

    if adata_filtered.shape[0] > 0:
        fig, ax = plt.subplots(figsize=(6, 5))

        gene_scores_dmso = []
        gene_scores = {}

        for condition in adata_filtered.obs.condition_concentraion.unique():
            if condition == "DMSO_0.1":
                gene_scores_dmso = adata_filtered[adata_filtered.obs['condition_concentraion'] == condition].obs[f'gene_score_{gene}'].values.tolist()
            else:
                gene_scores[condition] = adata_filtered[adata_filtered.obs['condition_concentraion'] == condition].obs[f'gene_score_{gene}'].values.tolist()

        conditions_sorted = ["DMSO_0.1"] + list(gene_scores.keys())
        data = [gene_scores_dmso] + [gene_scores.get(cond, []) for cond in gene_scores.keys()]

        positions = np.arange(len(conditions_sorted))
        width = 0.5

        bp = ax.boxplot(data, positions=positions, widths=width, patch_artist=True,
                        boxprops=dict(facecolor='C0', color='None', alpha=0.5),
                        whiskerprops=dict(color='C0'),
                        capprops=dict(color='C0'),
                        showfliers=False,
                        medianprops=dict(color='white'))

        for i, scores in enumerate(data):
            if scores:  # Check if scores is not empty
                ax.scatter([positions[i]] * len(scores), scores, color='C0', alpha=1.0, s=20)

        ax.set_xticks(positions)
        ax.set_xticklabels(conditions_sorted, rotation=45, ha='right')
        ax.set_title(f'Gene Scores - {gene}', fontsize=16)
        ax.set_ylabel("Gene Score", fontsize=12)

        plt.tight_layout()
        plt.show()
    else:
        print(f'No data to plot for the selected conditions: {conditions}')

def update_plots_receptors(adata, condition, gene_set, receptor, methods):
    plt.close('all')
    fig, axes = plt.subplots(1, len(methods), figsize=(14, 5))
    
    mask = adata.obs['sample_type'].isin(condition)
    adata_filtered = adata[mask]

    if adata_filtered.shape[0] > 0:
        for j, method in enumerate(methods):
            ax = axes[j] 
            expression = adata_filtered[:, receptor].X
            expression = np.array(expression).flatten()

            gene_score = adata_filtered.obs[f'{method}_gene_score_{gene_set}']
            gene_score = np.array(gene_score).flatten()

            sns.scatterplot(x=expression, y=gene_score, hue=adata_filtered.obs['sample_type'], palette=['red', 'green', 'blue'], ax=ax)
            ax.set_xlabel(f'Expression of {receptor}', fontsize=8)
            ax.set_ylabel('Gene Score', fontsize=8)
            ax.legend(title='Sample Type', fontsize=8)
            ax.set_title(method, fontsize=12)
        plt.tight_layout()
        fig.suptitle(f'{gene_set}', fontsize=16, y=1.1)
        plt.show()
    else:
        print(f'No data to plot for the selected condition: {condition}')

def update_plots_boxplot(adata, condition, gene, debug=False):
    plt.close('all')
    mask = adata.obs['sample_type'].isin(condition)
    adata_filtered = adata[mask]

    if adata_filtered.shape[0] > 0:
        fig, ax = plt.subplots(figsize=(6, 5))
        
        if debug: print(f'Condition: {condition}')
        gene_scores_dmso = adata_filtered[adata_filtered.obs['sample_type'] == condition[0]].obs[f'gene_score_{gene}'].values
        gene_scores_ag = adata_filtered[adata_filtered.obs['sample_type'] == condition[1]].obs[f'gene_score_{gene}'].values
        gene_scores_inh = adata_filtered[adata_filtered.obs['sample_type'] == condition[2]].obs[f'gene_score_{gene}'].values
        if debug: print(f'DM so: {gene_scores_dmso} \n AG: {gene_scores_ag} \n INH: {gene_scores_inh}')

        sns.boxplot(data=[gene_scores_dmso, gene_scores_ag, gene_scores_inh], notch=False,
                        boxprops=dict(alpha=0.5),
                        ax=ax)
        sns.stripplot(data=[gene_scores_dmso, gene_scores_ag, gene_scores_inh], 
            jitter=True, color=".3", linewidth=1, ax=ax)
        
        ax.set_xticks(range(3))
        ax.set_xticklabels([condition[0], condition[1], condition[2]], fontsize=12)
        ax.set_title(f'Gene Scores', fontsize=16)
        ax.set_ylabel("Gene Score", fontsize=12)
        plt.tight_layout()
        # plt.show()
        return fig
    else:
        print(f'No data to plot for the selected condition: {condition}')

def update_plots_boxplot_compare_sets(adata, condition, gene_set_1, gene_set_2):
    plt.close('all')
    mask = adata.obs['sample_type'].isin(condition)
    adata_filtered = adata[mask]

    if adata_filtered.shape[0] > 0:
        fig, ax = plt.subplots(1, 2, figsize=(12, 5))

        for i, gene_set in enumerate([gene_set_1, gene_set_2]):
            gene_scores_dmso = adata_filtered[adata_filtered.obs['sample_type'] == condition[0]].obs[f'gene_score_{gene_set}'].values
            gene_scores_ag = adata_filtered[adata_filtered.obs['sample_type'] == condition[1]].obs[f'gene_score_{gene_set}'].values
            gene_scores_inh = adata_filtered[adata_filtered.obs['sample_type'] == condition[2]].obs[f'gene_score_{gene_set}'].values

            # _, p_value_dmso_ag = stats.mannwhitneyu(gene_scores_dmso, gene_scores_ag)
            # _, p_value_dmso_inh = stats.mannwhitneyu(gene_scores_dmso, gene_scores_inh)

            sns.boxplot(data=[gene_scores_dmso, gene_scores_ag, gene_scores_inh], notch=False,
                        boxprops=dict(alpha=0.5),
                        ax=ax[i])
            sns.stripplot(data=[gene_scores_dmso, gene_scores_ag, gene_scores_inh], 
            jitter=True, color=".3", linewidth=1, ax=ax[i])
            ax[i].set_xticks(range(3))
            ax[i].set_xticklabels([condition[0], condition[1], condition[2]], fontsize=12)
            ax[i].set_title(f'Gene Scores - {gene_set}', fontsize=16)
            ax[i].set_ylabel("Gene Score", fontsize=12)

            # Add p-values as text below each subplot
            # ax[i].text(0.3, -0.2, f'p={p_value_dmso_ag:.3e}', ha='center', va='top', transform=ax[i].transAxes, fontsize=16)
            # ax[i].text(0.7, -0.2, f'p={p_value_dmso_inh:.3e}', ha='center', va='top', transform=ax[i].transAxes, fontsize=16)

        plt.tight_layout()
        plt.show()
    else:
        print(f'No data to plot for the selected condition: {condition}')

def update_plots_boxplot_Herring_Bhadur(adata, condition, gene_set, debug=False):
    plt.close('all')
    mask = adata.obs['sample_type'].isin(condition)
    adata_filtered = adata[mask]

    if adata_filtered.shape[0] > 0:
        studies = ['herring', 'bhadur']
        num_studies = len(studies)
        fig, axes = plt.subplots(2, num_studies, figsize=(6 * num_studies, 10))
        
        if num_studies == 1:
            axes = [axes]  # Convert single Axes object to a list
        
        if debug: print(f'Condition: {condition}')
        for i, study in enumerate(studies):
            gene_scores_dmso = adata_filtered[adata_filtered.obs['sample_type'] == condition[0]].obs[f'HT_gene_score_{gene_set}_{study}'].values
            gene_scores_ag = adata_filtered[adata_filtered.obs['sample_type'] == condition[1]].obs[f'HT_gene_score_{gene_set}_{study}'].values
            gene_scores_inh = adata_filtered[adata_filtered.obs['sample_type'] == condition[2]].obs[f'HT_gene_score_{gene_set}_{study}'].values
            if debug: print(f'DMSO: {gene_scores_dmso} \n AG: {gene_scores_ag} \n INH: {gene_scores_inh}')

            _, p_value_dmso_ag = stats.mannwhitneyu(gene_scores_dmso, gene_scores_ag)
            _, p_value_dmso_inh = stats.mannwhitneyu(gene_scores_dmso, gene_scores_inh)
            
            # Plot for Ag
            ax_ag = axes[0, i]
            sns.boxplot(data=[gene_scores_dmso, gene_scores_ag], notch=False,
                        boxprops=dict(alpha=0.5),
                        ax=ax_ag)
            sns.stripplot(data=[gene_scores_dmso, gene_scores_ag], 
                          jitter=True, color=".3", linewidth=1, ax=ax_ag)
            
            ax_ag.set_xticks(range(2))
            ax_ag.set_xticklabels([condition[0], condition[1]], fontsize=12)
            ax_ag.set_title(f'{study.capitalize()} - Gene Scores (Ag)', fontsize=16)
            ax_ag.set_ylabel("Gene Score", fontsize=12)
            ax_ag.text(0.7, -0.2, f'p={p_value_dmso_ag:.3e}', ha='center', va='top', transform=ax_ag.transAxes, fontsize=16)
            
            # Plot for Inh
            ax_inh = axes[1, i]
            sns.boxplot(data=[gene_scores_dmso, gene_scores_inh], notch=False,
                        boxprops=dict(alpha=0.5),
                        ax=ax_inh)
            sns.stripplot(data=[gene_scores_dmso, gene_scores_inh], 
                          jitter=True, color=".3", linewidth=1, ax=ax_inh)
            
            ax_inh.set_xticks(range(2))
            ax_inh.set_xticklabels([condition[0], condition[2]], fontsize=12)
            ax_inh.set_title(f'{study.capitalize()} - Gene Scores (Inh)', fontsize=16)
            ax_inh.set_ylabel("Gene Score", fontsize=12)
            ax_inh.text(0.7, -0.2, f'p={p_value_dmso_inh:.3e}', ha='center', va='top', transform=ax_inh.transAxes, fontsize=16)
        
        plt.tight_layout()
        plt.show()
    else:
        print(f'No data to plot for the selected condition: {condition}')

def update_plots_receptors_vs_score(adata, receptors, condition, gene_set):
    plt.close('all')
    
    fig, axes = plt.subplots(1, len(receptors), figsize=(len(receptors)*4, 5))
    mask = adata.obs['sample_type'].isin(condition)
    adata_filtered = adata[mask]

    if adata_filtered.shape[0] > 0:
        for j, gene in enumerate(receptors):
            ax = axes[j] 
            expression = adata_filtered[:, gene].X
            expression = np.array(expression).flatten()

            gene_score = adata_filtered.obs[f'gene_score_{gene_set}']
            gene_score = np.array(gene_score).flatten()

            sns.scatterplot(x=expression, y=gene_score, hue=adata_filtered.obs['sample_type'], palette=['red', 'green', 'blue'], ax=ax)
            ax.set_xlabel(f'Expression of {gene}', fontsize=12)
            ax.set_ylabel('Gene Score', fontsize=12)
            ax.legend(title='Sample Type', fontsize=12)
        plt.tight_layout()
        fig.suptitle(f'{gene_set} - gene set', fontsize=16, y=1.1)
        plt.show()
    else:
        print(f'No data to plot for the selected condition: {condition}')

def update_plots_gene_sets_comparison(adata, receptors, condition, gene_set_1, gene_set_2):
    plt.close('all')  
    fig, axes = plt.subplots(2, len(receptors), figsize=(14, 10))  
    
    for idx, gene_set in enumerate([gene_set_1, gene_set_2]):
        mask = adata.obs['sample_type'].isin(condition)
        adata_filtered = adata[mask]
        
        if adata_filtered.shape[0] > 0:
            for j, gene in enumerate(receptors):
                ax = axes[idx, j]
                expression = adata_filtered[:, gene].X
                expression = np.array(expression).flatten() 
                
                gene_score = adata_filtered.obs[f'gene_score_{gene_set}']
                gene_score = np.array(gene_score).flatten() 
                
                sns.scatterplot(x=expression, y=gene_score, hue=adata_filtered.obs['sample_type'], palette=['red', 'green', 'blue'], ax=ax)
                ax.set_xlabel(f'Expression of {gene}', fontsize=8)
                ax.set_ylabel(f'Gene Score', fontsize=8)
                ax.legend(title='Sample Type', fontsize=8)
            axes[idx, 0].set_title(f'{gene_set} - gene set', fontsize=12, pad=10)    
        else:
            print(f'No data to plot for the selected condition: {condition}')

    plt.tight_layout()
    plt.show()

def update_plots_regression(adata, receptors, condition, gene_set):
    plt.close('all')
    # print(f"Plotting for {condition} and {gene_set}")
    fig, axes = plt.subplots(1, len(receptors), figsize=(len(receptors) * 4, 4))
    
    mask = adata.obs['sample_type'].isin(condition)
    adata_filtered = adata[mask]

    condition_title = condition[1].split('_')[0]

    if adata_filtered.shape[0] > 0:
        results_per_gene = {}
        for j, gene in enumerate(receptors):
            ax = axes[j]
            
            expression = adata_filtered[:, gene].X
            expression = np.array(expression).flatten()
            
            gene_score = adata_filtered.obs[f'gene_score_{gene_set}']
            gene_score = np.array(gene_score).flatten()
 
            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(gene_score, expression)
            results_per_gene[gene] = {
                'slope': slope,
                'intercept': intercept,
                'r_value': r_value,
                'p_value': p_value,
                'std_err': std_err
            }
        
            sns.scatterplot(y=gene_score, x=expression, hue=adata_filtered.obs['sample_type'], palette=['red', 'green', 'blue'], ax=ax)
            if results_per_gene[gene]['p_value'] < 0.05:
                sns.regplot(y=gene_score, x=expression, ax=ax, scatter=False, color='red', label='Significant Fit')
            ax.set_xlabel(f'Expression of {gene}', fontsize=12)
            ax.set_ylabel('Gene Score', fontsize=12)
            ax.legend(title='Sample Type', fontsize=8)
            
        plt.tight_layout()
        fig.suptitle(f'gene set: {gene_set}, condition: {condition_title}', fontsize=16, y=1.1)
        plt.show()
    else:
        print(f'No data to plot for the selected condition: {condition}')

def update_plots_regression_matching(adata, gene_sets, condition):
    plt.close('all')
    
    num_gene_sets = len(gene_sets)
    fig, axes = plt.subplots(1, num_gene_sets, figsize=(num_gene_sets * 4, 4))
    
    mask = adata.obs['sample_type'].isin(condition)
    adata_filtered = adata[mask]

    if adata_filtered.shape[0] > 0:
        for i, (gene_set, receptor) in enumerate(gene_sets.items()):
            ax = axes[i] if num_gene_sets > 1 else axes
            
            expression = adata_filtered[:, receptor].X
            expression = np.array(expression).flatten()
            
            gene_score = adata_filtered.obs[f'gene_score_{gene_set}']
            gene_score = np.array(gene_score).flatten()
 
            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(gene_score, expression)
            
            sns.scatterplot(y=gene_score, x=expression, hue=adata_filtered.obs['sample_type'], palette=['red', 'green', 'blue'], ax=ax)
            if p_value < 0.05:
                sns.regplot(x=gene_score, y=expression, ax=ax, scatter=False, color='red', label='Significant Fit')
            ax.set_xlabel(f'{gene_set} Gene Score', fontsize=10)
            ax.set_ylabel(f'{receptor} Expression', fontsize=10)
            ax.set_title(f'{gene_set} - Gene Set \np-value: {p_value:.4f}, r: {r_value:.4f}', fontsize=12)
            ax.legend(title='Sample Type', fontsize=8)
            
        plt.tight_layout()
        fig.suptitle(f'Gene Set Scores vs Receptor Expression', fontsize=16, y=1.05)
        plt.show()
    else:
        print(f'No data to plot for the selected condition: {condition}')

def update_plots_exposure_score(adata, receptors, condition, gene_set):
    plt.close('all')
    
    mask = adata.obs['sample_type'].isin(condition)
    adata_filtered = adata[mask]

    if adata_filtered.shape[0] > 0:
        gene_scores_dmso = adata_filtered[adata_filtered.obs['sample_type'] == condition[0]].obs[f'gene_score_{gene_set}'].values
        gene_scores_ag = adata_filtered[adata_filtered.obs['sample_type'] == condition[1]].obs[f'gene_score_{gene_set}'].values
        gene_scores_inh = adata_filtered[adata_filtered.obs['sample_type'] == condition[2]].obs[f'gene_score_{gene_set}'].values

        print(f'GENE SET: {gene_set}')
        u_statistic, p_value_dmso_ag = stats.mannwhitneyu(gene_scores_dmso, gene_scores_ag)
        print(f'{condition[0]} vs {condition[1]}:')
        print(f'U Statistic: {u_statistic}, P-value: {p_value_dmso_ag}')

        u_statistic, p_value_dmso_inh = stats.mannwhitneyu(gene_scores_dmso, gene_scores_inh)
        print(f'{condition[0]} vs {condition[2]}:')
        print(f'U Statistic: {u_statistic}, P-value: {p_value_dmso_inh}')

        plt.figure(figsize=(6, 3))

        plt.subplot(1, 2, 1)
        # sns.boxplot(data=[gene_scores_dmso, gene_scores_ag], notch=True)
        sns.boxplot(data=[gene_scores_dmso, gene_scores_ag], notch=False,
                            boxprops=dict(alpha=0.5))
        sns.stripplot(data=[gene_scores_dmso, gene_scores_ag], 
            jitter=True, color=".3", linewidth=1)
        
        plt.xticks([0, 1], [condition[0], condition[1]])
        plt.title(f'Gene Scores: {condition[0]} vs {condition[1]}:')
        plt.ylabel("Gene Score")

        plt.subplot(1, 2, 2)
        # sns.boxplot(data=[gene_scores_dmso, gene_scores_inh], notch=True)
        sns.boxplot(data=[gene_scores_dmso, gene_scores_inh], notch=False,
                            boxprops=dict(alpha=0.5))
        sns.stripplot(data=[gene_scores_dmso, gene_scores_inh], 
            jitter=True, color=".3", linewidth=1)
        
        plt.xticks([0, 1], [condition[0], condition[2]])
        plt.title(f'Gene Scores: {condition[0]} vs {condition[2]}:')
        plt.ylabel("Gene Score")

        plt.tight_layout()
        plt.show()
    else:
        print(f'No data to plot for the selected condition: {condition}')

def plot_scores_for_conditions_old(adata, conditions, title_suffix, gene_set):
    plt.figure(figsize=(10, 6))
    
    scores_list = []
    labels_list = []

    for condition in conditions:
        scores = adata[adata.obs['sample_type'] == condition].obs[f'gene_score_{gene_set}'].values
        scores_list.append(scores)
        labels_list.append(condition)
        
    sns.boxplot(data=scores_list, notch=False,
                            boxprops=dict(alpha=0.5))
    sns.stripplot(data=scores_list, 
        jitter=True, color=".3", linewidth=1)
    
    plt.xticks(range(len(conditions)), labels=labels_list, rotation=45)
    plt.title(f'Gene Scores for {title_suffix}')
    plt.tight_layout()
    plt.show()

def update_plots_each_exposure_score(adata, conditions_set, gene_set):
    conditions_inh = [cond for cond in conditions_set if cond.endswith('_Inh')]
    conditions_ag = [cond for cond in conditions_set if cond.endswith('_Ag')]
    conditions_inh.append('DMSO')
    conditions_ag.append('DMSO')

    # Step 2: Plot scores for each category
    plot_scores_for_conditions_old(adata, conditions_inh, 'Inhibitors', gene_set)
    plot_scores_for_conditions_old(adata, conditions_ag, 'Agonists', gene_set)

def inflate_gene_expression_raw(adata, adjustment_factor=1.5, adjustment_percentage=10):
    adata_copy = adata.copy()
    
    # Get the total number of genes in the dataset
    total_genes = adata_copy.var_names.shape[0]
    
    # Calculate the number of genes to adjust (both inflate and deflate)
    num_genes_to_adjust = int(np.ceil(total_genes * (adjustment_percentage / 100)))
    if num_genes_to_adjust < 1:
        num_genes_to_adjust = 1
    
    # Select subsets of genes randomly for inflation and deflation
    genes_to_adjust = np.random.choice(adata_copy.var_names, size=num_genes_to_adjust*2, replace=False)
    genes_to_inflate = genes_to_adjust[:num_genes_to_adjust]
    genes_to_deflate = genes_to_adjust[num_genes_to_adjust:]
    
    # Inflate the expression of the selected genes
    n_inflated = 0
    for gene in genes_to_inflate:
        adata_copy[:, gene].X *= adjustment_factor
        n_inflated += 1
    
    # Deflate the expression of another set of selected genes by the reciprocal of the inflation factor
    n_deflated = 0
    deflation_factor = 1 / adjustment_factor
    for gene in genes_to_deflate:
        adata_copy[:, gene].X *= deflation_factor
        n_deflated += 1
            
    print(f'num_genes_inflated: {n_inflated}, num_genes_deflated: {n_deflated}, out of: {num_genes_to_adjust*2}/{total_genes}')
    
    return adata_copy

def inflate_gene_expression_log(adata, adjustment_factor=1.5, adjustment_percentage=10):
    adata_copy = adata.copy()
    
    # Since data is log-normalized, convert adjustment factors to log-space
    log_adjustment_factor = np.log(adjustment_factor)
    log_deflation_factor = -log_adjustment_factor  # This ensures a reciprocal effect in log-space
    
    # Get the total number of genes in the dataset
    total_genes = adata_copy.var_names.shape[0]
    
    # Calculate the number of genes to adjust (both inflate and deflate)
    num_genes_to_adjust = int(np.ceil(total_genes * (adjustment_percentage / 100)))
    if num_genes_to_adjust < 1:
        num_genes_to_adjust = 1
    
    # Select subsets of genes randomly for inflation and deflation
    genes_to_adjust = np.random.choice(adata_copy.var_names, size=num_genes_to_adjust*2, replace=False)
    genes_to_inflate = genes_to_adjust[:num_genes_to_adjust]
    genes_to_deflate = genes_to_adjust[num_genes_to_adjust:]
    
    # Apply adjustments in log-space for selected genes
    adata_copy[:, genes_to_inflate].X += log_adjustment_factor
    adata_copy[:, genes_to_deflate].X += log_deflation_factor
    
    print(f'num_genes_inflated: {num_genes_to_adjust}, num_genes_deflated: {num_genes_to_adjust}, out of: {num_genes_to_adjust*2}/{total_genes}')
    
    return adata_copy

def update_plots_genes_noise(adata, conditions_set, gene_sets, gene_set):
    print(gene_set)
    num_conditions = len(conditions_set)
    n_columns = 3
    num_rows = np.ceil(num_conditions / n_columns).astype(int)
    fig, axes = plt.subplots(num_rows, n_columns, figsize=(20, 3 * num_rows))
    axes = axes.flatten()

    valid_genes = list(set(gene_sets[gene_set]['targets']) & set(adata.var_names))

    # Calculate original data scores
    # sc.tl.score_genes(adata, gene_list=gene_sets[gene_set], score_name=f'gene_score_{gene_set}_original')
    sc.tl.score_genes(adata, gene_list=valid_genes, score_name=f'gene_score_{gene_set}_original')

    # Create a copy of the data and artificially inflate gene expression
    adata_modified = inflate_gene_expression_log(adata, adjustment_factor=1.5, adjustment_percentage=10)

    # Recalculate the score with modified data
    # sc.tl.score_genes(adata_modified, gene_list=gene_sets[gene_set], score_name=f'gene_score_{gene_set}_modified')
    sc.tl.score_genes(adata_modified, gene_list=valid_genes, score_name=f'gene_score_{gene_set}_modified')


    for i, condition in enumerate(conditions_set):
        # Get original and modified data scores for DMSO
        gene_scores_dmso_original = adata[adata.obs['sample_type'] == "DMSO"].obs[f'gene_score_{gene_set}_original'].values
        gene_scores_dmso_modified = adata_modified[adata_modified.obs['sample_type'] == "DMSO"].obs[f'gene_score_{gene_set}_modified'].values

        # Get original and modified data scores for the current condition
        gene_scores_condition_original = adata[adata.obs['sample_type'] == condition].obs[f'gene_score_{gene_set}_original'].values
        gene_scores_condition_modified = adata_modified[adata_modified.obs['sample_type'] == condition].obs[f'gene_score_{gene_set}_modified'].values

        # Plotting
        ax = axes[i]  # Access the correct subplot
        # sns.boxplot(data=[gene_scores_dmso_original, gene_scores_dmso_modified, gene_scores_condition_original, gene_scores_condition_modified], ax=ax)
        sns.boxplot(data=[gene_scores_dmso_original, gene_scores_dmso_modified, gene_scores_condition_original, gene_scores_condition_modified], notch=False,
                            boxprops=dict(alpha=0.5),
                            ax=ax)
        sns.stripplot(data=[gene_scores_dmso_original, gene_scores_dmso_modified, gene_scores_condition_original, gene_scores_condition_modified], 
            jitter=True, color=".3", linewidth=1, ax=ax)

        # Set tick locations explicitly
        ax.set_xticks(range(4))
        ax.set_xticklabels(['DMSO Orig', 'DMSO Mod', f'{condition} Orig', f'{condition} Mod'])

        if i % n_columns != 0:  # Only remove y-labels for non-first column plots
            ax.set_ylabel('')

    # Hide any unused axes
    for j in range(i + 1, num_rows * n_columns):
        if j < len(axes):
            fig.delaxes(axes[j])

    plt.tight_layout()
    plt.show()

def update_plots_boxplot_EDCs(adata, condition, gene_set, measure, method):

    conditions = ['DMSO', condition]
    # print(conditions)
    plt.close('all')
    mask = adata.obs['Condition'].isin(conditions)
    adata_filtered = adata[mask]
    if adata_filtered.shape[0] > 0:
        fig, ax = plt.subplots(figsize=(6, 5))

        gene_scores_dmso_ag, gene_scores_dmso_inh = 0, 0
        gene_scores_ag, gene_scores_inh = {}, {}

        for condition in adata_filtered.obs.condition_concentraion.unique():
            if condition == "DMSO_0.1":
                gene_scores_dmso_ag = adata_filtered[adata_filtered.obs['condition_concentraion'] == condition].obs[f'{method}_gene_score_{measure}_{gene_set}_ag'].values
                gene_scores_dmso_inh = adata_filtered[adata_filtered.obs['condition_concentraion'] == condition].obs[f'{method}_gene_score_{measure}_{gene_set}_inh'].values
            else:
                gene_scores_ag[condition] = adata_filtered[adata_filtered.obs['condition_concentraion'] == condition].obs[f'{method}_gene_score_{measure}_{gene_set}_ag'].values
                gene_scores_inh[condition] = adata_filtered[adata_filtered.obs['condition_concentraion'] == condition].obs[f'{method}_gene_score_{measure}_{gene_set}_inh'].values

        conditions_sorted = ["DMSO_0.1"] + list(gene_scores_ag.keys()) # sorted(gene_scores_ag.keys())
        data_ag = [gene_scores_dmso_ag] + [gene_scores_ag[cond] for cond in gene_scores_ag.keys()]
        data_inh = [gene_scores_dmso_inh] + [gene_scores_inh[cond] for cond in gene_scores_inh.keys()]

        positions = np.arange(len(conditions_sorted))
        width = 0.35

        bp_ag = ax.boxplot(data_ag, positions=positions - width/2, widths=width, patch_artist=True, boxprops=dict(facecolor='C0', color='C0'),
                   whiskerprops=dict(color='C0'), capprops=dict(color='C0'), flierprops=dict(color='C0', markeredgecolor='C0'), medianprops=dict(color='white'))
        
        bp_inh = ax.boxplot(data_inh, positions=positions + width/2, widths=width, patch_artist=True, boxprops=dict(facecolor='C1', color='C1'),
                   whiskerprops=dict(color='C1'), capprops=dict(color='C1'), flierprops=dict(color='C1', markeredgecolor='C1'), medianprops=dict(color='white'))

        ax.set_xticks(positions)
        ax.set_xticklabels(conditions_sorted, rotation=45, ha='right')
        ax.set_title(f'{method} - Gene Scores', fontsize=16)
        ax.set_ylabel("Gene Score", fontsize=12)

        ax.legend([bp_ag["boxes"][0], bp_inh["boxes"][0]], ['AG', 'INH'], loc='upper right')

        plt.tight_layout()
        plt.show()
    else:
        print(f'No data to plot for the selected conditions: {conditions}')

def update_plots_boxplot_EDCs_with_scatter(adata, condition, gene_set, measure, method):

    conditions = ['DMSO', condition]
    # print(conditions)
    plt.close('all')
    mask = adata.obs['Condition'].isin(conditions)
    adata_filtered = adata[mask]
    if adata_filtered.shape[0] > 0:
        fig, ax = plt.subplots(figsize=(6, 5))

        gene_scores_dmso_ag, gene_scores_dmso_inh = 0, 0
        gene_scores_ag, gene_scores_inh = {}, {}

        for condition in adata_filtered.obs.condition_concentraion.unique():
            if condition == "DMSO_0.1":
                gene_scores_dmso_ag = adata_filtered[adata_filtered.obs['condition_concentraion'] == condition].obs[f'{method}_gene_score_{measure}_{gene_set}_ag'].values
                gene_scores_dmso_inh = adata_filtered[adata_filtered.obs['condition_concentraion'] == condition].obs[f'{method}_gene_score_{measure}_{gene_set}_inh'].values
            else:
                gene_scores_ag[condition] = adata_filtered[adata_filtered.obs['condition_concentraion'] == condition].obs[f'{method}_gene_score_{measure}_{gene_set}_ag'].values
                gene_scores_inh[condition] = adata_filtered[adata_filtered.obs['condition_concentraion'] == condition].obs[f'{method}_gene_score_{measure}_{gene_set}_inh'].values

        conditions_sorted = ["DMSO_0.1"] + list(gene_scores_ag.keys()) # sorted(gene_scores_ag.keys())
        data_ag = [gene_scores_dmso_ag] + [gene_scores_ag[cond] for cond in gene_scores_ag.keys()]
        data_inh = [gene_scores_dmso_inh] + [gene_scores_inh[cond] for cond in gene_scores_inh.keys()]

        positions = np.arange(len(conditions_sorted))
        width = 0.35

        bp_ag = ax.boxplot(data_ag, positions=positions - width/2, widths=width, patch_artist=True, 
                           boxprops=dict(facecolor='C0', color='None', alpha=0.5),
                           whiskerprops=dict(color='C0'), 
                           capprops=dict(color='C0'), 
                           showfliers=False,
                           medianprops=dict(color='white'))
        
        bp_inh = ax.boxplot(data_inh, positions=positions + width/2, widths=width, patch_artist=True, 
                            boxprops=dict(facecolor='C1', color='None', alpha=0.5),
                            whiskerprops=dict(color='C1'), 
                            capprops=dict(color='C1'), 
                            showfliers=False,
                            medianprops=dict(color='white'))

        for i, (ag_scores, inh_scores) in enumerate(zip(data_ag, data_inh)):
            ax.scatter([positions[i] - width/2] * len(ag_scores), ag_scores, color='C0', alpha=1.0, s=20)
            ax.scatter([positions[i] + width/2] * len(inh_scores), inh_scores, color='C1', alpha=1.0, s=20)
            
        ax.set_xticks(positions)
        ax.set_xticklabels(conditions_sorted, rotation=45, ha='right')
        ax.set_title(f'{method} - Gene Scores', fontsize=16)
        ax.set_ylabel("Gene Score", fontsize=12)

        # ax.legend([bp_ag["boxes"][0], bp_inh["boxes"][0]], ['AG', 'INH'], loc='upper right')

        plt.tight_layout()
        plt.show()
    else:
        print(f'No data to plot for the selected conditions: {conditions}')

def update_plots_boxplot_ag_and_inh_sets_Velmeshev(adata, conditions, type, gene_set, measure, methods):
    plt.close('all')
    if adata.shape[0] > 0:
        num_methods = len(methods)
        fig, axes = plt.subplots(1, num_methods, figsize=(6 * num_methods, 5))

        for i, method in enumerate(methods):
            gene_scores = []
            for condition in conditions:
                if condition == 'Control':
                    gene_scores_control = adata[adata.obs['Og_diagnosis'] == condition].obs[f'{method}_gene_score_{measure}_{gene_set}_{type}'].values
                else:
                    gene_scores_asd = adata[adata.obs['Og_diagnosis'] == condition].obs[f'{method}_gene_score_{measure}_{gene_set}_{type}'].values
            
            ax = axes[i]
            # sns.boxplot(data=[gene_scores_control, gene_scores_asd], notch=True, ax=ax)
            sns.boxplot(data=[gene_scores_control, gene_scores_asd], notch=False,
                            boxprops=dict(alpha=0.5),
                            ax=ax)
            sns.stripplot(data=[gene_scores_control, gene_scores_asd], 
              jitter=True, color=".3", linewidth=1, ax=ax)
            
            ax.set_xticks(range(2))
            ax.set_xticklabels([conditions[0], conditions[1]], fontsize=12)
            ax.set_title(f'{method} - Gene Scores', fontsize=16)
            ax.set_ylabel("Gene Score", fontsize=12)
        
        plt.tight_layout()
        plt.show()
    else:
        print(f'No data to plot for the selected condition: {condition}')

def update_plots_boxplot_ag_and_inh_sets_Velmeshev2(adata, type, measure, method):
    plt.close('all')
    if adata.shape[0] > 0:
        clusters = adata.obs['Og_cluster'].unique()
        num_rows = 9
        num_cols = len(clusters) // num_rows
        if len(clusters) % num_rows != 0:
            num_cols += 1
        
        fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(20, 40))
        axes = axes.flatten()
        
        plt.rcParams.update({'font.size': 14})
        
        for i, cluster in enumerate(clusters):
            ax = axes[i]
            
            cluster_data = adata.obs[adata.obs['Og_cluster'] == cluster]

            pattern = f'{method}_gene_score_{measure}_*_{type}'
            gene_set_score_names = [col for col in cluster_data.columns if fnmatch.fnmatch(col, pattern)]
            
            # Melt the data using the filtered columns
            melted_data = pd.melt(cluster_data, id_vars=['Og_diagnosis'], value_vars=gene_set_score_names, var_name='Score', value_name='Value')
            
            melted_data['Gene_set'] = melted_data['Score'].apply(lambda x: x.split('_')[-2])
    
            sns.boxplot(x='Gene_set', y='Value', hue='Og_diagnosis', data=melted_data, ax=ax)
            
            ax.set_title(f'Cluster:{cluster}  Type:{type}', fontsize=16)
            ax.set_xlabel('Gene_set', fontsize=14)
            ax.set_ylabel('Value', fontsize=14)
            ax.tick_params(axis='x', rotation=45, labelsize=12)
            ax.tick_params(axis='y', labelsize=12)
            
            ax.legend(title='Diagnosis', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)
        
        for j in range(i+1, len(axes)):
            fig.delaxes(axes[j])
        
        plt.tight_layout()
        
        plt.show()