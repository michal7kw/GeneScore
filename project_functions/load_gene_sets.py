import csv
import os
import pandas as pd
import anndata as ad
import numpy as np


def remove_duplicates_preserve_order_GRNs(data_dict):
    result = {}
    for set_selected, set_data in data_dict.items():
        result[set_selected] = {}
        for cell_type_selected, cell_type_data in set_data.items():
            result[set_selected][cell_type_selected] = {}
            for scored_gene_selected, gene_data in cell_type_data.items():
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

def load_GRNs_gene_sets(root_dir):
    gene_sets = {}

    # List of gene sets to process
    gene_set_list = [
        "all_excitatory",
        "all_inhibitory",
        "all_excitatory_all_ages",
        "all_inhibitory_all_ages"
    ]

    # Load data for each gene set
    for gene_set in gene_set_list:
        path = os.path.join(root_dir, f"output_hg19_{gene_set}", "celloracle")
        gene_sets[gene_set] = pd.read_csv(os.path.join(path, 'scores_grn_all.csv'))

    # Sets Formatting
    gene_sets_dict = {}
    gene_sets_dict_cell_type_first = {}

    for key, value in gene_sets.items():
        gene_sets_dict[key] = {}
        gene_sets_dict_cell_type_first[key] = {}

        for _, row in value.iterrows():
            goi = row['source']
            target = row['target']
            score1 = float(row['score']) * float(row['coef_mean'])
            score2 = float(row['coef_mean'])
            source = row['celltype']

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

    print(gene_sets_dict_cell_type_first.keys())

    return gene_sets_dict, gene_sets_dict_cell_type_first

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
    
def load_GRNs_gene_sets_legacy(root_dir="/group/testa/michal.kubacki/herring"):
    # %%
    gene_sets = {}

    # %% [markdown]
    # ## GRNs - basic

    # %%
    gene_set = pd.read_csv(os.path.join(root_dir, 'data', 'endocrineReceptorsGRNs_L.tsv'), sep='\t')
    gene_sets["base_set"] = gene_set

    cell_type_mapping = {
        'L2-3_CUX2': ['Exc_DL', 'Exc_Mig'], 
        'L4_RORB': 'Exc_Mat',
        'L5-6_THEMIS': 'Exc_Mat',
        'L5-6_TLE4': 'Exc_Mat',
        'PN_dev': ['InCGE', 'InMGE', 'IP', 'oRG', 'PgG2M', 'PgS', 'vRG']
    }

    reverse_mapping = {}
    for k, v in cell_type_mapping.items():
        if isinstance(v, list):
            for item in v:
                reverse_mapping[item] = k
        else:
            reverse_mapping[v] = k

    # Map the 'celltype' column to the corresponding set 2 cell types
    gene_sets["base_set"]['celltype'] = gene_sets["base_set"]['celltype'].map(reverse_mapping)


    # %% [markdown]
    # ## GRNs - Bhadur

    # %%
    gene_set = "output_hg19_all_excitatory_base_grn"
    path = os.path.join(root_dir, gene_set, "celloracle")

    gene_sets[gene_set] = pd.read_csv(os.path.join(path, 'combined_table_data_GRN.csv'))

    # %% [markdown]
    # ## GRNs - Herring

    # %%
    gene_set = "output_hg19_all_excitatory"
    path = os.path.join(root_dir, gene_set, "celloracle")

    gene_sets[gene_set] = pd.read_csv(os.path.join(path, 'scores_grn_all.csv'))

    # %% [markdown]
    # # Sets Formating

    # %%
    gere_sets_dict = {}

    for key, value in gene_sets.items():
        gere_sets_dict[key] = {}
        for _, row in value.iterrows():
            goi = row['source']
            target = row['target']
            score = float(row['score'])
            source = row['celltype']

            if goi not in gere_sets_dict[key]:
                gere_sets_dict[key][goi] = {}

            if source not in gere_sets_dict[key][goi]:
                gere_sets_dict[key][goi][source] = {'targets': [], 'weights': []}

            gere_sets_dict[key][goi][source]['targets'].append(target)
            gere_sets_dict[key][goi][source]['weights'].append(score)

    # %%
    gene_sets_dict_cell_type_first = {}

    for key, value in gene_sets.items():
        gene_sets_dict_cell_type_first[key] = {}
        for _, row in value.iterrows():
            goi = row['source']
            target = row['target']
            score = float(row['score'])
            source = row['celltype']

            if source not in gene_sets_dict_cell_type_first[key]:
                gene_sets_dict_cell_type_first[key][source] = {}

            if goi not in gene_sets_dict_cell_type_first[key][source]:
                gene_sets_dict_cell_type_first[key][source][goi] = {'targets': [], 'weights': []}

            gene_sets_dict_cell_type_first[key][source][goi]['targets'].append(target)
            gene_sets_dict_cell_type_first[key][source][goi]['weights'].append(score)

    # %%
    print(gene_sets_dict_cell_type_first.keys())

    return gere_sets_dict, gene_sets_dict_cell_type_first
