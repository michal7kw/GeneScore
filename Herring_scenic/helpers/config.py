import os

refference = 0
genomes = ["hg19", "hg38"]
genome = genomes[refference]

root_dir = "/group/testa/michal.kubacki/herring"
# root_dir = "/scratch/michal.kubacki/herring"

def set_output_folders(genome, suffix):
    tmp_dir = os.path.join(root_dir, "tmp")
    in_dir = os.path.join(root_dir, "data")
    if(genome == "hg38"):
        out_dir = os.path.join(root_dir, f"output_hg38_{suffix}")
    else:
        out_dir = os.path.join(root_dir, f"output_hg19_{suffix}")

    print(f"root_dir: {root_dir}")
    print(f"out_dir: {out_dir}")
    print(f"in_dir: {in_dir}")
    print(f"tmp_dir: {tmp_dir}")

    if (genome == "hg38"):
        data_folder = "GSE168408_RAW_hg38_bgzip"
    else:
        data_folder = "GSE168408_RAW"

    os.makedirs(out_dir, exist_ok = True)
    return out_dir, in_dir, root_dir, tmp_dir, data_folder

def select_files(genome, selected_fragments = ['2y', '4y']):
    input_dir = os.path.join(root_dir, "data")

    if (genome == "hg38"):
        data_folder = "GSE168408_RAW_hg38_bgzip"
    else:
        data_folder = "GSE168408_RAW"
    
    filenames = """
    GSM5138510_RL2366_ga22_snATAC_fragments.tsv GSM5138526_RL2209_1y_snATAC_fragments.tsv GSM5138542_RL2372_14y_snATAC_fragments.tsv GSM5138512_RL2207_ga24_snATAC_fragments.tsv GSM5138529_RL1784_2y_snATAC_fragments.tsv GSM5138544_RL1785_16y_snATAC_fragments.tsv GSM5138515_RL2367_1m_snATAC_fragments.tsv GSM5138532_RL2210_4y_snATAC_fragments.tsv GSM5138548_RL2085_20y_snATAC_fragments.tsv GSM5138518_RL1914_3m_snATAC_fragments.tsv GSM5138534_RL2364_6y_snATAC_fragments.tsv GSM5138550_RL2369_25y_snATAC_fragments.tsv GSM5138521_RL2208_6m_snATAC_fragments.tsv GSM5138536_RL1994_8y_snATAC_fragments.tsv GSM5138552_RL2373_40y_snATAC_fragments.tsv GSM5138523_RL2371_10m_snATAC_fragments.tsv GSM5138539_RL2368_10y_snATAC_fragments.tsv
    """.strip().split()

    fragments_dict = {}
    for filename in filenames:
        sample_id = filename.split('_')[2]
        if (genome == "hg38"):
            fragments_dict[sample_id] = os.path.join(input_dir, data_folder, f"{filename}.bgz")
        else:
            fragments_dict[sample_id] = os.path.join(input_dir, data_folder, f"{filename}.gz")  

    print(f"All fragments: {fragments_dict}")
    print(f"Selected fragments: {selected_fragments}") 
    fragments_dict = {key: fragments_dict[key] for key in selected_fragments if key in fragments_dict}
    return fragments_dict


# for file in *fragments.tsv; do tabix -p bed "$file"; done

def select_files_rna(genome, selected_fragments = ['2y', '4y']):
    input_dir = os.path.join(root_dir, "data")

    if (genome == "hg38"):
        data_folder = "GSE168408_RAW_hg38_bgzip"
    else:
        data_folder = "GSE168408_RAW"
    
    filenames = """
    GSM5138509_RL2103_ga22_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138511_RL2107_ga24_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138513_RL2121_ga34_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138514_RL1777_1m_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138516_RL1612_2m_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138517_RL2100_3m_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138519_RL2104_4m_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138520_RL2108_6m_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138522_RL2122_10m_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138524_RL2125_1y_a_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138525_RL2105_1y_b_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138528_RL1613_2y_a_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138530_RL2129_3y_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138533_RL2106_6y_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138535_RL1614_8y_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138537_RL2110_10y_a_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138538_RL2126_10y_b_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138540_RL2127_12y_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138541_RL2130_14y_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138543_RL2102_16y_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138545_RL2131_17y_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138546_RL2123_20y_a_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138547_RL2128_20y_b_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138549_RL2132_25y_snRNAseq_filtered_feature_bc_matrix.h5 GSM5138551_RL2124_40y_snRNAseq_filtered_feature_bc_matrix.h5 GSM5640800_RL2290_12mth_organoid_snRNAseq_filtered_feature_bc_matrix.h5 GSM5640801_RL2340_9mth_organoid_snRNAseq_filtered_feature_bc_matrix.h5 GSM5640802_RL2432_5mth_organoid_snRNAseq_filtered_feature_bc_matrix.h5 GSM6499429_RL2986_1y_c_snRNAseq_filtered_feature_bc_matrix.h5 GSM6499430_RL2987_1y_d_snRNAseq_filtered_feature_bc_matrix.h5
    """.strip().split()

    fragments_dict = {}
    for filename in filenames:
        sample_id = filename.split('_')[2]
        if (genome == "hg38"):
            fragments_dict[sample_id] = os.path.join(input_dir, data_folder, f"{filename}.bgz")
        else:
            fragments_dict[sample_id] = os.path.join(input_dir, data_folder, f"{filename}.gz")  

    print(f"All fragments: {fragments_dict}")
    print(f"Selected fragments: {selected_fragments}") 
    fragments_dict = {key: fragments_dict[key] for key in selected_fragments if key in fragments_dict}
    return fragments_dict