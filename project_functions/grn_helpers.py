
import os

def set_custom_folders(root_dir, suffix):
    tmp_dir = os.path.join(root_dir, "celloracle", "tmp")

    in_dir = os.path.join(root_dir, "data")
    in_dir_from_scenic = os.path.join(root_dir, suffix)

    out_dir = os.path.join(root_dir, suffix, "celloracle")

    print(f"root_dir: {root_dir}")
    print(f"out_dir: {out_dir}")
    print(f"in_dir: {in_dir}")
    print(f"tmp_dir: {tmp_dir}")

    os.makedirs(out_dir, exist_ok = True)
    return out_dir, in_dir, root_dir, tmp_dir, in_dir_from_scenic

def set_output_folders(root_dir,suffix):
    tmp_dir = os.path.join(root_dir, "tmp")
    in_dir = os.path.join(root_dir, "data")
    out_dir = os.path.join(root_dir, suffix)

    print(f"root_dir: {root_dir}")
    print(f"out_dir: {out_dir}")
    print(f"in_dir: {in_dir}")
    print(f"tmp_dir: {tmp_dir}")

    data_folder = "GSE168408_RAW"

    os.makedirs(out_dir, exist_ok = True)
    return out_dir, in_dir, root_dir, tmp_dir, data_folder

def select_files(root_dir, selected_fragments = ['2y', '4y']):
    input_dir = os.path.join(root_dir, "data")
    data_folder = "GSE168408_RAW"
    filenames = """
    GSM5138510_RL2366_ga22_snATAC_fragments.tsv GSM5138526_RL2209_1y_snATAC_fragments.tsv GSM5138542_RL2372_14y_snATAC_fragments.tsv GSM5138512_RL2207_ga24_snATAC_fragments.tsv GSM5138529_RL1784_2y_snATAC_fragments.tsv GSM5138544_RL1785_16y_snATAC_fragments.tsv GSM5138515_RL2367_1m_snATAC_fragments.tsv GSM5138532_RL2210_4y_snATAC_fragments.tsv GSM5138548_RL2085_20y_snATAC_fragments.tsv GSM5138518_RL1914_3m_snATAC_fragments.tsv GSM5138534_RL2364_6y_snATAC_fragments.tsv GSM5138550_RL2369_25y_snATAC_fragments.tsv GSM5138521_RL2208_6m_snATAC_fragments.tsv GSM5138536_RL1994_8y_snATAC_fragments.tsv GSM5138552_RL2373_40y_snATAC_fragments.tsv GSM5138523_RL2371_10m_snATAC_fragments.tsv GSM5138539_RL2368_10y_snATAC_fragments.tsv
    """.strip().split()

    fragments_dict = {}
    for filename in filenames:
        sample_id = filename.split('_')[2]
        fragments_dict[sample_id] = os.path.join(input_dir, data_folder, f"{filename}.gz")  

    print(f"All fragments: {fragments_dict}")
    print(f"Selected fragments: {selected_fragments}") 
    fragments_dict = {key: fragments_dict[key] for key in selected_fragments if key in fragments_dict}
    return fragments_dict