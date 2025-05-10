
# %%
import os
import re
import sys

# Set working directory
# work_dir = '/home/michal.kubacki/Githubs/GeneScore/trimmed_GRN_derivation'
# work_dir = 'D:/Github/GeneScore/trimmed_GRN_derivation'
work_dir = '/mnt/d/Github/GeneScore/trimmed_GRN_derivation'

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

from grn_helpers import set_output_folders

# %%
# neurons_set = "L2-3_CUX2"
neurons_set = "all_ex"
# neurons_set = "all_ex_all_ages"
root_dir = os.getenv('BASE_PATH')

# %%
out_dir, in_dir, root_dir, tmp_dir, data_folder = set_output_folders(root_dir, neurons_set)

dir_input_files = {
    "all_ex" : ['L2-3_CUX2_consensus_regions.bed', 'L4_RORB_consensus_regions.bed', 'L5-6_THEMIS_consensus_regions.bed', 'L5-6_TLE4_consensus_regions.bed', 'PN_dev_consensus_regions.bed'],
    "all_ex_all_ages" : ['L2-3_CUX2_consensus_regions.bed', 'L4_RORB_consensus_regions.bed', 'L5-6_THEMIS_consensus_regions.bed', 'L5-6_TLE4_consensus_regions.bed', 'PN_dev_consensus_regions.bed'],
    "L2-3_CUX2" : ['L2-3_CUX2_consensus_regions.bed']
}

input_files = dir_input_files[neurons_set]

# %%
# Function to process a single file
def process_file(file_name, out_dir):
    file_path = os.path.join(out_dir, "cell_type_consensus_regions", file_name)
    peaks = []
    
    with open(file_path, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            if len(fields) >= 8:
                # chr_start_end = fields[7].replace(':', '_').replace('-', '_')
                chr_start_end = fields[6]
                peaks.append(f'{chr_start_end}')
    
    return peaks


def extract_cell_type(file_name):
    pattern = r'^(.+)_consensus_regions\.bed$'
    match = re.match(pattern, file_name)
    if match:
        return match.group(1)
    else:
        return None

# Process each file
for file_name in input_files:
    # Extract the cell type from the file name
    cell_type = extract_cell_type(file_name)

    
    # Construct the output file name
    output_file = f'{cell_type}_peaks.csv'
    output_path = os.path.join(out_dir, output_file)
    
    # Process the file
    peaks = process_file(file_name, out_dir)
    
    # Save the peaks to a new file
    with open(output_path, 'w') as file:
        file.write(' '.join(peaks))

print("Processing complete. Output files saved with '_peaks.csv' suffix.")
