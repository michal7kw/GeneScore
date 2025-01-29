#!/bin/bash
#SBATCH --job-name=cistarget
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michal.kubacki@external.fht.org
#SBATCH --partition=cpuq
#SBATCH --time=32:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --mem=65G
#SBATCH --output=./logs/cistarget_%j.log

# Source the conda profile script
source /home/michal.kubacki/new_miniconda/etc/profile.d/conda.sh
conda activate scenicplus

module load bedtools2/2.31.0

######## Select correct configuration ##################
ROOT_DIR="/group/testa/michal.kubacki/herring"
# ROOT_DIR="/scratch/michal.kubacki/herring"


DATABASE_PREFIX="all_cells"
# DATABASE_PREFIX="all_excitatory"
# DATABASE_PREFIX="all_inhibitory"
# DATABASE_PREFIX="all_excitatory_all_ages"
# DATABASE_PREFIX="all_inhibitory_all_ages"
######################################################## 

REGION_BED="${ROOT_DIR}/output_hg19_${DATABASE_PREFIX}/consensus_peak_calling/consensus_regions.bed"
GENOME_FASTA="${ROOT_DIR}/data/hg19.fa"

CHROMSIZES="${ROOT_DIR}/data/hg19.chrom.sizes"
SCRIPT_DIR="${ROOT_DIR}/cistarget_database/create_cisTarget_databases"

OUT_DIR="${ROOT_DIR}/output_hg19_${DATABASE_PREFIX}/cistarget_motif_dataset"
CBDIR="${ROOT_DIR}/cistarget_database/aertslab_motif_colleciton/v10nr_clust_public/singletons"

FASTA_FILE="${OUT_DIR}/hg19_${DATABASE_PREFIX}.with_1kb_bg_padding.fa"
MOTIF_LIST="${OUT_DIR}/motifs.txt"

mkdir -p "${OUT_DIR}"

"${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh" \
        "${GENOME_FASTA}" \
        "${CHROMSIZES}" \
        "${REGION_BED}" \
        "${FASTA_FILE}" \
        1000 \
        yes

head -n 2 "${FASTA_FILE}"

touch "${MOTIF_LIST}"
ls "${CBDIR}" > "${MOTIF_LIST}"

"${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f "${FASTA_FILE}" \
    -M "${CBDIR}" \
    -m "${MOTIF_LIST}" \
    -o "${OUT_DIR}" \
    --bgpadding 1000 \
    -t 128
