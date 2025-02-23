#!/bin/bash
#SBATCH --job-name=GRN_pipeline_master
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michal.kubacki@external.fht.org
#SBATCH --partition=cpuq
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --output=./logs/master_pipeline.log

# Master script to run the entire GRN analysis pipeline
# This script will execute all analysis steps in order using SLURM dependencies

# Load environment variables
set -a
source .env
set +a

# Change to work directory
cd $WORK_DIR

# Create logs directory if it doesn't exist
mkdir -p logs

# Function to submit job and return job ID
submit_job() {
    local script=$1
    local dependency=$2
    
    if [ -z "$dependency" ]; then
        # Submit job without dependency
        job_id=$(sbatch "run_${script}.sh" | awk '{print $4}')
    else
        # Submit job with dependency
        job_id=$(sbatch --dependency=afterok:${dependency} "run_${script}.sh" | awk '{print $4}')
    fi
    
    echo $job_id
}

# Array of scripts in order of execution
scripts=(
    "1_herring_scRNA_processing"
    "2_herring_atac_processing"
    "3_celltype_consensus_regions"
    "4_coaccess"
    "5_prepare_peak_data_for_celloracle"
    "6_integrate_tss_with_coaccess_regions"
    "7_scan_motifs_for_different_cell_types"
    "8_celloracle_scan_motifs"
    "9_create_scores_procedural"
)

# Submit first job
echo "Submitting ${scripts[0]}..."
previous_job_id=$(submit_job "${scripts[0]}")
echo "Submitted ${scripts[0]} (Job ID: $previous_job_id)"

# Submit remaining jobs with dependencies
for ((i=1; i<${#scripts[@]}; i++)); do
    script="${scripts[$i]}"
    echo "Submitting $script (dependent on Job ID: $previous_job_id)..."
    
    job_id=$(submit_job "$script" "$previous_job_id")
    
    if [ -z "$job_id" ]; then
        echo "Failed to submit $script"
        exit 1
    fi
    
    echo "Submitted $script (Job ID: $job_id)"
    previous_job_id=$job_id
done

echo "All jobs submitted successfully!"
echo "Monitor progress with: squeue -u \$USER"
echo "Check individual job logs in the logs directory" 