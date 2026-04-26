#!/bin/bash
#SBATCH --job-name=14_kagg
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=PWSU0516
#SBATCH --output=logs/14_kallisto_aggregate_%j.log
# ═══════════════════════════════════════════════════════════════════════════════
# 14_kallisto_aggregate.sh — Aggregate Kallisto counts and plot vs featureCounts
# ═══════════════════════════════════════════════════════════════════════════════

SCRIPT_DIR="/users/PWSU0516/mufakiransari/ebola_class_project/scripts"
source "/users/PWSU0516/mufakiransari/ebola_class_project/pipeline.config"
source "/users/PWSU0516/mufakiransari/ebola_class_project/scripts/utils.sh"

STEP_NAME="14_kallisto_aggregate"
print_job_info
timer_start

# Checkpoint check
if checkpoint_exists "${STEP_NAME}"; then
    log_info "Already completed. Skipping."
    exit 0
fi

# Load Conda
module load "${MOD_CONDA}"
ensure_dirs "${KALLISTO_DIR}"

log_info "Running kallisto_aggregate.py..."
conda run -n "${CONDA_ENV}" python "${SCRIPT_DIR}/kallisto_aggregate.py"

validate_files "${PROJECT_DIR}/kallisto_output/count_matrix_kallisto.csv"

checkpoint_set "${STEP_NAME}"
elapsed_time
