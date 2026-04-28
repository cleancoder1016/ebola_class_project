#!/bin/bash
#SBATCH --job-name=19_hyb_kagg
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --account=PWSU0516
#SBATCH --output=logs/19_hybrid_kallisto_aggregate_%j.log
# ═══════════════════════════════════════════════════════════════════════════════
# 19_hybrid_kallisto_aggregate.sh — Aggregate hybrid Kallisto counts and
#                                    compare against hybrid featureCounts
# ═══════════════════════════════════════════════════════════════════════════════

SCRIPT_DIR="/users/PWSU0516/mufakiransari/ebola_class_project/scripts"
source "/users/PWSU0516/mufakiransari/ebola_class_project/pipeline.config"
source "/users/PWSU0516/mufakiransari/ebola_class_project/scripts/utils.sh"

STEP_NAME="19_hybrid_kallisto_aggregate"
print_job_info
timer_start

# Checkpoint check
if checkpoint_exists "${STEP_NAME}"; then
    log_info "Already completed. Skipping."
    exit 0
fi

# Load Conda
module load "${MOD_CONDA}"
ensure_dirs "${HYBRID_KALLISTO_DIR}"

log_info "Running hybrid_kallisto_aggregate.py..."
conda run -n "${CONDA_ENV}" python "${SCRIPT_DIR}/hybrid_kallisto_aggregate.py"

validate_files "${HYBRID_KALLISTO_DIR}/count_matrix_hybrid_kallisto.csv"

checkpoint_set "${STEP_NAME}"
elapsed_time
