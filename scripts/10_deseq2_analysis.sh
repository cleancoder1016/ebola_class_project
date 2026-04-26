#!/bin/bash
#SBATCH --job-name=10_deseq2
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --account=PWSU0516
#SBATCH --output=logs/10_deseq2_%j.log
# ═══════════════════════════════════════════════════════════════════════════════
# 10_deseq2_analysis.sh — Run DESeq2 exploratory / DE analysis in R
# ═══════════════════════════════════════════════════════════════════════════════

SCRIPT_DIR="/users/PWSU0516/mufakiransari/ebola_class_project/scripts"
source "/users/PWSU0516/mufakiransari/ebola_class_project/pipeline.config"
source "/users/PWSU0516/mufakiransari/ebola_class_project/scripts/utils.sh"

STEP_NAME="10_deseq2"
print_job_info
timer_start

# ── Checkpoint check ────────────────────────────────────────────────────────
if checkpoint_exists "${STEP_NAME}"; then
    log_info "Already completed. Skipping."
    exit 0
fi

# ── Setup ────────────────────────────────────────────────────────────────────
module load "${MOD_GCC}"
module load "${MOD_R}"
ensure_dirs "${DESEQ2_DIR}"

# ── Validate input ─────────────────────────────────────────────────────────
COUNTS_FILE="${COUNTS_DIR}/gene_counts.txt"
validate_files "${COUNTS_FILE}"

# ── Check for sample metadata ──────────────────────────────────────────────
METADATA_FILE="${PROJECT_DIR}/sample_metadata.csv"
if [[ -f "${METADATA_FILE}" ]]; then
    log_info "Found sample metadata: ${METADATA_FILE}"
    METADATA_ARG="${METADATA_FILE}"
else
    log_info "No sample_metadata.csv found — running exploratory analysis only."
    log_info "To enable DE analysis, create ${METADATA_FILE} with columns: sample_id, condition"
    METADATA_ARG=""
fi

# ── Run R analysis ─────────────────────────────────────────────────────────
log_info "Running DESeq2 analysis..."
Rscript "${SCRIPT_DIR}/deseq2_analysis.R" \
    "${COUNTS_FILE}" \
    "${DESEQ2_DIR}" \
    ${METADATA_ARG:+"${METADATA_ARG}"} \
    2>&1 | tee "${DESEQ2_DIR}/deseq2.log"

# ── Validate outputs ───────────────────────────────────────────────────────
if [[ -f "${DESEQ2_DIR}/pca_plot.pdf" ]]; then
    log_info "DESeq2 analysis completed successfully."
    ls -la "${DESEQ2_DIR}/"
else
    log_warn "DESeq2 analysis may have issues — check ${DESEQ2_DIR}/deseq2.log"
fi

checkpoint_set "${STEP_NAME}"
elapsed_time
