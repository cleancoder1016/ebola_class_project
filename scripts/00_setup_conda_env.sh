#!/bin/bash
#SBATCH --job-name=00_conda_setup
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=PWSU0516
#SBATCH --output=logs/00_conda_%j.log
# ═══════════════════════════════════════════════════════════════════════════════
# 00_setup_conda_env.sh — One-time conda environment creation
# ═══════════════════════════════════════════════════════════════════════════════
# Creates a local conda environment with tools not available as HPC modules:
#   - subread   (provides featureCounts)
#   - multiqc   (aggregated QC reporting)
#   - gffread   (GFF ↔ GTF conversion)
# ═══════════════════════════════════════════════════════════════════════════════

SCRIPT_DIR="/users/PWSU0516/mufakiransari/ebola_class_project/scripts"
source "/users/PWSU0516/mufakiransari/ebola_class_project/pipeline.config"
source "/users/PWSU0516/mufakiransari/ebola_class_project/scripts/utils.sh"

STEP_NAME="00_conda_setup"
print_job_info
timer_start

# ── Check if environment already exists ──────────────────────────────────────
if checkpoint_exists "${STEP_NAME}"; then
    log_info "Conda environment '${CONDA_ENV}' already set up. Skipping."
    exit 0
fi

ensure_dirs "${LOG_DIR}"

# ── Load conda module ───────────────────────────────────────────────────────
log_info "Loading conda module: ${MOD_CONDA}"
module load "${MOD_CONDA}"

# ── Create the environment ──────────────────────────────────────────────────
if conda env list | grep -q "^${CONDA_ENV} "; then
    log_info "Conda environment '${CONDA_ENV}' already exists. Updating..."
    conda install -n "${CONDA_ENV}" -y -c bioconda -c conda-forge \
        subread multiqc gffread
else
    log_info "Creating conda environment '${CONDA_ENV}'..."
    conda create -n "${CONDA_ENV}" -y -c bioconda -c conda-forge \
        python=3.10 subread multiqc gffread
fi

# ── Verify installation ────────────────────────────────────────────────────
log_info "Verifying installed tools..."
conda run -n "${CONDA_ENV}" featureCounts -v 2>&1 | head -1
conda run -n "${CONDA_ENV}" multiqc --version 2>&1 | head -1
conda run -n "${CONDA_ENV}" gffread --version 2>&1 | head -1

log_info "Conda environment '${CONDA_ENV}' is ready."
checkpoint_set "${STEP_NAME}"
elapsed_time
