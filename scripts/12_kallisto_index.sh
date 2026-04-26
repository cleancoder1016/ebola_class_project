#!/bin/bash
#SBATCH --job-name=12_kindex
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=PWSU0516
#SBATCH --output=logs/12_kallisto_index_%j.log
# ═══════════════════════════════════════════════════════════════════════════════
# 12_kallisto_index.sh — Generate Transcriptome from FASTA/GTF and Index
# ═══════════════════════════════════════════════════════════════════════════════

SCRIPT_DIR="/users/PWSU0516/mufakiransari/ebola_class_project/scripts"
source "/users/PWSU0516/mufakiransari/ebola_class_project/pipeline.config"
source "/users/PWSU0516/mufakiransari/ebola_class_project/scripts/utils.sh"

export KALLISTO_INDEX="${REFERENCE_DIR}/ebola_transcriptome.idx"
export TRANSCRIPTOME_FASTA="${REFERENCE_DIR}/ebola_transcriptome.fasta"

STEP_NAME="12_kallisto_index"
print_job_info
timer_start

# Checkpoint check
if checkpoint_exists "${STEP_NAME}"; then
    log_info "Already completed. Skipping."
    exit 0
fi

# Load Conda
module load "${MOD_CONDA}"
ensure_dirs "${REFERENCE_DIR}"

log_info "Extracting transcriptome FASTA using gffread..."
conda run -n "${CONDA_ENV}" gffread -w "${TRANSCRIPTOME_FASTA}" -g "${REFERENCE_FASTA}" "${REFERENCE_GTF}"
validate_files "${TRANSCRIPTOME_FASTA}"

log_info "Building Kallisto index..."
conda run -n "${CONDA_ENV}" kallisto index -i "${KALLISTO_INDEX}" "${TRANSCRIPTOME_FASTA}"
validate_files "${KALLISTO_INDEX}"

checkpoint_set "${STEP_NAME}"
elapsed_time
