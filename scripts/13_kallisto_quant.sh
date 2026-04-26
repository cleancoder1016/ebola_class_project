#!/bin/bash
#SBATCH --job-name=13_kquant
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --account=PWSU0516
#SBATCH --array=1-50%10
#SBATCH --output=logs/13_kallisto_quant_%A_%a.log
# ═══════════════════════════════════════════════════════════════════════════════
# 13_kallisto_quant.sh — Pseudo-alignment and quantification with Kallisto
# ═══════════════════════════════════════════════════════════════════════════════

SCRIPT_DIR="/users/PWSU0516/mufakiransari/ebola_class_project/scripts"
source "/users/PWSU0516/mufakiransari/ebola_class_project/pipeline.config"
source "/users/PWSU0516/mufakiransari/ebola_class_project/scripts/utils.sh"

export KALLISTO_DIR="${PROJECT_DIR}/kallisto_output"
export KALLISTO_INDEX="${REFERENCE_DIR}/ebola_transcriptome.idx"
export TRANSCRIPTOME_FASTA="${REFERENCE_DIR}/ebola_transcriptome.fasta"

STEP_NAME="13_kallisto_quant"
SRR=$(get_accession ${SLURM_ARRAY_TASK_ID})

print_job_info
timer_start

# Checkpoint check
if checkpoint_exists "${STEP_NAME}" "${SRR}"; then
    log_info "Already completed ${SRR}. Skipping."
    exit 0
fi

# Load Conda
module load "${MOD_CONDA}"
ensure_dirs "${KALLISTO_DIR}/${SRR}"

# Define input files
FQ1="${TRIMMED_DIR}/${SRR}/${SRR}_1_paired.fastq.gz"
FQ2="${TRIMMED_DIR}/${SRR}/${SRR}_2_paired.fastq.gz"

validate_files "${FQ1}" "${FQ2}" "${KALLISTO_INDEX}"

# Run Kallisto Quantification
log_info "Running kallisto quant for ${SRR}..."
conda run -n "${CONDA_ENV}" kallisto quant \
    -i "${KALLISTO_INDEX}" \
    -o "${KALLISTO_DIR}/${SRR}" \
    -b 100 \
    -t "${THREADS}" \
    "${FQ1}" "${FQ2}" \
    2>&1 | tee "${KALLISTO_DIR}/${SRR}/kallisto.log"

validate_files "${KALLISTO_DIR}/${SRR}/abundance.tsv"

checkpoint_set "${STEP_NAME}" "${SRR}"
elapsed_time
