#!/bin/bash
#SBATCH --job-name=18_hyb_kquant
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --account=PWSU0516
#SBATCH --array=1-356%25
#SBATCH --output=logs/18_hybrid_kallisto_quant_%A_%a.log
# ═══════════════════════════════════════════════════════════════════════════════
# 18_hybrid_kallisto_quant.sh — Pseudo-alignment with Kallisto against hybrid
#                                (Macaque + Ebola) transcriptome
# ═══════════════════════════════════════════════════════════════════════════════

SCRIPT_DIR="/users/PWSU0516/mufakiransari/ebola_class_project/scripts"
source "/users/PWSU0516/mufakiransari/ebola_class_project/pipeline.config"
source "/users/PWSU0516/mufakiransari/ebola_class_project/scripts/utils.sh"

STEP_NAME="18_hybrid_kallisto_quant"
SRR=$(get_accession ${SLURM_ARRAY_TASK_ID})
CURRENT_SRR="$SRR"

print_job_info
timer_start

# Checkpoint check
if checkpoint_exists "${STEP_NAME}" "${SRR}"; then
    log_info "Already completed ${SRR}. Skipping."
    exit 0
fi

# Load Conda
module load "${MOD_CONDA}"
ensure_dirs "${HYBRID_KALLISTO_DIR}/${SRR}"

# Define input files
FQ1="${TRIMMED_DIR}/${SRR}/${SRR}_1_paired.fastq.gz"
FQ2="${TRIMMED_DIR}/${SRR}/${SRR}_2_paired.fastq.gz"

validate_files "${FQ1}" "${FQ2}" "${KALLISTO_HYBRID_INDEX}"

# Run Kallisto Quantification against hybrid transcriptome
log_info "Running kallisto quant (hybrid) for ${SRR}..."
conda run -n "${CONDA_ENV}" kallisto quant \
    -i "${KALLISTO_HYBRID_INDEX}" \
    -o "${HYBRID_KALLISTO_DIR}/${SRR}" \
    -b 100 \
    -t "${THREADS}" \
    "${FQ1}" "${FQ2}" \
    2>&1 | tee "${HYBRID_KALLISTO_DIR}/${SRR}/kallisto.log"

validate_files "${HYBRID_KALLISTO_DIR}/${SRR}/abundance.tsv"

checkpoint_set "${STEP_NAME}" "${SRR}"
elapsed_time
