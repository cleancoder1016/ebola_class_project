#!/bin/bash
#SBATCH --job-name=04_fastqc_trim
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --account=PWSU0516
#SBATCH --array=1-356%10
#SBATCH --output=logs/04_fastqc_trim_%A_%a.log
# ═══════════════════════════════════════════════════════════════════════════════
# 04_fastqc_trimmed.sh — Quality control on trimmed FASTQ reads
# ═══════════════════════════════════════════════════════════════════════════════

SCRIPT_DIR="/users/PWSU0516/mufakiransari/ebola_class_project/scripts"
source "/users/PWSU0516/mufakiransari/ebola_class_project/pipeline.config"
source "/users/PWSU0516/mufakiransari/ebola_class_project/scripts/utils.sh"

STEP_NAME="04_fastqc_trimmed"
SRR=$(get_accession)
CURRENT_SRR="$SRR"
print_job_info
timer_start

# ── Checkpoint check ────────────────────────────────────────────────────────
if checkpoint_exists "${STEP_NAME}" "$SRR"; then
    log_info "Already completed. Skipping."
    exit 0
fi

# ── Setup ────────────────────────────────────────────────────────────────────
module load "${MOD_FASTQC}"
ensure_dirs "${QC_TRIMMED_DIR}"

# ── Locate trimmed paired reads ─────────────────────────────────────────────
R1P="${TRIMMED_DIR}/${SRR}/${SRR}_1_paired.fastq.gz"
R2P="${TRIMMED_DIR}/${SRR}/${SRR}_2_paired.fastq.gz"
validate_files "$R1P" "$R2P"

# ── Run FastQC on trimmed reads ─────────────────────────────────────────────
log_info "Running FastQC on trimmed reads..."
fastqc "$R1P" "$R2P" \
    --outdir "${QC_TRIMMED_DIR}" \
    --threads "${SLURM_CPUS_PER_TASK:-4}" \
    --quiet

log_info "FastQC trimmed QC complete."
checkpoint_set "${STEP_NAME}" "$SRR"
elapsed_time
