#!/bin/bash
#SBATCH --job-name=01_sra2fq
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=48G
#SBATCH --account=PWSU0516
#SBATCH --array=1-356%10
#SBATCH --output=logs/01_sra2fq_%A_%a.log
# ═══════════════════════════════════════════════════════════════════════════════
# 01_sra_to_fastq.sh — Convert SRA files to FASTQ (with retry + checkpointing)
# ═══════════════════════════════════════════════════════════════════════════════

SCRIPT_DIR="/users/PWSU0516/mufakiransari/ebola_class_project/scripts"
source "/users/PWSU0516/mufakiransari/ebola_class_project/pipeline.config"
source "/users/PWSU0516/mufakiransari/ebola_class_project/scripts/utils.sh"

STEP_NAME="01_sra_to_fastq"
SRR=$(get_accession)
CURRENT_SRR="$SRR"
print_job_info
timer_start

# ── Checkpoint check ────────────────────────────────────────────────────────
if checkpoint_exists "${STEP_NAME}" "$SRR"; then
    log_info "Already converted. Skipping."
    exit 0
fi

# ── Setup ────────────────────────────────────────────────────────────────────
module load "${MOD_SRA}"
ensure_dirs "${LOG_DIR}" "${FASTQ_DIR}/${SRR}" "${SRA_DIR}"

# ── Prefetch with retry ─────────────────────────────────────────────────────
log_info "Prefetching SRA data..."
if [[ -f "${SRA_DIR}/${SRR}/${SRR}.sra" ]]; then
    log_info "SRA file already exists: ${SRA_DIR}/${SRR}/${SRR}.sra"
else
    retry 3 prefetch "$SRR" --max-size 100G -O "${SRA_DIR}/"
fi

# ── Convert to FASTQ ────────────────────────────────────────────────────────
log_info "Converting SRA → FASTQ with fasterq-dump..."
fasterq-dump "${SRA_DIR}/${SRR}/${SRR}.sra" \
    --split-3 \
    -e "${THREADS}" \
    -O "${FASTQ_DIR}/${SRR}" \
    --temp "${PROJECT_DIR}" \
    --force

# ── Validate outputs ────────────────────────────────────────────────────────
validate_files "${FASTQ_DIR}/${SRR}/${SRR}_1.fastq" "${FASTQ_DIR}/${SRR}/${SRR}_2.fastq"

# ── Report read counts ──────────────────────────────────────────────────────
R1_COUNT=$(fastq_read_count "${FASTQ_DIR}/${SRR}/${SRR}_1.fastq")
R2_COUNT=$(fastq_read_count "${FASTQ_DIR}/${SRR}/${SRR}_2.fastq")
log_info "R1 reads: ${R1_COUNT} | R2 reads: ${R2_COUNT}"

# ── Compress FASTQs ─────────────────────────────────────────────────────────
log_info "Compressing FASTQ files with gzip..."
gzip -f "${FASTQ_DIR}/${SRR}/${SRR}_1.fastq"
gzip -f "${FASTQ_DIR}/${SRR}/${SRR}_2.fastq"
validate_files "${FASTQ_DIR}/${SRR}/${SRR}_1.fastq.gz" "${FASTQ_DIR}/${SRR}/${SRR}_2.fastq.gz"

log_info "Conversion complete."
checkpoint_set "${STEP_NAME}" "$SRR"
elapsed_time
