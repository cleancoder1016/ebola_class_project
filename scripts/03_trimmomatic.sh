#!/bin/bash
#SBATCH --job-name=03_trimmomatic
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=16G
#SBATCH --account=PWSU0516
#SBATCH --array=1-356%25
#SBATCH --output=logs/03_trimmomatic_%A_%a.log
# ═══════════════════════════════════════════════════════════════════════════════
# 03_trimmomatic.sh — Adapter removal and quality trimming (paired-end)
# ═══════════════════════════════════════════════════════════════════════════════

SCRIPT_DIR="/users/PWSU0516/mufakiransari/ebola_class_project/scripts"
source "/users/PWSU0516/mufakiransari/ebola_class_project/pipeline.config"
source "/users/PWSU0516/mufakiransari/ebola_class_project/scripts/utils.sh"

STEP_NAME="03_trimmomatic"
SRR=$(get_accession)
CURRENT_SRR="$SRR"
print_job_info
timer_start

# ── Checkpoint check ────────────────────────────────────────────────────────
if checkpoint_exists "${STEP_NAME}" "$SRR"; then
    log_info "Already trimmed. Skipping."
    exit 0
fi

# ── Setup ────────────────────────────────────────────────────────────────────
module load "${MOD_TRIMMOMATIC}"
ensure_dirs "${TRIMMED_DIR}/${SRR}"

# ── Locate input files ──────────────────────────────────────────────────────
R1="${FASTQ_DIR}/${SRR}/${SRR}_1.fastq.gz"
R2="${FASTQ_DIR}/${SRR}/${SRR}_2.fastq.gz"
if [[ ! -f "$R1" ]]; then
    R1="${FASTQ_DIR}/${SRR}/${SRR}_1.fastq"
    R2="${FASTQ_DIR}/${SRR}/${SRR}_2.fastq"
fi
validate_files "$R1" "$R2"

# ── Define output files ─────────────────────────────────────────────────────
OUT_1P="${TRIMMED_DIR}/${SRR}/${SRR}_1_paired.fastq.gz"
OUT_1U="${TRIMMED_DIR}/${SRR}/${SRR}_1_unpaired.fastq.gz"
OUT_2P="${TRIMMED_DIR}/${SRR}/${SRR}_2_paired.fastq.gz"
OUT_2U="${TRIMMED_DIR}/${SRR}/${SRR}_2_unpaired.fastq.gz"

# ── Count input reads ───────────────────────────────────────────────────────
BEFORE_COUNT=$(fastq_read_count "$R1")
log_info "Input reads (R1): ${BEFORE_COUNT}"

# ── Run Trimmomatic ─────────────────────────────────────────────────────────
log_info "Running Trimmomatic PE..."
trimmomatic PE \
    -threads "${THREADS}" \
    -phred33 \
    "$R1" "$R2" \
    "$OUT_1P" "$OUT_1U" \
    "$OUT_2P" "$OUT_2U" \
    ILLUMINACLIP:"${TRIM_ILLUMINACLIP}" \
    LEADING:"${TRIM_LEADING}" \
    TRAILING:"${TRIM_TRAILING}" \
    SLIDINGWINDOW:"${TRIM_SLIDINGWINDOW}" \
    MINLEN:"${TRIM_MINLEN}" \
    2>&1 | tee "${TRIMMED_DIR}/${SRR}/${SRR}_trimmomatic.log"

# ── Validate outputs ────────────────────────────────────────────────────────
validate_files "$OUT_1P" "$OUT_2P"

# ── Report survival rate ────────────────────────────────────────────────────
AFTER_COUNT=$(fastq_read_count "$OUT_1P")
report_survival_rate "$BEFORE_COUNT" "$AFTER_COUNT" "Paired reads"

# ── Warn if survival rate is low ────────────────────────────────────────────
if (( BEFORE_COUNT > 0 )); then
    PCT=$(awk "BEGIN { printf \"%d\", ($AFTER_COUNT / $BEFORE_COUNT) * 100 }")
    if (( PCT < 70 )); then
        log_warn "Low survival rate (${PCT}%) — check adapter contamination or quality"
    fi
fi

log_info "Trimming complete."
checkpoint_set "${STEP_NAME}" "$SRR"
elapsed_time
