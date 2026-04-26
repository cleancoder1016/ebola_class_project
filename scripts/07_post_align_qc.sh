#!/bin/bash
#SBATCH --job-name=07_align_qc
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --account=PWSU0516
#SBATCH --array=1-50%10
#SBATCH --output=logs/07_align_qc_%A_%a.log
# ═══════════════════════════════════════════════════════════════════════════════
# 07_post_align_qc.sh — Comprehensive post-alignment quality control
# ═══════════════════════════════════════════════════════════════════════════════
# Per sample:
#   1. samtools flagstat   — alignment summary
#   2. samtools stats      — detailed metrics (insert size, GC, etc.)
#   3. samtools idxstats   — per-reference read counts
#   4. samtools coverage   — genome coverage breadth/depth
#   5. Picard CollectInsertSizeMetrics — insert size distribution
#   6. Write summary line to shared TSV
# ═══════════════════════════════════════════════════════════════════════════════

SCRIPT_DIR="/users/PWSU0516/mufakiransari/ebola_class_project/scripts"
source "/users/PWSU0516/mufakiransari/ebola_class_project/pipeline.config"
source "/users/PWSU0516/mufakiransari/ebola_class_project/scripts/utils.sh"

STEP_NAME="07_post_align_qc"
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
module load "${MOD_SAMTOOLS}"
module load "${MOD_PICARD}"
ensure_dirs "${QC_ALIGN_DIR}/${SRR}"

# ── Locate BAM ──────────────────────────────────────────────────────────────
BAM="${ALIGNED_DIR}/${SRR}/${SRR}.dedup.bam"
if [[ ! -f "$BAM" ]]; then
    BAM="${ALIGNED_DIR}/${SRR}/${SRR}.sorted.bam"
fi
validate_files "$BAM"

# ── 1. samtools flagstat ────────────────────────────────────────────────────
log_info "Running samtools flagstat..."
samtools flagstat -@ "${SLURM_CPUS_PER_TASK:-4}" "$BAM" \
    > "${QC_ALIGN_DIR}/${SRR}/${SRR}.flagstat.txt"

# ── 2. samtools stats ──────────────────────────────────────────────────────
log_info "Running samtools stats..."
samtools stats -@ "${SLURM_CPUS_PER_TASK:-4}" "$BAM" \
    > "${QC_ALIGN_DIR}/${SRR}/${SRR}.stats.txt"

# ── 3. samtools idxstats ────────────────────────────────────────────────────
log_info "Running samtools idxstats..."
samtools idxstats "$BAM" \
    > "${QC_ALIGN_DIR}/${SRR}/${SRR}.idxstats.txt"

# ── 4. samtools coverage ────────────────────────────────────────────────────
log_info "Calculating genome coverage..."
samtools coverage "$BAM" \
    > "${QC_ALIGN_DIR}/${SRR}/${SRR}.coverage.txt"

# ── 5. samtools depth (mean depth) ──────────────────────────────────────────
log_info "Calculating mean depth..."
MEAN_DEPTH=$(samtools depth -a "$BAM" | awk '{sum+=$3; n++} END {if(n>0) printf "%.1f", sum/n; else print "0"}')
log_info "Mean genome depth: ${MEAN_DEPTH}x"

# ── 6. Picard CollectInsertSizeMetrics ──────────────────────────────────────
log_info "Collecting insert size metrics..."
java -jar "${PICARD_JAR:-$(which picard 2>/dev/null || echo picard)}" CollectInsertSizeMetrics \
    -I "$BAM" \
    -O "${QC_ALIGN_DIR}/${SRR}/${SRR}.insert_size_metrics.txt" \
    -H "${QC_ALIGN_DIR}/${SRR}/${SRR}.insert_size_histogram.pdf" \
    --VALIDATION_STRINGENCY LENIENT \
    2>&1 || {
        log_warn "Picard InsertSizeMetrics failed (may be single-end or low reads)."
    }

# ── 7. BAM integrity check ─────────────────────────────────────────────────
log_info "Checking BAM integrity..."
if samtools quickcheck "$BAM"; then
    log_info "BAM integrity: OK ✓"
else
    log_warn "BAM integrity check FAILED for ${SRR}"
fi

# ── 8. Extract key metrics and write to shared summary ─────────────────────
SUMMARY_TSV="${QC_ALIGN_DIR}/alignment_summary.tsv"

# Parse flagstat
TOTAL_READS=$(grep "in total" "${QC_ALIGN_DIR}/${SRR}/${SRR}.flagstat.txt" | awk '{print $1}')
MAPPED_READS=$(grep "mapped (" "${QC_ALIGN_DIR}/${SRR}/${SRR}.flagstat.txt" | head -1 | awk '{print $1}')
MAPPED_PCT=$(grep "mapped (" "${QC_ALIGN_DIR}/${SRR}/${SRR}.flagstat.txt" | head -1 | grep -oP '\([\d.]+' | tr -d '(')
PROPERLY_PAIRED=$(grep "properly paired" "${QC_ALIGN_DIR}/${SRR}/${SRR}.flagstat.txt" | awk '{print $1}')

# Write header if file doesn't exist
if [[ ! -f "${SUMMARY_TSV}" ]]; then
    echo -e "Sample\tTotal_Reads\tMapped_Reads\tMapped_Pct\tProperly_Paired\tMean_Depth" \
        > "${SUMMARY_TSV}"
fi

# Append sample summary (use flock to prevent race conditions in array jobs)
(
    flock -x 200
    echo -e "${SRR}\t${TOTAL_READS}\t${MAPPED_READS}\t${MAPPED_PCT}\t${PROPERLY_PAIRED}\t${MEAN_DEPTH}"
) 200>>"${SUMMARY_TSV}"

log_info "Alignment QC summary: Total=${TOTAL_READS} Mapped=${MAPPED_READS} (${MAPPED_PCT}%) Depth=${MEAN_DEPTH}x"
checkpoint_set "${STEP_NAME}" "$SRR"
elapsed_time
