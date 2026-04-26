#!/bin/bash
#SBATCH --job-name=06_hisat2
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --account=PWSU0516
#SBATCH --array=1-50%10
#SBATCH --output=logs/06_hisat2_%A_%a.log
# ═══════════════════════════════════════════════════════════════════════════════
# 06_hisat2_align.sh — Align trimmed reads with HISAT2 + Picard deduplication
# ═══════════════════════════════════════════════════════════════════════════════
# Pipeline per sample:
#   1. HISAT2 alignment → SAM
#   2. samtools sort → sorted BAM (piped, no intermediate SAM)
#   3. Picard MarkDuplicates → deduplicated BAM
#   4. samtools index → BAI
# ═══════════════════════════════════════════════════════════════════════════════

SCRIPT_DIR="/users/PWSU0516/mufakiransari/ebola_class_project/scripts"
source "/users/PWSU0516/mufakiransari/ebola_class_project/pipeline.config"
source "/users/PWSU0516/mufakiransari/ebola_class_project/scripts/utils.sh"

STEP_NAME="06_hisat2_align"
SRR=$(get_accession)
CURRENT_SRR="$SRR"
print_job_info
timer_start

# ── Checkpoint check ────────────────────────────────────────────────────────
if checkpoint_exists "${STEP_NAME}" "$SRR"; then
    log_info "Already aligned. Skipping."
    exit 0
fi

# ── Setup ────────────────────────────────────────────────────────────────────
module load "${MOD_HISAT2}"
module load "${MOD_SAMTOOLS}"
module load "${MOD_PICARD}"
ensure_dirs "${ALIGNED_DIR}/${SRR}"

# ── Locate trimmed reads ────────────────────────────────────────────────────
R1="${TRIMMED_DIR}/${SRR}/${SRR}_1_paired.fastq.gz"
R2="${TRIMMED_DIR}/${SRR}/${SRR}_2_paired.fastq.gz"
validate_files "$R1" "$R2"

# ── Output files ────────────────────────────────────────────────────────────
SORTED_BAM="${ALIGNED_DIR}/${SRR}/${SRR}.sorted.bam"
DEDUP_BAM="${ALIGNED_DIR}/${SRR}/${SRR}.dedup.bam"
DEDUP_METRICS="${ALIGNED_DIR}/${SRR}/${SRR}.dedup_metrics.txt"
ALIGN_SUMMARY="${ALIGNED_DIR}/${SRR}/${SRR}.hisat2_summary.txt"

# ── Step 1+2: HISAT2 alignment → piped to samtools sort ────────────────────
log_info "Aligning with HISAT2 and sorting..."
hisat2 \
    -x "${HISAT2_INDEX_PREFIX}" \
    -1 "$R1" \
    -2 "$R2" \
    -p "${THREADS}" \
    ${HISAT2_PARAMS} \
    --rg-id "${SRR}" \
    --rg "SM:${SRR}" \
    --rg "PL:ILLUMINA" \
    --rg "LB:lib1" \
    --new-summary \
    --summary-file "${ALIGN_SUMMARY}" \
    2> "${ALIGNED_DIR}/${SRR}/${SRR}.hisat2.stderr" \
    | samtools sort \
        -@ "${THREADS}" \
        -m 2G \
        -O BAM \
        -o "${SORTED_BAM}" \
        -

validate_files "${SORTED_BAM}"

# ── Report alignment rate ───────────────────────────────────────────────────
if [[ -f "${ALIGN_SUMMARY}" ]]; then
    log_info "HISAT2 alignment summary:"
    cat "${ALIGN_SUMMARY}" | while IFS= read -r line; do
        log_info "  $line"
    done
fi

# ── Step 3: Picard MarkDuplicates ───────────────────────────────────────────
log_info "Marking duplicates with Picard..."
java -jar "${PICARD_JAR:-$(which picard 2>/dev/null || echo picard)}" MarkDuplicates \
    -I "${SORTED_BAM}" \
    -O "${DEDUP_BAM}" \
    -M "${DEDUP_METRICS}" \
    --REMOVE_DUPLICATES false \
    --CREATE_INDEX true \
    --VALIDATION_STRINGENCY LENIENT \
    --TMP_DIR "${ALIGNED_DIR}/${SRR}" \
    2>&1 | tee -a "${ALIGNED_DIR}/${SRR}/${SRR}.picard.log" || {
        # If Picard fails (e.g., low read count), fall back to using sorted BAM
        log_warn "Picard MarkDuplicates failed — using sorted BAM without deduplication."
        cp "${SORTED_BAM}" "${DEDUP_BAM}"
        samtools index "${DEDUP_BAM}"
    }

validate_files "${DEDUP_BAM}"

# ── Step 4: Index the final BAM ─────────────────────────────────────────────
if [[ ! -f "${DEDUP_BAM}.bai" && ! -f "${DEDUP_BAM%.*}.bai" ]]; then
    log_info "Indexing BAM..."
    samtools index "${DEDUP_BAM}"
fi

# ── Clean up intermediate sorted BAM ────────────────────────────────────────
if [[ -f "${DEDUP_BAM}" && "${DEDUP_BAM}" != "${SORTED_BAM}" ]]; then
    rm -f "${SORTED_BAM}"
    log_info "Removed intermediate sorted BAM."
fi

log_info "Alignment pipeline complete."
checkpoint_set "${STEP_NAME}" "$SRR"
elapsed_time
