#!/bin/bash
#SBATCH --job-name=08_fcount
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --account=PWSU0516
#SBATCH --output=logs/08_featurecounts_%j.log
# ═══════════════════════════════════════════════════════════════════════════════
# 08_featurecounts.sh — Gene-level read quantification across all samples
# ═══════════════════════════════════════════════════════════════════════════════
# Single job that processes ALL BAMs at once (featureCounts supports
# multi-BAM input and produces a single count matrix).
# ═══════════════════════════════════════════════════════════════════════════════

SCRIPT_DIR="/users/PWSU0516/mufakiransari/ebola_class_project/scripts"
source "/users/PWSU0516/mufakiransari/ebola_class_project/pipeline.config"
source "/users/PWSU0516/mufakiransari/ebola_class_project/scripts/utils.sh"

STEP_NAME="08_featurecounts"
print_job_info
timer_start

# ── Checkpoint check ────────────────────────────────────────────────────────
if checkpoint_exists "${STEP_NAME}"; then
    log_info "Already completed. Skipping."
    exit 0
fi

# ── Setup ────────────────────────────────────────────────────────────────────
module load "${MOD_CONDA}"
ensure_dirs "${COUNTS_DIR}"

# ── Build list of all BAM files ─────────────────────────────────────────────
module load "${MOD_SAMTOOLS}"
BAM_FILES=()
MISSING=0
SKIPPED=0

while IFS= read -r SRR; do
    [[ -z "$SRR" ]] && continue
    BAM="${ALIGNED_DIR}/${SRR}/${SRR}.dedup.bam"
    if [[ ! -f "$BAM" ]]; then
        BAM="${ALIGNED_DIR}/${SRR}/${SRR}.sorted.bam"
    fi
    if [[ -f "$BAM" ]]; then
        # Quick check: read just 1 paired-end read (instant vs counting all reads)
        HAS_PAIRED=$(samtools view -f 1 "$BAM" 2>/dev/null | head -1 | wc -l)
        if (( HAS_PAIRED > 0 )); then
            BAM_FILES+=("$BAM")
        else
            log_warn "BAM for ${SRR} has no paired-end reads, skipping."
            (( SKIPPED++ )) || true
        fi
    else
        log_warn "BAM not found for ${SRR}, skipping."
        (( MISSING++ )) || true
    fi
done < "${ACCESSION_LIST}"

log_info "Found ${#BAM_FILES[@]} valid BAM files (${MISSING} missing, ${SKIPPED} skipped - no PE reads)"

if (( ${#BAM_FILES[@]} == 0 )); then
    die "No BAM files found. Cannot run featureCounts."
fi

# ── Validate annotation ────────────────────────────────────────────────────
ANNOTATION="${REFERENCE_GTF}"
if [[ ! -f "${ANNOTATION}" || ! -s "${ANNOTATION}" ]]; then
    # Fall back to GFF
    ANNOTATION="${REFERENCE_GFF}"
    FC_FORMAT_FLAG="-F GFF"
    log_warn "GTF not found, falling back to GFF: ${ANNOTATION}"
else
    FC_FORMAT_FLAG="-F GTF"
fi
validate_files "${ANNOTATION}"

# ── Run featureCounts ───────────────────────────────────────────────────────
log_info "Running featureCounts on ${#BAM_FILES[@]} BAM files..."

conda run -n "${CONDA_ENV}" featureCounts \
    -a "${ANNOTATION}" \
    ${FC_FORMAT_FLAG} \
    -t "${FC_FEATURE_TYPE}" \
    -g "${FC_ATTRIBUTE}" \
    ${FC_EXTRA_PARAMS} \
    -T "${THREADS}" \
    -o "${COUNTS_DIR}/gene_counts.txt" \
    "${BAM_FILES[@]}" \
    2>&1 | tee "${COUNTS_DIR}/featurecounts.log"

# ── Validate output ────────────────────────────────────────────────────────
validate_files "${COUNTS_DIR}/gene_counts.txt" "${COUNTS_DIR}/gene_counts.txt.summary"

# ── Create clean count matrix ───────────────────────────────────────────────
# featureCounts output has 6 metadata columns before the counts.
# Create a cleaned version with just Geneid and sample columns.
log_info "Creating clean count matrix..."
cut -f1,7- "${COUNTS_DIR}/gene_counts.txt" | \
    tail -n +2 | \
    sed 's|[^\t]*/||g; s|\.dedup\.bam||g; s|\.sorted\.bam||g' \
    > "${COUNTS_DIR}/gene_counts_clean.txt"

# ── Report summary statistics ──────────────────────────────────────────────
log_info "featureCounts assignment summary:"
cat "${COUNTS_DIR}/gene_counts.txt.summary" | while IFS= read -r line; do
    log_info "  $line"
done

NUM_GENES=$(tail -n +3 "${COUNTS_DIR}/gene_counts.txt" | wc -l)
NUM_SAMPLES=$(head -2 "${COUNTS_DIR}/gene_counts.txt" | tail -1 | awk -F'\t' '{print NF - 6}')
log_info "Count matrix: ${NUM_GENES} genes × ${NUM_SAMPLES} samples"

checkpoint_set "${STEP_NAME}"
elapsed_time
