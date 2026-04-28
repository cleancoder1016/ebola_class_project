#!/bin/bash
#SBATCH --job-name=17_hybrid_fc
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=48G
#SBATCH --account=PWSU0516
#SBATCH --output=logs/17_hybrid_featurecounts_%j.log
# ═══════════════════════════════════════════════════════════════════════════════
# 17_hybrid_featurecounts.sh — Gene-level quantification on hybrid genome
# ═══════════════════════════════════════════════════════════════════════════════
# Processes all hybrid-aligned BAMs at once using the hybrid GTF annotation.
# Produces a count matrix with ~30,000+ macaque genes AND 7 Ebola genes.
# ═══════════════════════════════════════════════════════════════════════════════

SCRIPT_DIR="/users/PWSU0516/mufakiransari/ebola_class_project/scripts"
source "/users/PWSU0516/mufakiransari/ebola_class_project/pipeline.config"
source "/users/PWSU0516/mufakiransari/ebola_class_project/scripts/utils.sh"

STEP_NAME="17_hybrid_featurecounts"
print_job_info
timer_start

# ── Checkpoint check ────────────────────────────────────────────────────────
if checkpoint_exists "${STEP_NAME}"; then
    log_info "Already completed. Skipping."
    exit 0
fi

# ── Setup ────────────────────────────────────────────────────────────────────
module load "${MOD_CONDA}"
ensure_dirs "${HYBRID_COUNTS_DIR}"

# ── Build list of all hybrid BAM files ──────────────────────────────────────
BAM_FILES=()
MISSING=0

while IFS= read -r SRR; do
    [[ -z "$SRR" ]] && continue
    BAM="${HYBRID_ALIGNED_DIR}/${SRR}/${SRR}.dedup.bam"
    if [[ ! -f "$BAM" ]]; then
        BAM="${HYBRID_ALIGNED_DIR}/${SRR}/${SRR}.sorted.bam"
    fi
    if [[ -f "$BAM" ]]; then
        BAM_FILES+=("$BAM")
    else
        log_warn "Hybrid BAM not found for ${SRR}, skipping."
        (( MISSING++ ))
    fi
done < "${ACCESSION_LIST}"

log_info "Found ${#BAM_FILES[@]} hybrid BAM files (${MISSING} missing)"

if (( ${#BAM_FILES[@]} == 0 )); then
    die "No hybrid BAM files found. Cannot run featureCounts."
fi

# ── Validate annotation ────────────────────────────────────────────────────
ANNOTATION="${HYBRID_GTF}"
if [[ ! -f "${ANNOTATION}" || ! -s "${ANNOTATION}" ]]; then
    die "Hybrid GTF not found: ${ANNOTATION}"
fi
validate_files "${ANNOTATION}"

# ── Run featureCounts ───────────────────────────────────────────────────────
# For the hybrid genome, we use gene_id as the attribute since the macaque
# GTF uses standard Ensembl/NCBI gene_id identifiers
log_info "Running featureCounts on ${#BAM_FILES[@]} hybrid BAM files..."

conda run -n "${CONDA_ENV}" featureCounts \
    -a "${ANNOTATION}" \
    -F GTF \
    -t "exon" \
    -g "gene_id" \
    ${FC_EXTRA_PARAMS} \
    -T "${THREADS}" \
    -o "${HYBRID_COUNTS_DIR}/gene_counts.txt" \
    "${BAM_FILES[@]}" \
    2>&1 | tee "${HYBRID_COUNTS_DIR}/featurecounts.log"

# ── Validate output ────────────────────────────────────────────────────────
validate_files "${HYBRID_COUNTS_DIR}/gene_counts.txt" "${HYBRID_COUNTS_DIR}/gene_counts.txt.summary"

# ── Create clean count matrix ───────────────────────────────────────────────
# featureCounts output has 6 metadata columns before the counts.
# Create a cleaned version with just Geneid and sample columns.
log_info "Creating clean count matrix..."
cut -f1,7- "${HYBRID_COUNTS_DIR}/gene_counts.txt" | \
    tail -n +2 | \
    sed 's|[^\t]*/||g; s|\.dedup\.bam||g; s|\.sorted\.bam||g' \
    > "${HYBRID_COUNTS_DIR}/gene_counts_clean.txt"

# ── Report summary statistics ──────────────────────────────────────────────
log_info "Hybrid featureCounts assignment summary:"
cat "${HYBRID_COUNTS_DIR}/gene_counts.txt.summary" | while IFS= read -r line; do
    log_info "  $line"
done

NUM_GENES=$(tail -n +3 "${HYBRID_COUNTS_DIR}/gene_counts.txt" | wc -l)
NUM_SAMPLES=$(head -2 "${HYBRID_COUNTS_DIR}/gene_counts.txt" | tail -1 | awk -F'\t' '{print NF - 6}')
log_info "Hybrid count matrix: ${NUM_GENES} genes × ${NUM_SAMPLES} samples"

checkpoint_set "${STEP_NAME}"
elapsed_time
