#!/bin/bash
#SBATCH --job-name=11_multiqc
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --account=PWSU0516
#SBATCH --output=logs/11_multiqc_%j.log
# ═══════════════════════════════════════════════════════════════════════════════
# 11_multiqc_report.sh — Aggregate all QC reports into a single HTML report
# ═══════════════════════════════════════════════════════════════════════════════
# Aggregates outputs from:
#   - FastQC (raw + trimmed)
#   - Trimmomatic logs
#   - HISAT2 alignment summaries
#   - Picard MarkDuplicates metrics
#   - Picard InsertSizeMetrics
#   - samtools flagstat / stats / idxstats
#   - featureCounts assignment summary
#   - bcftools stats
# ═══════════════════════════════════════════════════════════════════════════════

SCRIPT_DIR="/users/PWSU0516/mufakiransari/ebola_class_project/scripts"
source "/users/PWSU0516/mufakiransari/ebola_class_project/pipeline.config"
source "/users/PWSU0516/mufakiransari/ebola_class_project/scripts/utils.sh"

STEP_NAME="11_multiqc"
print_job_info
timer_start

# ── Checkpoint check ────────────────────────────────────────────────────────
if checkpoint_exists "${STEP_NAME}"; then
    log_info "Already completed. Skipping."
    exit 0
fi

# ── Setup ────────────────────────────────────────────────────────────────────
module load "${MOD_CONDA}"
ensure_dirs "${MULTIQC_DIR}"

# ── Run MultiQC ────────────────────────────────────────────────────────────
log_info "Running MultiQC across all pipeline outputs..."
conda run -n "${CONDA_ENV}" multiqc \
    "${QC_RAW_DIR}" \
    "${QC_TRIMMED_DIR}" \
    "${TRIMMED_DIR}" \
    "${ALIGNED_DIR}" \
    "${QC_ALIGN_DIR}" \
    "${COUNTS_DIR}" \
    "${VARIANTS_DIR}" \
    "${LOG_DIR}" \
    --outdir "${MULTIQC_DIR}" \
    --filename "multiqc_report" \
    --title "Ebola RNA-seq Pipeline — QC Report (PRJNA938511)" \
    --comment "50 samples from the 2014 West African Ebola outbreak. Pipeline: SRA → FastQC → Trimmomatic → HISAT2 → Picard → featureCounts → bcftools." \
    --force \
    --no-data-dir \
    2>&1 | tee "${MULTIQC_DIR}/multiqc.log"

# ── Validate output ────────────────────────────────────────────────────────
if [[ -f "${MULTIQC_DIR}/multiqc_report.html" ]]; then
    log_info "MultiQC report generated successfully."
    log_info "Report: ${MULTIQC_DIR}/multiqc_report.html"
    log_info "Report size: $(du -h "${MULTIQC_DIR}/multiqc_report.html" | cut -f1)"
else
    log_warn "MultiQC report not found — check ${MULTIQC_DIR}/multiqc.log"
fi

# ── Print final pipeline summary ───────────────────────────────────────────
log_info ""
log_info "╔══════════════════════════════════════════════════════╗"
log_info "║       EBOLA RNA-SEQ PIPELINE COMPLETE               ║"
log_info "╚══════════════════════════════════════════════════════╝"
log_info ""
log_info "Key outputs:"
log_info "  📊 MultiQC Report:    ${MULTIQC_DIR}/multiqc_report.html"
log_info "  🔢 Count Matrix:      ${COUNTS_DIR}/gene_counts_clean.txt"
log_info "  📈 DESeq2 Results:    ${DESEQ2_DIR}/"
log_info "  🔬 Variant Calls:     ${VARIANTS_DIR}/"
log_info "  🎯 Aligned BAMs:      ${ALIGNED_DIR}/"
log_info "  📋 Alignment Summary: ${QC_ALIGN_DIR}/alignment_summary.tsv"
log_info "  📋 Variant Summary:   ${VARIANTS_DIR}/variant_summary.tsv"
log_info ""

checkpoint_set "${STEP_NAME}"
elapsed_time
