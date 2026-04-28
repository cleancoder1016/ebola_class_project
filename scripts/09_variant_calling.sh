#!/bin/bash
#SBATCH --job-name=09_variants
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --account=PWSU0516
#SBATCH --array=1-356%10
#SBATCH --output=logs/09_variants_%A_%a.log
# ═══════════════════════════════════════════════════════════════════════════════
# 09_variant_calling.sh — Per-sample variant calling with bcftools + SnpEff
# ═══════════════════════════════════════════════════════════════════════════════
# Pipeline:
#   1. bcftools mpileup  — pileup generation
#   2. bcftools call      — variant calling (multi-allelic)
#   3. bcftools filter    — quality filtering
#   4. bcftools norm      — normalize indels
#   5. SnpEff             — functional annotation
#   6. bcftools consensus — consensus FASTA sequence
#   7. bcftools stats     — variant summary
# ═══════════════════════════════════════════════════════════════════════════════

SCRIPT_DIR="/users/PWSU0516/mufakiransari/ebola_class_project/scripts"
source "/users/PWSU0516/mufakiransari/ebola_class_project/pipeline.config"
source "/users/PWSU0516/mufakiransari/ebola_class_project/scripts/utils.sh"

STEP_NAME="09_variant_calling"
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
module load "${MOD_BCFTOOLS}"
module load "${MOD_SAMTOOLS}"
module load "${MOD_HTSLIB}"
module load "${MOD_SNPEFF}"
ensure_dirs "${VARIANTS_DIR}/${SRR}"

# ── Locate BAM ──────────────────────────────────────────────────────────────
BAM="${ALIGNED_DIR}/${SRR}/${SRR}.dedup.bam"
if [[ ! -f "$BAM" ]]; then
    BAM="${ALIGNED_DIR}/${SRR}/${SRR}.sorted.bam"
fi
validate_files "$BAM" "${REFERENCE_FASTA}"

# ── Output files ────────────────────────────────────────────────────────────
RAW_VCF="${VARIANTS_DIR}/${SRR}/${SRR}.raw.vcf.gz"
FILTERED_VCF="${VARIANTS_DIR}/${SRR}/${SRR}.filtered.vcf.gz"
NORM_VCF="${VARIANTS_DIR}/${SRR}/${SRR}.norm.vcf.gz"
ANNOTATED_VCF="${VARIANTS_DIR}/${SRR}/${SRR}.annotated.vcf.gz"
CONSENSUS="${VARIANTS_DIR}/${SRR}/${SRR}.consensus.fasta"
STATS_FILE="${VARIANTS_DIR}/${SRR}/${SRR}.bcftools_stats.txt"

# ── Step 1+2: mpileup → call ───────────────────────────────────────────────
log_info "Calling variants with bcftools..."
bcftools mpileup \
    ${BCFTOOLS_MPILEUP_PARAMS} \
    --fasta-ref "${REFERENCE_FASTA}" \
    --annotate FORMAT/AD,FORMAT/DP \
    "$BAM" \
    | bcftools call \
        ${BCFTOOLS_CALL_PARAMS} \
        --ploidy 1 \
        -Oz \
        -o "${RAW_VCF}"

tabix -p vcf "${RAW_VCF}"

# ── Count raw variants ─────────────────────────────────────────────────────
RAW_COUNT=$(bcftools view -H "${RAW_VCF}" | wc -l)
log_info "Raw variants called: ${RAW_COUNT}"

# ── Step 3: Quality filtering ──────────────────────────────────────────────
log_info "Filtering variants (QUAL>${VCF_FILTER_QUAL}, DP>${VCF_FILTER_DP})..."
bcftools filter \
    -s "LowQual" \
    -e "QUAL<${VCF_FILTER_QUAL} || INFO/DP<${VCF_FILTER_DP}" \
    "${RAW_VCF}" \
    | bcftools view -f "PASS" -Oz -o "${FILTERED_VCF}"

tabix -p vcf "${FILTERED_VCF}"

FILTERED_COUNT=$(bcftools view -H "${FILTERED_VCF}" | wc -l)
log_info "Variants after filtering: ${FILTERED_COUNT}"

# ── Step 4: Normalize indels ───────────────────────────────────────────────
log_info "Normalizing variants..."
bcftools norm \
    -f "${REFERENCE_FASTA}" \
    -m -both \
    "${FILTERED_VCF}" \
    -Oz -o "${NORM_VCF}"

tabix -p vcf "${NORM_VCF}"

# ── Step 5: SnpEff annotation ──────────────────────────────────────────────
log_info "Annotating variants with SnpEff..."
# SnpEff may not have Ebola genome built-in — attempt annotation, skip on failure
if snpEff databases 2>/dev/null | grep -qi "ebola\|KJ660346"; then
    snpEff ann \
        -no-downstream -no-upstream -no-intergenic \
        "ebola" \
        "${NORM_VCF}" \
        2> "${VARIANTS_DIR}/${SRR}/${SRR}.snpeff.log" \
        | bgzip > "${ANNOTATED_VCF}"
    tabix -p vcf "${ANNOTATED_VCF}"
    log_info "SnpEff annotation complete."
else
    log_warn "SnpEff Ebola database not found — skipping annotation."
    cp "${NORM_VCF}" "${ANNOTATED_VCF}"
    tabix -p vcf "${ANNOTATED_VCF}"
fi

# ── Step 6: Generate consensus FASTA ───────────────────────────────────────
log_info "Generating consensus sequence..."
bcftools consensus \
    -f "${REFERENCE_FASTA}" \
    -s "${SRR}" \
    "${FILTERED_VCF}" \
    | sed "s/^>.*/>>${SRR}/" \
    > "${CONSENSUS}" 2>/dev/null || {
        # If consensus fails due to sample name, try without -s
        bcftools consensus \
            -f "${REFERENCE_FASTA}" \
            "${FILTERED_VCF}" \
            | sed "s/^>.*/>${SRR}/" \
            > "${CONSENSUS}"
    }

# ── Step 7: Variant statistics ─────────────────────────────────────────────
log_info "Generating variant statistics..."
bcftools stats "${FILTERED_VCF}" > "${STATS_FILE}"

# ── Write to shared summary ────────────────────────────────────────────────
VARIANT_SUMMARY_TSV="${VARIANTS_DIR}/variant_summary.tsv"
if [[ ! -f "${VARIANT_SUMMARY_TSV}" ]]; then
    echo -e "Sample\tRaw_Variants\tFiltered_Variants\tSNPs\tIndels" > "${VARIANT_SUMMARY_TSV}"
fi

SNPS=$(bcftools view -v snps -H "${FILTERED_VCF}" | wc -l)
INDELS=$(bcftools view -v indels -H "${FILTERED_VCF}" | wc -l)

(
    flock -x 200
    echo -e "${SRR}\t${RAW_COUNT}\t${FILTERED_COUNT}\t${SNPS}\t${INDELS}"
) 200>>"${VARIANT_SUMMARY_TSV}"

log_info "Variant summary: Raw=${RAW_COUNT} Filtered=${FILTERED_COUNT} SNPs=${SNPS} Indels=${INDELS}"
checkpoint_set "${STEP_NAME}" "$SRR"
elapsed_time
