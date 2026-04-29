#!/bin/bash
#SBATCH --job-name=05_hisat2_idx
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --account=PWSU0516
#SBATCH --output=logs/05_hisat2_index_%j.log
# ═══════════════════════════════════════════════════════════════════════════════
# 05_hisat2_index.sh — Download Ebola reference + build HISAT2 index
# ═══════════════════════════════════════════════════════════════════════════════
# Downloads the Ebola Makona reference genome (KJ660346.2) from NCBI,
# converts GFF → GTF, extracts splice sites, and builds the HISAT2 index.
# ═══════════════════════════════════════════════════════════════════════════════

SCRIPT_DIR="/users/PWSU0516/mufakiransari/ebola_class_project/scripts"
source "/users/PWSU0516/mufakiransari/ebola_class_project/pipeline.config"
source "/users/PWSU0516/mufakiransari/ebola_class_project/scripts/utils.sh"

STEP_NAME="05_hisat2_index"
print_job_info
timer_start

# ── Checkpoint check ────────────────────────────────────────────────────────
if checkpoint_exists "${STEP_NAME}"; then
    log_info "HISAT2 index already built. Skipping."
    exit 0
fi

# ── Setup ────────────────────────────────────────────────────────────────────
module load "${MOD_HISAT2}"
module load "${MOD_SAMTOOLS}"
module load "${MOD_CONDA}"
ensure_dirs "${REFERENCE_DIR}" "${REFERENCE_DIR}/hisat2_index"

# ── Download reference genome FASTA ─────────────────────────────────────────
if [[ -f "${REFERENCE_FASTA}" ]]; then
    log_info "Reference FASTA already exists."
else
    log_info "Downloading Ebola reference genome (${REFERENCE_ID})..."

    # Use NCBI Datasets CLI or efetch — fall back to direct URL
    NCBI_FASTA_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${REFERENCE_ID}&rettype=fasta&retmode=text"
    retry 3 curl -sS -o "${REFERENCE_FASTA}" "${NCBI_FASTA_URL}"
    validate_files "${REFERENCE_FASTA}"
    log_info "Reference FASTA downloaded: $(wc -c < "${REFERENCE_FASTA}") bytes"
fi

# ── Download GFF3 annotation ────────────────────────────────────────────────
if [[ -f "${REFERENCE_GFF}" ]]; then
    log_info "GFF annotation already exists."
else
    log_info "Downloading GFF3 annotation..."

    NCBI_GFF_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${REFERENCE_ID}&rettype=gff3&retmode=text"
    retry 3 curl -sS -o "${REFERENCE_GFF}" "${NCBI_GFF_URL}"
    validate_files "${REFERENCE_GFF}"
    log_info "GFF3 annotation downloaded."
fi

# ── Convert GFF3 → GTF ─────────────────────────────────────────────────────
if [[ -f "${REFERENCE_GTF}" ]]; then
    log_info "GTF annotation already exists."
else
    log_info "Converting GFF3 → GTF using gffread..."
    conda run -n "${CONDA_ENV}" gffread "${REFERENCE_GFF}" -T -o "${REFERENCE_GTF}"
    validate_files "${REFERENCE_GTF}"
    log_info "GTF conversion complete."
fi

# ── Create FASTA index ──────────────────────────────────────────────────────
if [[ ! -f "${REFERENCE_FASTA}.fai" ]]; then
    log_info "Creating samtools FASTA index..."
    samtools faidx "${REFERENCE_FASTA}"
fi

# ── Build HISAT2 index ──────────────────────────────────────────────────────
# NOTE: Ebola is a non-segmented negative-sense RNA virus with NO introns.
# Splice sites and exon annotations are not applicable and cause HISAT2-build
# to crash with 'Nongraph exception' on this tiny ~19kb genome.
# We build a simple (non-graph) index instead.
log_info "Building HISAT2 index (simple mode — viral genome, no introns)..."
hisat2-build -p "${SLURM_CPUS_PER_TASK:-4}" "${REFERENCE_FASTA}" "${HISAT2_INDEX_PREFIX}"

# ── Validate index files ────────────────────────────────────────────────────
INDEX_FILES=$(ls "${HISAT2_INDEX_PREFIX}"*.ht2 2>/dev/null | wc -l)
if (( INDEX_FILES == 0 )); then
    die "HISAT2 index build failed — no .ht2 files found"
fi
log_info "HISAT2 index built successfully (${INDEX_FILES} index files)"

log_info "Reference setup complete."
checkpoint_set "${STEP_NAME}"
elapsed_time
