#!/bin/bash
#SBATCH --job-name=15_hybrid_build
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --account=PWSU0516
#SBATCH --output=logs/15_hybrid_genome_build_%j.log
# ═══════════════════════════════════════════════════════════════════════════════
# 15_hybrid_genome_build.sh — Download macaque genome and build hybrid
#                              (Macaque Mmul_10 + Ebola KJ660346.2)
# ═══════════════════════════════════════════════════════════════════════════════
# Creates a combined reference for dual host+pathogen transcriptome analysis.
# Steps:
#   1. Download macaque genome FASTA (Mmul_10)
#   2. Download macaque GTF annotation
#   3. Concatenate macaque + Ebola FASTA → hybrid FASTA
#   4. Concatenate macaque + Ebola GTF  → hybrid GTF
#   5. Build HISAT2 index for hybrid genome
#   6. Extract transcriptome and build Kallisto index for hybrid
# ═══════════════════════════════════════════════════════════════════════════════

SCRIPT_DIR="/users/PWSU0516/mufakiransari/ebola_class_project/scripts"
source "/users/PWSU0516/mufakiransari/ebola_class_project/pipeline.config"
source "/users/PWSU0516/mufakiransari/ebola_class_project/scripts/utils.sh"

STEP_NAME="15_hybrid_genome_build"
print_job_info
timer_start

# ── Checkpoint check ────────────────────────────────────────────────────────
if checkpoint_exists "${STEP_NAME}"; then
    log_info "Hybrid genome already built. Skipping."
    exit 0
fi

# ── Setup ────────────────────────────────────────────────────────────────────
module load "${MOD_HISAT2}"
module load "${MOD_SAMTOOLS}"
module load "${MOD_CONDA}"
ensure_dirs "${HYBRID_GENOME_DIR}" "${HYBRID_GENOME_DIR}/hisat2_index"

# ── Temporary files for macaque downloads ────────────────────────────────────
MACAQUE_FASTA_GZ="${HYBRID_GENOME_DIR}/macaque_genomic.fna.gz"
MACAQUE_FASTA="${HYBRID_GENOME_DIR}/macaque_genomic.fna"
MACAQUE_GTF_GZ="${HYBRID_GENOME_DIR}/macaque_genomic.gtf.gz"
MACAQUE_GTF="${HYBRID_GENOME_DIR}/macaque_genomic.gtf"

# ── Step 1: Download macaque genome FASTA ───────────────────────────────────
if [[ -f "${MACAQUE_FASTA}" ]]; then
    log_info "Macaque FASTA already downloaded."
else
    log_info "Downloading macaque genome (Mmul_10, GCF_003339765.3)..."
    retry 3 curl -sS -L -o "${MACAQUE_FASTA_GZ}" "${MACAQUE_FASTA_URL}"
    validate_files "${MACAQUE_FASTA_GZ}"
    log_info "Decompressing macaque FASTA..."
    gunzip -f "${MACAQUE_FASTA_GZ}"
    validate_files "${MACAQUE_FASTA}"
    log_info "Macaque FASTA ready: $(wc -c < "${MACAQUE_FASTA}") bytes"
fi

# ── Step 2: Download macaque GTF annotation ─────────────────────────────────
if [[ -f "${MACAQUE_GTF}" ]]; then
    log_info "Macaque GTF already downloaded."
else
    log_info "Downloading macaque GTF annotation..."
    retry 3 curl -sS -L -o "${MACAQUE_GTF_GZ}" "${MACAQUE_GTF_URL}"
    validate_files "${MACAQUE_GTF_GZ}"
    log_info "Decompressing macaque GTF..."
    gunzip -f "${MACAQUE_GTF_GZ}"
    validate_files "${MACAQUE_GTF}"
    log_info "Macaque GTF ready: $(wc -l < "${MACAQUE_GTF}") lines"
fi

# ── Step 3: Concatenate FASTA (macaque + Ebola) ────────────────────────────
if [[ -f "${HYBRID_FASTA}" ]]; then
    log_info "Hybrid FASTA already exists."
else
    log_info "Building hybrid FASTA (macaque + Ebola)..."
    validate_files "${REFERENCE_FASTA}"
    cat "${MACAQUE_FASTA}" "${REFERENCE_FASTA}" > "${HYBRID_FASTA}"
    validate_files "${HYBRID_FASTA}"
    log_info "Hybrid FASTA created: $(grep -c '^>' "${HYBRID_FASTA}") sequences"
fi

# ── Step 4: Concatenate GTF (macaque + Ebola) ──────────────────────────────
if [[ -f "${HYBRID_GTF}" ]]; then
    log_info "Hybrid GTF already exists."
else
    log_info "Building hybrid GTF (macaque + Ebola)..."
    validate_files "${REFERENCE_GTF}"
    cat "${MACAQUE_GTF}" "${REFERENCE_GTF}" > "${HYBRID_GTF}"
    validate_files "${HYBRID_GTF}"
    log_info "Hybrid GTF created: $(wc -l < "${HYBRID_GTF}") lines"
fi

# ── Step 5: Build HISAT2 index for hybrid genome ───────────────────────────
# The macaque genome is ~3 GB; HISAT2-build needs ~64 GB RAM for large genomes.
log_info "Building HISAT2 index for hybrid genome (this will take ~1-2 hours)..."
hisat2-build \
    -p "${SLURM_CPUS_PER_TASK:-8}" \
    "${HYBRID_FASTA}" \
    "${HISAT2_HYBRID_INDEX_PREFIX}"

INDEX_FILES=$(ls "${HISAT2_HYBRID_INDEX_PREFIX}"*.ht2 2>/dev/null | wc -l)
if (( INDEX_FILES == 0 )); then
    die "HISAT2 hybrid index build failed — no .ht2 files found"
fi
log_info "HISAT2 hybrid index built successfully (${INDEX_FILES} index files)"

# ── Step 6: Build Kallisto index for hybrid transcriptome ──────────────────
log_info "Extracting hybrid transcriptome FASTA using gffread..."
conda run -n "${CONDA_ENV}" gffread -w "${HYBRID_TRANSCRIPTOME_FASTA}" \
    -g "${HYBRID_FASTA}" "${HYBRID_GTF}"
validate_files "${HYBRID_TRANSCRIPTOME_FASTA}"

log_info "Building Kallisto index for hybrid transcriptome..."
conda run -n "${CONDA_ENV}" kallisto index \
    -i "${KALLISTO_HYBRID_INDEX}" \
    "${HYBRID_TRANSCRIPTOME_FASTA}"
validate_files "${KALLISTO_HYBRID_INDEX}"

# ── Create FASTA index for hybrid genome ────────────────────────────────────
log_info "Creating samtools FASTA index for hybrid genome..."
samtools faidx "${HYBRID_FASTA}"

log_info "Hybrid genome build complete."
checkpoint_set "${STEP_NAME}"
elapsed_time
