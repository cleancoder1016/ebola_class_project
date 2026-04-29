#!/bin/bash
#SBATCH --job-name=02_fastqc_raw
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --account=PWSU0516
#SBATCH --array=1-356%50
#SBATCH --output=logs/02_fastqc_raw_%A_%a.log
# ═══════════════════════════════════════════════════════════════════════════════
# 02_fastqc_raw.sh — Quality control on raw FASTQ reads
# ═══════════════════════════════════════════════════════════════════════════════

SCRIPT_DIR="/users/PWSU0516/mufakiransari/ebola_class_project/scripts"
source "/users/PWSU0516/mufakiransari/ebola_class_project/pipeline.config"
source "/users/PWSU0516/mufakiransari/ebola_class_project/scripts/utils.sh"

STEP_NAME="02_fastqc_raw"
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
module load "${MOD_FASTQC}"
ensure_dirs "${QC_RAW_DIR}"

# ── Locate input files ──────────────────────────────────────────────────────
R1="${FASTQ_DIR}/${SRR}/${SRR}_1.fastq.gz"
R2="${FASTQ_DIR}/${SRR}/${SRR}_2.fastq.gz"

# Fall back to uncompressed if gzipped doesn't exist
if [[ ! -f "$R1" ]]; then
    R1="${FASTQ_DIR}/${SRR}/${SRR}_1.fastq"
    R2="${FASTQ_DIR}/${SRR}/${SRR}_2.fastq"
fi

validate_files "$R1" "$R2"

# ── Run FastQC ──────────────────────────────────────────────────────────────
log_info "Running FastQC on raw reads..."
fastqc "$R1" "$R2" \
    --outdir "${QC_RAW_DIR}" \
    --threads "${SLURM_CPUS_PER_TASK:-4}" \
    --quiet

# ── Validate outputs ────────────────────────────────────────────────────────
# FastQC names output based on input filename (without .gz)
log_info "FastQC raw QC complete."
checkpoint_set "${STEP_NAME}" "$SRR"
elapsed_time
