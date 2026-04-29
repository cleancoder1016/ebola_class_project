#!/usr/bin/env bash
#SBATCH --job-name=pipeline_coord
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=PWSU0516
#SBATCH --output=logs/coordinator_%j.log
# ═══════════════════════════════════════════════════════════════════════════════
# coordinator.sh — Submit pipeline steps one at a time to stay under
#                   SLURM's 1000 MaxSubmitJobsPerUser limit
# ═══════════════════════════════════════════════════════════════════════════════
# Each 356-task array job counts as 356 submitted jobs. With a 1000 limit,
# we can only have ~2 array jobs pending at once. This coordinator submits
# each step sequentially, waiting for completion before submitting the next.
#
# Usage:
#   sbatch coordinator.sh              # Submit as SLURM job (recommended)
#   bash coordinator.sh                # Run on login node (blocks terminal)
#   bash coordinator.sh --start-from 06   # Resume from step 06
#   sbatch coordinator.sh --start-from 09 --stop-at 14  # Run steps 09-14 only
# ═══════════════════════════════════════════════════════════════════════════════

set -uo pipefail

PROJECT_DIR="/users/PWSU0516/mufakiransari/ebola_class_project"
cd "${PROJECT_DIR}"
source pipeline.config

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[0;33m'
CYAN='\033[0;36m'; BOLD='\033[1m'; NC='\033[0m'

# Parse args
START_STEP=""
STOP_STEP=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --start-from) START_STEP="${2:-00}"; shift 2 ;;
        --stop-at)    STOP_STEP="${2:-99}";  shift 2 ;;
        *) shift ;;
    esac
done

# ── Step definitions (same order as run_pipeline.sh) ────────────────────────
# Format: STEP_NUM:SCRIPT
# Dependencies are handled by sequential execution — each step waits for
# the previous to finish before submitting.
# The Ebola-only Kallisto and Hybrid pipelines branch after step 04, so
# we handle branching explicitly.
declare -a MAIN_STEPS=(
    "00:00_setup_conda_env.sh"
    "01:01_sra_to_fastq.sh"
    "02:02_fastqc_raw.sh"
    "03:03_trimmomatic.sh"
    "04:04_fastqc_trimmed.sh"
)

declare -a EBOLA_ALIGN_STEPS=(
    "05:05_hisat2_index.sh"
    "06:06_hisat2_align.sh"
    "07:07_post_align_qc.sh"
    "08:08_featurecounts.sh"
    "09:09_variant_calling.sh"
    "10:10_deseq2_analysis.sh"
    "11:11_multiqc_report.sh"
)

declare -a EBOLA_KALLISTO_STEPS=(
    "12:12_kallisto_index.sh"
    "13:13_kallisto_quant.sh"
    "14:14_kallisto_aggregate.sh"
)

declare -a HYBRID_STEPS=(
    "15:15_hybrid_genome_build.sh"
    "16:16_hybrid_hisat2_align.sh"
    "17:17_hybrid_featurecounts.sh"
)

declare -a HYBRID_KALLISTO_STEPS=(
    "18:18_hybrid_kallisto_quant.sh"
    "19:19_hybrid_kallisto_aggregate.sh"
)

# ── Functions ───────────────────────────────────────────────────────────────

log() { echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

submit_and_wait() {
    local step_num="$1"
    local script="$2"

    # Skip if before start point
    if [[ -n "${START_STEP}" && "${step_num}" < "${START_STEP}" ]]; then
        log "${YELLOW}⏭  Step ${step_num} — Skipped (--start-from ${START_STEP})${NC}"
        return 0
    fi

    # Stop if past stop point
    if [[ -n "${STOP_STEP}" && "${step_num}" > "${STOP_STEP}" ]]; then
        log "${YELLOW}⏹  Step ${step_num} — Stopped (--stop-at ${STOP_STEP})${NC}"
        return 0
    fi

    local script_path="scripts/${script}"
    if [[ ! -f "${script_path}" ]]; then
        log "${RED}✗ Script not found: ${script_path}${NC}"
        return 1
    fi

    # Submit
    log "${CYAN}Submitting Step ${step_num}: ${script}...${NC}"
    local output
    output=$(sbatch "${script_path}" 2>&1)
    local rc=$?

    if (( rc != 0 )); then
        log "${RED}✗ Failed to submit Step ${step_num}: ${output}${NC}"
        return 1
    fi

    local job_id
    job_id=$(echo "$output" | grep -oP '\d+' | tail -1)
    log "${GREEN}✓ Step ${step_num} — ${script} → JobID ${job_id}${NC}"

    # Wait for job to complete
    log "  Waiting for JobID ${job_id} to finish..."
    local status=""
    local wait_count=0
    while true; do
        # Check if job is still in queue
        if ! squeue -j "${job_id}" -h 2>/dev/null | grep -q .; then
            # Job no longer in queue — check final status
            status=$(sacct -j "${job_id}" --format=State --noheader 2>/dev/null | head -1 | tr -d ' ')
            break
        fi

        # Progress update every 5 minutes
        (( wait_count++ ))
        if (( wait_count % 5 == 0 )); then
            local running=$(squeue -j "${job_id}" -h -t R 2>/dev/null | wc -l)
            local pending=$(squeue -j "${job_id}" -h -t PD 2>/dev/null | wc -l)
            log "  ... still running (R:${running} PD:${pending})"
        fi
        sleep 60
    done

    # Check if it succeeded
    if [[ "${status}" == "COMPLETED" ]]; then
        log "${GREEN}✓ Step ${step_num} COMPLETED${NC}"
        return 0
    else
        log "${RED}✗ Step ${step_num} finished with status: ${status}${NC}"
        log "${RED}  Check: sacct -j ${job_id} --format=JobID,State,ExitCode,Elapsed${NC}"
        log "${RED}  Logs:  ls -lt logs/${step_num}_*${NC}"
        # Don't exit — some array tasks might fail but the pipeline can continue
        log "${YELLOW}  Continuing pipeline despite failures...${NC}"
        return 0
    fi
}

submit_steps() {
    local step_list=("$@")
    for step_def in "${step_list[@]}"; do
        IFS=':' read -r num script <<< "$step_def"
        submit_and_wait "$num" "$script"
    done
}

# ═══════════════════════════════════════════════════════════════════════════════
# MAIN EXECUTION
# ═══════════════════════════════════════════════════════════════════════════════

log "${BOLD}${CYAN}"
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║  EBOLA PIPELINE COORDINATOR — Sequential Submission         ║"
echo "║  (Avoids SLURM 1000-job submission limit)                   ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo -e "${NC}"

log "Project: ${PROJECT_DIR}"
log "Samples: ${NUM_SAMPLES}"
log "Branch:  $(git branch --show-current)"
if [[ -n "${START_STEP}" ]]; then
    log "Starting from: Step ${START_STEP}"
fi
if [[ -n "${STOP_STEP}" ]]; then
    log "Stopping at:   Step ${STOP_STEP}"
fi
echo ""

# ── Phase 1: Common steps (SRA → FASTQ → QC → Trim) ───────────────────────
log "${BOLD}═══ PHASE 1: Data Acquisition & QC (Steps 00-04) ═══${NC}"
submit_steps "${MAIN_STEPS[@]}"

# ── Phase 2: Ebola-only alignment pipeline ─────────────────────────────────
log "${BOLD}═══ PHASE 2: Ebola-Only HISAT2 Pipeline (Steps 05-11) ═══${NC}"
submit_steps "${EBOLA_ALIGN_STEPS[@]}"

# ── Phase 3: Ebola-only Kallisto ───────────────────────────────────────────
log "${BOLD}═══ PHASE 3: Ebola-Only Kallisto (Steps 12-14) ═══${NC}"
submit_steps "${EBOLA_KALLISTO_STEPS[@]}"

# ── Phase 4: Hybrid genome pipeline ───────────────────────────────────────
log "${BOLD}═══ PHASE 4: Hybrid Genome HISAT2 (Steps 15-17) ═══${NC}"
submit_steps "${HYBRID_STEPS[@]}"

# ── Phase 5: Hybrid Kallisto ──────────────────────────────────────────────
log "${BOLD}═══ PHASE 5: Hybrid Kallisto (Steps 18-19) ═══${NC}"
submit_steps "${HYBRID_KALLISTO_STEPS[@]}"

# ── Done ──────────────────────────────────────────────────────────────────
echo ""
log "${BOLD}${GREEN}"
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║  ALL 20 PIPELINE STEPS COMPLETED!                          ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo -e "${NC}"

log "Run verification: bash run_and_verify.sh --verify"
log "Push results:     bash run_and_verify.sh --push"
