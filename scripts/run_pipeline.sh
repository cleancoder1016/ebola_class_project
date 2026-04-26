#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════════════
# run_pipeline.sh — Master orchestration script for the Ebola RNA-seq pipeline
# ═══════════════════════════════════════════════════════════════════════════════
#
# Usage:
#   bash scripts/run_pipeline.sh                  # Run entire pipeline
#   bash scripts/run_pipeline.sh --start-from 06  # Start from step 06
#   bash scripts/run_pipeline.sh --only 03        # Run only step 03
#   bash scripts/run_pipeline.sh --dry-run        # Show commands without submitting
#   bash scripts/run_pipeline.sh --status         # Check status of submitted jobs
#   bash scripts/run_pipeline.sh --reset 03       # Clear checkpoint for step 03
#
# ═══════════════════════════════════════════════════════════════════════════════

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../pipeline.config"

# ── Color codes ─────────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[0;33m'
BLUE='\033[0;34m'; CYAN='\033[0;36m'; BOLD='\033[1m'; NC='\033[0m'

# ── Pipeline step definitions ───────────────────────────────────────────────
# Format: STEP_NUM:SCRIPT_NAME:DEPENDENCY_TYPE
# Dependency types: "none", "prev" (afterok on previous), "multi:X,Y" (afterok on steps X and Y)
declare -a STEPS=(
    "00:00_setup_conda_env.sh:none"
    "01:01_sra_to_fastq.sh:prev"
    "02:02_fastqc_raw.sh:prev"
    "03:03_trimmomatic.sh:prev"
    "04:04_fastqc_trimmed.sh:prev"
    "05:05_hisat2_index.sh:dep:00"
    "06:06_hisat2_align.sh:multi:04,05"
    "07:07_post_align_qc.sh:prev"
    "08:08_featurecounts.sh:prev"
    "09:09_variant_calling.sh:dep:07"
    "10:10_deseq2_analysis.sh:dep:08"
    "11:11_multiqc_report.sh:multi:09,10"
    "12:12_kallisto_index.sh:dep:00"
    "13:13_kallisto_quant.sh:multi:04,12"
    "14:14_kallisto_aggregate.sh:dep:13"
)

# ── Parse arguments ─────────────────────────────────────────────────────────
DRY_RUN=false
START_FROM=""
ONLY_STEP=""
CHECK_STATUS=false
RESET_STEP=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)     DRY_RUN=true ;;
        --start-from)  START_FROM="$2"; shift ;;
        --only)        ONLY_STEP="$2"; shift ;;
        --status)      CHECK_STATUS=true ;;
        --reset)       RESET_STEP="$2"; shift ;;
        -h|--help)
            echo "Usage: bash $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --dry-run          Show commands without submitting"
            echo "  --start-from NN    Start from step NN (skip earlier steps)"
            echo "  --only NN          Run only step NN"
            echo "  --status           Check status of submitted pipeline jobs"
            echo "  --reset NN         Clear checkpoint for step NN (forces re-run)"
            echo "  -h, --help         Show this help"
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
    shift
done

# ── Job tracking file ──────────────────────────────────────────────────────
JOB_TRACKER="${LOG_DIR}/pipeline_jobs.txt"
mkdir -p "${LOG_DIR}"

# ── Status check mode ─────────────────────────────────────────────────────
if $CHECK_STATUS; then
    echo -e "${BOLD}${CYAN}═══ Pipeline Job Status ═══${NC}"
    if [[ -f "${JOB_TRACKER}" ]]; then
        echo ""
        while IFS=$'\t' read -r step_num job_id script_name timestamp; do
            status=$(sacct -j "$job_id" --format=State --noheader 2>/dev/null | head -1 | tr -d ' ' || echo "UNKNOWN")
            case "$status" in
                COMPLETED)  color="${GREEN}" ;;
                RUNNING)    color="${BLUE}" ;;
                PENDING)    color="${YELLOW}" ;;
                FAILED|CANCELLED|TIMEOUT) color="${RED}" ;;
                *)          color="${NC}" ;;
            esac
            printf "  Step %s %-30s JobID=%-10s ${color}%s${NC}\n" \
                "$step_num" "$script_name" "$job_id" "$status"
        done < "${JOB_TRACKER}"
        echo ""
    else
        echo "  No jobs submitted yet."
    fi
    exit 0
fi

# ── Reset mode ─────────────────────────────────────────────────────────────
if [[ -n "${RESET_STEP}" ]]; then
    CHECKPOINT_PATTERN="${CHECKPOINT_DIR}/*${RESET_STEP}*"
    echo -e "${YELLOW}Clearing checkpoints matching: ${CHECKPOINT_PATTERN}${NC}"
    rm -rf ${CHECKPOINT_PATTERN} 2>/dev/null || true
    echo -e "${GREEN}Checkpoints cleared.${NC}"
    exit 0
fi

# ── Pre-flight checks ─────────────────────────────────────────────────────
echo -e "${BOLD}${CYAN}"
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║           EBOLA RNA-SEQ PIPELINE — SUBMISSION               ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo -e "${NC}"

echo -e "${BOLD}Pre-flight checks:${NC}"

# Check accession list
if [[ ! -f "${ACCESSION_LIST}" ]]; then
    echo -e "  ${RED}✗ Accession list not found: ${ACCESSION_LIST}${NC}"
    exit 1
fi
SAMPLE_COUNT=$(grep -c . "${ACCESSION_LIST}" || true)
echo -e "  ${GREEN}✓${NC} Accession list: ${SAMPLE_COUNT} samples"

# Check scripts exist
MISSING_SCRIPTS=0
for step_def in "${STEPS[@]}"; do
    IFS=':' read -r num script deps <<< "$step_def"
    if [[ ! -f "${SCRIPT_DIR}/${script}" ]]; then
        echo -e "  ${RED}✗ Missing script: ${script}${NC}"
        ((MISSING_SCRIPTS++))
    fi
done
if (( MISSING_SCRIPTS > 0 )); then
    echo -e "  ${RED}${MISSING_SCRIPTS} script(s) missing. Aborting.${NC}"
    exit 1
fi
echo -e "  ${GREEN}✓${NC} All ${#STEPS[@]} pipeline scripts found"

# Check SRA files
SRA_COUNT=$(ls -d "${SRA_DIR}"/SRR* 2>/dev/null | wc -l)
echo -e "  ${GREEN}✓${NC} SRA files: ${SRA_COUNT} found"

# Check disk space (warn if < 10GB free)
FREE_GB=$(df --output=avail "${PROJECT_DIR}" | tail -1 | awk '{print int($1/1048576)}')
if (( FREE_GB < 10 )); then
    echo -e "  ${YELLOW}⚠ Low disk space: ${FREE_GB} GB free${NC}"
else
    echo -e "  ${GREEN}✓${NC} Disk space: ${FREE_GB} GB available"
fi

echo ""

# ── Submit jobs ────────────────────────────────────────────────────────────
declare -A JOB_IDS  # step_num → job_id

# Clear job tracker for new run
if ! $DRY_RUN; then
    > "${JOB_TRACKER}"
fi

submit_step() {
    local step_num="$1"
    local script="$2"
    local dep_spec="$3"

    # Skip logic
    if [[ -n "${START_FROM}" && "${step_num}" < "${START_FROM}" ]]; then
        echo -e "  ${YELLOW}⏭  Step ${step_num} — Skipped (--start-from ${START_FROM})${NC}"
        return 0
    fi
    if [[ -n "${ONLY_STEP}" && "${step_num}" != "${ONLY_STEP}" ]]; then
        return 0
    fi

    # Build dependency string
    local dep_flag=""
    case "$dep_spec" in
        none)
            dep_flag=""
            ;;
        prev)
            # Find the previous step's job ID
            local prev_num
            prev_num=$(printf "%02d" $(( 10#${step_num} - 1 )))
            if [[ -n "${JOB_IDS[$prev_num]:-}" ]]; then
                dep_flag="--dependency=afterok:${JOB_IDS[$prev_num]}"
            fi
            ;;
        multi:*)
            # Dependencies on multiple specific steps
            local dep_steps="${dep_spec#multi:}"
            local dep_ids=""
            IFS=',' read -ra DEPS <<< "$dep_steps"
            for d in "${DEPS[@]}"; do
                d_padded=$(printf "%02d" $((10#$d)))
                if [[ -n "${JOB_IDS[$d_padded]:-}" ]]; then
                    [[ -n "$dep_ids" ]] && dep_ids="${dep_ids}:"
                    dep_ids="${dep_ids}${JOB_IDS[$d_padded]}"
                fi
            done
            [[ -n "$dep_ids" ]] && dep_flag="--dependency=afterok:${dep_ids}"
            ;;
        dep:*)
            # Dependency on a single specific step
            local dep_step="${dep_spec#dep:}"
            dep_step=$(printf "%02d" $((10#$dep_step)))
            if [[ -n "${JOB_IDS[$dep_step]:-}" ]]; then
                dep_flag="--dependency=afterok:${JOB_IDS[$dep_step]}"
            fi
            ;;
    esac

    local cmd="sbatch ${dep_flag} ${SCRIPT_DIR}/${script}"

    if $DRY_RUN; then
        echo -e "  ${BLUE}[DRY-RUN]${NC} Step ${step_num}: ${cmd}"
    else
        local output
        output=$(eval "${cmd}" 2>&1)
        local job_id
        job_id=$(echo "$output" | grep -oP '\d+' | tail -1)

        if [[ -n "$job_id" ]]; then
            JOB_IDS[$step_num]="$job_id"
            echo -e "  ${GREEN}✓${NC} Step ${step_num} — ${script} → JobID ${BOLD}${job_id}${NC} ${dep_flag:+(${dep_flag})}"
            echo -e "${step_num}\t${job_id}\t${script}\t$(date '+%Y-%m-%d %H:%M:%S')" >> "${JOB_TRACKER}"
        else
            echo -e "  ${RED}✗${NC} Step ${step_num} — Failed to submit: ${output}"
            return 1
        fi
    fi
}

echo -e "${BOLD}Submitting pipeline jobs:${NC}"
echo ""

for step_def in "${STEPS[@]}"; do
    IFS=':' read -r num script deps_raw <<< "$step_def"
    # Rejoin deps (they may contain colons)
    deps="${step_def#*:*:}"
    submit_step "$num" "$script" "$deps"
done

echo ""
if $DRY_RUN; then
    echo -e "${YELLOW}${BOLD}DRY RUN — no jobs submitted.${NC}"
else
    echo -e "${GREEN}${BOLD}All jobs submitted successfully!${NC}"
    echo ""
    echo -e "Monitor progress:"
    echo -e "  ${CYAN}bash scripts/run_pipeline.sh --status${NC}"
    echo -e "  ${CYAN}squeue -u \$USER${NC}"
    echo -e "  ${CYAN}tail -f logs/*.log${NC}"
    echo ""
    echo -e "Job IDs saved to: ${JOB_TRACKER}"
fi
