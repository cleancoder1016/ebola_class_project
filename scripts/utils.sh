#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════════════
# utils.sh — Shared utility library for the Ebola RNA-seq pipeline
# ═══════════════════════════════════════════════════════════════════════════════
# Source this file after pipeline.config in every script:
#   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#   source "${SCRIPT_DIR}/../pipeline.config"
#   source "${SCRIPT_DIR}/utils.sh"
# ═══════════════════════════════════════════════════════════════════════════════

set -euo pipefail

# ── Color codes (disabled if not a terminal) ─────────────────────────────────
if [[ -t 1 ]]; then
    RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[0;33m'
    BLUE='\033[0;34m'; CYAN='\033[0;36m'; BOLD='\033[1m'; NC='\033[0m'
else
    RED=''; GREEN=''; YELLOW=''; BLUE=''; CYAN=''; BOLD=''; NC=''
fi

# ── Timestamped Logging ──────────────────────────────────────────────────────
# Usage: log_info "message"   log_warn "message"   log_error "message"
_log() {
    local level="$1"; shift
    local color="$1"; shift
    local timestamp
    timestamp="$(date '+%Y-%m-%d %H:%M:%S')"
    local prefix=""
    [[ -n "${STEP_NAME:-}" ]] && prefix="[${STEP_NAME}] "
    [[ -n "${CURRENT_SRR:-}" ]] && prefix="${prefix}[${CURRENT_SRR}] "
    echo -e "${color}${BOLD}[${level}]${NC} ${CYAN}${timestamp}${NC} ${prefix}$*"
}

log_info()  { _log "INFO"  "${GREEN}"  "$@"; }
log_warn()  { _log "WARN"  "${YELLOW}" "$@"; }
log_error() { _log "ERROR" "${RED}"    "$@"; }

# ── Die: log error and exit ──────────────────────────────────────────────────
die() {
    log_error "$@"
    exit 1
}

# ── Timer Functions ──────────────────────────────────────────────────────────
# Usage:
#   timer_start
#   ... do work ...
#   elapsed_time   # prints "Elapsed: 0h 5m 32s"
_TIMER_START=0
timer_start() { _TIMER_START=$(date +%s); }
elapsed_time() {
    local end diff h m s
    end=$(date +%s)
    diff=$(( end - _TIMER_START ))
    h=$(( diff / 3600 ))
    m=$(( (diff % 3600) / 60 ))
    s=$(( diff % 60 ))
    log_info "Elapsed: ${h}h ${m}m ${s}s"
}

# ── Checkpoint System ────────────────────────────────────────────────────────
# Each step+sample creates a marker file in $CHECKPOINT_DIR/<step>/<SRR>.done
# If the marker exists, the step is skipped on re-run.
#
# Usage:
#   if checkpoint_exists "03_trimmomatic" "$SRR"; then
#       log_info "Already completed, skipping."
#       exit 0
#   fi
#   ... do work ...
#   checkpoint_set "03_trimmomatic" "$SRR"

checkpoint_exists() {
    local step="$1" sample="${2:-global}"
    [[ -f "${CHECKPOINT_DIR}/${step}/${sample}.done" ]]
}

checkpoint_set() {
    local step="$1" sample="${2:-global}"
    mkdir -p "${CHECKPOINT_DIR}/${step}"
    date '+%Y-%m-%d %H:%M:%S' > "${CHECKPOINT_DIR}/${step}/${sample}.done"
    log_info "Checkpoint set: ${step}/${sample}"
}

# ── Get Accession by SLURM Array Task ID ─────────────────────────────────────
# Usage: SRR=$(get_accession)   — reads $SLURM_ARRAY_TASK_ID from environment
get_accession() {
    local task_id="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID not set}"
    local acc
    acc=$(sed -n "${task_id}p" "${ACCESSION_LIST}")
    [[ -z "$acc" ]] && die "No accession found at line ${task_id} in ${ACCESSION_LIST}"
    echo "$acc"
}

# ── File Validation ──────────────────────────────────────────────────────────
# Usage: validate_files "file1.fastq" "file2.fastq"
# Dies if any file is missing or empty.
validate_files() {
    local f
    for f in "$@"; do
        if [[ ! -f "$f" ]]; then
            die "Expected output file not found: $f"
        elif [[ ! -s "$f" ]]; then
            die "Output file is empty: $f"
        fi
    done
    log_info "Validated ${#} output file(s) ✓"
}

# ── Ensure Directories Exist ─────────────────────────────────────────────────
# Usage: ensure_dirs "$FASTQ_DIR" "$LOG_DIR" "$TRIMMED_DIR/$SRR"
ensure_dirs() {
    local d
    for d in "$@"; do
        mkdir -p "$d"
    done
}

# ── Count Lines in FASTQ ────────────────────────────────────────────────────
# Returns the number of reads (lines / 4) in a FASTQ file.
fastq_read_count() {
    local file="$1"
    local lines
    if [[ "$file" == *.gz ]]; then
        lines=$(zcat "$file" | wc -l)
    else
        lines=$(wc -l < "$file")
    fi
    echo $(( lines / 4 ))
}

# ── Report Read Survival Rate ───────────────────────────────────────────────
# Usage: report_survival_rate <before_count> <after_count> <label>
report_survival_rate() {
    local before="$1" after="$2" label="${3:-reads}"
    if (( before > 0 )); then
        local pct
        pct=$(awk "BEGIN { printf \"%.1f\", ($after / $before) * 100 }")
        log_info "${label}: ${after}/${before} (${pct}%) survived"
    else
        log_warn "Before count is 0, cannot compute survival rate"
    fi
}

# ── SLURM Job Header (for log context) ──────────────────────────────────────
print_job_info() {
    log_info "════════════════════════════════════════════════════════"
    log_info "Job ID:       ${SLURM_JOB_ID:-N/A}"
    log_info "Array Task:   ${SLURM_ARRAY_TASK_ID:-N/A}"
    log_info "Node:         ${SLURM_NODELIST:-N/A}"
    log_info "CPUs:         ${SLURM_CPUS_PER_TASK:-N/A}"
    log_info "Memory:       ${SLURM_MEM_PER_NODE:-N/A}"
    log_info "Working Dir:  $(pwd)"
    log_info "════════════════════════════════════════════════════════"
}

# ── Retry a command with exponential backoff ────────────────────────────────
# Usage: retry 3 prefetch "$SRR" --max-size 100G
retry() {
    local max_attempts="$1"; shift
    local attempt=1
    local wait_time=10
    while (( attempt <= max_attempts )); do
        log_info "Attempt ${attempt}/${max_attempts}: $*"
        if "$@"; then
            return 0
        fi
        log_warn "Attempt ${attempt} failed."
        if (( attempt < max_attempts )); then
            log_info "Retrying in ${wait_time}s..."
            sleep "$wait_time"
            wait_time=$(( wait_time * 2 ))
        fi
        (( attempt++ ))
    done
    die "All ${max_attempts} attempts failed for: $*"
}
