#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════════════
# run_and_verify.sh — Master script to run pipeline, verify, and push results
# ═══════════════════════════════════════════════════════════════════════════════
# Usage (on OSC):
#   bash run_and_verify.sh           # Full run: execute + verify + push
#   bash run_and_verify.sh --verify  # Only verify existing results
#   bash run_and_verify.sh --push    # Only push results to GitHub
# ═══════════════════════════════════════════════════════════════════════════════

set -euo pipefail

PROJECT_DIR="/users/PWSU0516/mufakiransari/ebola_class_project"
BRANCH="mufakir/fix-pipeline-356-hybrid"
EXPECTED_SAMPLES=356
EXPECTED_EBOLA_GENES=7

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[0;33m'
CYAN='\033[0;36m'; BOLD='\033[1m'; NC='\033[0m'

cd "${PROJECT_DIR}"

# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 1: PRE-FLIGHT CHECKS
# ═══════════════════════════════════════════════════════════════════════════════
preflight() {
    echo -e "${BOLD}${CYAN}"
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║        PRE-FLIGHT CHECKS                                   ║"
    echo "╚══════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"

    FAIL=0

    # Check branch
    CURRENT_BRANCH=$(git branch --show-current)
    if [[ "${CURRENT_BRANCH}" == "${BRANCH}" ]]; then
        echo -e "  ${GREEN}✓${NC} On branch: ${BRANCH}"
    else
        echo -e "  ${RED}✗${NC} Wrong branch: ${CURRENT_BRANCH} (expected ${BRANCH})"
        echo -e "  ${YELLOW}  Fix: git checkout ${BRANCH}${NC}"
        FAIL=1
    fi

    # Check accession list
    if [[ -f srrAccession.txt ]]; then
        COUNT=$(wc -l < srrAccession.txt)
        echo -e "  ${GREEN}✓${NC} Accession list: ${COUNT} samples"
    else
        echo -e "  ${RED}✗${NC} srrAccession.txt not found"
        FAIL=1
    fi

    # Check all scripts exist
    MISSING=0
    for s in scripts/{00,01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19}_*.sh; do
        if [[ ! -f "$s" ]]; then
            echo -e "  ${RED}✗${NC} Missing: $s"
            MISSING=1
        fi
    done
    if (( MISSING == 0 )); then
        echo -e "  ${GREEN}✓${NC} All 20 pipeline scripts found"
    else
        FAIL=1
    fi

    # Check Python scripts
    for py in scripts/kallisto_aggregate.py scripts/hybrid_kallisto_aggregate.py; do
        if [[ -f "$py" ]]; then
            echo -e "  ${GREEN}✓${NC} Found: $(basename $py)"
        else
            echo -e "  ${RED}✗${NC} Missing: $py"
            FAIL=1
        fi
    done

    # Check pipeline.config
    NS=$(grep "^NUM_SAMPLES=" pipeline.config | cut -d= -f2)
    if [[ "$NS" == "356" ]]; then
        echo -e "  ${GREEN}✓${NC} NUM_SAMPLES=${NS}"
    else
        echo -e "  ${RED}✗${NC} NUM_SAMPLES=${NS} (expected 356)"
        FAIL=1
    fi

    # Check disk space
    FREE_GB=$(df --output=avail "${PROJECT_DIR}" | tail -1 | awk '{print int($1/1048576)}')
    if (( FREE_GB < 100 )); then
        echo -e "  ${YELLOW}⚠${NC}  Low disk space: ${FREE_GB} GB free (recommend 800+ GB)"
    else
        echo -e "  ${GREEN}✓${NC} Disk space: ${FREE_GB} GB available"
    fi

    if (( FAIL > 0 )); then
        echo -e "\n  ${RED}${BOLD}Pre-flight checks FAILED. Fix issues above before running.${NC}"
        exit 1
    fi

    echo -e "\n  ${GREEN}${BOLD}All pre-flight checks passed!${NC}\n"
}

# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 2: RUN PIPELINE
# ═══════════════════════════════════════════════════════════════════════════════
run_pipeline() {
    echo -e "${BOLD}${CYAN}"
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║        RUNNING PIPELINE (20 steps, 356 samples)            ║"
    echo "╚══════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"

    # Clean old checkpoints
    echo -e "  ${YELLOW}Clearing old checkpoints...${NC}"
    rm -rf .checkpoints/
    mkdir -p logs

    # Submit pipeline
    echo -e "  ${CYAN}Submitting all 20 pipeline steps...${NC}\n"
    bash scripts/run_pipeline.sh

    echo -e "\n${BOLD}${GREEN}Pipeline submitted!${NC}"
    echo -e "${BOLD}Now waiting for jobs to complete...${NC}"
    echo ""
    echo -e "You can monitor progress in another terminal:"
    echo -e "  ${CYAN}bash scripts/run_pipeline.sh --status${NC}"
    echo -e "  ${CYAN}squeue -u \$USER${NC}"
    echo -e "  ${CYAN}bash run_and_verify.sh --monitor${NC}"
    echo ""
}

# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 2b: MONITOR (wait for completion)
# ═══════════════════════════════════════════════════════════════════════════════
monitor() {
    echo -e "${BOLD}${CYAN}"
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║        MONITORING PIPELINE PROGRESS                        ║"
    echo "╚══════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"

    while true; do
        RUNNING=$(squeue -u "$USER" -h 2>/dev/null | wc -l)
        SRA_DONE=$(ls sra_files/ 2>/dev/null | wc -l)
        FASTQ_DONE=$(ls fastq_outputs/ 2>/dev/null | wc -l)
        ALIGNED_DONE=$(ls aligned/ 2>/dev/null | wc -l)
        KALLISTO_DONE=$(ls kallisto_output/SRR*/abundance.tsv 2>/dev/null | wc -l)
        HYBRID_ALIGNED=$(ls hybrid_aligned/ 2>/dev/null | wc -l)
        HYBRID_KALL=$(ls hybrid_kallisto_output/SRR*/abundance.tsv 2>/dev/null | wc -l)

        echo -e "  [$(date '+%H:%M:%S')] Jobs running: ${BOLD}${RUNNING}${NC} | SRA: ${SRA_DONE}/356 | FASTQ: ${FASTQ_DONE}/356 | Aligned: ${ALIGNED_DONE}/356 | Kallisto: ${KALLISTO_DONE}/356 | HybridAlign: ${HYBRID_ALIGNED}/356 | HybridKall: ${HYBRID_KALL}/356"

        if (( RUNNING == 0 )); then
            echo -e "\n  ${GREEN}${BOLD}All SLURM jobs finished!${NC}"
            break
        fi
        sleep 60
    done
}

# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 3: VERIFY RESULTS
# ═══════════════════════════════════════════════════════════════════════════════
verify() {
    echo -e "${BOLD}${CYAN}"
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║        VERIFYING PIPELINE RESULTS                          ║"
    echo "╚══════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"

    PASS=0
    FAIL=0
    WARN=0

    check_pass() { echo -e "  ${GREEN}✓${NC} $1"; (( PASS++ )); }
    check_fail() { echo -e "  ${RED}✗${NC} $1"; (( FAIL++ )); }
    check_warn() { echo -e "  ${YELLOW}⚠${NC}  $1"; (( WARN++ )); }

    # ── 1. SRA Downloads ───────────────────────────────────────────────
    echo -e "\n${BOLD}  [1/8] SRA Downloads${NC}"
    SRA_COUNT=$(ls -d sra_files/SRR* 2>/dev/null | wc -l)
    if (( SRA_COUNT >= 300 )); then
        check_pass "SRA files: ${SRA_COUNT}/${EXPECTED_SAMPLES}"
    elif (( SRA_COUNT > 0 )); then
        check_warn "SRA files: ${SRA_COUNT}/${EXPECTED_SAMPLES} (some failed)"
    else
        check_fail "SRA files: 0 found"
    fi

    # ── 2. FASTQ Files ─────────────────────────────────────────────────
    echo -e "\n${BOLD}  [2/8] FASTQ Conversion${NC}"
    FASTQ_COUNT=$(ls -d fastq_outputs/SRR* 2>/dev/null | wc -l)
    if (( FASTQ_COUNT >= 300 )); then
        check_pass "FASTQ directories: ${FASTQ_COUNT}/${EXPECTED_SAMPLES}"
    elif (( FASTQ_COUNT > 0 )); then
        check_warn "FASTQ directories: ${FASTQ_COUNT}/${EXPECTED_SAMPLES}"
    else
        check_fail "FASTQ directories: 0 found"
    fi

    # ── 3. Ebola-only Alignment ────────────────────────────────────────
    echo -e "\n${BOLD}  [3/8] Ebola-Only Alignment (HISAT2)${NC}"
    ALIGNED_COUNT=$(ls -d aligned/SRR* 2>/dev/null | wc -l)
    if (( ALIGNED_COUNT >= 300 )); then
        check_pass "Aligned samples: ${ALIGNED_COUNT}/${EXPECTED_SAMPLES}"
    elif (( ALIGNED_COUNT > 0 )); then
        check_warn "Aligned samples: ${ALIGNED_COUNT}/${EXPECTED_SAMPLES}"
    else
        check_fail "Aligned samples: 0 found"
    fi

    # ── 4. Ebola-only featureCounts ────────────────────────────────────
    echo -e "\n${BOLD}  [4/8] Ebola-Only Count Matrix (featureCounts)${NC}"
    if [[ -f counts/gene_counts_clean.txt ]]; then
        GENE_ROWS=$(tail -n +2 counts/gene_counts_clean.txt | wc -l)
        SAMPLE_COLS=$(head -1 counts/gene_counts_clean.txt | tr '\t' '\n' | tail -n +2 | wc -l)
        check_pass "Count matrix exists: ${GENE_ROWS} genes × ${SAMPLE_COLS} samples"
        if (( GENE_ROWS == EXPECTED_EBOLA_GENES )); then
            check_pass "Gene count matches expected (${EXPECTED_EBOLA_GENES} Ebola genes)"
        else
            check_warn "Gene count: ${GENE_ROWS} (expected ${EXPECTED_EBOLA_GENES})"
        fi
        if (( SAMPLE_COLS >= 300 )); then
            check_pass "Sample count: ${SAMPLE_COLS}"
        else
            check_warn "Sample count: ${SAMPLE_COLS} (expected ~${EXPECTED_SAMPLES})"
        fi
    else
        check_fail "counts/gene_counts_clean.txt not found"
    fi

    # ── 5. Ebola-only Kallisto ─────────────────────────────────────────
    echo -e "\n${BOLD}  [5/8] Ebola-Only Kallisto${NC}"
    KALL_DONE=$(ls kallisto_output/SRR*/abundance.tsv 2>/dev/null | wc -l)
    if (( KALL_DONE >= 300 )); then
        check_pass "Kallisto quantified: ${KALL_DONE}/${EXPECTED_SAMPLES}"
    elif (( KALL_DONE > 0 )); then
        check_warn "Kallisto quantified: ${KALL_DONE}/${EXPECTED_SAMPLES}"
    else
        check_fail "Kallisto: 0 samples quantified"
    fi

    if [[ -f kallisto_output/count_matrix_kallisto.csv ]]; then
        KALL_GENES=$(tail -n +2 kallisto_output/count_matrix_kallisto.csv | wc -l)
        KALL_SAMPLES=$(head -1 kallisto_output/count_matrix_kallisto.csv | tr ',' '\n' | tail -n +2 | wc -l)
        check_pass "Kallisto matrix: ${KALL_GENES} genes × ${KALL_SAMPLES} samples"
    else
        check_fail "kallisto_output/count_matrix_kallisto.csv not found"
    fi

    if [[ -f report/figures/kallisto_vs_star_scatter.png ]]; then
        check_pass "Kallisto vs HISAT2 scatter plot generated"
    else
        check_warn "Scatter plot not found (report/figures/kallisto_vs_star_scatter.png)"
    fi

    # ── 6. Hybrid Genome Build ─────────────────────────────────────────
    echo -e "\n${BOLD}  [6/8] Hybrid Genome (Macaque + Ebola)${NC}"
    if [[ -f reference/hybrid/hybrid_macaque_ebola.fasta ]]; then
        HYBRID_SIZE=$(du -sh reference/hybrid/hybrid_macaque_ebola.fasta | cut -f1)
        check_pass "Hybrid FASTA exists: ${HYBRID_SIZE}"
    else
        check_fail "Hybrid FASTA not found"
    fi

    HYBRID_IDX=$(ls reference/hybrid/hisat2_index/hybrid*.ht2 2>/dev/null | wc -l)
    if (( HYBRID_IDX > 0 )); then
        check_pass "HISAT2 hybrid index: ${HYBRID_IDX} files"
    else
        check_fail "HISAT2 hybrid index not built"
    fi

    if [[ -f reference/hybrid/hybrid_transcriptome.idx ]]; then
        check_pass "Kallisto hybrid index exists"
    else
        check_fail "Kallisto hybrid index not built"
    fi

    # ── 7. Hybrid Alignment + Counts ───────────────────────────────────
    echo -e "\n${BOLD}  [7/8] Hybrid Alignment & Counts${NC}"
    HYBRID_ALIGNED=$(ls -d hybrid_aligned/SRR* 2>/dev/null | wc -l)
    if (( HYBRID_ALIGNED >= 300 )); then
        check_pass "Hybrid aligned: ${HYBRID_ALIGNED}/${EXPECTED_SAMPLES}"
    elif (( HYBRID_ALIGNED > 0 )); then
        check_warn "Hybrid aligned: ${HYBRID_ALIGNED}/${EXPECTED_SAMPLES}"
    else
        check_fail "Hybrid aligned: 0 found"
    fi

    if [[ -f hybrid_counts/gene_counts_clean.txt ]]; then
        HYB_GENES=$(tail -n +2 hybrid_counts/gene_counts_clean.txt | wc -l)
        HYB_SAMPLES=$(head -1 hybrid_counts/gene_counts_clean.txt | tr '\t' '\n' | tail -n +2 | wc -l)
        check_pass "Hybrid count matrix: ${HYB_GENES} genes × ${HYB_SAMPLES} samples"
        if (( HYB_GENES > 1000 )); then
            check_pass "Hybrid gene count > 1000 (includes macaque genes)"
        else
            check_warn "Hybrid gene count only ${HYB_GENES} (expected ~30,000+)"
        fi
    else
        check_fail "hybrid_counts/gene_counts_clean.txt not found"
    fi

    # ── 8. Hybrid Kallisto ─────────────────────────────────────────────
    echo -e "\n${BOLD}  [8/8] Hybrid Kallisto${NC}"
    HYB_KALL=$(ls hybrid_kallisto_output/SRR*/abundance.tsv 2>/dev/null | wc -l)
    if (( HYB_KALL >= 300 )); then
        check_pass "Hybrid Kallisto quantified: ${HYB_KALL}/${EXPECTED_SAMPLES}"
    elif (( HYB_KALL > 0 )); then
        check_warn "Hybrid Kallisto quantified: ${HYB_KALL}/${EXPECTED_SAMPLES}"
    else
        check_fail "Hybrid Kallisto: 0 samples"
    fi

    if [[ -f hybrid_kallisto_output/count_matrix_hybrid_kallisto.csv ]]; then
        HK_GENES=$(tail -n +2 hybrid_kallisto_output/count_matrix_hybrid_kallisto.csv | wc -l)
        HK_SAMPLES=$(head -1 hybrid_kallisto_output/count_matrix_hybrid_kallisto.csv | tr ',' '\n' | tail -n +2 | wc -l)
        check_pass "Hybrid Kallisto matrix: ${HK_GENES} genes × ${HK_SAMPLES} samples"
    else
        check_fail "hybrid_kallisto_output/count_matrix_hybrid_kallisto.csv not found"
    fi

    if [[ -f report/figures/hybrid_kallisto_vs_hisat2_scatter.png ]]; then
        check_pass "Hybrid scatter plot generated"
    else
        check_warn "Hybrid scatter plot not found"
    fi

    # ── Summary ────────────────────────────────────────────────────────
    echo -e "\n${BOLD}${CYAN}═══════════════════════════════════════════════════════════${NC}"
    echo -e "  ${GREEN}✓ Passed: ${PASS}${NC}  ${RED}✗ Failed: ${FAIL}${NC}  ${YELLOW}⚠ Warnings: ${WARN}${NC}"
    echo -e "${BOLD}${CYAN}═══════════════════════════════════════════════════════════${NC}\n"

    if (( FAIL > 0 )); then
        echo -e "${RED}${BOLD}VERIFICATION FAILED — ${FAIL} check(s) failed.${NC}"
        echo -e "Check logs/ directory for errors:"
        echo -e "  ${CYAN}ls -lt logs/ | head -20${NC}"
        echo -e "  ${CYAN}grep -l 'ERROR\\|FATAL\\|die' logs/*.log${NC}"
        return 1
    elif (( WARN > 0 )); then
        echo -e "${YELLOW}${BOLD}VERIFICATION PASSED WITH WARNINGS — review above.${NC}"
        return 0
    else
        echo -e "${GREEN}${BOLD}ALL CHECKS PASSED!${NC}"
        return 0
    fi
}

# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 4: PUSH RESULTS TO GITHUB
# ═══════════════════════════════════════════════════════════════════════════════
push_results() {
    echo -e "${BOLD}${CYAN}"
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║        PUSHING RESULTS TO GITHUB                           ║"
    echo "╚══════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"

    cd "${PROJECT_DIR}"

    # Make sure .gitignore is in place (data dirs excluded)
    if [[ ! -f .gitignore ]]; then
        echo -e "  ${RED}✗${NC} .gitignore missing! Data files could be pushed."
        exit 1
    fi

    # Stage only count matrices, figures, and summary files (NOT raw data)
    echo -e "  ${CYAN}Staging result files...${NC}"

    # Count matrices (small CSV/TSV files)
    git add -f counts/gene_counts_clean.txt 2>/dev/null || true
    git add -f counts/gene_counts.txt.summary 2>/dev/null || true
    git add -f kallisto_output/count_matrix_kallisto.csv 2>/dev/null || true
    git add -f kallisto_output/tpm_matrix_kallisto.csv 2>/dev/null || true
    git add -f hybrid_counts/gene_counts_clean.txt 2>/dev/null || true
    git add -f hybrid_counts/gene_counts.txt.summary 2>/dev/null || true
    git add -f hybrid_kallisto_output/count_matrix_hybrid_kallisto.csv 2>/dev/null || true
    git add -f hybrid_kallisto_output/tpm_matrix_hybrid_kallisto.csv 2>/dev/null || true

    # Figures
    git add -f report/figures/*.png 2>/dev/null || true

    # Variant summary
    git add -f variants/variant_summary.tsv 2>/dev/null || true

    # DESeq2 results (PDFs and CSVs only)
    git add -f deseq2_results/*.pdf 2>/dev/null || true
    git add -f deseq2_results/*.csv 2>/dev/null || true

    # Featurecounts logs
    git add -f counts/featurecounts.log 2>/dev/null || true
    git add -f hybrid_counts/featurecounts.log 2>/dev/null || true

    echo ""
    echo -e "  ${CYAN}Files staged:${NC}"
    git diff --cached --stat

    # Commit
    TIMESTAMP=$(date '+%Y-%m-%d %H:%M')
    echo ""
    echo -e "  ${CYAN}Committing...${NC}"
    git commit -m "results: pipeline output — 356 samples, Ebola-only + hybrid genome

Generated on OSC at ${TIMESTAMP}

Ebola-only results:
- featureCounts matrix: 7 Ebola genes × 356 samples
- Kallisto matrix: 7 Ebola genes × 356 samples
- Kallisto vs HISAT2 scatter plot

Hybrid (macaque + Ebola) results:
- featureCounts matrix: ~30,000+ genes × 356 samples
- Kallisto matrix: ~30,000+ genes × 356 samples
- Hybrid scatter plot

Also includes: variant summary, DESeq2 exploratory plots" || {
        echo -e "  ${YELLOW}Nothing to commit (results may already be pushed).${NC}"
        return 0
    }

    # Push
    echo ""
    echo -e "  ${CYAN}Pushing to origin/${BRANCH}...${NC}"
    git push origin "${BRANCH}"

    echo -e "\n  ${GREEN}${BOLD}Results pushed to GitHub!${NC}"
    echo -e "  PR link: https://github.com/cleancoder1016/ebola_class_project/pull/new/${BRANCH}"
}

# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════
case "${1:-full}" in
    --verify)
        verify
        ;;
    --push)
        push_results
        ;;
    --monitor)
        monitor
        ;;
    --preflight)
        preflight
        ;;
    full|"")
        preflight
        run_pipeline
        echo ""
        echo -e "${BOLD}${YELLOW}═══════════════════════════════════════════════════════════${NC}"
        echo -e "${BOLD}${YELLOW}  Pipeline submitted! Jobs are running on SLURM.${NC}"
        echo -e "${BOLD}${YELLOW}  This will take ~24-36 hours to complete.${NC}"
        echo -e "${BOLD}${YELLOW}═══════════════════════════════════════════════════════════${NC}"
        echo ""
        echo -e "When all jobs finish, run:"
        echo -e "  ${CYAN}bash run_and_verify.sh --verify${NC}    # Check results"
        echo -e "  ${CYAN}bash run_and_verify.sh --push${NC}      # Push to GitHub"
        echo ""
        echo -e "Or monitor now:"
        echo -e "  ${CYAN}bash run_and_verify.sh --monitor${NC}   # Watch until done"
        ;;
    *)
        echo "Usage: bash run_and_verify.sh [--verify|--push|--monitor|--preflight]"
        echo ""
        echo "  (no args)    Full run: preflight → submit pipeline"
        echo "  --verify     Verify all pipeline results"
        echo "  --push       Push result files to GitHub"
        echo "  --monitor    Watch SLURM queue until all jobs finish"
        echo "  --preflight  Only run pre-flight checks"
        exit 0
        ;;
esac
