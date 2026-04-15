#!/bin/bash
#SBATCH --job-name=one_srr_ebola
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=48G
#SBATCH --account=PWSU0516
#SBATCH --output=logs/%x_%j.log

set -euo pipefail

module purge
module load sratoolkit/3.0.2
module load star/2.7.11b
module load fastqc/0.12.1

PROJECT_DIR="${SLURM_SUBMIT_DIR}"
ACCESSION=$(sed -n '1p' "${PROJECT_DIR}/srrAccession.txt" | tr -d '\r[:space:]')

STAR_BIN="STAR"
FASTQC_BIN="fastqc"

GENOME_DIR="/fs/scratch/PWSU0516/siva/ebola_class_project/genome/macaque_index"

# STAR ReadsPerGene.out.tab:
# 2 = unstranded, 3 = forward-stranded, 4 = reverse-stranded
COUNT_COLUMN=4

STATUS_DIR="${PROJECT_DIR}/count_tsv"
STATUS_FILE="${STATUS_DIR}/sample_status.tsv"

if [ -z "${ACCESSION}" ]; then
    echo "Error: first line of srrAccession.txt is empty"
    exit 1
fi

mkdir -p \
    "${PROJECT_DIR}/logs" \
    "${PROJECT_DIR}/sra_files" \
    "${PROJECT_DIR}/fastq_outputs" \
    "${PROJECT_DIR}/fastqc_reports/${ACCESSION}" \
    "${PROJECT_DIR}/star_alignments/${ACCESSION}" \
    "${PROJECT_DIR}/count_tsv" \
    "${PROJECT_DIR}/tmp"

if [ ! -f "${STATUS_FILE}" ]; then
    printf "sample\tstatus\tmode\tmessage\tcount_tsv\n" > "${STATUS_FILE}"
fi

JOB_TMP="${TMPDIR:-${PROJECT_DIR}/tmp/${SLURM_JOB_ID}}"
mkdir -p "${JOB_TMP}"

cleanup() {
    [ -n "${ACCESSION:-}" ] || return 0

    rm -f "${PROJECT_DIR}/fastq_outputs/${ACCESSION}.fastq.gz" 2>/dev/null || true
    rm -f "${PROJECT_DIR}/fastq_outputs/${ACCESSION}_1.fastq.gz" 2>/dev/null || true
    rm -f "${PROJECT_DIR}/fastq_outputs/${ACCESSION}_2.fastq.gz" 2>/dev/null || true
    rm -f "${PROJECT_DIR}/fastq_outputs/${ACCESSION}_3.fastq.gz" 2>/dev/null || true
    rm -rf "${PROJECT_DIR}/sra_files/${ACCESSION}" 2>/dev/null || true
}

trap cleanup EXIT

log_status() {
    local sample="$1"
    local status="$2"
    local mode="$3"
    local message="$4"
    local count_tsv="$5"
    printf "%s\t%s\t%s\t%s\t%s\n" "$sample" "$status" "$mode" "$message" "$count_tsv" >> "${STATUS_FILE}"
}

echo "Processing accession: ${ACCESSION}"

if ! command -v "${STAR_BIN}" >/dev/null 2>&1; then
    log_status "${ACCESSION}" "FAILED" "NA" "STAR not found on PATH" "NA"
    echo "Error: STAR not found"
    exit 1
fi

if ! command -v "${FASTQC_BIN}" >/dev/null 2>&1; then
    log_status "${ACCESSION}" "FAILED" "NA" "FastQC not found on PATH" "NA"
    echo "Error: FastQC not found"
    exit 1
fi

prefetch "${ACCESSION}" --max-size 100G -O "${PROJECT_DIR}/sra_files/"

SRA_PATH="${PROJECT_DIR}/sra_files/${ACCESSION}/${ACCESSION}.sra"
if [ ! -f "${SRA_PATH}" ]; then
    log_status "${ACCESSION}" "FAILED" "NA" "SRA file not found after prefetch" "NA"
    echo "Error: SRA file not found: ${SRA_PATH}"
    exit 1
fi

fasterq-dump "${SRA_PATH}" \
    --split-3 \
    -e "${SLURM_CPUS_PER_TASK}" \
    -O "${PROJECT_DIR}/fastq_outputs/" \
    -t "${JOB_TMP}"

if command -v pigz >/dev/null 2>&1; then
    find "${PROJECT_DIR}/fastq_outputs" -maxdepth 1 -type f -name "${ACCESSION}*.fastq" -print0 | \
        xargs -0 -r -n 1 -P "${SLURM_CPUS_PER_TASK}" pigz
else
    find "${PROJECT_DIR}/fastq_outputs" -maxdepth 1 -type f -name "${ACCESSION}*.fastq" -print0 | \
        xargs -0 -r -n 1 gzip
fi

R1="${PROJECT_DIR}/fastq_outputs/${ACCESSION}_1.fastq.gz"
R2="${PROJECT_DIR}/fastq_outputs/${ACCESSION}_2.fastq.gz"
R3="${PROJECT_DIR}/fastq_outputs/${ACCESSION}_3.fastq.gz"
SE="${PROJECT_DIR}/fastq_outputs/${ACCESSION}.fastq.gz"

COUNT_TSV="${PROJECT_DIR}/count_tsv/${ACCESSION}.counts.tsv"
COUNTS_FILE="${PROJECT_DIR}/star_alignments/${ACCESSION}/${ACCESSION}.ReadsPerGene.out.tab"
BAM_FILE="${PROJECT_DIR}/star_alignments/${ACCESSION}/${ACCESSION}.Aligned.sortedByCoord.out.bam"

MODE="NA"

if [ -f "${R1}" ] && [ -f "${R2}" ]; then
    MODE="PAIRED"

    "${FASTQC_BIN}" -t "${SLURM_CPUS_PER_TASK}" \
        -o "${PROJECT_DIR}/fastqc_reports/${ACCESSION}" \
        "${R1}" "${R2}"

    "${STAR_BIN}" \
        --runThreadN "${SLURM_CPUS_PER_TASK}" \
        --genomeDir "${GENOME_DIR}" \
        --readFilesIn "${R1}" "${R2}" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${PROJECT_DIR}/star_alignments/${ACCESSION}/${ACCESSION}." \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --seedSearchStartLmax 1

elif [ -f "${SE}" ] && [ ! -f "${R1}" ] && [ ! -f "${R2}" ] && [ ! -f "${R3}" ]; then
    MODE="SINGLE"

    "${FASTQC_BIN}" -t "${SLURM_CPUS_PER_TASK}" \
        -o "${PROJECT_DIR}/fastqc_reports/${ACCESSION}" \
        "${SE}"

    "${STAR_BIN}" \
        --runThreadN "${SLURM_CPUS_PER_TASK}" \
        --genomeDir "${GENOME_DIR}" \
        --readFilesIn "${SE}" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${PROJECT_DIR}/star_alignments/${ACCESSION}/${ACCESSION}." \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --seedSearchStartLmax 1

else
    MODE="UNSUPPORTED"
    MESSAGE="Unsupported FASTQ pattern; expected single file or _1/_2 pair"
    echo "${MESSAGE}"
    ls -1 "${PROJECT_DIR}/fastq_outputs/${ACCESSION}"*.fastq.gz 2>/dev/null || true
    log_status "${ACCESSION}" "SKIPPED" "${MODE}" "${MESSAGE}" "NA"
    exit 0
fi

if [ ! -f "${COUNTS_FILE}" ]; then
    log_status "${ACCESSION}" "FAILED" "${MODE}" "STAR count file missing" "NA"
    echo "Error: STAR count file missing"
    exit 1
fi

awk -v col="${COUNT_COLUMN}" 'BEGIN{OFS="\t"} NR>4 {print $1, $col}' "${COUNTS_FILE}" > "${COUNT_TSV}"

if [ ! -s "${COUNT_TSV}" ]; then
    log_status "${ACCESSION}" "FAILED" "${MODE}" "Cleaned count TSV is empty" "NA"
    echo "Error: cleaned count TSV is empty"
    exit 1
fi

if [ -f "${BAM_FILE}" ]; then
    rm -f "${PROJECT_DIR}/fastq_outputs/${ACCESSION}*.fastq.gz"
    rm -rf "${PROJECT_DIR}/sra_files/${ACCESSION}"

    rm -f "${BAM_FILE}"
    log_status "${ACCESSION}" "SUCCESS" "${MODE}" "Completed successfully" "${COUNT_TSV}"
    echo "Done: ${ACCESSION}"
else
    log_status "${ACCESSION}" "FAILED" "${MODE}" "BAM file missing after STAR" "NA"
    echo "Error: BAM file missing"
    exit 1
fi
