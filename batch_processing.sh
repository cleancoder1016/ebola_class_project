#!/bin/bash
#SBATCH --job-name=ebola_batch
#SBATCH --time=16:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --account=PWSU0516
#SBATCH --array=2-356%5 
#SBATCH --output=logs/%x_%A_%a.log

set -euo pipefail

module purge
module load sratoolkit/3.0.2
module load star/2.7.11b
module load fastqc/0.12.1
module load subread/2.0.8 # featurecounts

PROJECT_DIR="${SLURM_SUBMIT_DIR}"

ACCESSION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PROJECT_DIR}/srrAccession.txt" | tr -d '\r[:space:]')

STAR_BIN="STAR"
FASTQC_BIN="fastqc"
FASTP_BIN="${PROJECT_DIR}/bin/fastp"
GENOME_DIR="/fs/scratch/PWSU0516/siva/ebola_class_project/genome/hybrid_index"
HYBRID_GTF="/fs/scratch/PWSU0516/siva/ebola_class_project/genome/hybrid_annotation.gtf"

COUNT_COLUMN=4

STATUS_DIR="${PROJECT_DIR}/count_tsv"
STATUS_FILE="${STATUS_DIR}/sample_status.tsv"

if [ -z "${ACCESSION}" ]; then
    echo "Error: line ${SLURM_ARRAY_TASK_ID} of srrAccession.txt is empty"
    exit 1
fi

mkdir -p \
    "${PROJECT_DIR}/logs" \
    "${PROJECT_DIR}/sra_files" \
    "${PROJECT_DIR}/fastq_outputs" \
    "${PROJECT_DIR}/trimmed_fastq" \
    "${PROJECT_DIR}/fastqc_reports/pre_trim/${ACCESSION}" \
    "${PROJECT_DIR}/fastqc_reports/post_trim/${ACCESSION}" \
    "${PROJECT_DIR}/fastp_reports/${ACCESSION}" \
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

    rm -f "${PROJECT_DIR}/fastq_outputs/${ACCESSION}*.fastq" 2>/dev/null || true
    rm -f "${PROJECT_DIR}/fastq_outputs/${ACCESSION}*.fastq.gz" 2>/dev/null || true
    rm -rf "${PROJECT_DIR}/sra_files/${ACCESSION}" 2>/dev/null || true
    rm -rf "${PROJECT_DIR}/star_alignments/${ACCESSION}/_STARtmp" 2>/dev/null || true
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

R1_RAW="${PROJECT_DIR}/fastq_outputs/${ACCESSION}_1.fastq.gz"

if [ -f "${PROJECT_DIR}/fastq_outputs/${ACCESSION}_3.fastq.gz" ]; then
    R2_RAW="${PROJECT_DIR}/fastq_outputs/${ACCESSION}_3.fastq.gz"
    echo "Detected paired reads as _1 and _3"
elif [ -f "${PROJECT_DIR}/fastq_outputs/${ACCESSION}_2.fastq.gz" ]; then
    R2_RAW="${PROJECT_DIR}/fastq_outputs/${ACCESSION}_2.fastq.gz"
    echo "Detected paired reads as _1 and _2"
else
    R2_RAW="MISSING"
fi

SE_RAW="${PROJECT_DIR}/fastq_outputs/${ACCESSION}.fastq.gz"

# trimmed outputs
R1_TRIM="${PROJECT_DIR}/trimmed_fastq/${ACCESSION}_1_trim.fastq.gz"
R2_TRIM="${PROJECT_DIR}/trimmed_fastq/${ACCESSION}_2_trim.fastq.gz"

COUNT_TSV="${PROJECT_DIR}/count_tsv/${ACCESSION}.counts.tsv"
COUNTS_FILE="${PROJECT_DIR}/star_alignments/${ACCESSION}/${ACCESSION}.ReadsPerGene.out.tab"
BAM_FILE="${PROJECT_DIR}/star_alignments/${ACCESSION}/${ACCESSION}.Aligned.sortedByCoord.out.bam"

MODE="NA"

if [ -f "${R1_RAW}" ] && [ -f "${R2_RAW}" ]; then
    MODE="PAIRED"

    # pre-trim fastqc
    "${FASTQC_BIN}" -t "${SLURM_CPUS_PER_TASK}" \
        -o "${PROJECT_DIR}/fastqc_reports/pre_trim/${ACCESSION}" \
        "${R1_RAW}" "${R2_RAW}"

    # trimming
    "${FASTP_BIN}" --thread "${SLURM_CPUS_PER_TASK}" \
        -i "${R1_RAW}" -I "${R2_RAW}" \
        -o "${R1_TRIM}" -O "${R2_TRIM}" \
        --umi --umi_loc=read1 --umi_len=12 \
        --html "${PROJECT_DIR}/fastp_reports/${ACCESSION}/${ACCESSION}.fastp.html" \
        --json "${PROJECT_DIR}/fastp_reports/${ACCESSION}/${ACCESSION}.fastp.json"

    # post-trim QC
    "${FASTQC_BIN}" -t "${SLURM_CPUS_PER_TASK}" \
        -o "${PROJECT_DIR}/fastqc_reports/post_trim/${ACCESSION}" \
        "${R1_TRIM}" "${R2_TRIM}"


    "${STAR_BIN}" \
        --runThreadN "${SLURM_CPUS_PER_TASK}" \
        --genomeDir "${GENOME_DIR}" \
        --readFilesIn "${R1_TRIM}" "${R2_TRIM}" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${PROJECT_DIR}/star_alignments/${ACCESSION}/${ACCESSION}." \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM NM MD \
        --quantMode GeneCounts \
        --seedSearchStartLmax 1

elif [ -f "${SE_RAW}" ] && [ ! -f "${R1_RAW}" ]; then
    MODE="SINGLE"

    log_status "${ACCESSION}" "SKIPPED" "${MODE}" "Script optimized for Paired-End" "NA"
    exit 0

else
    MODE="UNSUPPORTED"
    MESSAGE="Unsupported FASTQ pattern; expected single file or _1/_3 pair"
    echo "${MESSAGE}"
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
