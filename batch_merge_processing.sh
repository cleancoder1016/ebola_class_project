#!/bin/bash
#SBATCH --job-name=ebola_merge
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --account=PWSU0516
#SBATCH --output=logs/%x_%j.log


set -euo pipefail

module purge
module load sratoolkit/3.0.2
module load star/2.7.11b
module load fastqc/0.12.1
module load subread/2.0.8

PROJECT_DIR="${SLURM_SUBMIT_DIR}"

STAR_BIN="STAR"
FASTQC_BIN="fastqc"
FASTP_BIN="${PROJECT_DIR}/bin/fastp"
GENOME_DIR="/fs/scratch/PWSU0516/siva/ebola_class_project/genome/hybrid_index"
HYBRID_GTF="/fs/scratch/PWSU0516/siva/ebola_class_project/genome/hybrid_annotation.gtf"

COUNT_COLUMN=4

STATUS_DIR="${PROJECT_DIR}/count_tsv"
STATUS_FILE="${STATUS_DIR}/sample_status.tsv"

if [ -z "${SAMPLE_NAME:-}" ]; then
    echo "Error: SAMPLE_NAME not set"
    exit 1
fi
if [ -z "${SRRS:-}" ]; then
    echo "Error: SRRS not set"
    exit 1
fi

read -r -a SRR_LIST <<< "${SRRS}"
echo "Processing merged sample: ${SAMPLE_NAME}"
echo "Lanes (${#SRR_LIST[@]}): ${SRR_LIST[*]}"

FASTQ_OUT="${PROJECT_DIR}/fastq_outputs/${SAMPLE_NAME}"
TRIM_OUT="${PROJECT_DIR}/trimmed_fastq/${SAMPLE_NAME}"

mkdir -p \
    "${PROJECT_DIR}/logs" \
    "${PROJECT_DIR}/sra_files" \
    "${FASTQ_OUT}" \
    "${TRIM_OUT}" \
    "${PROJECT_DIR}/fastqc_reports/pre_trim/${SAMPLE_NAME}" \
    "${PROJECT_DIR}/fastqc_reports/post_trim/${SAMPLE_NAME}" \
    "${PROJECT_DIR}/fastp_reports/${SAMPLE_NAME}" \
    "${PROJECT_DIR}/star_alignments/${SAMPLE_NAME}" \
    "${PROJECT_DIR}/count_tsv" \
    "${PROJECT_DIR}/tmp"

if [ ! -f "${STATUS_FILE}" ]; then
    printf "sample\tstatus\tmode\tmessage\tcount_tsv\n" > "${STATUS_FILE}"
fi

JOB_TMP="${TMPDIR:-${PROJECT_DIR}/tmp/${SLURM_JOB_ID}}"
mkdir -p "${JOB_TMP}"

log_status() {
    printf "%s\t%s\t%s\t%s\t%s\n" "$1" "$2" "$3" "$4" "$5" >> "${STATUS_FILE}"
}

cleanup() {
    rm -rf "${FASTQ_OUT}" "${TRIM_OUT}" "${JOB_TMP}" 2>/dev/null || true
    for srr in "${SRR_LIST[@]}"; do
        rm -rf "${PROJECT_DIR}/sra_files/${srr}" 2>/dev/null || true
    done
    rm -rf "${PROJECT_DIR}/star_alignments/${SAMPLE_NAME}/_STARtmp" 2>/dev/null || true
}
trap cleanup EXIT

if ! command -v "${STAR_BIN}" >/dev/null 2>&1; then
    log_status "${SAMPLE_NAME}" "FAILED" "NA" "STAR not found on PATH" "NA"
    exit 1
fi
if ! command -v "${FASTQC_BIN}" >/dev/null 2>&1; then
    log_status "${SAMPLE_NAME}" "FAILED" "NA" "FastQC not found on PATH" "NA"
    exit 1
fi
if [ ! -x "${FASTP_BIN}" ]; then
    log_status "${SAMPLE_NAME}" "FAILED" "NA" "fastp not found at ${FASTP_BIN}" "NA"
    exit 1
fi

echo "Fetching ${#SRR_LIST[@]} lanes"
for SRR in "${SRR_LIST[@]}"; do
    echo "  prefetch ${SRR}"
    prefetch "${SRR}" --max-size 100G -O "${PROJECT_DIR}/sra_files/"

    SRA_PATH="${PROJECT_DIR}/sra_files/${SRR}/${SRR}.sra"
    if [ ! -f "${SRA_PATH}" ]; then
        log_status "${SAMPLE_NAME}" "FAILED" "NA" "SRA file missing for ${SRR}" "NA"
        echo "Error: ${SRA_PATH} not found"
        exit 1
    fi

    echo "  fasterq-dump ${SRR}"
    fasterq-dump "${SRA_PATH}" \
        --split-3 \
        -e "${SLURM_CPUS_PER_TASK}" \
        -O "${FASTQ_OUT}/" \
        -t "${JOB_TMP}"

    # Compress immediately to save disk
    if command -v pigz >/dev/null 2>&1; then
        find "${FASTQ_OUT}" -maxdepth 1 -type f -name "${SRR}*.fastq" -print0 | \
            xargs -0 -r -n 1 -P "${SLURM_CPUS_PER_TASK}" pigz
    else
        find "${FASTQ_OUT}" -maxdepth 1 -type f -name "${SRR}*.fastq" -print0 | \
            xargs -0 -r -n 1 gzip
    fi

    rm -rf "${PROJECT_DIR}/sra_files/${SRR}"
done

echo "Detecting layout and concatenating lanes"
FIRST="${SRR_LIST[0]}"

R1_MERGED="${FASTQ_OUT}/${SAMPLE_NAME}_1.fastq.gz"
R2_MERGED="${FASTQ_OUT}/${SAMPLE_NAME}_2.fastq.gz"

if [ -f "${FASTQ_OUT}/${FIRST}_3.fastq.gz" ]; then
    PAIR_SUFFIX="3"
    echo "  Detected paired pattern: _1 and _3"
elif [ -f "${FASTQ_OUT}/${FIRST}_2.fastq.gz" ]; then
    PAIR_SUFFIX="2"
    echo "  Detected paired pattern: _1 and _2"
else
    log_status "${SAMPLE_NAME}" "SKIPPED" "SINGLE" "Script optimized for Paired-End" "NA"
    echo "Single-end detected — skipping (not supported)"
    exit 0
fi

R1_LANES=()
R2_LANES=()
for SRR in "${SRR_LIST[@]}"; do
    f1="${FASTQ_OUT}/${SRR}_1.fastq.gz"
    f2="${FASTQ_OUT}/${SRR}_${PAIR_SUFFIX}.fastq.gz"
    if [ ! -f "$f1" ] || [ ! -f "$f2" ]; then
        log_status "${SAMPLE_NAME}" "FAILED" "PAIRED" "Missing FASTQ for lane ${SRR}" "NA"
        echo "Error: expected ${f1} and ${f2}"
        exit 1
    fi
    R1_LANES+=("$f1")
    R2_LANES+=("$f2")
done

echo "  Concatenating R1 lanes -> ${R1_MERGED}"
cat "${R1_LANES[@]}" > "${R1_MERGED}"
echo "  Concatenating R2 lanes -> ${R2_MERGED}"
cat "${R2_LANES[@]}" > "${R2_MERGED}"

rm -f "${R1_LANES[@]}" "${R2_LANES[@]}"

echo "Pre-trim FastQC"
"${FASTQC_BIN}" -t "${SLURM_CPUS_PER_TASK}" \
    -o "${PROJECT_DIR}/fastqc_reports/pre_trim/${SAMPLE_NAME}" \
    "${R1_MERGED}" "${R2_MERGED}"

echo "fastp trimming"
R1_TRIM="${TRIM_OUT}/${SAMPLE_NAME}_1_trim.fastq.gz"
R2_TRIM="${TRIM_OUT}/${SAMPLE_NAME}_2_trim.fastq.gz"

"${FASTP_BIN}" --thread "${SLURM_CPUS_PER_TASK}" \
    -i "${R1_MERGED}" -I "${R2_MERGED}" \
    -o "${R1_TRIM}" -O "${R2_TRIM}" \
    --umi --umi_loc=read1 --umi_len=12 \
    --html "${PROJECT_DIR}/fastp_reports/${SAMPLE_NAME}/${SAMPLE_NAME}.fastp.html" \
    --json "${PROJECT_DIR}/fastp_reports/${SAMPLE_NAME}/${SAMPLE_NAME}.fastp.json"

rm -f "${R1_MERGED}" "${R2_MERGED}"

echo "Post-trim FastQC"
"${FASTQC_BIN}" -t "${SLURM_CPUS_PER_TASK}" \
    -o "${PROJECT_DIR}/fastqc_reports/post_trim/${SAMPLE_NAME}" \
    "${R1_TRIM}" "${R2_TRIM}"

echo "STAR alignment"
"${STAR_BIN}" \
    --runThreadN "${SLURM_CPUS_PER_TASK}" \
    --genomeDir "${GENOME_DIR}" \
    --readFilesIn "${R1_TRIM}" "${R2_TRIM}" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${PROJECT_DIR}/star_alignments/${SAMPLE_NAME}/${SAMPLE_NAME}." \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS nM NM MD \
    --quantMode GeneCounts \
    --seedSearchStartLmax 1

echo "Extracting counts"
COUNTS_FILE="${PROJECT_DIR}/star_alignments/${SAMPLE_NAME}/${SAMPLE_NAME}.ReadsPerGene.out.tab"
BAM_FILE="${PROJECT_DIR}/star_alignments/${SAMPLE_NAME}/${SAMPLE_NAME}.Aligned.sortedByCoord.out.bam"
COUNT_TSV="${PROJECT_DIR}/count_tsv/${SAMPLE_NAME}.counts.tsv"

if [ ! -f "${COUNTS_FILE}" ]; then
    log_status "${SAMPLE_NAME}" "FAILED" "PAIRED" "STAR count file missing" "NA"
    echo "Error: STAR count file missing: ${COUNTS_FILE}"
    exit 1
fi

awk -v col="${COUNT_COLUMN}" 'BEGIN{OFS="\t"} NR>4 {print $1, $col}' \
    "${COUNTS_FILE}" > "${COUNT_TSV}"

if [ ! -s "${COUNT_TSV}" ]; then
    log_status "${SAMPLE_NAME}" "FAILED" "PAIRED" "Count TSV is empty" "NA"
    echo "Error: count TSV is empty"
    exit 1
fi

if [ -f "${BAM_FILE}" ]; then
    rm -f "${BAM_FILE}"
    log_status "${SAMPLE_NAME}" "SUCCESS" "PAIRED" \
        "Completed (${#SRR_LIST[@]} lanes merged)" "${COUNT_TSV}"
    echo "Done: ${SAMPLE_NAME} -> ${COUNT_TSV}"
else
    log_status "${SAMPLE_NAME}" "FAILED" "PAIRED" "BAM file missing after STAR" "NA"
    echo "Error: BAM file missing"
    exit 1
fi
