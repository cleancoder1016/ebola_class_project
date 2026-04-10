#!/bin/bash
#SBATCH --job-name=sra_array
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=48G
#SBATCH --account=PWSU0516
#SBATCH --array=1-891%10 
#SBATCH --output=logs/sra_%A_%a.log

# load sra
module load sratoolkit/3.0.2

# make the dirs
mkdir -p logs fastq_outputs sra_files

# define sra variable
ACCESSION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" srrAccession.txt)
mkdir -p ./fastq_outputs/$ACCESSION

# prefetch accession  instead of a download , prefetch resumes if connection drops
echo "Ok now Task ID $SLURM_ARRAY_TASK_ID prefetcing sra data from  $ACCESSION"
prefetch $ACCESSION --max-size 100G -O ./sra_files/

# convert to fastq
echo "converting sra to fastq..."
fasterq-dump ./sra_files/$ACCESSION/$ACCESSION.sra --split-files -e 10 -O ./fastq_outputs/$ACCESSION
echo "Results stored at ./fastq_outputs/."
echo "Conversion Completed..."

