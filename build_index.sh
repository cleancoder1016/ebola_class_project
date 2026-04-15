#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --time=05:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=48G
#SBATCH --account=PWSU0516

module load STAR/2.7.11b

mkdir -p macaque_index

STAR --runMode genomeGenerate \
     --runThreadN 10 \
     --genomeDir ./macaque_index \
     --genomeFastaFiles Macaca_mulatta.Mmul_10.dna.toplevel.fa \
     --sjdbGTFfile Macaca_mulatta.Mmul_10.111.gtf \
     --sjdbOverhang 99
