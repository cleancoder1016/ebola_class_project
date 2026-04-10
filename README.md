# Ebola RNA-seq Pipeline
An HPC pipeline for analyzing 50 SRR runs from the 2014 outbreak (PRJNA257197). Automates SRA extraction, QC, alignment, and variant calling. Built for graduate capstones as part of the coursework. In this, we are using HPC clusters from Ohio Super Computer, the clusters we used for this project are cardinal and ascend clusters.


## Overview
**Dataset:** PRJNA257197 - 2014 West African Ebola Outbreak
**Scale:** For simplicity, collected only 50 SRR runs (441 MB SRA / ~ 2.1 GB FASTQ)
**Workflow:**
```mermaid
flowchart TD
    A[srrAccession.txt\n50 SRR IDs] --> B[SLURM Array Job\n1-50 % 10 concurrent]
    B --> C[prefetch\nusing SRA Toolkit 3.0.2]
    C --> D{Download\nSuccessful?}
    D -- No --> C
    D -- Yes --> E[fasterq-dump\n--split-files -e 10]
    E --> F[Paired-end FASTQs\nR1 / R2 per accession]
    F --> G[FastQC\nRaw QC Report]
    G --> H[Trimmomatic\nAdapter & Quality Trimming]
    H --> I[FastQC\nPost-trim QC Report]
    I --> J[STAR / HISAT2\nAlignment to Reference Genome]
    J --> K[BAM Files\nSorted & Indexed]
    K --> L[featureCounts\nRead Quantification]
    L --> M[Count Matrix\nGenes × Samples]
    style A fill:#4a90d9,color:#fff
    style O fill:#27ae60,color:#fff
    style D fill:#e67e22,color:#fff
```
