# Ebola RNA-seq Pipeline

An end-to-end HPC pipeline for analyzing 356 SRA runs from the 2014 West African Ebola Outbreak (PRJNA938511). Built for the Ohio Supercomputer Center (OSC) Ascend cluster with SLURM array job parallelization, checkpointing, and dual-quantification (HISAT2 + Kallisto).

## Overview

| Property | Value |
|---|---|
| **Dataset** | PRJNA938511 — 2014 West African Ebola Outbreak |
| **Scale** | 356 SRR runs → 344 processed (~0.83 TB raw data) |
| **Reference** | KJ660346.2 (Ebola virus Makona) |
| **HPC** | Ohio Supercomputer Center (Ascend) |
| **Tools** | Trimmomatic, HISAT2, featureCounts, Kallisto, bcftools, DESeq2 |

## Pipeline Results Summary

| Metric | Value |
|---|---|
| Samples downloaded | 344 / 356 (96.6%) |
| featureCounts samples | 276 (paired-end validated) |
| Kallisto samples | 344 |
| Variant calls | 344 samples |
| Ebola genes quantified | 7 (NP, VP35, VP40, GP, VP30, VP24, L) |
| Total Ebola read pairs | 15,131,200 |
| Samples with detectable virus | 177 / 276 (64.1%) |
| Kallisto vs featureCounts (per-sample r) | 0.9978 |
| Generated figures | 13 publication plots |

## Pipeline Architecture

```mermaid
flowchart TD
    subgraph S1["Stage 1 — Data Acquisition"]
        A[srrAccession.txt<br/>356 SRR IDs] --> B[01 SRA → FASTQ<br/>SLURM Array · prefetch + fasterq-dump]
    end

    subgraph S2["Stage 2 — QC & Trimming"]
        B --> C[02 FastQC Raw<br/>SLURM Array]
        C --> D[03 Trimmomatic<br/>SLURM Array · PE adapter + quality trim]
        D --> E[04 FastQC Trimmed<br/>SLURM Array]
    end

    subgraph S3["Stage 3 — HISAT2 Alignment"]
        F[05 HISAT2 Index<br/>Download ref + index] --> G
        E --> G[06 HISAT2 Align<br/>SLURM Array · align + dedup]
        G --> H[07 Post-Align QC<br/>flagstat + depths]
        H --> I[08 featureCounts<br/>276 paired-end samples]
    end

    subgraph S6["Stage 4 — Kallisto Pseudo-alignment"]
        K1[12 Kallisto Index<br/>Build transcriptome k-mers] --> K2
        E --> K2[13 Kallisto Quant<br/>SLURM Array · 344 samples]
        K2 --> K3[14 Kallisto Aggregate<br/>Python → count matrix]
    end

    subgraph S4["Stage 5 — Variant Calling"]
        G -.-> J[09 Variant Calling<br/>SLURM Array · bcftools · 344 samples]
    end

    subgraph S5["Stage 6 — Analysis & Reporting"]
        I --> K[10 DESeq2<br/>PCA + Heatmaps]
        K3 --> K
        J -.-> L[11 MultiQC<br/>Aggregated QC report]
        K -.-> L
    end

    style S1 fill:#1a1a2e,stroke:#e94560,color:#fff
    style S2 fill:#16213e,stroke:#0f3460,color:#fff
    style S3 fill:#1a1a2e,stroke:#533483,color:#fff
    style S4 fill:#350000,stroke:#e94560,color:#fff
    style S6 fill:#0f3460,stroke:#e94560,color:#fff
    style S5 fill:#1a1a2e,stroke:#533483,color:#fff
```

## Quick Start

### 1. Setup
```bash
sbatch scripts/00_setup_conda_env.sh
```

### 2. Run Full Pipeline
The coordinator submits steps sequentially to stay under SLURM's 1,000-job limit:
```bash
# Run all steps
sbatch coordinator.sh

# Resume from a specific step
sbatch coordinator.sh --start-from 08

# Run a subset of steps
sbatch coordinator.sh --start-from 09 --stop-at 14
```

### 3. Monitor
```bash
squeue -u $USER
tail -f logs/coordinator_<JOBID>.log
```

### 4. Generate Figures
```bash
python3 scripts/generate_plots.py     # 6 core figures
python3 scripts/generate_plots_v2.py  # 7 additional figures
```

## Output Directories

| Directory | Contents |
|---|---|
| `counts/` | featureCounts gene count matrix (7 genes × 276 samples) |
| `kallisto_output/` | Kallisto count + TPM matrices (7 genes × 344 samples) |
| `variants/` | Per-sample VCF files, consensus FASTA, variant summary |
| `deseq2_results/` | PCA plot, sample distance heatmap, expression distribution |
| `report/figures/` | 13 publication-quality PNG figures |
| `logs/` | Per-step SLURM logs |
| `.checkpoints/` | Checkpoint files for resume support |

## Configuration

All pipeline parameters (module versions, trimming settings, HISAT2 options, filter thresholds) are defined in `pipeline.config`. The pipeline uses checkpoint files in `.checkpoints/` to skip completed samples on re-runs.

## Contributors

- **Mufakir Ansari** — Kallisto pseudo-alignment pipeline, pipeline orchestration, cross-validation analysis, figure generation
