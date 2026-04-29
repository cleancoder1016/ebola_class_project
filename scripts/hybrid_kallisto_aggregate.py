#!/usr/bin/env python3
"""
hybrid_kallisto_aggregate.py

Aggregates Kallisto abundance.tsv files from the hybrid (macaque + Ebola)
transcriptome run. Maps transcript-level counts to gene-level, correlates
against the hybrid featureCounts matrix, and generates a comparative scatter plot.
"""
import os
import sys
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

PROJECT_DIR = "/users/PWSU0516/mufakiransari/ebola_class_project"
HYBRID_KALLISTO_DIR = os.path.join(PROJECT_DIR, "hybrid_kallisto_output")
HYBRID_COUNTS_DIR = os.path.join(PROJECT_DIR, "hybrid_counts")
REPORT_FIG_DIR = os.path.join(PROJECT_DIR, "report", "figures")
os.makedirs(REPORT_FIG_DIR, exist_ok=True)

def strip_prefix(tid):
    """Strip NCBI transcript prefixes to get clean gene identifiers."""
    for prefix in ["rna-", "gene-", "transcript-"]:
        if tid.startswith(prefix):
            return tid[len(prefix):]
    return tid

def main():
    print("Aggregating hybrid Kallisto results...")
    samples = glob.glob(os.path.join(HYBRID_KALLISTO_DIR, "SRR*"))
    all_counts = {}
    all_tpm = {}
    
    for sample_dir in samples:
        srr = os.path.basename(sample_dir)
        tsv_file = os.path.join(sample_dir, "abundance.tsv")
        if not os.path.exists(tsv_file):
            continue
        
        df = pd.read_csv(tsv_file, sep="\t")
        df["gene_name"] = df["target_id"].apply(strip_prefix)
        
        # Group by gene_name and sum (e.g. for multiple transcripts per gene)
        gene_df = df.groupby("gene_name").sum(numeric_only=True)
        
        all_counts[srr] = gene_df["est_counts"]
        all_tpm[srr] = gene_df["tpm"]

    if not all_counts:
        print("No hybrid Kallisto abundance files found.")
        sys.exit(1)

    counts_df = pd.DataFrame(all_counts).fillna(0)
    tpm_df = pd.DataFrame(all_tpm).fillna(0)
    
    counts_df.to_csv(os.path.join(HYBRID_KALLISTO_DIR, "count_matrix_hybrid_kallisto.csv"))
    tpm_df.to_csv(os.path.join(HYBRID_KALLISTO_DIR, "tpm_matrix_hybrid_kallisto.csv"))
    print(f"Aggregated {counts_df.shape[1]} samples. Matrix shape: {counts_df.shape}")
    print(f"  Genes: {counts_df.shape[0]}")
    print(f"  Samples: {counts_df.shape[1]}")

    # Verify hybrid featureCounts matrix exists
    fc_clean = os.path.join(HYBRID_COUNTS_DIR, "gene_counts_clean.txt")
    if not os.path.exists(fc_clean):
        print(f"Warning: hybrid featureCounts matrix not found at {fc_clean}")
        print("Skipping correlation analysis. Kallisto matrices saved successfully.")
        sys.exit(0)

    fc_df = pd.read_csv(fc_clean, sep="\t").set_index("Geneid")
    
    # Identify common samples and genes
    common_cols = list(set(counts_df.columns).intersection(set(fc_df.columns)))
    common_genes = list(set(counts_df.index).intersection(set(fc_df.index)))
    
    print(f"Intersecting {len(common_cols)} samples and {len(common_genes)} genes for correlation.")
    
    if len(common_cols) == 0 or len(common_genes) == 0:
        print("Warning: No overlapping samples/genes between Kallisto and featureCounts.")
        print("This may indicate different gene ID formats. Check GTF annotation consistency.")
        sys.exit(0)
    
    kallisto_subset = counts_df.loc[common_genes, common_cols].astype(float)
    fc_subset = fc_df.loc[common_genes, common_cols].astype(float)
    
    # Flatten matrices for scatter plot
    df_plot = pd.DataFrame({
        "Kallisto_Estimated_Counts": kallisto_subset.values.flatten(),
        "HISAT2_featureCounts": fc_subset.values.flatten()
    })
    
    # Remove zeros for log-log plot to avoid inf
    df_plot = df_plot[(df_plot["Kallisto_Estimated_Counts"] > 0) & (df_plot["HISAT2_featureCounts"] > 0)]
    
    corr = df_plot.corr(method="spearman").iloc[0, 1]
    
    # Plotting
    plt.figure(figsize=(8, 8))
    sns.scatterplot(data=df_plot, x="HISAT2_featureCounts", y="Kallisto_Estimated_Counts", 
                    alpha=0.3, edgecolor=None, color='#2c7fb8', s=5)
    
    # Add x=y reference line
    min_val = df_plot.min().min()
    max_val = df_plot.max().max()
    plt.plot([min_val, max_val], [min_val, max_val], 'r--', label='x = y')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('HISAT2 + featureCounts (log scale)')
    plt.ylabel('Kallisto Estimated Counts (log scale)')
    plt.title(f'Hybrid Genome: Kallisto vs HISAT2 Quantification\nSpearman Rho = {corr:.3f}')
    plt.legend()
    plt.grid(True, which="both", ls="-", alpha=0.2)
    
    plot_path = os.path.join(REPORT_FIG_DIR, "hybrid_kallisto_vs_hisat2_scatter.png")
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    print(f"Saved scatter plot to {plot_path}")

if __name__ == "__main__":
    main()
