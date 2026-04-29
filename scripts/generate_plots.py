#!/usr/bin/env python3
"""
generate_plots.py — Publication-quality figures for Ebola RNA-seq pipeline
═══════════════════════════════════════════════════════════════════════════
Generates:
  1. Per-gene total expression bar chart
  2. Per-sample library size distribution
  3. Kallisto vs featureCounts gene-level correlation
  4. Viral load heatmap (top 50 samples)
  5. Per-gene expression boxplots across samples
  6. Sample viral load ranked bar chart
"""

import os
import sys
import csv
import numpy as np

# Try to import matplotlib; install if missing
try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend for HPC
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from matplotlib.colors import LinearSegmentedColormap
except ImportError:
    print("[ERROR] matplotlib not found. Install with: pip install matplotlib")
    sys.exit(1)

# ── Configuration ──────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)

COUNTS_FILE = os.path.join(PROJECT_DIR, "counts_git", "gene_counts_clean.txt")
KALLISTO_FILE = os.path.join(PROJECT_DIR, "kallisto_git", "count_matrix_kallisto.csv")
OUTPUT_DIR = os.path.join(PROJECT_DIR, "report", "figures")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Ebola gene order (5' → 3' on genome)
EBOLA_GENE_ORDER = ["NP", "VP35", "VP40", "GP", "VP30", "VP24", "L"]
EBOLA_COLORS = {
    "NP": "#E63946", "VP35": "#457B9D", "VP40": "#2A9D8F",
    "GP": "#E9C46A", "VP30": "#F4A261", "VP24": "#264653", "L": "#6A4C93"
}


def load_featurecounts(filepath):
    """Load featureCounts clean matrix (TSV: Geneid, sample1, sample2, ...)"""
    genes = {}
    samples = []
    with open(filepath, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        samples = header[1:]  # First column is Geneid
        for row in reader:
            gene = row[0]
            counts = [float(x) for x in row[1:]]
            genes[gene] = counts
    return genes, samples


def load_kallisto(filepath):
    """Load Kallisto count matrix (CSV: gene_name, sample1, sample2, ...)"""
    genes = {}
    samples = []
    with open(filepath, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        samples = header[1:]
        for row in reader:
            gene = row[0]
            counts = [float(x) for x in row[1:]]
            genes[gene] = counts
    return genes, samples


def plot_gene_expression_bar(fc_genes, output_dir):
    """Plot 1: Total expression per Ebola gene (bar chart)"""
    print("[INFO] Generating per-gene expression bar chart...")

    genes_ordered = [g for g in EBOLA_GENE_ORDER if g in fc_genes]
    totals = [sum(fc_genes[g]) for g in genes_ordered]
    colors = [EBOLA_COLORS.get(g, "#888888") for g in genes_ordered]

    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(genes_ordered, totals, color=colors, edgecolor='white', linewidth=1.2)

    # Add value labels on bars
    for bar, total in zip(bars, totals):
        ax.text(bar.get_x() + bar.get_width()/2., bar.get_height(),
                f'{total:,.0f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

    ax.set_xlabel("Ebola Gene", fontsize=13, fontweight='bold')
    ax.set_ylabel("Total Read Pairs (featureCounts)", fontsize=13, fontweight='bold')
    ax.set_title("Total Expression per Ebola Gene\n(276 samples, featureCounts)", fontsize=15, fontweight='bold')
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x/1e6:.1f}M'))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "gene_expression_total.png"), dpi=200)
    plt.close()
    print("  ✓ gene_expression_total.png")


def plot_sample_library_sizes(fc_genes, samples, output_dir):
    """Plot 2: Per-sample total Ebola reads (ranked bar chart)"""
    print("[INFO] Generating sample library size plot...")

    # Calculate total Ebola reads per sample
    n_samples = len(samples)
    lib_sizes = [0.0] * n_samples
    for gene, counts in fc_genes.items():
        for i, c in enumerate(counts):
            lib_sizes[i] += c

    # Sort by library size
    sorted_pairs = sorted(zip(samples, lib_sizes), key=lambda x: x[1], reverse=True)
    sorted_samples, sorted_sizes = zip(*sorted_pairs)

    # Only plot top 50 for readability
    n_show = min(50, len(sorted_samples))

    fig, ax = plt.subplots(figsize=(14, 7))

    # Color gradient: high = red, low = blue
    max_size = max(sorted_sizes[:n_show])
    colors = [plt.cm.YlOrRd(s / max_size) for s in sorted_sizes[:n_show]]

    ax.bar(range(n_show), sorted_sizes[:n_show], color=colors, edgecolor='none')
    ax.set_xticks(range(n_show))
    ax.set_xticklabels([s.replace("SRR38105", "") for s in sorted_samples[:n_show]],
                       rotation=90, fontsize=7)
    ax.set_xlabel("Sample (SRR38105xxx)", fontsize=12, fontweight='bold')
    ax.set_ylabel("Total Ebola Read Pairs", fontsize=12, fontweight='bold')
    ax.set_title(f"Top {n_show} Samples by Ebola Viral Load\n(featureCounts)", fontsize=14, fontweight='bold')
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x/1e6:.1f}M' if x >= 1e6 else f'{x/1e3:.0f}K'))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "sample_viral_load_ranked.png"), dpi=200)
    plt.close()
    print("  ✓ sample_viral_load_ranked.png")


def plot_gene_boxplots(fc_genes, output_dir):
    """Plot 3: Per-gene expression boxplots across all samples"""
    print("[INFO] Generating per-gene expression boxplots...")

    genes_ordered = [g for g in EBOLA_GENE_ORDER if g in fc_genes]

    # Log2(count + 1) transform
    data = [np.log2(np.array(fc_genes[g]) + 1) for g in genes_ordered]
    colors = [EBOLA_COLORS.get(g, "#888888") for g in genes_ordered]

    fig, ax = plt.subplots(figsize=(10, 6))
    bp = ax.boxplot(data, labels=genes_ordered, patch_artist=True, showfliers=False,
                    medianprops=dict(color='black', linewidth=2))

    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_xlabel("Ebola Gene", fontsize=13, fontweight='bold')
    ax.set_ylabel("log₂(count + 1)", fontsize=13, fontweight='bold')
    ax.set_title("Per-Gene Expression Distribution Across Samples\n(276 samples, featureCounts)", fontsize=14, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "gene_expression_boxplots.png"), dpi=200)
    plt.close()
    print("  ✓ gene_expression_boxplots.png")


def plot_kallisto_vs_featurecounts(fc_genes, fc_samples, kal_genes, kal_samples, output_dir):
    """Plot 4: Kallisto vs featureCounts per-gene total correlation"""
    print("[INFO] Generating Kallisto vs featureCounts correlation...")

    common_genes = [g for g in EBOLA_GENE_ORDER if g in fc_genes and g in kal_genes]
    if len(common_genes) == 0:
        print("  [WARN] No common genes found, skipping.")
        return

    fc_totals = [sum(fc_genes[g]) for g in common_genes]
    kal_totals = [sum(kal_genes[g]) for g in common_genes]
    colors = [EBOLA_COLORS.get(g, "#888888") for g in common_genes]

    fig, ax = plt.subplots(figsize=(8, 8))

    for i, gene in enumerate(common_genes):
        ax.scatter(fc_totals[i], kal_totals[i], c=colors[i], s=120, edgecolors='white',
                   linewidth=1.5, zorder=3)
        ax.annotate(gene, (fc_totals[i], kal_totals[i]),
                    xytext=(8, 8), textcoords='offset points',
                    fontsize=11, fontweight='bold')

    # Add diagonal reference line
    max_val = max(max(fc_totals), max(kal_totals)) * 1.1
    ax.plot([0, max_val], [0, max_val], '--', color='grey', alpha=0.5, linewidth=1)

    # Correlation
    if len(fc_totals) > 1:
        corr = np.corrcoef(fc_totals, kal_totals)[0, 1]
        ax.text(0.05, 0.95, f'r = {corr:.4f}', transform=ax.transAxes,
                fontsize=13, fontweight='bold', verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.set_xlabel("featureCounts (total read pairs)", fontsize=13, fontweight='bold')
    ax.set_ylabel("Kallisto (estimated counts)", fontsize=13, fontweight='bold')
    ax.set_title("featureCounts vs Kallisto\nPer-Gene Total Expression", fontsize=14, fontweight='bold')
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x/1e6:.1f}M'))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x/1e6:.1f}M'))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_aspect('equal', adjustable='box')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "kallisto_vs_featurecounts.png"), dpi=200)
    plt.close()
    print("  ✓ kallisto_vs_featurecounts.png")


def plot_viral_load_heatmap(fc_genes, samples, output_dir):
    """Plot 5: Heatmap of gene expression for top 40 samples by viral load"""
    print("[INFO] Generating viral load heatmap...")

    genes_ordered = [g for g in EBOLA_GENE_ORDER if g in fc_genes]
    n_samples = len(samples)

    # Calculate total viral load per sample
    lib_sizes = [0.0] * n_samples
    for gene in genes_ordered:
        for i, c in enumerate(fc_genes[gene]):
            lib_sizes[i] += c

    # Get top 40 samples
    n_show = min(40, n_samples)
    sorted_indices = sorted(range(n_samples), key=lambda i: lib_sizes[i], reverse=True)[:n_show]

    # Build matrix: genes x samples (log2 + 1)
    matrix = []
    for gene in genes_ordered:
        row = [np.log2(fc_genes[gene][i] + 1) for i in sorted_indices]
        matrix.append(row)
    matrix = np.array(matrix)

    sample_labels = [samples[i].replace("SRR38105", "") for i in sorted_indices]

    fig, ax = plt.subplots(figsize=(16, 5))

    # Custom colormap: dark blue → white → dark red
    cmap = LinearSegmentedColormap.from_list("ebola",
        ["#1a1a2e", "#16213e", "#0f3460", "#53a8b6", "#f5f5dc", "#e9c46a", "#f4a261", "#e76f51", "#9b2226"])

    im = ax.imshow(matrix, aspect='auto', cmap=cmap, interpolation='nearest')

    ax.set_yticks(range(len(genes_ordered)))
    ax.set_yticklabels(genes_ordered, fontsize=11, fontweight='bold')
    ax.set_xticks(range(len(sample_labels)))
    ax.set_xticklabels(sample_labels, rotation=90, fontsize=7)
    ax.set_xlabel("Sample (SRR38105xxx, ranked by total viral load)", fontsize=11, fontweight='bold')
    ax.set_title(f"Ebola Gene Expression Heatmap — Top {n_show} Samples\nlog₂(count + 1)", fontsize=14, fontweight='bold')

    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("log₂(count + 1)", fontsize=11)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "viral_load_heatmap.png"), dpi=200)
    plt.close()
    print("  ✓ viral_load_heatmap.png")


def plot_summary_stats(fc_genes, samples, output_dir):
    """Plot 6: Summary statistics text panel"""
    print("[INFO] Generating summary statistics panel...")

    n_samples = len(samples)
    genes_ordered = [g for g in EBOLA_GENE_ORDER if g in fc_genes]
    total_reads = sum(sum(fc_genes[g]) for g in genes_ordered)
    lib_sizes = [sum(fc_genes[g][i] for g in genes_ordered) for i in range(n_samples)]

    nonzero_samples = sum(1 for s in lib_sizes if s > 100)
    median_lib = sorted(lib_sizes)[n_samples // 2]
    max_lib = max(lib_sizes)
    max_sample = samples[lib_sizes.index(max_lib)]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.axis('off')

    stats_text = f"""
    ╔══════════════════════════════════════════════════════════════╗
    ║           EBOLA RNA-SEQ PIPELINE — RESULTS SUMMARY          ║
    ╠══════════════════════════════════════════════════════════════╣
    ║                                                              ║
    ║   Total Samples Processed:     {n_samples:>6,}                         ║
    ║   Samples with Viral Reads:    {nonzero_samples:>6,}  ({100*nonzero_samples/n_samples:.1f}%)              ║
    ║   Total Ebola Read Pairs:      {total_reads:>12,}                ║
    ║   Median Reads/Sample:         {median_lib:>12,}                ║
    ║   Max Reads/Sample:            {max_lib:>12,}                ║
    ║   Top Sample:                  {max_sample:<20}           ║
    ║   Ebola Genes Quantified:      {len(genes_ordered):>6}                         ║
    ║                                                              ║
    ║   Pipeline:  SRA → FastQC → Trimmomatic → HISAT2 →          ║
    ║              Picard → featureCounts + Kallisto → DESeq2      ║
    ║                                                              ║
    ╚══════════════════════════════════════════════════════════════╝
    """

    ax.text(0.5, 0.5, stats_text, transform=ax.transAxes,
            fontsize=11, fontfamily='monospace',
            verticalalignment='center', horizontalalignment='center',
            bbox=dict(boxstyle='round', facecolor='#f0f0f0', alpha=0.8))

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "pipeline_summary.png"), dpi=200)
    plt.close()
    print("  ✓ pipeline_summary.png")


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("═══════════════════════════════════════════════════════════")
    print("  Ebola Pipeline — Generating Publication Figures")
    print("═══════════════════════════════════════════════════════════\n")

    # Load featureCounts
    if not os.path.exists(COUNTS_FILE):
        print(f"[ERROR] featureCounts file not found: {COUNTS_FILE}")
        sys.exit(1)
    fc_genes, fc_samples = load_featurecounts(COUNTS_FILE)
    print(f"[INFO] featureCounts: {len(fc_genes)} genes × {len(fc_samples)} samples")

    # Load Kallisto
    kal_genes, kal_samples = None, None
    if os.path.exists(KALLISTO_FILE):
        kal_genes, kal_samples = load_kallisto(KALLISTO_FILE)
        print(f"[INFO] Kallisto: {len(kal_genes)} genes × {len(kal_samples)} samples")
    else:
        print(f"[WARN] Kallisto file not found: {KALLISTO_FILE}")

    print(f"\n[INFO] Output directory: {OUTPUT_DIR}\n")

    # Generate all plots
    plot_gene_expression_bar(fc_genes, OUTPUT_DIR)
    plot_sample_library_sizes(fc_genes, fc_samples, OUTPUT_DIR)
    plot_gene_boxplots(fc_genes, OUTPUT_DIR)
    plot_viral_load_heatmap(fc_genes, fc_samples, OUTPUT_DIR)
    plot_summary_stats(fc_genes, fc_samples, OUTPUT_DIR)

    if kal_genes:
        plot_kallisto_vs_featurecounts(fc_genes, fc_samples, kal_genes, kal_samples, OUTPUT_DIR)

    print(f"\n═══════════════════════════════════════════════════════════")
    print(f"  Done! {len(os.listdir(OUTPUT_DIR))} figures saved to:")
    print(f"  {OUTPUT_DIR}")
    print(f"═══════════════════════════════════════════════════════════")
