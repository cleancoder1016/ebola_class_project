#!/usr/bin/env python3
"""
generate_plots_v2.py — Additional publication figures for Ebola RNA-seq
Generates plots from featureCounts summary, variant calls, and Kallisto data.
"""
import os, sys, csv, glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(PROJECT, "report", "figures")
os.makedirs(OUT, exist_ok=True)

GENE_ORDER = ["NP","VP35","VP40","GP","VP30","VP24","L"]
COLORS = {"NP":"#E63946","VP35":"#457B9D","VP40":"#2A9D8F",
           "GP":"#E9C46A","VP30":"#F4A261","VP24":"#264653","L":"#6A4C93"}

def load_tsv(path, skip_comment=False):
    rows = []
    with open(path) as f:
        for line in f:
            if skip_comment and line.startswith('#'): continue
            rows.append(line.strip().split('\t'))
    return rows

def load_counts(path):
    data = load_tsv(path)
    header = data[0]; samples = header[1:]
    genes = {}
    for row in data[1:]:
        genes[row[0]] = [float(x) for x in row[1:]]
    return genes, samples

def load_kallisto_csv(path):
    genes, samples = {}, []
    with open(path) as f:
        reader = csv.reader(f)
        header = next(reader); samples = header[1:]
        for row in reader:
            genes[row[0]] = [float(x) for x in row[1:]]
    return genes, samples

# ═══════════════════════════════════════════════════════════════
# PLOT 1: Gene Proportion Stacked Bar (top 30 samples)
# ═══════════════════════════════════════════════════════════════
def plot_gene_proportions(fc_genes, samples):
    print("[1/7] Gene proportion stacked bar...")
    n = len(samples)
    totals = [sum(fc_genes[g][i] for g in GENE_ORDER if g in fc_genes) for i in range(n)]
    
    # Top 30 by viral load
    top_idx = sorted(range(n), key=lambda i: totals[i], reverse=True)[:30]
    top_idx = [i for i in top_idx if totals[i] > 0]
    
    fig, ax = plt.subplots(figsize=(14, 7))
    bottom = np.zeros(len(top_idx))
    
    for gene in GENE_ORDER:
        if gene not in fc_genes: continue
        vals = np.array([fc_genes[gene][i]/totals[i]*100 if totals[i]>0 else 0 for i in top_idx])
        ax.bar(range(len(top_idx)), vals, bottom=bottom, label=gene,
               color=COLORS.get(gene,"#888"), edgecolor='white', linewidth=0.5)
        bottom += vals
    
    ax.set_xticks(range(len(top_idx)))
    ax.set_xticklabels([samples[i].replace("SRR38105","") for i in top_idx], rotation=90, fontsize=7)
    ax.set_ylabel("% of Total Reads", fontsize=12, fontweight='bold')
    ax.set_xlabel("Sample (ranked by total viral load)", fontsize=12, fontweight='bold')
    ax.set_title("Gene Composition per Sample — Top 30 by Viral Load", fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9)
    ax.set_ylim(0, 100)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "gene_proportions_stacked.png"), dpi=200)
    plt.close()
    print("  ✓ gene_proportions_stacked.png")

# ═══════════════════════════════════════════════════════════════
# PLOT 2: Viral Load Histogram (all samples)
# ═══════════════════════════════════════════════════════════════
def plot_viral_load_histogram(fc_genes, samples):
    print("[2/7] Viral load histogram...")
    n = len(samples)
    totals = [sum(fc_genes[g][i] for g in GENE_ORDER if g in fc_genes) for i in range(n)]
    
    # Log10 transform (add 1 to handle zeros)
    log_totals = [np.log10(t+1) for t in totals]
    zero_count = sum(1 for t in totals if t == 0)
    nonzero = sum(1 for t in totals if t > 0)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left: histogram of log10 viral load
    ax1.hist(log_totals, bins=40, color='#E63946', edgecolor='white', alpha=0.8)
    ax1.set_xlabel("log₁₀(Total Ebola Reads + 1)", fontsize=12, fontweight='bold')
    ax1.set_ylabel("Number of Samples", fontsize=12, fontweight='bold')
    ax1.set_title("Distribution of Viral Load", fontsize=14, fontweight='bold')
    ax1.axvline(np.median(log_totals), color='navy', linestyle='--', linewidth=2, label=f'Median')
    ax1.legend(fontsize=10)
    ax1.spines['top'].set_visible(False); ax1.spines['right'].set_visible(False)
    
    # Right: pie chart of zero vs nonzero
    ax2.pie([nonzero, zero_count], labels=[f'Detected\n({nonzero})', f'Not Detected\n({zero_count})'],
            colors=['#E63946','#d3d3d3'], autopct='%1.1f%%', startangle=90,
            textprops={'fontsize':12, 'fontweight':'bold'}, wedgeprops={'edgecolor':'white','linewidth':2})
    ax2.set_title("Ebola Detection Rate", fontsize=14, fontweight='bold')
    
    plt.suptitle(f"Viral Load Distribution — {n} Samples", fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "viral_load_histogram.png"), dpi=200, bbox_inches='tight')
    plt.close()
    print("  ✓ viral_load_histogram.png")

# ═══════════════════════════════════════════════════════════════
# PLOT 3: Sample Clustering Dendrogram
# ═══════════════════════════════════════════════════════════════
def plot_dendrogram(fc_genes, samples):
    print("[3/7] Sample clustering dendrogram...")
    genes = [g for g in GENE_ORDER if g in fc_genes]
    n = len(samples)
    
    # Build matrix: samples x genes (log2+1)
    matrix = np.zeros((n, len(genes)))
    for j, g in enumerate(genes):
        for i in range(n):
            matrix[i, j] = np.log2(fc_genes[g][i] + 1)
    
    # Filter to samples with >100 total reads for meaningful clustering
    totals = matrix.sum(axis=1)
    mask = totals > 5  # log2 scale, ~32 reads minimum
    if mask.sum() < 5:
        print("  [WARN] Too few samples with reads for dendrogram, skipping")
        return
    
    sub_matrix = matrix[mask]
    sub_labels = [s.replace("SRR38105","") for s, m in zip(samples, mask) if m]
    
    # Only show top 60 for readability
    if len(sub_labels) > 60:
        top_idx = np.argsort(sub_matrix.sum(axis=1))[-60:]
        sub_matrix = sub_matrix[top_idx]
        sub_labels = [sub_labels[i] for i in top_idx]
    
    Z = linkage(sub_matrix, method='ward')
    
    fig, ax = plt.subplots(figsize=(16, 6))
    dendrogram(Z, labels=sub_labels, leaf_rotation=90, leaf_font_size=7,
               color_threshold=Z[-3,2] if len(Z) > 3 else 0, ax=ax)
    ax.set_ylabel("Ward Distance", fontsize=12, fontweight='bold')
    ax.set_title("Sample Clustering by Ebola Gene Expression\n(Ward's method, log₂ counts)", fontsize=14, fontweight='bold')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "sample_dendrogram.png"), dpi=200)
    plt.close()
    print("  ✓ sample_dendrogram.png")

# ═══════════════════════════════════════════════════════════════
# PLOT 4: Variant Summary (SNPs vs Indels)
# ═══════════════════════════════════════════════════════════════
def plot_variant_summary():
    print("[4/7] Variant summary...")
    var_file = os.path.join(PROJECT, "variants_git", "variant_summary.tsv")
    if not os.path.exists(var_file):
        print("  [WARN] variant_summary.tsv not found, skipping")
        return
    
    data = load_tsv(var_file)
    header = data[0]
    samples, snps, indels = [], [], []
    for row in data[1:]:
        if len(row) < 5: continue
        try:
            samples.append(row[0])
            snps.append(int(row[3]))
            indels.append(int(row[4]))
        except: continue
    
    if not samples:
        print("  [WARN] No variant data, skipping")
        return
    
    # Sort by total variants
    total = [s+i for s,i in zip(snps, indels)]
    order = sorted(range(len(samples)), key=lambda i: total[i], reverse=True)
    n_show = min(50, len(order))
    order = order[:n_show]
    
    fig, ax = plt.subplots(figsize=(14, 7))
    x = range(n_show)
    ax.bar(x, [snps[i] for i in order], label='SNPs', color='#E63946', edgecolor='white')
    ax.bar(x, [indels[i] for i in order], bottom=[snps[i] for i in order],
           label='Indels', color='#457B9D', edgecolor='white')
    ax.set_xticks(x)
    ax.set_xticklabels([samples[i].replace("SRR38105","") for i in order], rotation=90, fontsize=7)
    ax.set_xlabel("Sample", fontsize=12, fontweight='bold')
    ax.set_ylabel("Number of Variants", fontsize=12, fontweight='bold')
    ax.set_title(f"Variant Calls per Sample — Top {n_show}\n(bcftools, filtered)", fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "variant_summary.png"), dpi=200)
    plt.close()
    print("  ✓ variant_summary.png")

# ═══════════════════════════════════════════════════════════════
# PLOT 5: featureCounts Assignment Summary
# ═══════════════════════════════════════════════════════════════
def plot_assignment_summary():
    print("[5/7] featureCounts assignment summary...")
    summary_file = os.path.join(PROJECT, "counts_git", "gene_counts.txt.summary")
    if not os.path.exists(summary_file):
        print("  [WARN] gene_counts.txt.summary not found, skipping")
        return
    
    data = load_tsv(summary_file)
    categories = {}
    for row in data[1:]:  # skip header
        cat = row[0]
        total = sum(int(x) for x in row[1:])
        if total > 0:
            categories[cat] = total
    
    # Group small categories
    total_all = sum(categories.values())
    major = {k:v for k,v in categories.items() if v/total_all > 0.01}
    other = sum(v for k,v in categories.items() if v/total_all <= 0.01)
    if other > 0:
        major['Other'] = other
    
    labels = [k.replace("Unassigned_","") for k in major.keys()]
    sizes = list(major.values())
    colors_pie = ['#2A9D8F','#E63946','#457B9D','#E9C46A','#F4A261','#264653','#6A4C93','#d3d3d3']
    
    fig, ax = plt.subplots(figsize=(10, 8))
    wedges, texts, autotexts = ax.pie(sizes, labels=labels, autopct='%1.1f%%',
        colors=colors_pie[:len(sizes)], startangle=90,
        textprops={'fontsize':10}, wedgeprops={'edgecolor':'white','linewidth':2})
    for t in autotexts: t.set_fontweight('bold')
    ax.set_title("featureCounts Read Assignment Summary\n(all samples combined)", fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "assignment_summary_pie.png"), dpi=200)
    plt.close()
    print("  ✓ assignment_summary_pie.png")

# ═══════════════════════════════════════════════════════════════
# PLOT 6: Per-sample Kallisto vs featureCounts scatter
# ═══════════════════════════════════════════════════════════════
def plot_sample_correlation(fc_genes, fc_samples, kal_genes, kal_samples):
    print("[6/7] Per-sample Kallisto vs featureCounts...")
    
    # Total reads per sample for each method
    common = set(fc_samples) & set(kal_samples)
    if len(common) < 5:
        print(f"  [WARN] Only {len(common)} common samples, skipping")
        return
    
    fc_totals, kal_totals = [], []
    for s in common:
        fi = fc_samples.index(s)
        ki = kal_samples.index(s)
        fc_totals.append(sum(fc_genes[g][fi] for g in GENE_ORDER if g in fc_genes))
        kal_totals.append(sum(kal_genes[g][ki] for g in GENE_ORDER if g in kal_genes))
    
    fc_log = [np.log10(x+1) for x in fc_totals]
    kal_log = [np.log10(x+1) for x in kal_totals]
    
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(fc_log, kal_log, alpha=0.4, s=20, c='#E63946', edgecolors='white', linewidth=0.3)
    
    max_val = max(max(fc_log), max(kal_log)) * 1.05
    ax.plot([0, max_val], [0, max_val], '--', color='grey', alpha=0.5)
    
    corr = np.corrcoef(fc_log, kal_log)[0,1]
    ax.text(0.05, 0.95, f'r = {corr:.4f}\nn = {len(common)} samples',
            transform=ax.transAxes, fontsize=12, fontweight='bold', va='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    ax.set_xlabel("featureCounts — log₁₀(total reads + 1)", fontsize=12, fontweight='bold')
    ax.set_ylabel("Kallisto — log₁₀(total reads + 1)", fontsize=12, fontweight='bold')
    ax.set_title("Per-Sample Correlation\nfeatureCounts vs Kallisto", fontsize=14, fontweight='bold')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "sample_correlation_scatter.png"), dpi=200)
    plt.close()
    print("  ✓ sample_correlation_scatter.png")

# ═══════════════════════════════════════════════════════════════
# PLOT 7: Gene-Gene Correlation Matrix
# ═══════════════════════════════════════════════════════════════
def plot_gene_correlation(fc_genes):
    print("[7/7] Gene-gene correlation matrix...")
    genes = [g for g in GENE_ORDER if g in fc_genes]
    n_genes = len(genes)
    
    # Log2 transform
    data = {g: np.log2(np.array(fc_genes[g]) + 1) for g in genes}
    
    corr_matrix = np.zeros((n_genes, n_genes))
    for i, g1 in enumerate(genes):
        for j, g2 in enumerate(genes):
            corr_matrix[i,j] = np.corrcoef(data[g1], data[g2])[0,1]
    
    fig, ax = plt.subplots(figsize=(8, 7))
    cmap = LinearSegmentedColormap.from_list("corr", ["#264653","#2A9D8F","#E9C46A","#E76F51","#9B2226"])
    im = ax.imshow(corr_matrix, cmap=cmap, vmin=0, vmax=1, aspect='equal')
    
    ax.set_xticks(range(n_genes)); ax.set_yticks(range(n_genes))
    ax.set_xticklabels(genes, fontsize=12, fontweight='bold')
    ax.set_yticklabels(genes, fontsize=12, fontweight='bold')
    
    # Add correlation values
    for i in range(n_genes):
        for j in range(n_genes):
            color = 'white' if corr_matrix[i,j] > 0.7 else 'black'
            ax.text(j, i, f'{corr_matrix[i,j]:.2f}', ha='center', va='center',
                    fontsize=11, fontweight='bold', color=color)
    
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Pearson Correlation", fontsize=11)
    ax.set_title("Gene-Gene Expression Correlation\n(log₂ counts, 276 samples)", fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "gene_correlation_matrix.png"), dpi=200)
    plt.close()
    print("  ✓ gene_correlation_matrix.png")

# ═══════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("═" * 60)
    print("  Ebola Pipeline — Additional Publication Figures (v2)")
    print("═" * 60 + "\n")
    
    fc_file = os.path.join(PROJECT, "counts_git", "gene_counts_clean.txt")
    kal_file = os.path.join(PROJECT, "kallisto_git", "count_matrix_kallisto.csv")
    
    fc_genes, fc_samples = load_counts(fc_file)
    print(f"[INFO] featureCounts: {len(fc_genes)} genes × {len(fc_samples)} samples")
    
    kal_genes, kal_samples = None, None
    if os.path.exists(kal_file):
        kal_genes, kal_samples = load_kallisto_csv(kal_file)
        print(f"[INFO] Kallisto: {len(kal_genes)} genes × {len(kal_samples)} samples\n")
    
    plot_gene_proportions(fc_genes, fc_samples)
    plot_viral_load_histogram(fc_genes, fc_samples)
    plot_dendrogram(fc_genes, fc_samples)
    plot_variant_summary()
    plot_assignment_summary()
    if kal_genes:
        plot_sample_correlation(fc_genes, fc_samples, kal_genes, kal_samples)
    plot_gene_correlation(fc_genes)
    
    total = len([f for f in os.listdir(OUT) if f.endswith('.png')])
    print(f"\n{'═'*60}")
    print(f"  Done! {total} total figures in {OUT}")
    print(f"{'═'*60}")
