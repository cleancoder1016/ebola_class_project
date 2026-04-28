#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════════════════════
# deseq2_analysis.R — Exploratory & differential expression analysis
# ═══════════════════════════════════════════════════════════════════════════════
# Reads the featureCounts count matrix and performs:
#   - Pre-filtering (remove low-count genes)
#   - DESeq2 normalization & variance stabilizing transformation
#   - PCA plot
#   - Sample-to-sample distance heatmap
#   - Top variable genes heatmap
#   - Gene expression summary statistics
#
# If sample_metadata.csv is present, also performs differential expression.
# ═══════════════════════════════════════════════════════════════════════════════

# ── Parse arguments ─────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript deseq2_analysis.R <counts_file> <output_dir> [metadata_file]")
}

counts_file   <- args[1]
output_dir    <- args[2]
metadata_file <- if (length(args) >= 3) args[3] else NULL

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("═══════════════════════════════════════════════════════════\n")
cat("DESeq2 RNA-seq Analysis — Ebola Outbreak Dataset\n")
cat("Count matrix:", counts_file, "\n")
cat("Output dir:  ", output_dir, "\n")
cat("Metadata:    ", ifelse(is.null(metadata_file), "None (exploratory only)", metadata_file), "\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# ── Load libraries ──────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    for (pkg in c("DESeq2", "pheatmap", "RColorBrewer", "ggplot2", "ggrepel")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            if (pkg %in% c("DESeq2")) {
                BiocManager::install(pkg, ask = FALSE, update = FALSE)
            } else {
                install.packages(pkg, repos = "https://cloud.r-project.org")
            }
        }
    }
    library(DESeq2)
    library(pheatmap)
    library(RColorBrewer)
    library(ggplot2)
    library(ggrepel)
})

# ── Load count matrix ──────────────────────────────────────────────────────
cat("[INFO] Loading count matrix...\n")
counts_raw <- read.delim(counts_file, comment.char = "#", row.names = 1)

# Remove featureCounts metadata columns (Chr, Start, End, Strand, Length)
# These are the first 5 columns after Geneid in the raw featureCounts output
if (all(c("Chr", "Start", "End", "Strand", "Length") %in% colnames(counts_raw))) {
    counts_matrix <- counts_raw[, !(colnames(counts_raw) %in% c("Chr", "Start", "End", "Strand", "Length"))]
} else {
    counts_matrix <- counts_raw
}

# Clean column names (remove path prefixes, .bam suffixes)
colnames(counts_matrix) <- gsub(".*\\/", "", colnames(counts_matrix))
colnames(counts_matrix) <- gsub("\\.dedup\\.bam$|\\.sorted\\.bam$|\\.bam$", "", colnames(counts_matrix))

cat("[INFO] Raw count matrix dimensions:", nrow(counts_matrix), "genes x", ncol(counts_matrix), "samples\n")

# ── Pre-filtering ──────────────────────────────────────────────────────────
# Remove genes with fewer than 10 total reads across all samples
keep <- rowSums(counts_matrix) >= 10
counts_filtered <- counts_matrix[keep, ]
cat("[INFO] After filtering (>=10 total reads):", nrow(counts_filtered), "genes retained\n")

if (nrow(counts_filtered) == 0) {
    cat("[WARN] No genes passed filtering. Exiting.\n")
    quit(save = "no", status = 0)
}

# ── Create DESeq2 dataset ─────────────────────────────────────────────────
# Load metadata if available, otherwise create dummy
if (!is.null(metadata_file) && file.exists(metadata_file)) {
    cat("[INFO] Loading sample metadata from:", metadata_file, "\n")
    metadata <- read.csv(metadata_file, row.names = 1)
    # Ensure sample order matches
    common_samples <- intersect(colnames(counts_filtered), rownames(metadata))
    counts_filtered <- counts_filtered[, common_samples]
    metadata <- metadata[common_samples, , drop = FALSE]
    design_formula <- ~ condition  # Assumes a 'condition' column
} else {
    cat("[INFO] No metadata — performing exploratory analysis only.\n")
    metadata <- data.frame(
        row.names = colnames(counts_filtered),
        condition = rep("sample", ncol(counts_filtered))
    )
    design_formula <- ~ 1
}

dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(counts_filtered)),
    colData = metadata,
    design = design_formula
)

# ── Size factor estimation (handle genes with many zeros) ──────────────────
# Many Ebola samples have very low viral load → zeros in count matrix.
# Use type='poscounts' to handle this.
cat("[INFO] Estimating size factors (type=poscounts for zero-inflated data)...\n")
dds <- tryCatch(
    estimateSizeFactors(dds, type = "poscounts"),
    error = function(e) {
        cat("[WARN] poscounts failed, using manual size factors.\n")
        sizeFactors(dds) <- rep(1, ncol(dds))
        dds
    }
)

# ── Variance stabilizing transformation ────────────────────────────────────
# NOTE: Ebola has only 7 genes, which is fewer than vst()'s default nsub.
cat("[INFO] Running variance stabilizing transformation...\n")

# Robust fallback chain for transformation
vsd <- tryCatch({
    cat("[INFO] Trying varianceStabilizingTransformation...\n")
    varianceStabilizingTransformation(dds, blind = TRUE)
}, error = function(e1) {
    cat("[WARN] VST failed:", e1$message, "\n")
    tryCatch({
        cat("[INFO] Trying normTransform...\n")
        normTransform(dds)
    }, error = function(e2) {
        cat("[WARN] normTransform also failed:", e2$message, "\n")
        cat("[INFO] Falling back to manual log2(count+1) transformation.\n")
        # Create a manual log2-transformed DESeqTransform object
        log2_mat <- log2(counts(dds, normalized = FALSE) + 1)
        se <- SummarizedExperiment(
            assays = list(log2_mat),
            colData = colData(dds)
        )
        DESeqTransform(se)
    })
})

# ═══════════════════════════════════════════════════════════════════════════
# PLOT 1: PCA Plot
# ═══════════════════════════════════════════════════════════════════════════
cat("[INFO] Generating PCA plot...\n")
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"), 1)

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text_repel(size = 2.5, max.overlaps = 20) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    ggtitle("PCA — Ebola RNA-seq Samples") +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
    )

ggsave(file.path(output_dir, "pca_plot.pdf"), pca_plot, width = 10, height = 8)
ggsave(file.path(output_dir, "pca_plot.png"), pca_plot, width = 10, height = 8, dpi = 150)

# ═══════════════════════════════════════════════════════════════════════════
# PLOT 2: Sample-to-Sample Distance Heatmap
# ═══════════════════════════════════════════════════════════════════════════
cat("[INFO] Generating sample distance heatmap...\n")
sample_dists <- dist(t(assay(vsd)))
dist_matrix <- as.matrix(sample_dists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pdf(file.path(output_dir, "sample_distance_heatmap.pdf"), width = 12, height = 10)
pheatmap(dist_matrix,
    clustering_distance_rows = sample_dists,
    clustering_distance_cols = sample_dists,
    col = colors,
    main = "Sample-to-Sample Euclidean Distance",
    fontsize_row = 7,
    fontsize_col = 7
)
dev.off()

# ═══════════════════════════════════════════════════════════════════════════
# PLOT 3: Top Variable Genes Heatmap
# ═══════════════════════════════════════════════════════════════════════════
cat("[INFO] Generating top variable genes heatmap...\n")
top_n <- min(30, nrow(counts_filtered))
top_var_genes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), top_n)

pdf(file.path(output_dir, "top_variable_genes_heatmap.pdf"), width = 14, height = 8)
pheatmap(assay(vsd)[top_var_genes, ],
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    fontsize_col = 7,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    main = paste("Top", top_n, "Most Variable Genes (z-scores)")
)
dev.off()

# ═══════════════════════════════════════════════════════════════════════════
# PLOT 4: Per-Gene Expression Distribution
# ═══════════════════════════════════════════════════════════════════════════
cat("[INFO] Generating expression distribution plot...\n")
log_counts <- log2(counts_filtered + 1)
gene_means <- rowMeans(log_counts)

expr_df <- data.frame(mean_log2_expr = gene_means)
expr_plot <- ggplot(expr_df, aes(x = mean_log2_expr)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.8) +
    xlab("Mean log2(count + 1)") +
    ylab("Number of Genes") +
    ggtitle("Gene Expression Distribution") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "expression_distribution.pdf"), expr_plot, width = 8, height = 6)

# ═══════════════════════════════════════════════════════════════════════════
# Differential Expression (only if metadata with 'condition' column)
# ═══════════════════════════════════════════════════════════════════════════
if (!is.null(metadata_file) && file.exists(metadata_file) && "condition" %in% colnames(metadata)) {
    if (length(unique(metadata$condition)) >= 2) {
        cat("[INFO] Running differential expression analysis...\n")
        dds_de <- DESeq(dds)
        res <- results(dds_de, alpha = 0.05)
        res_ordered <- res[order(res$padj), ]

        # Save results table
        write.csv(as.data.frame(res_ordered),
            file = file.path(output_dir, "de_results.csv"),
            row.names = TRUE
        )

        sig_genes <- sum(res_ordered$padj < 0.05, na.rm = TRUE)
        cat("[INFO] Significant genes (padj < 0.05):", sig_genes, "\n")

        # MA Plot
        pdf(file.path(output_dir, "ma_plot.pdf"), width = 8, height = 6)
        plotMA(res, main = "MA Plot — Differential Expression", ylim = c(-5, 5))
        dev.off()

        # Volcano Plot
        vol_data <- as.data.frame(res)
        vol_data$significant <- ifelse(!is.na(vol_data$padj) & vol_data$padj < 0.05, "Significant", "Not Sig.")

        volcano <- ggplot(vol_data, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
            geom_point(alpha = 0.6, size = 1.5) +
            scale_color_manual(values = c("grey60", "firebrick")) +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
            xlab("Log2 Fold Change") +
            ylab("-Log10 P-value") +
            ggtitle("Volcano Plot") +
            theme_minimal(base_size = 12) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))

        ggsave(file.path(output_dir, "volcano_plot.pdf"), volcano, width = 8, height = 6)
    } else {
        cat("[WARN] Only one condition level found — skipping DE analysis.\n")
    }
} else {
    cat("[INFO] No metadata with 'condition' column — skipping DE analysis.\n")
}

# ═══════════════════════════════════════════════════════════════════════════
# Summary Statistics
# ═══════════════════════════════════════════════════════════════════════════
cat("[INFO] Writing summary statistics...\n")
summary_stats <- data.frame(
    Metric = c("Total genes (raw)", "Genes after filtering", "Total samples",
               "Total assigned reads", "Mean reads per sample", "Median reads per sample"),
    Value = c(
        nrow(counts_matrix),
        nrow(counts_filtered),
        ncol(counts_filtered),
        sum(counts_filtered),
        round(mean(colSums(counts_filtered))),
        round(median(colSums(counts_filtered)))
    )
)
write.csv(summary_stats, file.path(output_dir, "summary_statistics.csv"), row.names = FALSE)

# Per-sample library sizes
lib_sizes <- data.frame(
    Sample = colnames(counts_filtered),
    Library_Size = colSums(counts_filtered)
)
lib_sizes <- lib_sizes[order(lib_sizes$Library_Size), ]
write.csv(lib_sizes, file.path(output_dir, "library_sizes.csv"), row.names = FALSE)

cat("\n═══════════════════════════════════════════════════════════\n")
cat("DESeq2 analysis complete.\n")
cat("Output directory:", output_dir, "\n")
cat("═══════════════════════════════════════════════════════════\n")
