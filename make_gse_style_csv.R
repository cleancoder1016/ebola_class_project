library(data.table)

mat_file  <- "count_matrix.tsv"
meta_file <- "../prjna938511_metadata.txt"
out_file  <- "count_matrix_gse_style.csv"

# ── Load count matrix ─────────────────────────────────────────────────────────
message("Reading count matrix ...")
mat <- fread(mat_file, sep = "\t", header = TRUE, data.table = FALSE)
rownames(mat) <- mat[[1]]
mat <- mat[, -1, drop = FALSE]          # drop gene-ID column; keep counts only
message("  Dimensions: ", nrow(mat), " genes x ", ncol(mat), " samples")

# ── Build SRR → stripped library-name map ────────────────────────────────────
meta <- read.table(meta_file, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, quote = "")
meta$lib_short <- sub("^GSM[0-9]+_", "", meta$library_name)  # strip GSM##### prefix

srr2lib <- setNames(meta$lib_short, meta$run_accession)

missing <- setdiff(colnames(mat), names(srr2lib))
if (length(missing))
  warning(length(missing), " SRR(s) not found in metadata: ",
          paste(missing, collapse = ", "))

# ── Rename columns ────────────────────────────────────────────────────────────
colnames(mat) <- srr2lib[colnames(mat)]     # NA for any SRR not in metadata
message("  Columns renamed to library names.")

# ── Replace gene IDs with natural numbers ────────────────────────────────────
rownames(mat) <- seq_len(nrow(mat))
message("  Row names replaced with 1 .. ", nrow(mat))

# ── Write CSV ─────────────────────────────────────────────────────────────────
message("Writing ", out_file, " ...")
write.csv(mat, file = out_file, row.names = TRUE, quote = FALSE)
message("Done. Rows: ", nrow(mat), "  Cols: ", ncol(mat))
