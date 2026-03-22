# scran Normalization for P7 Audit
# Run this in RStudio after prep_normalization_inputs.py
# Reads raw count matrices, runs scran pooling normalization, saves normalized data

SCRIPT_DIR <- normalizePath(dirname(sys.frame(1)$ofile), mustWork = FALSE)
if (is.na(SCRIPT_DIR) || SCRIPT_DIR == "") SCRIPT_DIR <- getwd()
P7_DIR <- file.path(dirname(dirname(SCRIPT_DIR)), "projects", "p7_normalization")
setwd(P7_DIR)

library(scran)
library(scuttle)
library(SingleCellExperiment)

# ── Config ──
INPUT_DIR <- "r_inputs"
OUTPUT_DIR <- "scran_results"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

input_files <- list.files(INPUT_DIR, pattern = "\\.tsv\\.gz$", full.names = TRUE)
input_files <- input_files[!grepl("_labels\\.tsv", input_files)]

if (length(input_files) == 0) {
  stop("No input files found. Run prep_normalization_inputs.py first.")
}

cat("Found", length(input_files), "datasets to process\n\n")

for (f in input_files) {
  dataset_name <- gsub("\\.tsv\\.gz$", "", basename(f))
  out_file <- file.path(OUTPUT_DIR, paste0(dataset_name, "_scran.tsv.gz"))

  if (file.exists(out_file)) {
    cat("Skipping", dataset_name, "(already done)\n")
    next
  }

  cat("=== Processing", dataset_name, "===\n")

  # Read count matrix (genes x cells)
  mat <- read.table(f, header = TRUE, row.names = 1, sep = "\t",
                    check.names = FALSE)
  mat <- as.matrix(mat)
  cat("  Matrix:", nrow(mat), "genes x", ncol(mat), "cells\n")

  # Create SingleCellExperiment
  sce <- SingleCellExperiment(assays = list(counts = mat))

  # Run scran normalization
  cat("  Running scran (quickCluster + computeSumFactors + logNormCounts)...\n")
  t0 <- Sys.time()

  tryCatch({
    clust <- quickCluster(sce, min.size = 50)
    sce <- computeSumFactors(sce, clusters = clust)
    sce <- logNormCounts(sce)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat("  Done in", round(elapsed, 1), "seconds\n")

    # Extract normalized data
    norm_data <- as.matrix(logcounts(sce))
    cat("  Normalized:", nrow(norm_data), "genes x", ncol(norm_data), "cells\n")

    # Save as gzipped TSV
    gz <- gzfile(out_file, "w")
    write.table(norm_data, gz, sep = "\t", quote = FALSE)
    close(gz)

    # Save runtime
    rt_file <- file.path(OUTPUT_DIR, paste0(dataset_name, "_runtime.txt"))
    writeLines(as.character(round(elapsed, 2)), rt_file)

    cat("  Saved:", out_file, "\n\n")
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    cat("  Skipping this dataset\n\n")
  })

  # Clean up memory
  rm(sce, mat)
  if (exists("norm_data")) rm(norm_data)
  gc()
}

cat("=== All done ===\n")
cat("Results in:", OUTPUT_DIR, "\n")
