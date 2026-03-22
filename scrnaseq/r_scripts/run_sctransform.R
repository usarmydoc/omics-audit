# SCTransform Normalization for P7 Audit
# Run this in RStudio after prep_normalization_inputs.py
# Reads raw count matrices, runs SCTransform, saves normalized data

SCRIPT_DIR <- normalizePath(dirname(sys.frame(1)$ofile), mustWork = FALSE)
if (is.na(SCRIPT_DIR) || SCRIPT_DIR == "") SCRIPT_DIR <- getwd()
P7_DIR <- file.path(dirname(dirname(SCRIPT_DIR)), "projects", "p7_normalization")
setwd(P7_DIR)

library(Seurat)
library(Matrix)

# ── Config ──
INPUT_DIR <- "r_inputs"
OUTPUT_DIR <- "sctransform_results"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

input_files <- list.files(INPUT_DIR, pattern = "\\.tsv\\.gz$", full.names = TRUE)
input_files <- input_files[!grepl("_labels\\.tsv", input_files)]

if (length(input_files) == 0) {
  stop("No input files found. Run prep_normalization_inputs.py first.")
}

cat("Found", length(input_files), "datasets to process\n\n")

for (f in input_files) {
  dataset_name <- gsub("\\.tsv\\.gz$", "", basename(f))
  out_file <- file.path(OUTPUT_DIR, paste0(dataset_name, "_sctransform.tsv.gz"))

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

  # Create Seurat object
  options(future.globals.maxSize = 2 * 1024^3)
  sobj <- CreateSeuratObject(counts = Matrix(mat, sparse = TRUE))

  # Run SCTransform
  cat("  Running SCTransform...\n")
  t0 <- Sys.time()
  sobj <- SCTransform(sobj, verbose = FALSE)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat("  Done in", round(elapsed, 1), "seconds\n")

  # Extract normalized data (SCT data layer)
  norm_data <- as.matrix(GetAssayData(sobj, assay = "SCT", layer = "data"))
  cat("  Normalized:", nrow(norm_data), "genes x", ncol(norm_data), "cells\n")

  # Save as gzipped TSV
  gz <- gzfile(out_file, "w")
  write.table(norm_data, gz, sep = "\t", quote = FALSE)
  close(gz)

  # Also save runtime
  rt_file <- file.path(OUTPUT_DIR, paste0(dataset_name, "_runtime.txt"))
  writeLines(as.character(round(elapsed, 2)), rt_file)

  cat("  Saved:", out_file, "\n\n")

  # Clean up memory
  rm(sobj, mat, norm_data)
  gc()
}

cat("=== All done ===\n")
cat("Results in:", OUTPUT_DIR, "\n")
