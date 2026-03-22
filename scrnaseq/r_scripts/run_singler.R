# SingleR Annotation for P5 Audit
# Run this in RStudio. It will:
# 1. Query the same datasets P5 uses (from saved dataset list)
# 2. Run SingleR on each dataset
# 3. Save predictions as TSV for P5 to load

SCRIPT_DIR <- normalizePath(dirname(sys.frame(1)$ofile), mustWork = FALSE)
if (is.na(SCRIPT_DIR) || SCRIPT_DIR == "") SCRIPT_DIR <- getwd()
P5_DIR <- file.path(dirname(dirname(SCRIPT_DIR)), "projects", "p5_annotation_tools")
setwd(P5_DIR)

library(SingleR)
library(celldex)
library(Matrix)

# ── Config ──
OUTPUT_DIR <- file.path(P5_DIR, "singler_results")
DATASET_DIR <- file.path(P5_DIR, "singler_inputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Load references (downloads once, cached after) ──
cat("Loading HumanPrimaryCellAtlasData...\n")
ref_human <- HumanPrimaryCellAtlasData()
cat("Loading MouseRNAseqData...\n")
ref_mouse <- MouseRNAseqData()

# ── Process each dataset ──
input_files <- list.files(DATASET_DIR, pattern = "\\.tsv\\.gz$", full.names = TRUE)

if (length(input_files) == 0) {
  cat("No input files found in", DATASET_DIR, "\n")
  cat("Run the Python prep script first to export matrices.\n")
  stop("No input files")
}

for (f in input_files) {
  dataset_name <- gsub("\\.tsv\\.gz$", "", basename(f))
  out_file <- file.path(OUTPUT_DIR, paste0(dataset_name, "_singler.tsv"))

  if (file.exists(out_file)) {
    cat("Skipping", dataset_name, "(already done)\n")
    next
  }

  cat("\n=== Processing", dataset_name, "===\n")

  # Read count matrix (genes x cells, tab-separated, gzipped)
  mat <- read.table(f, header = TRUE, row.names = 1, sep = "\t",
                    check.names = FALSE)
  mat <- as.matrix(mat)

  cat("  Matrix:", nrow(mat), "genes x", ncol(mat), "cells\n")

  # Pick reference based on organism (encoded in filename)
  if (grepl("mouse|mus_musculus", dataset_name, ignore.case = TRUE)) {
    ref <- ref_mouse
    labels <- ref$label.main
  } else {
    ref <- ref_human
    labels <- ref$label.main
  }

  # Run SingleR
  cat("  Running SingleR...\n")
  t0 <- Sys.time()
  results <- SingleR(test = mat, ref = ref, labels = labels)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat("  Done in", round(elapsed, 1), "seconds\n")

  # Save predictions
  preds <- data.frame(
    cell_id = colnames(mat),
    singler_label = results$labels,
    stringsAsFactors = FALSE
  )
  write.table(preds, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  Saved:", out_file, "\n")
}

cat("\n=== All done ===\n")
cat("Results in:", OUTPUT_DIR, "\n")
