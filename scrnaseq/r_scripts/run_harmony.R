# Harmony Batch Correction for P3 Audit
# Run this in RStudio after prep_harmony_inputs.py
# Reads PCA embeddings + donor labels, runs Harmony, saves corrected embeddings
# Then clusters with Leiden and computes ARI vs cell_type ground truth

SCRIPT_DIR <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
PROJECT_DIR <- file.path(SCRIPT_DIR, "..", "projects", "p3_batch_correction")
setwd(PROJECT_DIR)

library(harmony)

# ── Config ──
INPUT_DIR <- "harmony_inputs"
OUTPUT_DIR <- "harmony_results"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Find all PCA files
pca_files <- list.files(INPUT_DIR, pattern = "_pca\\.tsv$", full.names = TRUE)

if (length(pca_files) == 0) {
  stop("No input files found. Run prep_harmony_inputs.py first.")
}

cat("Found", length(pca_files), "datasets to process\n\n")

for (pca_file in pca_files) {
  dataset_name <- gsub("_pca\\.tsv$", "", basename(pca_file))
  meta_file <- file.path(INPUT_DIR, paste0(dataset_name, "_meta.tsv"))
  out_file <- file.path(OUTPUT_DIR, paste0(dataset_name, "_harmony.tsv"))

  if (file.exists(out_file)) {
    cat("Skipping", dataset_name, "(already done)\n")
    next
  }

  if (!file.exists(meta_file)) {
    cat("Skipping", dataset_name, "(no metadata file)\n")
    next
  }

  cat("=== Processing", dataset_name, "===\n")

  # Read PCA embedding and metadata
  pca <- as.matrix(read.table(pca_file, header = TRUE, sep = "\t"))
  meta <- read.table(meta_file, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, comment.char = "")

  cat("  PCA:", nrow(pca), "cells x", ncol(pca), "PCs,",
      length(unique(meta$donor_id)), "donors\n")

  # Run Harmony
  cat("  Running Harmony...\n")
  t0 <- Sys.time()

  tryCatch({
    harmony_out <- RunHarmony(
      pca,
      meta,
      "donor_id",
      verbose = FALSE
    )

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat("  Done in", round(elapsed, 1), "seconds\n")

    # Extract corrected embedding
    corrected <- as.matrix(harmony_out)

    # Save corrected PCA
    write.table(corrected, out_file, sep = "\t", row.names = FALSE,
                col.names = paste0("PC", 1:ncol(corrected)), quote = FALSE)

    # Save runtime
    rt_file <- file.path(OUTPUT_DIR, paste0(dataset_name, "_runtime.txt"))
    writeLines(as.character(round(elapsed, 2)), rt_file)

    cat("  Saved:", out_file, "\n\n")

  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n\n")
  })

  rm(pca, meta)
  if (exists("harmony_out")) rm(harmony_out)
  if (exists("corrected")) rm(corrected)
  gc()
}

cat("=== All done ===\n")
cat("Results in:", OUTPUT_DIR, "\n")
