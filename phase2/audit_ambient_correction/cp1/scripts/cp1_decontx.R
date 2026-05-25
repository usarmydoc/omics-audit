#!/usr/bin/env Rscript
# CP1 Deliverable A — DecontX per-gene contamination
# Usage: Rscript cp1_decontx.R <filt_dir> <out_prefix>
#   filt_dir: STARsolo Gene/filtered CellRanger-format dir
# decontX() run with defaults (no external background) per locked scope.
# Per-gene contamination = (orig_total - decontX_total) / orig_total on cells.

suppressPackageStartupMessages({
  library(celda); library(DropletUtils); library(Matrix)
  library(SingleCellExperiment); library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE)
filt_dir <- args[1]; out_prefix <- args[2]
set.seed(42)

cat("[DecontX] reading filtered matrix\n")
sce <- read10xCounts(filt_dir, col.names = TRUE)
gene_id  <- rowData(sce)$ID
gene_sym <- rowData(sce)$Symbol
cat(sprintf("[DecontX] %d genes x %d cells\n", nrow(sce), ncol(sce)))

# decontX needs >1 cell and removes all-zero genes internally; run with defaults
sce <- decontX(sce)
contam_cell <- colData(sce)$decontX_contamination      # per-cell contamination
adj <- decontXcounts(sce)                               # corrected counts

orig_total <- Matrix::rowSums(counts(sce))
corr_total <- Matrix::rowSums(adj)
df <- data.frame(gene_id = gene_id, gene_symbol = gene_sym,
                 orig_total = as.numeric(orig_total),
                 corr_total = as.numeric(corr_total),
                 stringsAsFactors = FALSE)
write.table(df, paste0(out_prefix, "_pergene.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

summ <- list(tool = "DecontX", version = as.character(packageVersion("celda")),
             n_cells = ncol(sce), n_genes = nrow(sce),
             contam_cell_mean = mean(contam_cell, na.rm = TRUE),
             contam_cell_median = median(contam_cell, na.rm = TRUE),
             contam_cell_min = min(contam_cell, na.rm = TRUE),
             contam_cell_max = max(contam_cell, na.rm = TRUE),
             global_frac_removed = 1 - sum(corr_total) / sum(orig_total))
writeLines(toJSON(summ, auto_unbox = TRUE, pretty = TRUE),
           paste0(out_prefix, "_summary.json"))
cat(sprintf("[DecontX] DONE contam_cell_mean=%.4f global_frac_removed=%.4f\n",
            mean(contam_cell, na.rm = TRUE), summ$global_frac_removed))
