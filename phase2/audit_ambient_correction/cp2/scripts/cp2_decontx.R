#!/usr/bin/env Rscript
# CP2 Deliverable B — DecontX under one ordering.
# Usage: Rscript cp2_decontx.R <filt_dir> <out_prefix> <O1|O2>
#   O1 (correct->QC): decontX on ALL filtered cells -> corrected -> C2 QC on corrected.
#   O2 (QC->correct): C2 QC on ORIGINAL filtered -> survivors -> decontX on survivors.
suppressPackageStartupMessages({
  library(celda); library(DropletUtils); library(Matrix); library(SingleCellExperiment); library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE)
filt_dir <- args[1]; out_prefix <- args[2]; ordering <- args[3]
set.seed(42)
SCRIPT_DIR <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(SCRIPT_DIR, "cp2_qc.R"))

sce <- read10xCounts(filt_dir, col.names = TRUE)
gene_id <- rowData(sce)$ID; gene_sym <- rowData(sce)$Symbol
cnt_all <- counts(sce); rownames(cnt_all) <- gene_id
n_in <- ncol(cnt_all)

if (ordering == "O2") {
  keep <- c2_qc_mask(cnt_all, gene_sym)
  sce_run <- sce[, keep]
  cat(sprintf("[DecontX O2] QC original: %d -> %d cells\n", n_in, ncol(sce_run)))
  sce_run <- decontX(sce_run)
  adj <- decontXcounts(sce_run); orig_sub <- counts(sce_run); mito_p95 <- attr(keep, "mito_p95")
  survivors <- colnames(sce_run)
} else {
  sce_run <- decontX(sce)
  adj_all <- decontXcounts(sce_run)
  keep <- c2_qc_mask(adj_all, gene_sym)              # QC on CORRECTED
  adj <- adj_all[, keep, drop = FALSE]; orig_sub <- counts(sce)[, keep, drop = FALSE]
  survivors <- colnames(adj); mito_p95 <- attr(keep, "mito_p95")
  cat(sprintf("[DecontX O1] correct all %d, QC corrected -> %d cells\n", n_in, length(survivors)))
}

data.frame(gene_id = gene_id, gene_symbol = gene_sym,
           orig_total = as.numeric(Matrix::rowSums(orig_sub)),
           corr_total = as.numeric(Matrix::rowSums(adj))) |>
  write.table(paste0(out_prefix, "_pergene.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
data.frame(barcode = survivors,
           orig_total = as.numeric(Matrix::colSums(orig_sub)),
           corr_total = as.numeric(Matrix::colSums(adj))) |>
  write.table(paste0(out_prefix, "_percell.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(survivors, paste0(out_prefix, "_survivors.txt"))
writeLines(toJSON(list(tool = "DecontX", ordering = ordering, n_in = n_in,
  n_survivors = length(survivors), retention = length(survivors) / n_in, mito_p95 = mito_p95,
  global_frac_removed = 1 - sum(Matrix::rowSums(adj)) / sum(Matrix::rowSums(orig_sub))),
  auto_unbox = TRUE, pretty = TRUE), paste0(out_prefix, "_summary.json"))
cat(sprintf("[DecontX %s] DONE survivors=%d retention=%.3f\n", ordering, length(survivors), length(survivors)/n_in))
