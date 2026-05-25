#!/usr/bin/env Rscript
# CP2 Deliverable B — SoupX under one ordering.
# Usage: Rscript cp2_soupx.R <raw_dir> <filt_dir> <out_prefix> <O1|O2>
#   O1 (correct->QC): SoupX on ALL filtered cells -> corrected -> C2 QC on corrected.
#   O2 (QC->correct): C2 QC on ORIGINAL filtered -> survivors -> SoupX on survivors.
# Correction itself = CP1 defaults (autoEstCont). Clusters via minimal Seurat.
suppressPackageStartupMessages({
  library(SoupX); library(DropletUtils); library(Matrix); library(Seurat); library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE)
raw_dir <- args[1]; filt_dir <- args[2]; out_prefix <- args[3]; ordering <- args[4]
set.seed(42)
SCRIPT_DIR <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(SCRIPT_DIR, "cp2_qc.R"))
try(RhpcBLASctl::blas_set_num_threads(8), silent = TRUE)  # OpenBLAS threads for PCA/ScaleData

cluster_cells <- function(mat) {
  so <- CreateSeuratObject(counts = mat)
  so <- NormalizeData(so, verbose = FALSE); so <- FindVariableFeatures(so, verbose = FALSE)
  so <- ScaleData(so, verbose = FALSE);     so <- RunPCA(so, npcs = 30, verbose = FALSE)
  so <- FindNeighbors(so, dims = 1:30, verbose = FALSE)
  so <- FindClusters(so, resolution = 0.5, verbose = FALSE)
  setNames(as.character(Idents(so)), colnames(so))[colnames(mat)]
}
run_soupx <- function(tod, toc, clusters) {
  sc <- SoupChannel(tod, toc); sc <- setClusters(sc, clusters)
  sc <- autoEstCont(sc, doPlot = FALSE)
  list(adj = adjustCounts(sc), rho = sc$metaData$rho)
}

cat(sprintf("[SoupX %s] reading\n", ordering))
tod <- counts(read10xCounts(raw_dir, col.names = TRUE))
filt_sce <- read10xCounts(filt_dir, col.names = TRUE)
toc_all <- counts(filt_sce); gene_id <- rowData(filt_sce)$ID; gene_sym <- rowData(filt_sce)$Symbol
rownames(tod) <- gene_id; rownames(toc_all) <- gene_id
n_in <- ncol(toc_all)

if (ordering == "O2") {
  keep <- c2_qc_mask(toc_all, gene_sym)              # QC on ORIGINAL counts
  toc <- toc_all[, keep, drop = FALSE]
  cat(sprintf("[SoupX O2] QC original: %d -> %d cells (mito_p95=%.2f)\n", n_in, ncol(toc), attr(keep, "mito_p95")))
  clusters <- cluster_cells(toc)
  res <- run_soupx(tod, toc, clusters); adj <- res$adj
  survivors <- colnames(adj); orig_sub <- toc; mito_p95 <- attr(keep, "mito_p95")
} else {                                              # O1: correct all, then QC corrected
  clusters <- cluster_cells(toc_all)
  res <- run_soupx(tod, toc_all, clusters); adj_all <- res$adj
  keep <- c2_qc_mask(adj_all, gene_sym)              # QC on CORRECTED counts
  adj <- adj_all[, keep, drop = FALSE]; orig_sub <- toc_all[, keep, drop = FALSE]
  survivors <- colnames(adj); mito_p95 <- attr(keep, "mito_p95")
  cat(sprintf("[SoupX O1] correct all %d, QC corrected -> %d cells (mito_p95=%.2f)\n", n_in, length(survivors), mito_p95))
}

# outputs (over survivors)
data.frame(gene_id = gene_id, gene_symbol = gene_sym,
           orig_total = as.numeric(Matrix::rowSums(orig_sub)),
           corr_total = as.numeric(Matrix::rowSums(adj))) |>
  write.table(paste0(out_prefix, "_pergene.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
data.frame(barcode = survivors,
           orig_total = as.numeric(Matrix::colSums(orig_sub)),
           corr_total = as.numeric(Matrix::colSums(adj))) |>
  write.table(paste0(out_prefix, "_percell.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(survivors, paste0(out_prefix, "_survivors.txt"))
writeLines(toJSON(list(tool = "SoupX", ordering = ordering, n_in = n_in,
  n_survivors = length(survivors), retention = length(survivors) / n_in,
  rho_mean = mean(res$rho), mito_p95 = mito_p95,
  global_frac_removed = 1 - sum(Matrix::rowSums(adj)) / sum(Matrix::rowSums(orig_sub))),
  auto_unbox = TRUE, pretty = TRUE), paste0(out_prefix, "_summary.json"))
cat(sprintf("[SoupX %s] DONE survivors=%d retention=%.3f\n", ordering, length(survivors), length(survivors)/n_in))
