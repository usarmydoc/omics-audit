#!/usr/bin/env Rscript
# CP3 — generate a corrected cell-by-gene matrix for one (tool, ordering).
# Usage: Rscript cp3_genmatrix.R <SoupX|DecontX> <raw_dir> <filt_dir> <out_dir> <O1|O2>
#   O1 (correct->QC): correct ALL filtered cells; cp6 QC happens downstream in pipeline.
#   O2 (QC->correct): cp6 QC on ORIGINAL filtered (min_genes>=200 & pct_mt<20), then
#                     correct survivors. (cp6 QC has no total_counts floor, unlike C2.)
# Writes <out_dir>/{matrix.mtx (cells x genes), barcodes.txt, ensembl.txt}.
suppressPackageStartupMessages({
  library(SoupX); library(celda); library(DropletUtils); library(Matrix)
  library(SingleCellExperiment); library(Seurat)
})
a <- commandArgs(trailingOnly = TRUE)
tool <- a[1]; raw_dir <- a[2]; filt_dir <- a[3]; out_dir <- a[4]; ordering <- a[5]
set.seed(42); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
try(RhpcBLASctl::blas_set_num_threads(8), silent = TRUE)

cp6_qc_keep <- function(counts, gene_sym) {           # cp6 QC: min_genes 200 & pct_mt < 20
  n_genes <- Matrix::colSums(counts > 0)
  mito <- grepl("^[Mm][Tt][-.]", gene_sym)
  tot <- Matrix::colSums(counts)
  pct_mt <- ifelse(tot > 0, 100 * Matrix::colSums(counts[mito, , drop = FALSE]) / tot, 0)
  (n_genes >= 200) & (pct_mt < 20)
}
cluster_cells <- function(m) {
  so <- CreateSeuratObject(counts = m); so <- NormalizeData(so, verbose = FALSE)
  so <- FindVariableFeatures(so, verbose = FALSE); so <- ScaleData(so, verbose = FALSE)
  so <- RunPCA(so, npcs = 30, verbose = FALSE); so <- FindNeighbors(so, dims = 1:30, verbose = FALSE)
  so <- FindClusters(so, resolution = 0.5, verbose = FALSE)
  setNames(as.character(Idents(so)), colnames(so))[colnames(m)]
}

filt <- read10xCounts(filt_dir, col.names = TRUE)
gene_id <- rowData(filt)$ID; gene_sym <- rowData(filt)$Symbol
toc_all <- counts(filt); rownames(toc_all) <- gene_id

if (ordering == "O2") {
  keep <- cp6_qc_keep(toc_all, gene_sym); toc <- toc_all[, keep, drop = FALSE]
  cat(sprintf("[%s O2] cp6-QC %d -> %d cells\n", tool, ncol(toc_all), ncol(toc)))
} else { toc <- toc_all; cat(sprintf("[%s O1] correcting all %d cells\n", tool, ncol(toc))) }

if (tool == "SoupX") {
  tod <- counts(read10xCounts(raw_dir, col.names = TRUE)); rownames(tod) <- gene_id
  sc <- SoupChannel(tod, toc); sc <- setClusters(sc, cluster_cells(toc))
  sc <- autoEstCont(sc, doPlot = FALSE); adj <- adjustCounts(sc, roundToInt = TRUE)
} else if (tool == "DecontX") {
  sce <- SingleCellExperiment(assays = list(counts = toc)); sce <- decontX(sce)
  adj <- decontXcounts(sce); adj <- round(adj)             # integer counts for downstream
} else stop("bad tool")

Matrix::writeMM(t(as(adj, "CsparseMatrix")), file.path(out_dir, "matrix.mtx"))  # cells x genes
writeLines(colnames(adj), file.path(out_dir, "barcodes.txt"))
writeLines(gene_id, file.path(out_dir, "ensembl.txt"))
cat(sprintf("[%s %s] DONE wrote %d cells x %d genes\n", tool, ordering, ncol(adj), nrow(adj)))
