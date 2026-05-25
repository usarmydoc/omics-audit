#!/usr/bin/env Rscript
# CP1 Deliverable A — SoupX per-gene contamination
# Usage: Rscript cp1_soupx.R <raw_dir> <filt_dir> <out_prefix>
#   raw_dir/filt_dir: STARsolo Gene/{raw,filtered} CellRanger-format dirs
#   out_prefix: writes <prefix>_pergene.tsv and <prefix>_summary.json
#
# Per-gene contamination defined uniformly across all 3 tools as the fraction
# of each gene's counts removed as ambient on the tool's native cell set:
#   contam_g = (orig_total_g - corrected_total_g) / orig_total_g
# SoupX requires clusters for autoEstCont(); we compute them with a standard
# minimal Seurat pipeline (clusters are a required input, not a tunable param).

suppressPackageStartupMessages({
  library(SoupX); library(DropletUtils); library(Matrix)
  library(Seurat); library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE)
raw_dir <- args[1]; filt_dir <- args[2]; out_prefix <- args[3]
set.seed(42)

cat("[SoupX] reading matrices\n")
tod_sce <- read10xCounts(raw_dir,  col.names = TRUE)   # raw (table of droplets)
toc_sce <- read10xCounts(filt_dir, col.names = TRUE)   # filtered (table of counts = cells)
tod <- counts(tod_sce); toc <- counts(toc_sce)
gene_id  <- rowData(toc_sce)$ID
gene_sym <- rowData(toc_sce)$Symbol
rownames(tod) <- rowData(tod_sce)$ID; rownames(toc) <- gene_id

# SoupChannel requires tod genes ⊇ toc genes, same ordering on shared genes.
stopifnot(identical(rownames(tod), rownames(toc)))
cat(sprintf("[SoupX] tod %d x %d ; toc %d x %d\n", nrow(tod), ncol(tod), nrow(toc), ncol(toc)))

# --- clusters for autoEstCont via standard minimal Seurat pipeline ---
cat("[SoupX] clustering filtered cells (Seurat, defaults)\n")
so <- CreateSeuratObject(counts = toc)
so <- NormalizeData(so, verbose = FALSE)
so <- FindVariableFeatures(so, verbose = FALSE)
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, npcs = 30, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:30, verbose = FALSE)
so <- FindClusters(so, resolution = 0.5, verbose = FALSE)
clusters <- setNames(as.character(Idents(so)), colnames(so))
clusters <- clusters[colnames(toc)]
cat(sprintf("[SoupX] %d clusters\n", length(unique(clusters))))

# --- SoupX ---
sc <- SoupChannel(tod, toc)
sc <- setClusters(sc, clusters)
sc <- autoEstCont(sc, doPlot = FALSE)          # default contamination estimation
rho <- sc$metaData$rho
adj <- adjustCounts(sc)                          # corrected counts (genes x cells)

# --- per-gene contamination on filtered cell set ---
orig_total <- Matrix::rowSums(toc)
corr_total <- Matrix::rowSums(adj)
df <- data.frame(gene_id = gene_id, gene_symbol = gene_sym,
                 orig_total = as.numeric(orig_total),
                 corr_total = as.numeric(corr_total),
                 stringsAsFactors = FALSE)
write.table(df, paste0(out_prefix, "_pergene.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

summ <- list(tool = "SoupX", version = as.character(packageVersion("SoupX")),
             n_cells = ncol(toc), n_genes = nrow(toc),
             n_clusters = length(unique(clusters)),
             rho_mean = mean(rho), rho_median = median(rho),
             rho_min = min(rho), rho_max = max(rho),
             global_frac_removed = 1 - sum(corr_total) / sum(orig_total))
writeLines(toJSON(summ, auto_unbox = TRUE, pretty = TRUE),
           paste0(out_prefix, "_summary.json"))
cat(sprintf("[SoupX] DONE rho_mean=%.4f global_frac_removed=%.4f\n",
            mean(rho), summ$global_frac_removed))
