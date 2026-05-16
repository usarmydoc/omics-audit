#!/usr/bin/env Rscript
# Phase 2a A6 — NEBULA + muscat-dream DE on per-cell count matrix.
#
# Reads:
#   counts.mtx.gz     — gene x cell sparse counts
#   meta.tsv          — cell metadata: cell_id, donor_id, group {disease/normal}
#   gene_names.tsv    — gene symbols matching count rows
# Writes:
#   <out_dir>/<dataset>__nebula.tsv
#   <out_dir>/<dataset>__muscatpb.tsv     (muscat::pbDS — pseudobulk DESeq2 reference)
#   <out_dir>/<dataset>__muscatmm.tsv     (muscat::mmDS — mixed model)
#
# Per-gene unified schema: gene_id, log2FC, pvalue, padj, tool, comparison_id

suppressPackageStartupMessages({
  library(Matrix)
  library(nebula)
  library(muscat)
  library(SingleCellExperiment)
  library(scater)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) stop("Usage: Rscript a6_nebula_muscat.R <counts.mtx.gz> <meta.tsv> <genes.tsv> <out_dir> <dataset_label>")
counts_path <- args[1]
meta_path   <- args[2]
genes_path  <- args[3]
out_dir     <- args[4]
ds_label    <- args[5]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("[load] reading", counts_path, "\n")
counts <- readMM(gzfile(counts_path))
genes <- read.table(genes_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)$V1
meta <- read.table(meta_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rownames(counts) <- genes
colnames(counts) <- meta$cell_id
cat("[load] ", nrow(counts), " genes x ", ncol(counts), " cells; donors=",
    length(unique(meta$donor_id)), ", groups=",
    paste(unique(meta$group), collapse = "/"), "\n", sep = "")

# Drop very low-expression genes
keep <- Matrix::rowSums(counts > 0) >= 10
counts <- counts[keep, ]
genes <- genes[keep]
cat("[filter] ", nrow(counts), " genes after expression filter\n", sep = "")

# De-duplicate gene names. Census top-K selection can produce duplicated symbols
# (paralog annotations sharing names); muscat::mmDS sets rownames(SCE)=gene_id
# and crashes with "duplicate 'row.names' are not allowed" otherwise.
n_dups <- sum(duplicated(genes))
if (n_dups > 0) {
  cat("[filter] de-duplicating ", n_dups, " gene names via make.unique\n", sep = "")
  genes <- make.unique(genes)
  rownames(counts) <- genes
}

group <- factor(meta$group, levels = c("normal", "disease"))
donor <- factor(meta$donor_id)
comp_id <- paste0("Census_", ds_label, "_DiseaseVsNormal")

write_unified <- function(tool, df) {
  out <- data.frame(
    gene_id = df$gene_id,
    log2FC  = df$log2FC,
    pvalue  = df$pvalue,
    padj    = df$padj,
    tool    = tool,
    comparison_id = comp_id,
    stringsAsFactors = FALSE
  )
  p <- file.path(out_dir, paste0(ds_label, "__", tool, ".tsv"))
  write.table(out, file = p, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("[write] ", tool, ": ", nrow(out), " genes -> ", p, "\n", sep = "")
}

# ---- NEBULA ----
cat("\n[nebula] starting\n")
t0 <- Sys.time()
tryCatch({
  design <- model.matrix(~ group)
  offset_vec <- log(Matrix::colSums(counts) + 1)
  id_vec <- as.numeric(donor)
  # group_cell reorders cells by id; returns NULL when already grouped
  data_g <- group_cell(count = counts, id = id_vec, pred = design, offset = offset_vec)
  if (is.null(data_g)) {
    fit <- nebula(count = counts, id = id_vec, pred = design,
                  offset = offset_vec, method = "LN", verbose = FALSE, ncore = 4)
  } else {
    fit <- nebula(count = data_g$count, id = data_g$id, pred = data_g$pred,
                  offset = data_g$offset, method = "LN", verbose = FALSE, ncore = 4)
  }
  res <- fit$summary
  # Pick the GROUP coefficient explicitly (not intercept). Bug fix 2026-05-16.
  group_logfc_cols <- setdiff(grep("^logFC_", colnames(res), value = TRUE), "logFC_(Intercept)")
  group_pval_cols  <- setdiff(grep("^p_",     colnames(res), value = TRUE), "p_(Intercept)")
  if (length(group_logfc_cols) == 0 || length(group_pval_cols) == 0) {
    stop("nebula returned no non-intercept coefficient columns; cols=",
         paste(colnames(res), collapse = ","))
  }
  coef_col <- group_logfc_cols[1]
  pval_col <- group_pval_cols[1]
  cat("[nebula] using coef=", coef_col, " pval=", pval_col, "\n", sep = "")
  res_df <- data.frame(
    gene_id = res$gene,
    log2FC  = res[[coef_col]] / log(2),  # NEBULA returns natural log; convert
    pvalue  = res[[pval_col]],
    padj    = p.adjust(res[[pval_col]], method = "BH"),
    stringsAsFactors = FALSE
  )
  write_unified("nebula", res_df)
  cat("[nebula] done in ", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n", sep = "")
}, error = function(e) {
  cat("[nebula] FAILED: ", conditionMessage(e), "\n", sep = "")
})

# ---- muscat: convert to SCE first ----
cat("\n[muscat] building SCE\n")
sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = DataFrame(sample_id = donor, group_id = group, cluster_id = factor("all"))
)
sce <- prepSCE(sce, kid = "cluster_id", gid = "group_id", sid = "sample_id", drop = TRUE)

# ---- muscat pbDS (pseudobulk DESeq2 — Phase 1 P4 method, reference) ----
# Bug fix 2026-05-16: muscat 1.24's aggregateData doesn't auto-propagate
# group_id into pb's colData, but pbDS asserts pbs[["group_id"]]. The pb
# colData typically only contains the aggregation grouping vars (cluster_id),
# while sample_id and group_id live in metadata(pb)$experiment_info.
# Diagnostic print + map experiment_info sample_id→group_id back into pb's
# colData using colnames(pb) which match sample_id by muscat convention.
cat("\n[muscat_pb] pseudobulk DESeq2 via muscat\n")
t0 <- Sys.time()
tryCatch({
  pb <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))
  cat("[muscat_pb] pb dims:", nrow(pb), "x", ncol(pb), "\n")
  cat("[muscat_pb] pb colData cols:", paste(colnames(colData(pb)), collapse = ","), "\n")
  cat("[muscat_pb] pb colnames head:", paste(head(colnames(pb), 3), collapse = ","), "\n")
  ei <- metadata(pb)$experiment_info
  cat("[muscat_pb] experiment_info cols:", paste(colnames(ei), collapse = ","), "\n")
  if (is.null(pb$group_id) && !is.null(ei) && !is.null(ei$group_id) && !is.null(ei$sample_id)) {
    # muscat convention: colnames(pb) are sample_ids (per-cluster aggregation).
    pb_samples <- colnames(pb)
    new_groups <- as.character(ei$group_id[match(pb_samples, as.character(ei$sample_id))])
    cat("[muscat_pb] mapping", length(pb_samples), "samples ->", sum(!is.na(new_groups)), "groups matched\n")
    if (sum(!is.na(new_groups)) == length(pb_samples)) {
      pb$group_id <- factor(new_groups, levels = levels(ei$group_id))
    } else {
      stop("could not map all pb columns to experiment_info group_id")
    }
  }
  res_pb <- pbDS(pb, method = "DESeq2", verbose = FALSE)
  tab <- res_pb$table$disease$all
  res_df <- data.frame(
    gene_id = tab$gene,
    log2FC  = tab$logFC,
    pvalue  = tab$p_val,
    padj    = tab$p_adj.loc,
    stringsAsFactors = FALSE
  )
  write_unified("muscat_pb_DESeq2", res_df)
  cat("[muscat_pb] done in ", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n", sep = "")
}, error = function(e) {
  cat("[muscat_pb] FAILED: ", conditionMessage(e), "\n", sep = "")
})

# ---- muscat mmDS (mixed-model DREAM) ----
# Bug fix 2026-05-16: muscat 1.24 removed n_threads; uses BPPARAM. Output column
# names are p_val/p_adj.loc (muscat schema), not P.Value/adj.P.Val (limma schema).
cat("\n[muscat_mm] mixed-model dream\n")
t0 <- Sys.time()
tryCatch({
  res_mm <- mmDS(sce, method = "dream", verbose = FALSE,
                 BPPARAM = MulticoreParam(6, progressbar = FALSE))
  tab <- res_mm$all
  res_df <- data.frame(
    gene_id = tab$gene,
    log2FC  = tab$logFC,
    pvalue  = tab$p_val,
    padj    = tab$p_adj.loc,
    stringsAsFactors = FALSE
  )
  write_unified("muscat_mm_dream", res_df)
  cat("[muscat_mm] done in ", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n", sep = "")
}, error = function(e) {
  cat("[muscat_mm] FAILED: ", conditionMessage(e), "\n", sep = "")
})

cat("\n[done] ", ds_label, "\n", sep = "")
