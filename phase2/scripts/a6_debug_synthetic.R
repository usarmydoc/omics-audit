#!/usr/bin/env Rscript
# Phase 2a A6 debug — confirm hypothesis on tiny synthetic data.
# Tests: (1) is the NEBULA bug really intercept-vs-group? (2) what does mmDS actually error with?

suppressPackageStartupMessages({
  library(Matrix)
  library(nebula)
  library(muscat)
  library(SingleCellExperiment)
})

set.seed(42)

# --- build small synthetic dataset ---
# 400 genes x 600 cells; 6 donors (3 disease + 3 normal); 100 cells per donor
n_genes <- 400
n_donors <- 6
cells_per_donor <- 100
n_cells <- n_donors * cells_per_donor

donor_ids <- rep(paste0("D", 1:n_donors), each = cells_per_donor)
group_ids <- rep(c("normal", "normal", "normal", "disease", "disease", "disease"), each = cells_per_donor)

# Baseline expression: NB with mean varying by gene
gene_means <- rgamma(n_genes, shape = 2, scale = 1)  # baseline log-mean
# Make 40 genes truly differentially expressed (10%)
de_genes <- sample.int(n_genes, 40)
effect <- rep(0, n_genes)
effect[de_genes] <- sample(c(-1, 1), 40, TRUE) * runif(40, 1.0, 2.0)  # log2 effect

# Library size variation per cell
lib <- exp(rnorm(n_cells, log(5000), 0.3))

mat <- matrix(0L, n_genes, n_cells)
for (j in seq_len(n_cells)) {
  group_eff <- if (group_ids[j] == "disease") effect else rep(0, n_genes)
  # donor random effect
  donor_eff <- rnorm(n_genes, 0, 0.1)
  mu <- exp(log(gene_means) + group_eff * log(2) + donor_eff) * (lib[j] / 5000)
  mat[, j] <- rnbinom(n_genes, mu = mu, size = 2)
}
rownames(mat) <- paste0("g", seq_len(n_genes))
colnames(mat) <- paste0("c", seq_len(n_cells))
counts <- as(mat, "CsparseMatrix")
cat(sprintf("[synth] %d genes x %d cells, %d donors, %d true DE genes\n",
            n_genes, n_cells, n_donors, length(de_genes)))

group <- factor(group_ids, levels = c("normal", "disease"))
donor <- factor(donor_ids)

# ---- NEBULA — show ALL columns of result ----
cat("\n========== NEBULA ==========\n")
design <- model.matrix(~ group)
offset_vec <- log(Matrix::colSums(counts) + 1)
id_vec <- as.numeric(donor)

# group_cell returns NULL when data is already sorted by id; handle both
data_g <- group_cell(count = counts, id = id_vec,
                     pred = design, offset = offset_vec)
if (is.null(data_g)) {
  cat("[nebula] group_cell returned NULL (cells already grouped) — using inputs as-is\n")
  fit <- nebula(count = counts, id = id_vec, pred = design,
                offset = offset_vec, method = "LN", verbose = FALSE)
} else {
  cat("[nebula] group_cell reordered cells — using grouped outputs\n")
  fit <- nebula(count = data_g$count, id = data_g$id, pred = data_g$pred,
                offset = data_g$offset, method = "LN", verbose = FALSE)
}
cat("summary columns:", paste(colnames(fit$summary), collapse = ", "), "\n")
cat("first 3 rows:\n")
print(head(fit$summary, 3))

# What the OLD code does (BUG)
old_coef_col <- grep("^logFC_", colnames(fit$summary), value = TRUE)[1]
old_pval_col <- grep("^p_", colnames(fit$summary), value = TRUE)[1]
cat(sprintf("\nOLD code picked: logFC=%s, pval=%s\n", old_coef_col, old_pval_col))
old_padj <- p.adjust(fit$summary[[old_pval_col]], method = "BH")
cat(sprintf("OLD: n_sig=%d/%d, p=0 count=%d, |log2FC|>5 count=%d\n",
            sum(old_padj < 0.05, na.rm = TRUE), nrow(fit$summary),
            sum(fit$summary[[old_pval_col]] == 0, na.rm = TRUE),
            sum(abs(fit$summary[[old_coef_col]] / log(2)) > 5, na.rm = TRUE)))

# What the FIXED code should do — pick the GROUP coefficient, not intercept
group_logfc_cols <- grep("^logFC_", colnames(fit$summary), value = TRUE)
group_pval_cols <- grep("^p_", colnames(fit$summary), value = TRUE)
# Drop intercept
group_logfc_cols <- setdiff(group_logfc_cols, "logFC_(Intercept)")
group_pval_cols <- setdiff(group_pval_cols, "p_(Intercept)")
new_coef_col <- group_logfc_cols[1]
new_pval_col <- group_pval_cols[1]
cat(sprintf("\nFIXED code picks: logFC=%s, pval=%s\n", new_coef_col, new_pval_col))
new_padj <- p.adjust(fit$summary[[new_pval_col]], method = "BH")
n_sig_new <- sum(new_padj < 0.05, na.rm = TRUE)
cat(sprintf("FIXED: n_sig=%d/%d (true=%d)\n", n_sig_new, nrow(fit$summary), length(de_genes)))

# Confirm: do the top-hit genes overlap the true DE set?
top_idx <- order(fit$summary[[new_pval_col]])[1:40]
top_genes <- fit$summary$gene[top_idx]
true_genes <- paste0("g", de_genes)
overlap <- length(intersect(top_genes, true_genes))
cat(sprintf("FIXED: top-40 NEBULA hits overlap with true DE: %d/40 (recovery rate)\n", overlap))

# ---- muscat mmDS — drill into the failure mode ----
cat("\n========== muscat mmDS ==========\n")
sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = DataFrame(sample_id = donor, group_id = group, cluster_id = factor("all"))
)
sce <- prepSCE(sce, kid = "cluster_id", gid = "group_id", sid = "sample_id", drop = TRUE)
cat("SCE built:", nrow(sce), "x", ncol(sce), "\n")
cat("metadata$experiment_info:\n")
print(metadata(sce)$experiment_info)

cat("\n[mmDS attempt 1] method='dream', BPPARAM=SerialParam (n_threads removed in muscat 1.24)\n")
mm_result <- tryCatch({
  suppressPackageStartupMessages(library(BiocParallel))
  res_mm <- mmDS(sce, method = "dream", verbose = TRUE,
                 BPPARAM = MulticoreParam(2, progressbar = FALSE))
  cat("[mmDS] SUCCESS\n")
  cat("class:", class(res_mm), "\n")
  cat("names:", paste(names(res_mm), collapse = ", "), "\n")
  if (!is.null(res_mm$all)) {
    cat("res_mm$all head:\n")
    print(head(res_mm$all, 3))
  }
  res_mm
}, error = function(e) {
  cat("[mmDS] ERROR class:", paste(class(e), collapse = ", "), "\n")
  cat("[mmDS] ERROR message:", conditionMessage(e), "\n")
  if (!is.null(e$call)) cat("[mmDS] ERROR call:", deparse(e$call), "\n")
  NULL
})

if (!is.null(mm_result)) {
  cat("\n[mmDS post] examining output structure\n")
  cat("class:", paste(class(mm_result), collapse = ", "), "\n")
  cat("names:", paste(names(mm_result), collapse = ", "), "\n")
  if (!is.null(mm_result$table)) {
    cat("table contrasts:", paste(names(mm_result$table), collapse = ", "), "\n")
    tab <- mm_result$table[[1]][[1]]
    cat("first table head:\n")
    print(head(tab, 3))
    cat("table columns:", paste(colnames(tab), collapse = ", "), "\n")
    n_sig_mm <- sum(p.adjust(tab$p_val, method = "BH") < 0.05, na.rm = TRUE)
    cat(sprintf("mmDS dream sig (BH<0.05): %d/%d (true=40)\n", n_sig_mm, nrow(tab)))
  }
}

cat("\n[done]\n")
