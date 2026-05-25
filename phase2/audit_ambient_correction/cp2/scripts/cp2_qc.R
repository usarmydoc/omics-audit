# CP2 shared C2 QC (fixed-floor, from QC-MAD audit). Sourced by SoupX/DecontX scripts.
# Remove cells with: n_genes < 200  OR  total_counts < 500  OR  mito% > 95th percentile.
# 95th percentile computed on the matrix being filtered (corrected for O1, original for O2).
# Returns logical keep-mask over columns (cells). gene_symbols used for mito detection
# (human MT-, mouse mt-), identical convention to CP1.
c2_qc_mask <- function(counts, gene_symbols) {
  total_counts <- Matrix::colSums(counts)
  n_genes      <- Matrix::colSums(counts > 0)
  mito_idx     <- grepl("^[Mm][Tt]-", gene_symbols)
  mito_counts  <- if (any(mito_idx)) Matrix::colSums(counts[mito_idx, , drop = FALSE]) else rep(0, ncol(counts))
  mito_pct     <- ifelse(total_counts > 0, 100 * mito_counts / total_counts, 0)
  p95          <- as.numeric(stats::quantile(mito_pct, 0.95, na.rm = TRUE))
  keep <- (n_genes >= 200) & (total_counts >= 500) & (mito_pct <= p95)
  attr(keep, "mito_p95") <- p95
  attr(keep, "n_mito_genes") <- sum(mito_idx)
  keep
}
