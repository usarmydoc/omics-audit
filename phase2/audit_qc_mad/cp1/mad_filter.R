#!/usr/bin/env Rscript
# CP1 MAD filtering via scuttle::isOutlier (Germain 2020 / Heumos method).
# Usage: Rscript mad_filter.R <qc_metrics.tsv> <out.tsv>
# log=TRUE on n_genes + total_counts (lower-tail); raw pct_mt (upper-tail).
# A cell is a MAD outlier if flagged on ANY of the three metrics.
suppressMessages(library(scuttle))
args <- commandArgs(trailingOnly = TRUE)
df <- read.delim(args[1], check.names = FALSE)
out <- args[2]

mad_set <- function(k) {
  lo_genes  <- isOutlier(df$n_genes_by_counts, nmads = k, type = "lower", log = TRUE)
  lo_counts <- isOutlier(df$total_counts,      nmads = k, type = "lower", log = TRUE)
  hi_mito   <- isOutlier(df$pct_counts_mt,     nmads = k, type = "higher", log = FALSE)
  lo_genes | lo_counts | hi_mito           # TRUE = outlier (filtered)
}

res <- data.frame(cell = df$cell,
                  mad3_outlier = mad_set(3),
                  mad5_outlier = mad_set(5))
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%s: %d cells; MAD3 filtered %d, MAD5 filtered %d\n",
            basename(args[1]), nrow(df), sum(res$mad3_outlier), sum(res$mad5_outlier)))
