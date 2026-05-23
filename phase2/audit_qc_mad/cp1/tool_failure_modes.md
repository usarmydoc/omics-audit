# Audit QC-MAD CP1 — tool failure modes

No method failed on any dataset. All 8 datasets produced 4-method outputs.

## MAD log-transform note

scuttle::isOutlier(log=TRUE) on n_genes_by_counts + total_counts (lower tail); pct_counts_mt raw (upper tail). nmads 3 and 5.
