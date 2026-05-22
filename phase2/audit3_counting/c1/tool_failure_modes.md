# Audit 3 C1 — Tool failure modes

_Generated: 2026-05-18T20:34:29_

## Per-dataset notes

All 9 datasets produced 4-of-4 tool outputs.

## Per-tool cell count compared to whitelist

Tools report different barcode universes. STARsolo's raw output includes the entire 10x whitelist (~6.8M barcodes). alevin-fry and kb_count report only barcodes that received reads. See `n_common_cells` in per_dataset_metrics.tsv for the intersection.

## kb_count smaller output sizes

CP3 surfaced that kb_count produced systematically smaller output files than
STAR_default, STAR_cr_mimic, and alevin_fry. Investigation during CP4 resolved
this as a barcode-universe difference, not a count or detection issue.

kb_count's default barcode whitelist filtering produces a smaller cell-by-gene
matrix because fewer barcodes pass kb_count's filtering, but the cells that
are retained have similar per-gene counts to other tools. kb_count detects
more genes per retained cell, not fewer. Smaller output size therefore
reflects narrower barcode retention, not data loss.

Empirical support (CP4 per_dataset_metrics.tsv): across all 9 datasets,
kb_count's gene detection rate (0.42–0.69) equals or exceeds its partner
tools (0.39–0.68). The smaller on-disk file is the `counts_unfiltered`
matrix holding only barcodes that received reads (~240K–1.4M) versus
STARsolo's raw output spanning the full whitelist (6.8M barcodes).

Implication for C2 (cell barcode calling agreement): kb_count's barcode
calling is more restrictive at default settings than the other three tools.
This should appear in C2's Jaccard metrics on called cell barcode sets.

## USA-mode suffix bug in alevin-fry

During CP4 analysis, an inconsistency surfaced in gene identifier handling
between alevin-fry outputs and the other three tools. alevin-fry's USA-mode
appends suffixes to gene IDs (e.g., ENSG00000123456-U for unspliced reads,
ENSG00000123456-S for spliced, ENSG00000123456-A for ambiguous) to
distinguish read categories. Comparison across tools requires collapsing
these suffixes back to the base gene ID before computing per-gene agreement
metrics (summing the relevant buckets — S+A for STARsolo-comparable
abundance — and dropping U).

The initial CP4 run did not collapse USA-mode suffixes, producing apparent
tool disagreement that was actually a gene-ID matching artifact (the
`strip_version()` helper collapsed the -U/-A suffix along with the Ensembl
version, leaving only the last-inserted -A "ambiguous" bucket — roughly 10%
of total reads — to stand in for alevin-fry's per-gene count). Bootstrap
CIs were tight on this buggy data, which is itself a methodological finding:
tight CIs reflect sampling variance, not correctness.

After collapsing USA-mode suffixes (summing S+A per gene), cross-tool
Spearman correlation rose from ~0.52 (pbmc_1k_v3; ~0.58 median across
datasets in the buggy pass) to ~0.96, with overlapping CIs across all
tool pairs. This is the headline finding of C1: tools converge on per-gene
counts when gene ID conventions are correctly harmonized.

Methodological implication: future audits should verify gene ID convention
consistency before computing per-gene comparisons. Tight bootstrap CIs on
data with systematic configuration errors can produce misleading high-
confidence findings. Buggy outputs preserved at
`superseded_2026-05-18_buggy_usa_strip/`; corrected analysis in the current
c1/ outputs.
