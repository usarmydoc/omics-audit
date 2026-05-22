# Audit 3 C1 — Tool failure modes

_Generated: 2026-05-18T19:51:04_

## Per-dataset notes

All 9 datasets produced 4-of-4 tool outputs.

## Per-tool cell count compared to whitelist

Tools report different barcode universes. STARsolo's raw output includes the entire 10x whitelist (~6.8M barcodes). alevin-fry and kb_count report only barcodes that received reads. See `n_common_cells` in per_dataset_metrics.tsv for the intersection.
