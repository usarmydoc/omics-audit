# Superseded CP4 outputs — buggy USA-suffix stripping (2026-05-18)

These outputs were produced by the first CP4 run and contained a bug in
`strip_version()`: `g.split(".")[0]` stripped both the Ensembl `.version`
AND any trailing `-U`/`-A` USA-mode suffix from alevin-fry's splici output.

This caused alevin-fry's three per-gene rows (`ENSG.version`, `ENSG.version-U`,
`ENSG.version-A` = Spliced/Unspliced/Ambiguous) to collapse to a single dict
key, of which only the last-inserted (the `-A` Ambiguous bucket) survived.

The pair-metrics functions then compared alevin-fry's "ambiguous-only" reads
(~10% of total) against STARsolo's full reads (~100%), producing artifactual
rho_gene_total ~0.58, log2_median = -7.96 on mitochondrial genes (the -A
bucket has ~0 MT reads), and the spurious "alevin-fry is the outlier" headline.

Per AUDIT_STANDARDS §1.5, files retained here for reproducibility.
Replaced 2026-05-18 by run with corrected USA collapse: S + A summed per
gene, U dropped (per alevin-fry conventions for abundance-style counting).

Buggy alevin-fry vs STARsolo: rho ≈ 0.52
Fixed alevin-fry vs STARsolo:  rho ≈ 0.95 (pbmc_1k_v3 verified)

See also `c1/verify_1779151834/` (single-dataset verification run).
