# CP3 — pipeline behavior, sanity checks, anomalies

_Deliverable C, 2026-05-25. 2 datasets × 6 conditions = 12 pipeline runs (5 corrected + reused baseline per dataset). All exit 0; no failures._

## Baseline sanity (matches Audit 3 C3)
The no-correction baseline is C3's `star_default` output, run through the same cp6 pipeline. Reference values reproduced exactly (intestine postQC 21,807 / 21 clusters / 14 labels; pbmc_5k 4,520 / 16 / 14) — confirmed before use as the comparison reference.

## CellBender cell loss on intestine (real, not a failure)
CellBender's intestine postQC dropped to 16,700 (from 21,807 baseline) because its ~35% ambient removal lowers per-cell counts below the cp6 QC floor (min_genes 200). This is the "correct→QC is stricter" effect (CP2) manifesting biologically, not a pipeline error. SoupX (floors, ~no count change) keeps 21,807; DecontX ~20,600.

## CellTypist sanity per condition
Label counts are sensible and stable across conditions (intestine 13–15 labels, pbmc 14–18) — no silent annotation degradation (e.g. collapse to 1 label or explosion). Re-annotated cells concentrate in related subtypes (T/NK compartment on PBMC), not random reassignment.

## §3.5 tight-CI check
ARI bootstrap CIs are narrow (±0.01) due to large cell N. This is sampling precision, not correctness. The underlying corrected matrices were validated in CP1 (gene-ID integrity, plausible magnitudes); the PBMC-vs-intestine effect bands (0.85–0.90 vs 0.50–0.70) are widely separated, so the contrast is not a CI artifact.

## Determinism / reproducibility
All steps seeded (scanpy random_state=0, R set.seed(42), Leiden igraph deterministic). CellBender correction reused from CP1 (ordering-invariant). scDblFinder/scry run as Rscript subprocesses (no rpy2), matching cp6.

## QC-definition note (deliberate)
CP3 uses cp6 QC (min_genes 200, max_mito 20%), not CP2's C2 fixed-floor, so the no-correction baseline (built with cp6 QC in C3) is directly comparable. The O1/O2 ordering axis is preserved as **whether each tool's correction was estimated on all filtered cells (O1) or cp6-QC-passed survivors (O2)** — the CP2 axis that mattered for DecontX. Documented in findings.

## Run robustness
Sequential pipelines (RAM-bound), matrix-gen SoupX∥DecontX concurrent, resumable (skip on existing summary). Total ~19 min. No halt; fails=0.
