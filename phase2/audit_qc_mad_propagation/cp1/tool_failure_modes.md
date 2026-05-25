# CP1 — failure modes, sanity checks, anomalies

_QC-MAD-Propagation CP1, 2026-05-25. 12 conditions, 0 pipeline failures._

## Census re-pull: cell-identity bug caught by sanity check
QC-MAD's `qc_metrics/<tissue>.tsv` `cell` column is the **0-based positional index** into the dataset_id-filtered, soma_joinid-ordered full pull (a leftover RangeIndex after QC-MAD's `rng.choice` subsample) — **NOT** the global Census `soma_joinid`. The first pull used `obs_coords=cell` directly and silently grabbed the wrong 50k cells (right count, wrong identity; median n_genes 1021 vs expected 2152). The CP0/CP1 sanity check (recompute median QC metrics vs QC-MAD saved) caught it. Fix: ordered obs read → map QC-MAD positions to global soma_joinids → `obs_coords` those. Verified exact: liver 2152, small_intestine 639, blood 1517 median n_genes all match QC-MAD to precision.

## QC method reproduction: exact
All 4 methods reproduced QC-MAD's saved per-method cell counts EXACTLY across all 3 datasets (e.g. small_intestine C2=38700, MAD5=50000; liver MAD3=40965; etc.). Confirms the propagation builds on QC-MAD's identical QC decisions.

## Pipeline: QC step skipped by design (not a failure)
The QC method IS the cell filter, so the downstream pipeline SKIPS its own cell QC (min_genes/max_mito) and runs from gene-filter (min_cells=3) onward on the method's survivors. Otherwise identical to CP6/CP3. Census `var.feature_name` supplies symbols (no STARsolo mapping). 12/12 runs exit 0.

## §3.5 tight-CI check
ARI bootstrap CIs are narrow (±0.005) due to large cell N. Sanity: the underlying QC cell sets reproduce QC-MAD exactly and metrics are within tissue expectations; the modest ARI (0.80–0.91) is precisely estimated, not a CI artifact. Cross-method intersections are non-trivial (38.7k–47k shared cells), so comparisons are informative.

## Orchestration incidents (operator, not data)
- **Thread oversubscription.** First parallel attempt ran 3 conditions concurrently with UNCAPPED threads (scanpy/BLAS/numba each grabbed all 24 cores) → load 56, thrashing, forced a manual kill. Fixed: cap OMP/OPENBLAS/MKL/NUMBA_NUM_THREADS=7, MAXJOBS=3 → 3×7=21 cores, load controlled (~7–21), RAM ~30/60 GB. (See feedback memory.)
- **Redundant bootstrap.** The compare script bootstrapped twice (inline per-pair CIs + a redundant pooled per_stratum pass) → ~36k ARI calls, long single-threaded tail. Slimmed: bootstrap once (inline), derive `per_stratum_bootstrap.tsv` from those CIs.
Neither affected results — both are runtime/operator issues, fixed.

## Determinism
scanpy random_state=0, Leiden igraph deterministic, R set.seed(42), CellBender N/A. Reproducible.
