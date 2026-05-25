# CP2 — Tool-design constraints, sanity checks, anomalies

_Deliverable B, 2026-05-24. 9 datasets × 3 tools × 2 orderings. All runs exit 0; fails=0; no halt._

## CellBender ordering-invariance — a tool-design constraint (not a failure)
CellBender's variational inference learns ambient from the **raw droplet distribution**; it cannot be refit on a QC-passed cell subset (no empty droplets → no ambient model). Therefore its correction is **ordering-invariant by construction**: O1 (correct→QC) and O2 (QC→correct) reuse the **same** CP1 corrected matrix and differ **only in the QC basis** (QC metrics computed on corrected counts vs original counts). Re-running CellBender per ordering would yield (stochastically near-)identical output at large GPU cost; reuse is the methodologically correct choice, documented in findings. SoupX and DecontX estimate from the cell-by-gene matrix and were fully re-run for both orderings.

## Sanity check: CellBender O1 vs O2 cell sets must differ (non-trivially)
Because the corrected matrix is identical across CellBender orderings, the only thing that can differ is the QC-survivor cell set (QC on corrected vs original counts). **If O1 and O2 survivor sets were literally identical, the QC-basis switch would be broken.** Check result: **survivor sets differ (non-identical) on all 9 datasets** — O2-only cells range 22 (pbmc_1k) to 2,347 (lung). The QC basis switch operates correctly. CellBender's `mito_p95` also differs by basis (e.g. pbmc_1k: O1 corrected 98.7% vs O2 original 95.3%), confirming distinct metric computation.

## Cell-retention asymmetry (real effect, not a bug)
O2 (QC-on-original) retains more cells than O1 (QC-on-corrected) on every dataset, because correction lowers counts below the C2 floor (200 genes / 500 counts). Largest on high-ambient tissue (CellBender lung O2-only=2,347; kidney 823). Expected from the mechanism; reported as a finding, not flagged as anomaly.

## DecontX ordering-sensitivity (real, tool-design)
DecontX's per-gene contamination changes between orderings (pooled Spearman v2 0.813 / v3 0.855; kidney as low as 0.377) because its mixture model refits on the QC-passed subpopulation. This is genuine tool behavior, not a configuration error — DecontX estimates contamination from the cell matrix, so the input cell set changes the estimate. Magnitudes remained plausible across orderings (no degenerate 0%/100% outputs); no §3.5 stop triggered.

## Gene-ID integrity
All tools emit the identical STARsolo Ensembl ID universe across both orderings (C1 USA-mode check). CellBender symbols sourced from STARsolo `features.tsv` (canonical map) as in CP1. No mismatch.

## QC retention rates
C2 retention 0.92–1.00 across all (dataset × tool × ordering) — sensible, never 0% or 100%. No degenerate QC behavior.

## Run robustness
Sequential smallest→largest, SoupX∥DecontX∥both-orderings concurrent (4 R procs), CellBender QC-only reuse. Resumable (per (tool,ordering) summary sentinel). Halt-on-2-consecutive-dataset-failures armed; not triggered. Peak RAM ~9–28 GB / 60 GB. Total ~13 min.
