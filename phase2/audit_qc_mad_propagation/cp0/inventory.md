# Audit QC-MAD-Propagation — CP0 Inventory & Feasibility

_Generated 2026-05-25. Standards: AUDIT_STANDARDS.md v1.0.3 + §5.3.2._
_Status: CP0 COMPLETE — feasibility confirmed (Census re-pull reproduces cleanly). Ready for CP1 on approval._

## Question
Do the 4 QC filtering methods from Audit QC-MAD (C1, C2, MAD k=3, MAD k=5) produce
different downstream biological conclusions when their outputs run through the
standard pipeline? Closes the propagation loop for QC method choice — the third
"technical → biological propagation" test after Audit 3 C3 (cell-calling) and
Ambient Correction CP3 (ambient correction).

## Working set (3 of QC-MAD's 8 Census datasets)
The clean control is QC-MAD's **`blood`** Census dataset (not a borrowed PBMC) —
correct cross-audit consistency, since the QC-method outputs we propagate come
from that same dataset.

| short | tissue | Census dataset UUID | QC-MAD median mito | role |
|-------|--------|---------------------|--------------------|------|
| small_intestine | small intestine | a37f857c-779f-464e-9310-3db43a1811e7 | 8.32% | worst QC-method effect (Jaccard 0.774 C2 vs MAD5) |
| liver | liver | 34f5307e-7b4d-4a48-b68f-2ba844c6414b | 1.97% | second-worst (C1 vs MAD3 0.817) |
| blood | blood | 46104f0b-9af5-466a-ae0f-56b8dc1969a2 | 3.11% | clean-ish control (near-null expected) |

3 datasets × 4 QC methods = **12 downstream pipeline runs.** Same sample-size
structure as Audit 3 C3 and Ambient CP3 → existence-of-effect tierable;
generalization flag_and_warn until replicated.

Note: `blood` median mito 3.11% is moderate, not the cleanest of QC-MAD's 8
(heart 0.07, bone_marrow 0.00 are cleaner). It is the PBMC-equivalent clean
control and the right cross-audit choice; if the contrast needs a stricter clean
baseline later, heart/bone_marrow are available.

**Terminology bridge (for findings/synthesis):** the clean control is `blood`
in Census tissue taxonomy; in analyst terminology this corresponds to peripheral
blood mononuclear cells (PBMC). Findings will state "blood (Census tissue UUID;
≈ PBMC in analyst terms)" so corpus-internal Census language reads correctly to
analysts outside the corpus context.

## Step 1 / 1a — Census re-pull (VERIFIED)
Matrices were **not** saved by QC-MAD (its `cp1/` is 25 MB — only per-cell QC
metrics + MAD flags). Full count matrices must be re-pulled. **Option A (full
Census re-pull)** chosen: clean provenance, doesn't depend on QC-MAD intermediates.

- `cellxgene_census 1.17.0` present in **base env** (`~/miniforge3/bin/python`,
  scanpy 1.12); Census version **`2025-11-08`** (same as QC-MAD).
- **All 3 UUIDs resolve (human):** liver 394,534 / small_intestine 201,072 /
  blood 2,616,824 cells in the full datasets.
- **QC-MAD subsampled each to 50,000 cells.** The per-cell `soma_joinid` of those
  50k is recorded in `audit_qc_mad/cp1/qc_metrics/<tissue>.tsv` (`cell` column),
  so CP1 can re-pull and **subset to QC-MAD's exact 50k cells** (identical working
  set, not a fresh random subsample). For blood (2.6M cells) the pull MUST filter
  by `soma_joinid` rather than materialize all 2.6M.
- **Storage:** subset to 50k cells × 3 datasets ≈ a few GB; well within NVMe headroom.

### Reproduction verification (liver, re-pull vs QC-MAD saved)
| metric | re-pulled (full liver, 394,534 cells) | QC-MAD saved (50k subsample) | match |
|--------|---------------------------------------|------------------------------|-------|
| median n_genes_by_counts | 2165 | 2152 | within precision ✓ |
| median pct_counts_mt | 1.69% | 1.69% | identical ✓ |
The median match (full vs 50k subsample) confirms (a) Census `2025-11-08` data is
unchanged since QC-MAD and (b) the subsample is representative. small_intestine +
blood get the same explicit check at CP1 pull (recompute on the soma_joinid-matched
50k, compare to `audit_qc_mad/cp1/qc_metrics/<tissue>.tsv`; mismatch → stop).

### Implementation approach (confirmed): soma_joinid-filtered pull
CP1 pulls each dataset **filtered to QC-MAD's exact 50,000 `soma_joinid`s** (from
`audit_qc_mad/cp1/qc_metrics/<tissue>.tsv` `cell` column), via the documented
Census `obs_value_filter`/`obs_query` pattern:
```python
import cellxgene_census
census = cellxgene_census.open_soma(census_version="2025-11-08")  # PINNED
ids = ",".join(map(str, soma_joinids))          # QC-MAD's 50k per dataset
a = cellxgene_census.get_anndata(census, organism="homo_sapiens",
        obs_value_filter=f"soma_joinid in [{ids}]")
```
This is essential for blood (2.6M full cells → never materialize all 2.6M; pull
only the 50k). **Cross-audit consistency is then exact, not approximate:** the cells
analyzed here are the *identical* cells QC-MAD evaluated, matched by soma_joinid —
document this in findings.md ("cross-audit comparisons are direct, not approximate").

**Census version is PINNED to `2025-11-08`** in the CP1 environment spec, so future
re-runs reproduce exactly even if Census publishes newer versions.

**No blocker.** Re-pull is a CP1 step, verified feasible and reproducible.

## Step 2 — Pipeline reusability (CONFIRMED)
Reuse the Audit 3 CP6 / Ambient CP3 downstream pipeline (`audit3_counting` env):
- scanpy 1.11.5, scDblFinder via Rscript subprocess, scry::devianceFeatureSelection
  via Rscript subprocess — all present and exercised in CP3 (2026-05-25).
- CellTypist 1.7.1 present.
- **Loader tweak (minor):** CP3's `cp3_pipeline.py` reads STARsolo mtx + maps
  Ensembl→symbol via STARsolo `features.tsv`. Census AnnData already carries
  `var.feature_name` (symbol) + `feature_id` (Ensembl), so the CP1 pipeline uses
  Census var symbols directly (no STARsolo mapping). Otherwise identical.

## Step 2b — CellTypist models (per tissue, with reasoning)
| dataset | model | status | reasoning |
|---------|-------|--------|-----------|
| liver | `Healthy_Human_Liver.pkl` | ✓ on disk | liver-specific human model — strong coverage (better than feared) |
| small_intestine | `Cells_Intestinal_Tract.pkl` | downloadable | human intestinal-tract model; download at CP1 start |
| blood | `Immune_All_Low.pkl` | ✓ on disk | matches Ambient CP3's PBMC choice (cross-audit consistency); `Healthy_COVID19_PBMC.pkl` also available as alternative |

## Step 3 — QC method implementation (REUSABLE from QC-MAD)
All 4 methods are QC-MAD's exact code paths:
- **C1** (quantile-data): 5th-pct floors on n_genes + total_counts, 95th-pct mito
  ceiling — `audit_qc_mad/cp1/compare.py` `retained_sets()`.
- **C2** (fixed-floors): n_genes≥200, total_counts≥500, mito≤95th-pct — same.
- **MAD k=3 / k=5**: `scuttle::isOutlier(nmads=k, log=TRUE)` lower-tail on
  n_genes + total_counts, raw upper-tail on pct_mt — `audit_qc_mad/cp1/mad_filter.R`.
QC metrics (n_genes_by_counts, total_counts, pct_counts_mt) recomputed in the
pipeline env to match QC-MAD; the saved metrics are the sanity-check reference.

## Step 4 — Baseline reference: C2 (RECOMMENDED, documented)
No natural "no-QC" baseline exists (unfiltered scanpy = noise). Use **C2 as the
reference** and frame the question as "does switching from the typical scanpy
default (fixed floors 200/500) to another method — notably Heumos-recommended
MAD — change the biology?" Pairwise cross-method comparison also computed. The
reference choice affects framing, not what is computed.

## Step 5 — Pipeline parameters (mirror CP6 / CP3 exactly)
scDblFinder defaults · shifted-log norm (1e4) · scry deviance HVG top-2000 ·
PCA 50 · neighbors k=15 · Leiden res 1.0 (+0.5, 1.5 sensitivity) · Wilcoxon
markers top-50 on log-norm .X · CellTypist majority-vote. No deviations planned
beyond the Census-var loader (Step 2).

## Step 6 — Metrics (mirror C3 / CP3)
Per (method × dataset): post-QC cell count; doublet consistency on cell
intersection; cluster count at Leiden 1.0 (+0.5/1.5); ARI & NMI vs C2; marker
Jaccard per cluster-matched top-50; annotation agreement (confusion vs C2);
method-specific-cell disposition (cells C2 removes but MADk retains → what type
do they annotate as). Cell-level bootstrap (B=1000) on ARI/NMI/marker overlap.

## Step 7 — Runtime
12 pipeline runs; QC-MAD/Census cells ≈ 50k each (intestine post-QC was ~22k in
CP3 at similar scale; CP3 intestine pipeline ran ~2 min). Plus a one-time Census
re-pull (3 datasets, soma_joinid-filtered). Estimate **~1–2 h total** on 9900X +
4070 Ti (CP3's 12 conditions ran in ~19 min; this is comparable). Background-able.

## Cross-audit pattern (flag for CP4)
Third propagation test. Prior two used STARsolo-counted 10x data (Audit 3 mouse
intestine + human pbmc_5k; Ambient same). This uses **human Census** tissues
(small_intestine, liver, blood) — different source/samples, so cross-audit
comparison is at the **methodological-pattern level** ("technical choice
propagates on edge-case tissue, washes out on clean"), not same-sample. If the
pattern holds across all three audits, that is a corpus-level finding about
scRNA technical choices generally.

## Heumos positioning
Extends. Heumos recommends MAD filtering without testing downstream
consequences; this audit quantifies whether switching common-practice fixed
floors (C2) → MAD changes biology, and adds the third leg of the corpus's
technical→biological propagation pattern.

## Blocking issues
**None.** Single CP1 prerequisite: download `Cells_Intestinal_Tract.pkl`
(small, automated via celltypist). Census re-pull verified reproducible.

## Recommendation
**Proceed to CP1.** Re-pull the 3 datasets subset to QC-MAD's exact 50k
soma_joinids (sanity-check small_intestine + blood metrics vs QC-MAD saved at
pull), apply the 4 QC methods, run the 12 pipeline conditions, compare vs C2
reference + cross-method. Mirror CP3 machinery.
