# Audit QC-MAD (low-quality cell filtering: MAD vs quantile) — Status

Canonical drive: `/mnt/nvme1/omics-audit/`. Lock: `phase2/repro.lock`.
Extends Phase 1 p1 (mito-threshold QC). Heumos 2023 positioning: extends.

## AUDIT STATUS: CLOSED (2026-05-23)

- **Checkpoints:** CP0–CP4 all complete (work window 2026-05-23; extends Phase 1 p1).
- **Datasets:** 8 Census (subset of Phase 1 p1, census 2025-11-08), 8 tissues
  incl. blood/bone_marrow (PBMC-like) + high-mito GI/pancreas, human + mouse.
- **Methods:** 4 (C1 pure-quantile, C2 fixed-floor, MAD k=3, MAD k=5).
- **Rules contributed:** 2 (both conditional, §5.3.2-tiered).
- **Standards contributed:** none (used existing §5.3.2 amendment).
- **Lock:** 13 entries, 0 drift, 0 missing (`FINAL_LOCK_VERIFICATION.md`).
- **Synthesis:** `AUDIT_QC_MAD_SYNTHESIS.md`.

**Staged-rules framing:** 6 rules total now sit across the corpus
(4 from Audit 3 in `audit3_counting/draft_rules/`, 2 from QC-MAD in
`audit_qc_mad/draft_rules/`) as **reference documentation** for scRNA analyses
(see `REFERENCE_USE.md`). **BioOrchestrator integration is not planned** per the
tabling decision — the rules function as reference material for manual
application, not engine input. (Supersedes earlier "awaiting batched BO update"
framing.)

## Checkpoints

- **CP0** — inventory + pre-flight. COMPLETE. Surfaced: Audit 3 matrices cleaned
  (`processed/` gone); resolved via Census re-pull (A1). pipeComp absent →
  `scuttle::isOutlier` (the method pipeComp wraps).
- **CP1** — 4-method filtering comparison. COMPLETE (commit 46b0c27). 8 Census
  datasets (subset of Phase 1 p1, census 2025-11-08), mito 0.07–10.83%,
  human+mouse. 6 pairwise Jaccards + disagreement concentration + bootstrap.
- **CP2** — synthesis + findings + §5.3.2 tiers. COMPLETE. `cp2/findings.md`.
- **CP3** — rule drafting. COMPLETE 2026-05-23. 2 rules drafted, both pass
  validator --strict-steps (0 errors, 0 warnings); hash-registered.
- **CP4** — final synthesis + corpus integration. NOT STARTED.

## Findings (CP2)

Headline: **4 QC filtering methods produce largely equivalent cell sets**
(pair Jaccard 0.90–0.97, all CIs overlap). **CONDITIONAL tier** (§5.3.2
equivalence) — strong equivalence + documented feature-dependence
(gene-count distribution drives the small_intestine exception).

- C2 (typical scanpy defaults) ≈ MAD3 (Heumos rec): Jaccard 0.969 — conditional.
- C1 (pure quantile) is the distinct method (descriptive).
- Disagreement driver = gene-count distribution, NOT tissue/mito (conditional).
- MAD k=3 aggressive vs k=5 permissive (descriptive).

Heumos relationship: **extends** — MAD is empirically distinct from pure
quantiles but very close to common fixed-floor defaults; the practical
magnitude of switching is smaller than Heumos's qualitative framing suggests.

## Rule candidates (CP3, decision pending)

1. `scrna_qc_filtering_method_equivalence` (conditional, info) — reassurance.
2. `scrna_qc_low_gene_dataset_caution` (conditional, warn) — actionable on
   low-gene-cell datasets.
Recommendation: draft both. Sub-findings 2 + 4 stay descriptive.

## Standards / corpus updates

- Heumos index: QC filtering row → `covered_by_audit_qc_mad` / extends.
- Coverage summary: QC filtering moved to "covered"; ambient RNA stays open.
- DEFERRED: MAD-vs-quantile gap marked RESOLVED; ambient RNA correction remains.
