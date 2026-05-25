# Audit QC-MAD-Propagation — STATUS

_Created 2026-05-25. Standards: AUDIT_STANDARDS.md v1.0.3 + §5.3.2._

## AUDIT STATUS: CLOSED (2026-05-25)
CP0–CP2 complete. **1 rule** (`scrna_qc_method_choice_modest_propagation`,
hard_default/info, validates `--strict-steps`). **Standards contributions:** none
(used existing §5.3.2 + `scrnaseq_qc_filtering` step). Lock CLEAN (44 ok, 0 drift,
0 missing). Cross-audit pattern surfaced in prose (not yet a rule). Synthesis:
`AUDIT_QC_MAD_PROPAGATION_SYNTHESIS.md`.

## Question
Do Audit QC-MAD's 4 QC filtering methods (C1 quantile-data, C2 fixed-floors,
MAD k=3, MAD k=5) produce different downstream biology when run through the
standard pipeline? Third propagation test (after Audit 3 C3 + Ambient CP3).

## Scope (locked at CP0)
3 QC-MAD Census datasets (small_intestine, liver, blood) × 4 QC methods = 12
downstream pipeline runs. Reference = C2 (typical scanpy default). Pipeline
mirrors CP6/CP3 exactly. Census re-pull (option A) subset to QC-MAD's exact 50k
soma_joinids per dataset.

## Checkpoints
- **CP0 — inventory & feasibility** — ✓ COMPLETE (2026-05-25). `cp0/inventory.md`.
  Census re-pull verified reproducible (liver metrics match QC-MAD to precision);
  3 UUIDs resolve; QC methods + pipeline reusable; models on disk (liver, blood)
  or downloadable (intestine); baseline = C2. **No blockers.** Awaiting CP1 approval.
- **CP1 — propagation runs** — ✓ COMPLETE (2026-05-25). 3 datasets (soma_joinid-
  matched to QC-MAD, verified exact) × 4 QC methods × cp6 pipeline = 12 runs, 0
  failures. **Finding: QC-method choice is a modest, tissue-INDEPENDENT downstream
  effect (ARI 0.80–0.91 vs C2 everywhere; annotation 88–98%) — the edge-case
  amplification of Audit 3 C3 / Ambient CP3 does NOT replicate.** small_intestine
  (ARI 0.84–0.87) ≈ blood control (0.80–0.88). 42 outputs registered. Awaiting CP2.
- **CP2 — synthesis + rule + closeout** — ✓ COMPLETE (2026-05-25). Synthesis,
  1 rule (validated), FINAL_LOCK_VERIFICATION (44 ok/0 drift), DEFERRED/AUDIT_INDEX/
  Heumos-index/coverage/REFERENCE_USE updated, cross-audit pattern in prose.
  **Audit CLOSED.** (CP3/CP4 folded into CP2 — findings clear, light closeout.)

## Open items for CP1
- Download `Cells_Intestinal_Tract.pkl` (human gut CellTypist model).
- Census-var loader tweak to CP3 pipeline (use Census feature_name symbols).
- Sanity-check small_intestine + blood re-pull metrics vs QC-MAD saved.
