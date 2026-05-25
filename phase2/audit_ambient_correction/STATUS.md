# Audit — Ambient RNA Correction — STATUS

_Created 2026-05-24. Standards: AUDIT_STANDARDS.md v1.0.3 + §5.3.2._

## AUDIT STATUS: CLOSED (2026-05-25)
CP0–CP4 all complete. **4 rules contributed** (`draft_rules/`, all validate
under `--strict-steps`). **Standards contributions:** none to
AUDIT_STANDARDS.md (used existing §5.3.2); one registry extension
`scrnaseq_ambient_correction` (no schema change; PENDING_AMENDMENTS.md).
Lock verification CLEAN (216 ok, 0 drift, 0 missing). Synthesis:
`AUDIT_AMBIENT_CORRECTION_SYNTHESIS.md`.

## Scope (locked at CP0)
3 tools (SoupX, CellBender, DecontX) × 3 deliverables:
- **A** — per-gene/per-cell contamination estimate comparison, 3 tools × 9 datasets (mirror Audit 3 C1).
- **B** — ordering analysis: correction-vs-QC, 2 sensible orderings × 3 tools × 9 datasets (4→2 collapse, see CP0 inventory Step 5).
- **C** — biological propagation on intestine (stress) + PBMC (control), mirror Audit 3 C3/CP6; does correction reduce C3's high-ambient permissiveness divergence?

Working set = the 9 Audit 3 datasets (cross-audit consistency). Upstream counting tool fixed (recommend STARsolo) to avoid confounding with Audit 3's counting variation.

## Checkpoints
- **CP0 — inventory & feasibility** — ✓ COMPLETE (2026-05-24). `cp0/inventory.md`. Matrices ✓ (no re-counting), R tools ✓ (SoupX 1.6.2, celda/decontX 1.26.0, DropletUtils 1.30.0), CellBender ✓ (env `ambient_cb`, GPU smoke test PASSED end-to-end with documented checkpoint patch — see `environment/cellbender_env.md`). Scope decisions confirmed. **Awaiting user go-ahead for CP1.**
- **CP1 — Deliverable A** — ✓ COMPLETE (2026-05-24). `cp1/`. 9 datasets × 3 tools, all exit 0. **Finding: tools are NOT equivalent** for per-gene contamination (cross-tool Spearman 0.39–0.57 mean; fails §5.3.2 equivalence). SoupX↔CellBender agree most (0.57); DecontX is the outlier. **SoupX floors at ~1% and fails to detect high ambient** (intestine 1.5% vs CellBender 35%, DecontX 21%); agreement *degrades* on high-ambient intestine (lowest of all 9). 60 outputs hash-registered in repro.lock (path-keyed via dge_native). Awaiting CP2 approval.
- **CP2 — Deliverable B (ordering)** — ✓ COMPLETE (2026-05-24). `cp2/`. 9 datasets × 3 tools × 2 orderings (O1 correct→QC, O2 QC→correct), C2 fixed-floor QC, all exit 0. **Finding: ordering effects are tool-specific and secondary to tool choice.** SoupX ordering-invariant (ρ~0.99); CellBender correction invariant by design (reused CP1, QC-only rerun — sanity check confirmed O1/O2 cell sets differ); **DecontX ordering-sensitive** (ρ 0.81 v2 / 0.86 v3, kidney 0.38). CP1 cross-tool divergence persists within both orderings (tool choice dominates ordering). Universal effect: **correct→QC is stricter** (lowers counts below QC floor → drops more cells; CellBender lung 2,347 O2-only), scaling with ambient burden. 115 outputs hash-registered. Runtime ~13 min (CellBender reuse + concurrency). Awaiting CP3 approval.
- **CP3 — Deliverable C (biological propagation)** — ✓ COMPLETE (2026-05-25). `cp3/`. intestine (high-ambient) + pbmc_5k (clean control), 6 conditions each (nocorr baseline reused from C3 + SoupX/DecontX O1/O2 + CellBender) through cp6 pipeline; all exit 0, ~19 min. **Finding: CP1/CP2 technical differences propagate to biology, scaled by ambient burden.** PBMC near-null (ARI vs baseline 0.85–0.90); **intestine substantially altered (ARI 0.50–0.70)**. Tool choice propagates (intestine cross-tool ARI 0.46–0.60; CellBender most disruptive, ARI 0.50, loses 23% of cells to correct→QC). DecontX ordering propagates on intestine (O1 0.57 vs O2 0.66), negligible on PBMC. Existence-of-effect hard_default; generalization flag_and_warn (n=1+1). 36 outputs registered.
- **CP4 — synthesis + rule drafting + closeout** — ✓ COMPLETE (2026-05-25). 4 rules drafted + validated (`--strict-steps`, 0 errors); `AUDIT_AMBIENT_CORRECTION_SYNTHESIS.md`; `FINAL_LOCK_VERIFICATION.md` (216 ok / 0 drift); registry extension `scrnaseq_ambient_correction`; DEFERRED/AUDIT_INDEX/Heumos-index/coverage/REFERENCE_USE updated. **Audit CLOSED.**

Actual: CP0–CP4 in ~2 days; analytical runs far under estimate (CP1 ~2.5 h, CP2 ~13 min, CP3 ~19 min) via CellBender reuse + concurrency.

## Scope decisions — CONFIRMED (user, 2026-05-24)
1. ✓ Upstream counting tool fixed = STARsolo `star_default` (kb as optional sensitivity).
2. ✓ Ordering design = O1 (correct→QC) / O2 (QC→correct) across all 3 tools; pre-cell-calling orderings excluded for SoupX/DecontX.
3. ✓ Deliverable C reuses existing Audit 3 C3 uncorrected baseline.
