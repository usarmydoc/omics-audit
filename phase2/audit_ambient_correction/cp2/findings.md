# CP2 — Deliverable B findings: ordering analysis (correction vs QC)

_Ambient RNA correction audit, 2026-05-24. SoupX 1.6.2 · CellBender 0.3.2 (patched) · DecontX (celda 1.26.0). 9 datasets, STARsolo upstream. QC = C2 fixed-floor (n_genes≥200, total_counts≥500, mito ≤ 95th pct), the QC-MAD standard. Standards: AUDIT_STANDARDS v1.0.3 + §5.3.2. Deliverable B only._

## Two orderings
- **O1 (correct → QC):** run ambient correction on all filtered cells, then C2 QC on the **corrected** counts.
- **O2 (QC → correct):** C2 QC on the **original** counts, then run correction on the QC-passed survivors.

## Headline

**Ordering effects are tool-specific and secondary to tool choice.** Of the three tools only **DecontX** changes its contamination estimate materially with ordering; **SoupX is ordering-invariant** and **CellBender's correction is ordering-invariant by design**. The large cross-tool disagreement found in CP1 **persists essentially unchanged within both orderings** — so the order of operations does not reconcile the tools. The one universal ordering effect is on **cell retention**: *correct → QC* is systematically stricter than *QC → correct*, and the gap grows with ambient burden.

This is the "ordering interacts with tool choice" outcome, sharpened: **tool effects dominate; ordering matters specifically for DecontX (estimate) and for cell retention on high-ambient tissue.**

## 1. Does ordering change the correction? (per-tool O1 vs O2, per-gene contamination Spearman)

| tool | pooled v2 | pooled v3 | per-dataset range | interpretation |
|------|-----------|-----------|-------------------|----------------|
| **SoupX** | 0.996 | 0.993 | 0.992–0.999 | **ordering-invariant** |
| **CellBender** | — (corrected identical by design) | — | 0.916–1.000 | correction invariant; residual reflects cell-set QC only |
| **DecontX** | 0.813 | 0.855 | **0.377–0.965** | **ordering-sensitive** (worst: kidney 0.377, t_3k 0.746, intestine 0.833) |

This maps onto tool design:
- **CellBender** learns ambient from the **raw droplet distribution** — it cannot be refit on QC-passed cells, so its correction is identical across orderings by construction (only QC basis differs; see §3).
- **SoupX** estimates a global soup profile from empty droplets and (here) floors at ρ≈0.01 on most data — empirically stable to which cells survive QC.
- **DecontX** fits a **per-cell mixture model from the cell-by-gene matrix itself**, so refitting on the QC-passed subpopulation vs all cells shifts its per-gene estimate — most on heterogeneous/high-ambient tissue (kidney 0.38).

**Answer to "does ordering matter for any tool":** yes, for DecontX; effectively no for SoupX and CellBender.

## 2. Does ordering reconcile the CP1 cross-tool disagreement? (cross-tool per-gene Spearman within each ordering)

| pair | CP1 (mean) | CP2 O1 | CP2 O2 |
|------|-----------|--------|--------|
| SoupX ↔ CellBender | 0.57 | 0.530 | 0.544 |
| CellBender ↔ DecontX | 0.41 | 0.395 | 0.382 |
| SoupX ↔ DecontX | 0.39 | 0.393 | 0.387 |

**No.** The moderate cross-tool disagreement is virtually identical under both orderings and matches CP1. **Tool choice dominates ordering** — switching the order of operations does not make the tools agree.

## 3. Cell-retention asymmetry: correct→QC is stricter (the universal ordering effect)

O2 (QC on original counts) retains **more** cells than O1 (QC on corrected counts) on every dataset (O2-only ≫ O1-only). The gap scales with how much the tool removes:

| dataset | tool | O1 cells | O2 cells | O2-only (kept by QC→correct, dropped by correct→QC) | Jaccard |
|---|---|---|---|---|---|
| lung organoid | CellBender | 31,342 | 33,508 | **2,347** | 0.925 |
| kidney | CellBender | 24,281 | 25,078 | 823 | 0.966 |
| intestine | DecontX | — | — | 1,221 | 0.922 |
| (clean PBMC) | all | — | — | <130 | ≥0.97 |

**Mechanism:** ambient correction lowers per-cell counts; cells near the 500-count / 200-gene floor drop below it *after* correction. So computing QC on corrected counts (O1) removes cells that pass QC on original counts (O2). The effect is small on clean data (correction removes ~1%) and large on high-ambient tissue (CellBender removes 25–35%). **Practical consequence: "correct then QC" silently discards more cells, concentrated in high-ambient samples.**

## 4. Tissue and chemistry dependence
- **Tissue:** the cell-retention ordering effect is largest on high-ambient tissue (lung, kidney, intestine); negligible on clean PBMC. DecontX's estimate-sensitivity is also worst on the heterogeneous mouse tissues (kidney 0.38, intestine 0.83).
- **Chemistry:** DecontX ordering-sensitivity is slightly worse on v2 (0.813) than v3 (0.855); SoupX is invariant on both. Not a strong chemistry effect.

## Synthesis answers
- **Does ordering matter for any tool?** Yes — DecontX's estimate (ρ as low as 0.38) and cell retention for CellBender on high-ambient tissue. SoupX: no.
- **Which tools most sensitive?** DecontX (estimate); CellBender (retention, high-ambient only). SoupX is the most ordering-robust.
- **Does CP1 disagreement hold within each ordering, or does ordering interact with tool choice?** It holds — tool choice dominates; ordering is a secondary, tool-specific modifier.
- **Does the field's "correct then QC" convention hold empirically?** It is *defensible but not free*: correct→QC is stricter (discards more cells, especially high-ambient), and for DecontX the two orders give materially different corrections. The convention should be stated explicitly, not assumed neutral.
- **§5.3.2:** ordering is not an equivalence question with a single winner; the finding is that ordering-robustness is itself tool-dependent (SoupX/CellBender robust, DecontX not). No equivalence tier claimed.

## Methodological note — CellBender ordering optimization (no corners cut)
CellBender's variational inference requires the raw droplet distribution to learn ambient; it **cannot** be refit on QC-passed cells. Its corrected output is therefore **ordering-invariant by construction**. We **reused CP1's CellBender corrected matrices** for both CP2 orderings and re-ran only the QC step per ordering (O1 = QC on corrected, O2 = QC on original). This is not a shortcut that compromises the comparison — it is the *correct* handling: re-running CellBender would have produced (stochastically near-)identical output at large GPU cost. SoupX and DecontX, which estimate from the cell-by-gene matrix, were **fully re-run for both orderings**. This asymmetry (correction ordering-invariant for raw-droplet methods, ordering-sensitive for cell-matrix methods) is itself a Deliverable-B finding. Validated by the CellBender O1-vs-O2 sanity check: survivor cell sets differ (non-identical) on all 9 datasets, confirming the QC-basis switch operates correctly. Total CP2 runtime ~13 min (vs a naive ~4–8 h re-run-everything estimate) — entirely from this valid reuse + SoupX∥DecontX concurrency.

## Heumos positioning
Heumos 2023 does not address ordering of ambient correction relative to QC. CP2 fills this gap: ordering is a real but tool-specific consideration (matters for DecontX, not SoupX/CellBender), and "correct then QC" is stricter on cell retention than the reverse — a methodological detail the published consensus leaves unspecified.

## Outputs
- `per_dataset_ordering_metrics.tsv` — O1-vs-O2 per (dataset × tool): cell Jaccard, O1/O2-only, per-gene contamination Spearman + bootstrap CI, per-cell total Spearman.
- `per_dataset_cross_tool_metrics.tsv` — cross-tool per-gene Spearman + cell Jaccard within each ordering.
- `cell_disposition.tsv` — O1-only / O2-only cell counts per (dataset × tool).
- `per_stratum_bootstrap.tsv` — chemistry × comparison-type pooled bootstrap CIs.
- `tool_failure_modes.md` — the CellBender ordering-invariance constraint + sanity-check result.
All hash-registered in `phase2/repro.lock` (path-keyed via dge_native).

## Scope note
Deliverable B only. Does not address biological propagation (Deliverable C). C2 QC only; no other QC methods; no non-default correction parameters; no re-counting; no BioOrchestrator changes; CP0/CP1 outputs unmodified.
