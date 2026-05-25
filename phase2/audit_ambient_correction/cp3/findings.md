# CP3 — Deliverable C findings: biological propagation

_Ambient RNA correction audit, 2026-05-25. SoupX 1.6.2 · CellBender 0.3.2 · DecontX (celda 1.26.0). 2 datasets (intestine = high-ambient stress; pbmc_5k = clean control), 6 conditions each (no-correction baseline + SoupX/DecontX O1/O2 + CellBender). Downstream = Audit 3 CP6 pipeline, identical params. Baseline reused from C3 star_default. Standards: v1.0.3 + §5.3.2._

## Headline

**The technical differences from CP1 (tool disagreement) and CP2 (DecontX ordering sensitivity) propagate into biological conclusions — and the magnitude is governed by ambient burden.** On the clean PBMC control, all five corrected conditions preserve the baseline biology well (ARI 0.85–0.90 vs no-correction). On high-ambient intestine, correction **substantially alters** clustering and annotation (ARI 0.50–0.70 vs baseline), the tools **diverge from each other** (pairwise ARI 0.46–0.66), and **DecontX's ordering sensitivity propagates** (ARI O1 0.57 vs O2 0.66). CellBender is the most disruptive on intestine. **The framing assumption holds: ambient correction matters most where ambient is highest.**

## 1. Does correction change biology, and is it ambient-dependent? (vs-baseline)

| | PBMC (clean) | intestine (high-ambient) |
|---|---|---|
| ARI vs baseline (range) | **0.85 – 0.90** | **0.50 – 0.70** |
| annotation agreement | 0.86 – 0.93 | 0.54 – 0.78 |
| marker Jaccard (cluster-matched top-50) | 0.76 – 0.86 | 0.48 – 0.76 |

The effect is ~2× larger on intestine across every metric. PBMC is not a perfect null (ARI ~0.87, i.e. correction does nudge clustering of clean data), but the **contrast is large and consistent** — exactly the ambient-burden-driven pattern C3 established. Bootstrap CIs are tight (±0.01; §3.5 precision-not-correctness) and the PBMC vs intestine effect bands do not overlap.

### Per-condition vs baseline (intestine)
| condition | postQC cells | clusters (r1.0) | ARI [95% CI] | annot. agree | marker Jacc |
|---|---|---|---|---|---|
| no-correction | 21,807 | 21 | — | 1.00 | 1.00 |
| SoupX O1 | 21,807 | 25 | 0.663 [0.654, 0.671] | 0.769 | 0.704 |
| SoupX O2 | 21,807 | 24 | 0.698 [0.689, 0.706] | 0.784 | 0.761 |
| CellBender | **16,700** | 24 | **0.501** [0.493, 0.510] | **0.542** | 0.522 |
| DecontX O1 | 20,640 | 25 | 0.566 [0.559, 0.574] | 0.710 | 0.476 |
| DecontX O2 | 20,593 | 25 | 0.661 [0.652, 0.669] | 0.719 | 0.522 |

All corrected conditions add clusters (24–25 vs 21) on intestine — correction sharpens/splits structure.

## 2. Does CP1's tool disagreement propagate? (cross-tool, matched conditions)

| pair | intestine ARI | PBMC ARI |
|---|---|---|
| SoupX ↔ CellBender | 0.495 | 0.886 |
| SoupX ↔ DecontX (O1) | 0.595 | 0.832 |
| CellBender ↔ DecontX (O1) | **0.460** | 0.776 |

**Yes.** On intestine the tools produce materially different clusterings (pairwise ARI 0.46–0.60) — CellBender↔DecontX disagree most (0.46), mirroring CP1's per-gene divergence. On clean PBMC the tools largely agree (0.78–0.89). **Tool choice changes the biological answer, specifically on high-ambient tissue.**

## 3. Does CP2's DecontX ordering sensitivity propagate? (O1 vs O2)

| tool | intestine O1-vs-O2 ARI | PBMC O1-vs-O2 ARI |
|---|---|---|
| **DecontX** | **0.639** | 0.818 |
| SoupX | 0.660 | 0.906 |

**Yes, and it is tissue-dependent.** DecontX's ordering choice changes its intestine clustering (O1-vs-O2 ARI 0.64; and vs baseline O1 0.566 vs O2 0.661 — a 0.10 ARI gap, with O1 marker Jaccard 0.476 vs O2 0.522). On clean PBMC the ordering effect is negligible (ARI 0.82). SoupX is less ordering-sensitive than DecontX on intestine, consistent with CP2. **CP2's "DecontX is the ordering-sensitive tool" finding carries through to biology, on high-ambient tissue.**

## 4. Does the "correct→QC is stricter" effect (CP2) change biology?
**Yes, via cell loss.** CellBender removes ~35% of counts on intestine; computing QC on those corrected counts drops the dataset from 21,807 → **16,700 cells** (−23%), versus SoupX (21,807, no loss) and DecontX (~20,600). The most aggressive correction, combined with QC-on-corrected, discards the most cells — and produces the lowest ARI/annotation agreement (0.50/0.54). On clean PBMC no such cell loss occurs. The ordering/QC interaction is biologically consequential only under high ambient.

## 5. Where do contested cells go?
Cells re-annotated relative to baseline concentrate in **closely related subtypes**, not random reassignment. PBMC: reassignments are within the T/NK compartment (e.g. CellBender re-labels 471 baseline "Tcm/Naive helper T" cells; SoupX/DecontX similar, 264–274) — correction shifts fine T-subtype boundaries. (Full table: `contested_cell_disposition.tsv`.) This indicates correction perturbs sub-clustering granularity rather than gross cell-type identity.

## Synthesis answers
- **Tool choice (CP1) → biology?** Yes, on intestine (cross-tool ARI 0.46–0.60); minimal on clean PBMC.
- **DecontX ordering (CP2) → biology?** Yes, on intestine (O1-vs-O2 ARI 0.64; 0.10 ARI gap vs baseline); negligible on PBMC.
- **"correct→QC stricter" → biology?** Yes — CellBender loses 23% of intestine cells, lowest concordance.
- **Does correction matter on intestine?** Substantially (ARI to baseline as low as 0.50).
- **Is PBMC null?** Near-null relative to intestine (ARI 0.85–0.90) — confirms the effect is ambient-burden-driven.

## §5.3.2 / sample-size tiering
- **Existence-of-effect** (ambient correction alters biology, scaled by ambient burden; tool & ordering choices propagate on high-ambient tissue): supported by a clear, large PBMC-vs-intestine contrast with tight CIs → **hard_default** for existence.
- **Generalization** (exact magnitudes, which tool is most disruptive in general): n=1 high-ambient + n=1 control, same structure as Audit 3 CP6 → **flag_and_warn** until replicated on additional high-ambient tissues.

## Heumos positioning
Heumos 2023 recommends ambient correction as an initial step without quantifying downstream biological consequences of tool/ordering choice. CP3 shows the recommendation's qualitative framing has **substantive downstream effects on high-ambient tissue**: tool choice and (for DecontX) ordering change clusters, markers, and annotations; on clean tissue the effect is modest. The choice is not biologically neutral where it matters most.

## Outputs
- `all_per_condition_metrics.tsv` — per (dataset × condition): postQC, clusters, ARI/NMI vs baseline + bootstrap CI, annotation agreement, marker Jaccard.
- `cross_condition_comparison.tsv` — pairwise ARI/NMI/annotation across all conditions.
- `contested_cell_disposition.tsv` — baseline labels of re-annotated cells per condition.
- `<sub>/<condition>.{h5ad,obs.tsv,markers.tsv,summary.json}` — full pipeline outputs per condition.
All hash-registered in `phase2/repro.lock` (path-keyed via dge_native).

## Scope note
Deliverable C only. intestine + pbmc_5k only; baseline reused from C3; cp6 QC params (min_genes 200, max_mito 20) used throughout for C3 comparability (not CP2's C2 fixed-floor — the ordering axis is preserved as whether each tool's correction was fit on all cells vs QC-passed cells). No CP4 synthesis/rule-drafting performed; CP0/CP1/CP2 outputs unmodified.
