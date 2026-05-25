# CP1 — Deliverable A findings: per-gene contamination estimate comparison

_Ambient RNA correction audit, 2026-05-24. SoupX 1.6.2 · CellBender 0.3.2 (patched) · DecontX (celda 1.26.0). 9 datasets, STARsolo upstream (locked), per-gene contamination = (orig−corrected)/orig on each tool's native cell set. Standards: AUDIT_STANDARDS v1.0.3 + §5.3.2. Deliverable A only — ordering (B) and biological propagation (C) not covered here._

## Headline

**The three ambient-correction tools do NOT agree on per-gene contamination.** Cross-tool per-gene Spearman is **moderate at best and never high** (pairwise means 0.39–0.57; all 27 dataset×pair values in 0.24–0.69), far below an equivalence finding. They disagree on **both magnitude** (how much ambient) **and ranking** (which genes are contaminated). This is the opposite of Audit 3's counting-tool convergence (ρ≈0.95+): **the choice of ambient-correction tool materially changes the result.**

## Cross-tool agreement (per-gene Spearman)

| Pair | mean | min | max |
|------|------|-----|-----|
| **SoupX ↔ CellBender** | **0.574** | 0.491 | 0.691 |
| CellBender ↔ DecontX | 0.411 | 0.292 | 0.509 |
| SoupX ↔ DecontX | 0.391 | 0.237 | 0.530 |

SoupX and CellBender track each other best; **DecontX is the consistent outlier** (lowest agreement with both). No pair reaches even ρ=0.7 on average.

### By chemistry (pooled bootstrap, B=1000 gene-level, tight CIs)

| chemistry | SoupX↔CellBender | CellBender↔DecontX | SoupX↔DecontX |
|---|---|---|---|
| v2 (2 datasets, 54,382 genes) | 0.606 [0.600, 0.612] | 0.426 [0.419, 0.432] | 0.253 [0.245, 0.261] |
| v3 (7 datasets, 216,080 genes) | 0.572 [0.569, 0.576] | 0.505 [0.502, 0.509] | 0.374 [0.370, 0.378] |

**§3.5 caveat (CP4 amendment):** the bootstrap CIs are extremely narrow (±0.003–0.008) because pooled gene N is huge. This is **sampling precision, not correctness** — it means we are very *confident* the agreement is *only moderate* (~0.57 best), not that the tools agree well. The tightness must not be read as concordance.

## Magnitude divergence (global fraction removed)

| dataset | SoupX | CellBender | DecontX |
|---|---|---|---|
| PBMC 1k/5k/10k, neuron, kidney | **0.010** | 0.011–0.016 | 0.026–0.078 |
| t_3k_v2 / pbmc_4k_v2 | 0.085 / 0.056 | 0.018 | 0.029 / 0.102 |
| lung organoid | **0.010** | 0.047 | **0.257** |
| intestine (high-ambient) | **0.015** | **0.352** | **0.211** |

Three distinct behaviors:
- **CellBender** — low on clean data (~1–2%), but **scales with true ambient burden** (lung 4.7%, intestine 35%). Tracks the biology.
- **DecontX** — moderate-to-aggressive everywhere (2.6–25.7%); most aggressive on clean/organoid data.
- **SoupX** — **floors at ρ≈0.01 on 6/9 datasets and fails to detect high ambient** (intestine 1.5% vs CellBender's 35%). Detects *moderate* ambient (v2: 5.6–8.5%) but not *extreme* ambient.

## High-ambient case (intestine) — agreement DEGRADES, it does not improve

A natural hypothesis is that more ambient = more signal = better cross-tool agreement. **The data show the opposite.** Intestine (the C3 high-ambient stress dataset) has the **lowest agreement of all 9 datasets**:

| pair | intestine ρ | mean contam A | mean contam B |
|---|---|---|---|
| SoupX↔CellBender | 0.491 | 0.014 | 0.220 |
| CellBender↔DecontX | 0.292 | 0.220 | 0.179 |
| SoupX↔DecontX | 0.237 | 0.014 | 0.175 |

Mechanism: at extreme ambient, **SoupX floors (1.4%) while CellBender (22%) and DecontX (18%) detect heavy contamination** — a >10× magnitude gap that drags rank-agreement down. The moderate "agreement holds at v2 ambient" seen at intermediate burden (v2 SoupX↔CellBender 0.61) **does not extend to extreme ambient** because SoupX stops responding. Agreement vs ambient is **non-monotonic**: acceptable at low–moderate burden, collapsing at high burden.

## Stratified detail
- **Mitochondrial genes:** SoupX↔CellBender agree well on mito in clean data (ρ 0.72–0.87) but mito agreement **collapses or goes negative at high ambient** (lung −0.07, intestine −0.12). Any pair with DecontX agrees poorly on mito throughout (kidney CellBender↔DecontX mito ρ = −0.35). Mito-gene contamination ranking is tool-specific.
- **High-expression genes:** SoupX↔CellBender agree best here (ρ up to 0.76); DecontX pairs stay low (≤0.25).
- **Low-expression genes:** agreement is ~0 for all pairs; CellBender zeros low-expression ambient, making the bottom-decile Spearman undefined for its pairs.

## §5.3.2 tier reasoning
The audit question is equivalence-shaped ("do the tools measure the same contamination?"), so §5.3.2 criteria apply. **The finding fails equivalence at every tier:** equivalence hard_default requires cross-tool Spearman ρ ≥ 0.90 median with all per-pair lower bounds ≥ 0.80; observed medians are 0.39–0.57 with maxima ≤ 0.69. **Conclusion: the three tools are NOT equivalent for per-gene contamination estimation** — this is a *divergence* finding, not an equivalence finding. Tool choice is consequential, and the divergence is structured (SoupX conservative/floor-prone, CellBender ambient-tracking, DecontX aggressive).

## Heumos 2023 positioning
Heumos names **SoupX and CellBender** as ambient-correction options without recommending between them, and does not cover **DecontX**. This audit extends Heumos by:
1. **Quantifying SoupX vs CellBender agreement** (Heumos does not measure it): moderate (ρ≈0.57) and **breaking down on high-ambient tissue**, where SoupX under-detects and CellBender does not.
2. **Adding DecontX** to the comparison (Heumos omits it): it is the most divergent of the three (lowest rank agreement, most aggressive magnitude).

The practical implication — does this magnitude/ranking divergence propagate to clustering, markers, and DE? — is **Deliverable C** (not run here).

## Outputs (this checkpoint)
- `per_gene_contamination.tsv` — long format, per-gene contamination per (dataset × tool).
- `per_dataset_metrics.tsv` — Spearman, log2-ratio, directional, mito/non-mito/expr-decile strata, bootstrap CI per (dataset × pair).
- `per_stratum_bootstrap.tsv` — chemistry × pair pooled bootstrap CIs.
- `tool_failure_modes.md` — gene-ID integrity, SoupX floor behavior, CellBender symbol/zeroing, DecontX mito, cell-set caveat.
All hash-registered in `phase2/repro.lock`.

## Scope note
Deliverable A only. Does **not** address ordering (B) or biological propagation (C). No ordering analysis, no non-default parameters, no re-counting, no BioOrchestrator changes were performed.
