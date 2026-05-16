# Phase 2a — Findings

All compute local desktop. All outputs hash-registered in
`/mnt/nvme1/omics-audit/phase2/repro.lock`. No new modalities introduced —
extensions on existing scRNA-seq + bulk audit work.

**Provenance categories** for each finding:

- `as_original`: rule matches prior audit finding exactly
- `refines_prior`: original claim correct but needs nuance
- `contradicts_prior`: specific claim doesn't hold under broader testing
- `extends_prior`: new finding not in original audit

---

## TOP-LEVEL FINDING 1 — Tool concordance: agreement on biology vs disagreement on ranking

**Provenance:** `extends_prior` (Phase 2 main + Phase 2a confirmation)

The sharpest, most defensible framing in the entire Phase 2 corpus:

| Source | Pair | Direction agreement | log2FC Pearson | Top-100 Jaccard |
|---|---|---|---|---|
| TCGA bulk | DESeq2_vs_edgeR | **1.000** | 0.951 | 0.310 |
| TCGA bulk | DESeq2_vs_limma | **1.000** | 0.903 | 0.222 |
| TCGA bulk | edgeR_vs_limma | **1.000** | 0.858 | 0.198 |
| GTEx-pair | DESeq2_vs_edgeR | **1.000** | 0.987 | 0.211 |
| GTEx-pair | DESeq2_vs_limma | **1.000** | 0.971 | 0.145 |
| GTEx-pair | edgeR_vs_limma | **1.000** | 0.961 | 0.150 |
| Census scRNA | pseudobulk_vs_wilcoxon | **1.000** | 0.507 | 0.101 |

**Tools agree on the biology — direction of effect is unanimous, log2FC
magnitude is highly correlated (0.86-0.99 for bulk). Disagreement is
specifically about ranking by padj.** Top-N Jaccard is a brittle proxy
for "tool concordance" — it captures only one facet (ranking
agreement) and conceals the more important agreement on direction and
magnitude.

For Census pseudobulk vs cellwise Wilcoxon, log2FC Pearson is lower
(0.51) because the two methods are fundamentally different (per-donor
aggregate vs per-cell tests), not just different tools applied to the
same data. That 0.51 is biology, not noise.

**Implication for BioOrchestrator:** any rule that reports "tools
disagree at Jaccard X" without simultaneously reporting direction
agreement + log2FC correlation is misleading. See
`phase2/draft_rules/tool_concordance_reporting.yaml`.

---

## TOP-LEVEL FINDING 2 — Clustering selection metric matters as much as resolution itself

**Provenance:** `extends_prior` — P9 established that resolution should
be tuned per dataset but didn't address selection-metric choice.

**A2 result across 15 P9 datasets, 8 resolutions, 5 metrics:**

| Metric | Median optimal res | Mean | Distribution |
|---|---|---|---|
| **ARI** | **0.10** | 0.24 | 9/15 datasets pick 0.1 — VERY low |
| **V-measure** | 0.50 | 0.67 | spread 0.1-2.0 |
| Silhouette | 0.80 | 0.72 | spread 0.1-1.5 |
| Homogeneity | 2.0 (max) | 1.97 | **14/15 pick max — monotonically increasing in resolution** |
| Completeness | 0.1 (min) | 0.11 | **14/15 pick min — monotonically decreasing in resolution** |

**Key sub-findings:**

1. **Homogeneity and completeness are monotonic in resolution** —
   14/15 datasets pick the extreme value of the resolution sweep for
   each. They are **useless as standalone selection criteria**. Only
   meaningful as components of V-measure.

2. **ARI systematically under-recommends resolution.** On the 9/15
   datasets where ARI and V-measure disagree, **V-measure picks 5.1×
   higher resolution on average** (mean 0.90 vs 0.18).

3. **No two metrics agree more than 60% of the time** (ARI/completeness:
   9/15; ARI/V-measure: 6/15; ARI/silhouette: 3/15; V-measure/silhouette: 5/15).

**Implication for BioOrchestrator:** the choice of evaluation metric
shifts the recommended resolution by ~5× on disagreeing datasets. A
pipeline using ARI will under-cluster; one using V-measure will cluster
closer to reference granularity. Recommended rule: report multiple
metrics + explicit under/over-clustering tolerance, or use V-measure
as the primary selector. Never use homogeneity or completeness alone.
See `phase2/draft_rules/clustering_metric_selection.yaml`.

Outputs:
- `phase2a/p9_regen/{dataset_short}/` — per-dataset regen with per-cell
  labels, per-cluster markers at optimal resolution, HVG list
- `phase2a/p9_clustering_metrics.tsv` — 120-row flat metric table
- `phase2a/p9_optimal_by_metric.tsv` — 15-row optimal-by-metric table

---

## TOP-LEVEL FINDING 3 — B4's 80% DEG change is method-specific (SVA preprocessing)

**Provenance:** `refines_prior` (B4 audit finding) — directly refines
the prior B4 conclusion.

**B4's headline finding** was that batch correction method choice
altered ~80% of DEG lists vs uncorrected. **Phase 2a A7 finds this is
preprocessing-dependent for SVA but robust for ComBat and limma.**

| Method | Default deg_stability (change rate) | Native deg_stability (change rate) | |delta| |
|---|---|---|---|
| ComBat | 0.299 (70% change) | 0.276 (72% change) | **0.023 — finding holds** |
| limma | 0.276 (72% change) | 0.227 (77% change) | **0.049 — finding holds** |
| **SVA** | 0.283 (72% change) | **0.567 (43% change)** | **0.285 — finding does NOT hold** |

**Mechanistic interpretation:** SVA on log-CPM is more aggressive at
removing surrogate variable signal than svaseq on raw counts. The
~80% DEG change attributable to "SVA" was really attributable to the
SVA-on-log-CPM combination; the count-native svaseq variant preserves
much more of the original DEG list.

**Default-vs-native inter-method Jaccard** (do the two preprocessing
variants of one method agree?):

| Method | top-100 Jaccard between default and native |
|---|---|
| ComBat (default) vs ComBat-Seq (native) | 0.281 |
| limma RBE (default) vs voom (native) | 0.492 |
| SVA (default) vs svaseq (native) | 0.358 |

All three preprocessing pairs disagree substantially (Jaccard <0.5),
showing preprocessing choice is itself a meaningful determinant of DEG
list — even when the underlying method family is held constant.

**Implication for BioOrchestrator:** the rule about SVA should not pick
sides between SVA and svaseq, but should flag when SVA is invoked
without explicit preprocessing documentation. See
`phase2/draft_rules/sva_preprocessing_sensitivity.yaml`.

Outputs:
- `phase2a/b4_native/{cancer}/{cancer}__{variant}.tsv` — 16 cancers ×
  up to 7 variants each
- `phase2a/b4_native_compare.tsv` — pairwise Jaccard table

---

## A1 — P4 top-K sensitivity

**Provenance:** `refines_prior` — qualitative claim of pseudobulk-vs-Wilcoxon
disagreement holds; specific 0.094 number is K=100-specific.

21 P4 datasets analyzed.

| Metric | Mean | Median | Range |
|---|---|---|---|
| Jaccard top-50 | 0.078 | 0.053 | 0.00 – 0.33 |
| Jaccard top-100 | 0.103 | 0.064 | 0.00 – 0.42 |
| Jaccard top-200 | 0.143 | 0.102 | 0.01 – 0.55 |
| Jaccard padj<0.05 | 0.240 | 0.219 | 0.00 – 0.57 |
| Jaccard padj<0.01 | 0.165 | 0.106 | 0.00 – 0.58 |

**K-stability verdict: NOT stable.** Jaccard rises monotonically with K
(0.08 → 0.10 → 0.14 from K=50 to K=200, 1.8×). FDR-thresholded overlap
is higher than top-K. Published 0.094 figure is K=100-specific. The
qualitative claim that tools disagree holds; the specific magnitude
depends on K.

Output: `phase2a/p4_topk_sensitivity.tsv` (21 rows).

---

## A3 — Per-tissue mito threshold — reframed to quantile-on-data

**Provenance:** `refines_prior` — P1's "tissue-specific" claim correct,
but the deliverable is not a per-tissue table.

**Original deliverable plan:** per-tissue threshold table from 60 P1
dataset-tissue rows.

**What the data shows:** No tissue has n ≥ 3 datasets in P1. All 42
tissues flagged `insufficient_data`. Empirical p5/p50/p95/p99 emitted
for all 42 (some with n=1, useless); no `recommended_threshold` value
issued. Biology spot-check on n=2 tissues mostly matches expectations
(gut/liver/pancreas HIGH; brain/blood LOW) but kidney p95=2.2% is an
obvious one-dataset outlier — expected HIGH for metabolic tissue.

**Reframed deliverable: quantile-on-data rule.** The BioOrchestrator
rule should not depend on per-tissue tables. It should compute the
threshold from the working dataset:

> *Compute the 95th percentile of mitochondrial fraction in the working
> dataset; use that as the cell-filtering threshold. If the dataset has
> < 100 cells (insufficient for stable percentile estimation), fall back
> to the conservative 20% default.*

This avoids the tissue-table dependency entirely and adapts to the
data the user actually has. See
`phase2/draft_rules/mito_threshold_quantile.yaml`.

The original `phase2a/p1_per_tissue_thresholds.tsv` is retained as
documentation of why per-tissue tables are not currently deployable —
not deleted, marked clearly as `insufficient_data` per row.

---

## A4 — Resolution vs annotation granularity

**Provenance:** `as_original` — confirms P9's sweep-required conclusion.

15 P9 datasets (not 20 — P9 actually has 15 unique datasets with
duplicate runs handled by max-ARI).

| Test | Statistic | p-value | Notes |
|---|---|---|---|
| Pearson linear | r = 0.378, R²=0.14 | 0.165 | not significant — no linear fit |
| Spearman rank | ρ = 0.589 | **0.021** | significant rank correlation |
| Log-log slope | 0.921, R²=0.33 | 0.025 | proportional-ish, 67% variance unexplained |

**Honest verdict: WEAK linear, moderate rank-based.** There IS a
relationship (more cell types → higher optimal resolution in rank
order, slope ~1 in log space), but noisy enough that you can't
reliably *predict* the optimal value from n_celltypes alone. Sweep
remains the only safe recommendation. n_celltypes can be a starting
heuristic but not a substitute for the sweep.

Outputs: `phase2a/p9_resolution_vs_granularity.tsv`,
`phase2a/p9_resolution_vs_granularity_summary.json`.

---

## A5 — Demuxlet doublet ground truth

**Provenance:** `as_original` — confirms paper's synthetic finding
generalizes to natural ground truth.

**Citation note:** Phase 2a spec cited GSE96583 as "Kang 2018
HEK293T/Jurkat mixture". GSE96583 is actually the Kang 2018 PBMC
pooled-donor demuxlet study (8 donors, ~14K cells, pools A/B/C).
Doublets are still heterotypic by construction (cross-donor PBMC
mixing), which is the property the spec wanted to test, so the
finding stands.

Results on 12,282 cells with demuxlet calls (11,726 singlets + 556
doublets; ambs excluded):

| Tool | AUROC | F1 (default thr) | Precision | Recall |
|---|---|---|---|---|
| Scrublet | 0.741 | 0.035 | 0.833 | 0.018 |
| **scDblFinder** | **0.834** | **0.464** | 0.347 | 0.700 |

**scDblFinder outperforms Scrublet on natural heterotypic ground
truth.** AUROC delta ~0.09. Scrublet's auto-threshold is too strict
(default recall 1.8%) — high precision comes at almost no recall.
scDblFinder catches 70% of true doublets at default threshold.

Confirms the paper's synthetic heterotypic finding generalizes to
natural heterotypic ground truth.

Outputs:
- `phase2a/p2_demuxlet_benchmark.tsv` (2-row summary)
- `phase2a/p2_demuxlet_per_cell.tsv.gz` (per-cell scores)

---

## A6 — muscat-dream inflation: 26-dataset verification

**Provenance status:** `refines_prior` — §1.9a's center-of-distribution
claim holds, but the bimodal extremes span more tissues and organisms
than §1.9a stated.

**Run:** 2026-05-16, 26 datasets (20 H. sapiens + 6 Mus musculus),
235.6 min wall, ~9 min/dataset. Six per-dataset bug fixes documented in
`phase2/scripts/a6_nebula_muscat.R` and CLOSEOUT.md A6 section.

### Aggregate metrics (n=26)

| Metric | Value |
|---|---|
| muscat-dream inflation median | 0.78× (§1.9a 5-dataset claim: 0.85×) |
| muscat-dream inflation range | 0.00× to 57.83× — strongly bimodal |
| n datasets > 2× inflation | 7 of 26 (27%) |
| n datasets > 10× inflation | 4 of 26 (15%) |
| top-100 Jaccard mm vs pb (median) | 0.176 — 21/26 < 0.30 |
| log2FC Spearman mm vs pb (median) | 0.841 (direction/magnitude agree on shared hits) |

### HIGH inflation set (>50% genes called significant by mm)

| dataset | % sig | tissue | cell type |
|---|---:|---|---|
| c7775e88 | 96.7% | blood (COVID-19) | naive CD4 T cells |
| b617ee1b | 88.8% | breast (multi-cancer) | T cells |
| 16023185 | 86.9% | colon adenocarcinoma | stem cells |
| **a19d1667** | 77.3% | skin (keloid/scleroderma) | **skin fibroblast (§1.9a case ✓)** |
| a48343a2 | 59.5% | musculature/neuromuscular | type IIb muscle cell |

### LOW activity set (<10% genes sig) — 13 datasets

Notably includes **b2dd6bc9 (MOUSE skin fibroblast at 0.03%)** — directly
contradicting any tissue-specific generalization. Same tissue, different
organism, opposite behavior. The inflation is dataset-driven, not
tissue-driven.

### Verdict framing

The disagreement between muscat-dream and pseudobulk is in **which genes
pass FDR** (Jaccard 0.176 — they pick mostly different top hits), not
in the **direction or magnitude** of the effect on the shared subset
(log2FC Spearman 0.841 — when both tools call a gene, they agree).
This means muscat-dream isn't producing wrong biology — it's producing
wrong significance calls due to pseudoreplication inflation when the
random-effect variance is poorly estimated.

### Rule encoded

`phase2/draft_rules/muscat_dream_inflation_warn.yaml` —
`scrna_muscat_dream_cell_level_de_inflation_warn` —
severity=warn, confidence_tier=conditional. Recommendation: report
muscat-dream side-by-side with pseudobulk; flag dataset-level
inflation by joint criterion (mm/pb ratio > 2× AND top-100 Jaccard
< 0.30). Validates clean against schema v1.0.2.

### Outputs (all hash-registered)

- Per-dataset: `phase2a/a6_mixed_model/<short>/<short>__<tool>.tsv`
  (78 TSVs = 26 datasets × 3 tools)
- Aggregate: `phase2a/a6_analysis_full.tsv`
  (sha256 `5988c3c136a815f9...`)
- Run log: `phase2a/logs/a6_rerun_dedup_2026-05-16.log`

---

## Scope expansion note — P4 datasets: 21 → 26

Phase 2 P4 regen succeeded on **26 datasets** for Wilcoxon vs the
original audit's 21. The 5 additional datasets succeeded with
Wilcoxon but failed pseudobulk's ≥3 donors per group requirement
(Phase 2 regen criteria identical to Phase 1).

The 5 additional Census datasets (short IDs): `5d6d404a, 7b6bab5a,
d3cb449b, ebc2e1ff, f14bc322`. These represent newly-eligible Census
content within version `2025-11-08` that wasn't included when Phase 1
ran in March 2026 (likely due to Census data version internal updates
or pseudobulk eligibility shifts).

**For comparison against original 21-dataset findings**, use the
21-dataset subset (intersected against
`phase1_march2026_superseded/p4_pseudobulk.tsv` dataset_id list).
**For BioOrchestrator rule input**, 26-dataset Wilcoxon results are
appropriate where the rule doesn't compare to the original 21-dataset
audit specifically. The `audit1_inputs.tsv` flags which datasets are
in the original 21 vs the additional 5.

---

## Draft rule YAML inventory

All rules in `phase2/draft_rules/`. All include
`prior_audit_relationship` metadata per Issue 7 of the closeout.

| Rule | Source | Provenance | Status |
|---|---|---|---|
| `bioc_version_sensitivity_v2.yaml` | Step 5b | extends_prior | drafted; two-tiered |
| `tool_concordance_reporting.yaml` | Tool concordance finding | extends_prior | drafted |
| `clustering_metric_selection.yaml` | A2 | extends_prior | drafted |
| `sva_preprocessing_sensitivity.yaml` | A7 | refines_prior (B4) | drafted |
| `mito_threshold_quantile.yaml` | A3 (reframed) | refines_prior (P1) | drafted |
| `deg_tool_version_sensitivity.yaml` | earlier draft | (superseded by bioc_version_sensitivity_v2) | retained for history |

A6 rule deferred until A6 completes per Issue 6.

---

## Outstanding

- A6 result (NEBULA + muscat on all 21) — overnight unattended; verifies
  §1.9a framing per closeout Issue 6
- Supervised review of all 5 draft rule YAMLs before BioOrchestrator
  integration
