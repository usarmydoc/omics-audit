# CP1 — QC-MAD-Propagation findings: does QC filtering-method choice propagate to biology?

_2026-05-25. 3 Census datasets (blood/liver/small_intestine, soma_joinid-matched to QC-MAD's exact cells) × 4 QC methods (C1 quantile-data, C2 fixed-floors [reference], MAD k=3, MAD k=5) × the Audit 3 CP6 / Ambient CP3 downstream pipeline = 12 runs. Standards: v1.0.3 + §5.3.2._

## Headline

**QC filtering-method choice is a modest, tissue-INDEPENDENT downstream variance source — the edge-case-amplification pattern from the prior two propagation audits does NOT replicate.** Across all 3 datasets, switching from the typical fixed-floor default (C2) to any other method (C1, MAD k=3/5) gives cluster-assignment **ARI 0.80–0.91 vs C2** and **annotation agreement 88–98%**, with the high-disagreement edge-case tissue (small_intestine) **no worse** than the clean control (blood). This contrasts sharply with cell-calling (Audit 3 C3) and ambient correction (Ambient CP3), which propagated *strongly on high-ambient/edge-case tissue* (ARI 0.50–0.70) and washed out on clean tissue.

## Cross-method vs C2 reference (Leiden r1.0)

| dataset | comparison | n∩ | ARI [95% CI] | NMI | marker Jaccard | annot. agree |
|---------|-----------|----|----|----|----|----|
| blood (control) | C1 vs C2 | 44,538 | 0.798 [0.793, 0.803] | 0.903 | 0.837 | 94.3% |
| | MAD3 vs C2 | 47,008 | 0.879 [0.874, 0.883] | 0.919 | 0.867 | 93.6% |
| | MAD5 vs C2 | 47,421 | 0.836 [0.831, 0.841] | 0.898 | 0.844 | 94.3% |
| liver | C1 vs C2 | 44,627 | 0.908 [0.904, 0.912] | 0.965 | 0.868 | 98.1% |
| | MAD3 vs C2 | 39,941 | 0.878 [0.873, 0.882] | 0.956 | 0.790 | 98.2% |
| | MAD5 vs C2 | 44,042 | 0.911 [0.907, 0.915] | 0.973 | 0.868 | 98.3% |
| small_intestine (edge) | C1 vs C2 | 38,700 | 0.867 [0.862, 0.871] | 0.915 | 0.760 | 89.5% |
| | MAD3 vs C2 | 38,700 | 0.872 [0.868, 0.877] | 0.914 | 0.743 | 87.8% |
| | MAD5 vs C2 | 38,700 | 0.839 [0.834, 0.844] | 0.903 | 0.742 | 88.0% |

Bootstrap CIs (B=1000, cell-level) are tight (±0.005); the modest effect is precise, not noisy.

## Key observations

1. **No edge-case amplification.** small_intestine ARI (0.84–0.87) ≈ blood (0.80–0.88); liver is actually the *highest* agreement (0.88–0.91). The tissue where QC-MAD found the largest *cell-set* disagreement (small_intestine, C2 vs MAD5 Jaccard 0.774) does **not** show the largest *biological* disagreement. The effect is roughly uniform across tissues.

2. **The control isn't near-null either.** Even on blood, QC-method choice gives ARI ~0.80–0.88 (not ~0.95+). QC method perturbs clustering *a little, everywhere* — a low-grade uniform effect, not an on/off edge-case effect.

3. **Divergent cells don't reorganize the rest.** Where MAD5 keeps **11,300 more** (low-gene) cells than C2 on small_intestine (MAD3 +10,534; C1 +6,099), the **shared** cells still co-cluster at ARI 0.84. The extra borderline-quality cells are added without massively restructuring the high-quality core (`method_specific_cell_disposition.tsv`).

4. **Annotation is robust.** Cell-type label *counts* are nearly identical across methods (blood 19 all; liver 10–11; small_intestine 23–24), and per-cell annotation agreement is 88–98%. QC method choice rarely changes the cell-type call.

## Cross-audit pattern (the third propagation test)
| audit | technical choice | edge-case tissue | clean tissue |
|-------|------------------|------------------|--------------|
| Audit 3 C3 | cell-calling | strong (ARI ↓ to ~0.5) | near-null |
| Ambient CP3 | ambient correction | strong (ARI 0.50–0.70) | near-null (0.85–0.90) |
| **QC-MAD-prop CP1** | **QC filtering method** | **modest (ARI 0.84–0.87)** | **modest (ARI 0.80–0.88)** |

**The pattern BREAKS for QC method.** Cell-calling and ambient correction reshape *which/what* counts exist and scale with ambient burden, so they amplify on edge-case tissue. QC filtering method only decides *which borderline low-quality cells to drop* — a few % of low-information cells near thresholds, whose removal/retention nudges cluster boundaries uniformly but isn't amplified by tissue complexity. This is the spec's anticipated "if it breaks, that's also informative" outcome and is a genuine corpus-level contrast: **not all scRNA technical choices propagate the same way.**

## §5.3.2 / sample-size tiering
Existence-of-a-modest-effect (QC method perturbs clustering ~uniformly, ARI 0.80–0.91, tight CIs): supported, but it is a *non-amplification* finding across n=3 datasets (1 edge + 1 mid + 1 control). The notable result — that QC method does NOT show the edge-case amplification — is the headline; generalization flag_and_warn (n=3, same structure as C3/CP3). To be set in CP2 synthesis.

## Reproduction integrity
QC-MAD's exact cells reproduced via Census re-pull (position→soma_joinid mapping; median QC metrics matched QC-MAD to precision: liver 2152, small_intestine 639, blood 1517). All 4 QC methods reproduced QC-MAD's saved cell counts EXACTLY (every C1/C2/MAD3/MAD5 count, all 3 datasets). 12/12 pipeline runs succeeded (0 failures).

## Heumos positioning
Heumos recommends MAD filtering (Germain 2020) over fixed thresholds without testing downstream consequences. QC-MAD showed the methods give equivalent cell sets in most cases; CP1 shows that even where they differ most (small_intestine), the **biological** consequence is modest and not tissue-amplified. So the MAD-vs-default choice is **methodological hygiene that largely washes out downstream** — extends Heumos by quantifying that the choice is biologically low-stakes, unlike cell-calling/ambient correction.

## Outputs
`per_comparison_metrics.tsv` (per dataset×pair: ARI r1.0/0.5/1.5 + CI, NMI, marker Jaccard, annotation agreement), `per_stratum_bootstrap.tsv`, `method_specific_cell_disposition.tsv`, `per_dataset_pipeline_outputs/*.h5ad`, `per_condition/<ds>/<method>.{obs,markers,summary}`. Hash-registered in repro.lock.

## Scope
CP1 only. 3 datasets; C2 reference; cp6 pipeline params; no CP2 synthesis/rules. Does not modify QC-MAD or other prior audits.
