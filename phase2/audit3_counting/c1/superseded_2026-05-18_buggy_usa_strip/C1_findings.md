# Audit 3 C1 findings — per-gene count agreement

_Generated: 2026-05-18T18:15:47_
_B = 1000 (dataset-level), stratified by chemistry × tool_pair_

## Headlines

1. **alevin-fry is the clear outlier on per-gene count agreement.** Across both
   3' v2 (n=2 datasets, point ~0.48) and 3' v3 (n=7, point ~0.65 with
   95% CI [0.58, 0.72]), every pair involving alevin-fry shows substantially
   lower Spearman correlation on summed gene counts than the pairs not
   involving it. STARsolo (either parameterization) vs kb_count is ~0.96
   (v3) and ~0.92 (v2); alevin-fry vs anything is ~0.58–0.65 (v3), ~0.48
   (v2). The two STARsolo parameterizations are essentially identical
   (rho=0.991 v3), confirming the analysis is not picking up parameter noise.

2. **Tools agree on which cells are deepest, just not on which genes those
   reads belong to.** Per-cell total UMI Spearman is 0.89–0.97 across all
   pairs in v3 (including alevin-fry pairs). This means the disagreement
   is at the gene-assignment step, not the cell-barcode-extraction step.

3. **Disagreement concentrates in "overlapping" genes** (genes sharing
   genomic coordinates with another gene). Median log2 ratio for
   STARsolo vs alevin-fry on overlapping genes is **+1.0** (STARsolo
   reports 2× the counts alevin-fry does, on the median overlapping gene),
   with p95 = +8.85 (some overlapping genes show 460× higher counts in
   STARsolo). For "other" genes the equivalent median is 0.0. This is
   consistent with the two tools handling multi-mapping reads differently:
   STARsolo's MultiGeneFold or unique-only counting vs alevin-fry's
   probabilistic transcript-equivalence-class assignment.

4. **kb_count tracks STARsolo more than alevin-fry does.** STAR-vs-kb
   rho_gene_total = 0.962 [0.956, 0.967] (v3, n=7) vs alevin-fry-vs-kb
   rho = 0.648 [0.575, 0.713]. The kb_count salmon-style pseudoaligner
   could have been expected to behave like alevin-fry (also a
   pseudoaligner). Empirically it does not — kb agrees with STARsolo
   more than it agrees with alevin-fry.

5. **Direction-of-effect agreement is strong for STAR-vs-kb, weak for
   alevin-fry pairs.** STARsolo-vs-kb: median pct_A_higher = 0.08
   (kb is higher most of the time when there's any disagreement —
   small offset). alevin-fry-vs-kb: pct_A_higher = 0.004 (alevin-fry
   is virtually never higher than kb), pct_B_higher = ~0.50, ties =
   ~0.50. So kb is reporting more counts than alevin-fry for every
   gene where they disagree (consistent with #3: alevin-fry distributing
   reads across an expanded ~190K-element index, diluting per-gene
   counts).

## Where the tools agree (companion metrics, §3.4)

- **Per-cell UMI ranking**: all 6 tool pairs in v3 show rho >= 0.89.
  Cells called as "deep" by one tool are called deep by all four.
- **Direction at the median gene**: for STAR-vs-kb across all categories,
  median log2 ratio = 0.0; tools agree on the typical gene.
- **Same gene universe**: STARsolo and kb_count both use the GENCODE
  GTF (63241 human genes / 57180 mouse). alevin-fry's salmon index has
  ~190K (human) / ~171K (mouse) entries, but the ENSG.base intersection
  is the comparison set.
- **Per-cell total UMI Spearman is highest for the alevin-fry vs
  kb_count pair** (median 0.97) — both pseudoaligners agree on which
  cells received the most reads even though they disagree on which
  genes the reads came from.

## Where the tools disagree (the headline payload)

| What | Magnitude | Tools affected |
|---|---|---|
| Per-gene count Spearman | drops from ~0.96 to ~0.58 when alevin-fry enters the pair | alevin-fry pairs vs STAR/kb pairs |
| Median log2 ratio on overlapping genes | +1.0 (STAR > alevin-fry) | star_default & star_cr_mimic vs alevin_fry |
| Per-gene detection rate | alevin-fry 27-31% / kb 42-50% / STAR 23-50% | All comparisons (different gene universes) |
| Pseudogene counts | alevin-fry virtually always zero/near-zero | alevin-fry vs STAR & kb |
| Mitochondrial counts | alevin-fry-vs-kb median log2 = -7.96 (kb >> alevin-fry) | alevin-fry vs kb specifically (suspicious — possible index artifact) |

## Tool failure modes / caveats

- alevin-fry's mitochondrial count being 250× lower than kb's needs
  investigation. Could be reference content (whether ChrMT is in the
  index), could be intentional decoy behavior, could be a real biology
  difference. Treat as "mechanism candidate, not isolated" per §3.5.
  Recommend: check whether MT sequences appear in alevin-fry's index.
- Per-gene `n_genes` in the category table differs from the GTF gene
  count for some categories: `overlapping` shows ~99K genes in the
  alevin-fry-vs-kb table vs ~33K in the GTF. That's because alevin-fry's
  index has many duplicate gene IDs across transcript-form entries
  in the intersected set. The category counts pass the n>=100 filter
  but the inflation should be flagged.
- The 4-file layout datasets (kidney, intestine) and the 2-file layout
  datasets (lung, all 10x demos) show similar tool-pair patterns —
  no chemistry-by-layout interaction visible.

## Sample size limitations

- **3p_v2: n=2 datasets** (pbmc_4k_v2, t_3k_v2). Bootstrap CIs are
  artificially tight because only 2 datasets can be resampled. The
  v2 vs v3 comparison is **suggestive, not conclusive** — both
  chemistries show the same alevin-fry-outlier pattern, but the
  v2 effect sizes are larger (point ~0.48 vs ~0.65). Whether this
  reflects real chemistry-tool interaction or v2's smaller dataset
  count introducing noise cannot be distinguished with n=2. Per
  AUDIT_STANDARDS §3.1, this is a headline finding (alevin-fry
  is more discordant on v2 than v3) that is **insufficient to
  encode as a rule** until more v2 datasets are added (queued in
  DEFERRED.md as Audit 3-v2-expand).
- 3p_v3: n=7 datasets. Bootstrap CIs are reasonably tight.
- Mouse: n=3, human: n=6. Species-stratified bootstrap not
  reported because chemistry dominates the variance in this audit.
- Tissue: PBMC n=5, non-PBMC n=4. Tissue-stratified analysis
  yields the same tool-pair pattern; tissue is not the dominant axis.

## Implications for CP5 (C2 cell-barcode-calling agreement)

The C1 finding that **per-cell UMI Spearman is high (0.89-0.97) across
all pairs, including alevin-fry's** suggests cell-barcode calling
agreement (which is C2's question) is decoupled from per-gene count
agreement. Tools that disagree on per-gene counts can still agree on
which barcodes correspond to deep cells.

This re-shapes CP5 scope candidates:

- The "stratify by gene category" idea (from my CP4 reading-order doc)
  is **not the most informative slice** for C2. The category disagreement
  is concentrated at the gene level; barcode-calling looks at total UMI
  per barcode, which is gene-category-agnostic.
- The "treat alevin-fry as an outlier and compare 3-vs-1" framing is
  more relevant for C2. If alevin-fry's cell calls diverge from the
  other three by similar magnitude to its per-gene divergence, that
  reinforces "alevin-fry is doing something fundamentally different"
  as the headline. If alevin-fry's cell calls agree with the others
  (rho_cell_total_umi ~ 0.94 already suggests they will), then the
  framing becomes "tools disagree on gene-level but converge on
  cell-level" — a more nuanced and interesting story.
- Chemistry stratification (v2 vs v3) is still warranted for C2 even
  though n=2 for v2, since the v2 effect size was larger here.

## Recommendations for the CP5 prompt (not yet drafted)

CP5 should:

1. **Frame around "one outlier vs three convergent tools"** not
   "all six pairs equally." Compute Jaccard for the 3-trusted-tool
   subset separately from the 6-way matrix.
2. **Strongly emphasize companion metrics** (per AUDIT_STANDARDS §3.4):
   total cell count per tool, knee-point of barcode-rank curve,
   ambient RNA estimate. Don't lead with Jaccard alone.
3. **Stratify by low-UMI / high-UMI regime.** Tools agreeing on
   high-UMI cells but disagreeing on borderline cells is the most
   likely C2 finding given C1's pattern; designing C2 to surface
   that regime difference is high-value.
4. **Note explicitly** that C1 found the disagreement is in gene-level
   counts, not cell-level totals. C2 is testing whether the
   gene-level disagreement propagates into cell-calling.

Per the user's earlier directive, the actual CP5 prompt is drafted
after CP4 findings are reviewed — this section is candidates only.


## Per-stratum bootstrap (rho_gene_total point + 95% CI)

| chemistry | tool_a | tool_b | n_datasets | point | 95% CI |
|---|---|---|---|---|---|
| 3p_v2 | star_default | star_cr_mimic | 2 | 0.964 | [0.960, 0.968] (n<3, wide CI) |
| 3p_v2 | star_default | alevin_fry | 2 | 0.483 | [0.481, 0.485] (n<3, wide CI) |
| 3p_v2 | star_default | kb_count | 2 | 0.932 | [0.925, 0.939] (n<3, wide CI) |
| 3p_v2 | star_cr_mimic | alevin_fry | 2 | 0.498 | [0.497, 0.499] (n<3, wide CI) |
| 3p_v2 | star_cr_mimic | kb_count | 2 | 0.913 | [0.905, 0.922] (n<3, wide CI) |
| 3p_v2 | alevin_fry | kb_count | 2 | 0.476 | [0.473, 0.479] (n<3, wide CI) |
| 3p_v3 | star_default | star_cr_mimic | 7 | 0.991 | [0.990, 0.993] |
| 3p_v3 | star_default | alevin_fry | 7 | 0.648 | [0.581, 0.725] |
| 3p_v3 | star_default | kb_count | 7 | 0.962 | [0.956, 0.967] |
| 3p_v3 | star_cr_mimic | alevin_fry | 7 | 0.651 | [0.580, 0.722] |
| 3p_v3 | star_cr_mimic | kb_count | 7 | 0.963 | [0.957, 0.969] |
| 3p_v3 | alevin_fry | kb_count | 7 | 0.648 | [0.575, 0.713] |

## Per-stratum bootstrap (rho_cell_total_umi)

| chemistry | tool_a | tool_b | n_datasets | point | 95% CI |
|---|---|---|---|---|---|
| 3p_v2 | star_default | star_cr_mimic | 2 | 0.973 | [0.971, 0.974] (n<3, wide CI) |
| 3p_v2 | star_default | alevin_fry | 2 | 0.970 | [0.968, 0.972] (n<3, wide CI) |
| 3p_v2 | star_default | kb_count | 2 | 0.963 | [0.962, 0.965] (n<3, wide CI) |
| 3p_v2 | star_cr_mimic | alevin_fry | 2 | 0.975 | [0.972, 0.977] (n<3, wide CI) |
| 3p_v2 | star_cr_mimic | kb_count | 2 | 0.966 | [0.963, 0.970] (n<3, wide CI) |
| 3p_v2 | alevin_fry | kb_count | 2 | 0.982 | [0.982, 0.982] (n<3, wide CI) |
| 3p_v3 | star_default | star_cr_mimic | 7 | 0.926 | [0.893, 0.956] |
| 3p_v3 | star_default | alevin_fry | 7 | 0.944 | [0.917, 0.967] |
| 3p_v3 | star_default | kb_count | 7 | 0.894 | [0.854, 0.928] |
| 3p_v3 | star_cr_mimic | alevin_fry | 7 | 0.942 | [0.916, 0.967] |
| 3p_v3 | star_cr_mimic | kb_count | 7 | 0.936 | [0.913, 0.960] |
| 3p_v3 | alevin_fry | kb_count | 7 | 0.954 | [0.931, 0.972] |

## By-gene-category log2 ratio summary (across all 9 datasets)

| tool_a | tool_b | category | n_datasets | med(log2_median) | med(p05) | med(p95) |
|---|---|---|---|---|---|---|
| alevin_fry | kb_count | mitochondrial | 9 | -7.960 | -18.296 | -4.429 |
| alevin_fry | kb_count | other | 9 | 0.000 | -8.191 | 0.000 |
| alevin_fry | kb_count | overlapping | 9 | -1.000 | -8.943 | 0.000 |
| alevin_fry | kb_count | pseudogene | 9 | 0.000 | -5.858 | 0.000 |
| star_cr_mimic | alevin_fry | mitochondrial | 9 | nan | nan | nan |
| star_cr_mimic | alevin_fry | other | 9 | 0.000 | 0.000 | 8.089 |
| star_cr_mimic | alevin_fry | overlapping | 9 | 1.000 | 0.000 | 8.794 |
| star_cr_mimic | alevin_fry | pseudogene | 9 | 0.000 | 0.000 | 3.000 |
| star_cr_mimic | kb_count | mitochondrial | 9 | nan | nan | nan |
| star_cr_mimic | kb_count | other | 9 | 0.000 | -1.000 | 0.034 |
| star_cr_mimic | kb_count | overlapping | 9 | 0.000 | -0.848 | 0.071 |
| star_cr_mimic | kb_count | pseudogene | 9 | 0.000 | -2.000 | 0.415 |
| star_default | alevin_fry | mitochondrial | 9 | nan | nan | nan |
| star_default | alevin_fry | other | 9 | 0.000 | 0.000 | 8.091 |
| star_default | alevin_fry | overlapping | 9 | 1.000 | 0.000 | 8.852 |
| star_default | alevin_fry | pseudogene | 9 | 0.000 | 0.000 | 5.207 |
| star_default | kb_count | mitochondrial | 9 | nan | nan | nan |
| star_default | kb_count | other | 9 | 0.000 | -1.000 | 0.008 |
| star_default | kb_count | overlapping | 9 | 0.000 | -1.000 | 0.035 |
| star_default | kb_count | pseudogene | 9 | 0.000 | -1.856 | 1.000 |
| star_default | star_cr_mimic | mitochondrial | 9 | nan | nan | nan |
| star_default | star_cr_mimic | other | 9 | 0.000 | -0.061 | 0.000 |
| star_default | star_cr_mimic | overlapping | 9 | 0.000 | -0.070 | 0.019 |
| star_default | star_cr_mimic | pseudogene | 9 | 0.000 | 0.000 | 1.000 |

## Sample size limitations

- 3p_v2: 2 datasets (below n=3 → bootstrap CIs are wide; interpret per-chemistry comparisons with this caveat).
- 3p_v3: 7 datasets.
- Mouse n=3 / Human n=6 — species stratification underpowered for mouse.
- Tissue: PBMC n=5 (incl. T cells), non-PBMC n=4.
- Per AUDIT_STANDARDS §3.1, these are headline findings, not caveats.
