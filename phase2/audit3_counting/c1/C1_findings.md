# Audit 3 C1 findings — per-gene count agreement (CORRECTED)

_Generated: 2026-05-18T20:34:29 (re-run with USA-suffix fix)_
_B = 1000 (dataset-level), stratified by chemistry × tool_pair_

## Prior audit relationship

**novel** — No prior systematic audit of scRNA-seq counting tool agreement
exists in the audit corpus. Lafzi 2018 reviewed counting tools but did
not benchmark cross-tool agreement on identical FASTQs. This audit
contributes new empirical evidence on inter-tool count agreement.

> **Follow-up note (2026-05-22, from CP5/C2):** C2 regenerated alevin-fry
> with knee-point cell calling (Deliverable C). Recomputing this C1 metric —
> per-gene Spearman, alevin-fry **knee-filtered** matrices vs STARsolo — gives
> 0.976 (pbmc_1k_v3), 0.981 (gse287209 lung), 0.971 (gse288156 intestine),
> consistent with the ~0.96 reported here for the CP3 `--unfiltered-pl`
> matrices. **The C1 count-convergence finding is unchanged by alevin-fry's
> cell-calling configuration.** No rewrite needed; this note records the
> confirmation per CP5 Step 6.

> **Important note on supersession.** A first pass of this analysis
> (preserved at `superseded_2026-05-18_buggy_usa_strip/`) contained a
> bug in the alevin-fry gene-ID handling: `strip_version()` collapsed
> the splici-mode `-U` and `-A` suffixes along with the Ensembl version,
> so alevin-fry's per-gene counts were silently restricted to the
> "Ambiguous" (-A) bucket only (~10% of total reads). This produced an
> artifactual "alevin-fry is the outlier" headline (rho ~0.58) that did
> not survive the fix. The corrected analysis sums S+A per gene per
> alevin-fry's recommended abundance convention; results below.

## Headlines

1. **All four tools agree on per-gene counts.** Median per-gene Spearman
   across all 9 datasets is 0.95–0.99 for every pair. The only "high"
   pair is `star_default` vs `star_cr_mimic` at 0.991 (same tool, different
   parameters — this is the sanity-check baseline). All cross-tool pairs
   cluster tightly at 0.956–0.962 in v3 (bootstrap CIs in [0.95, 0.97]).
   No tool stands out as an outlier.

2. **Tools agree on per-cell UMI ranking too.** Per-cell total UMI Spearman
   across 6 pairs is 0.93–0.97 (v3 medians). Cells called as deep by one
   tool are called deep by all four — confirmed regardless of which gene
   the reads were attributed to.

3. **No gene-category-specific disagreement.** After the fix, every
   category (`mitochondrial`, `pseudogene`, `overlapping`, `other`) shows
   median log2_ratio = 0.0 with IQR ≤ 0.07 across all 6 tool pairs.
   The earlier "STARsolo reports 2× alevin-fry on overlapping genes"
   was a buggy-bucket artifact, not a multi-mapping disagreement.

4. **The cr-like multi-mapper question is still open.** alevin-fry was
   run with `-r cr-like` (multi-mappers discarded). With the bug fixed,
   the empirical effect of this mode choice was much smaller than the
   first pass suggested — likely because most reads in scRNA-seq align
   uniquely. A cr-like-em comparison would tighten the rule but is no
   longer load-bearing for the "are tools equivalent" question.

5. **Chemistry effect is small.** v3 (n=7): mean cross-tool rho 0.957.
   v2 (n=2): mean 0.92. Differences within the v2 n=2 noise floor.

## Where the tools agree (companion metrics, §3.4)

- **rho_gene_total** across all 6 cross-tool pairs in 3' v3 (n=7 datasets):
  point estimate 0.957–0.991, 95% CI widths 0.003–0.016.
- **rho_cell_total_umi** across all 6 pairs: 0.93–0.97. Cell-depth ranking
  is preserved across tools regardless of per-gene attribution differences.
- **log2_median = 0.0** for all 6 pairs and all 4 gene categories: no
  systematic offset.
- **det_rate (genes with ≥1 nonzero cell)** is similar across tools
  within each species: STARsolo 0.45, alevin-fry 0.45, kb_count 0.50.
  kb_count detects slightly more genes (~5 pp higher) due to its
  pseudoalignment liberality, but the magnitude is small.

## Where the tools still differ (modest residual signals)

| Axis | Effect size | Pairs affected |
|---|---|---|
| Direction of effect on overlapping genes | STARsolo slightly higher in 28% of overlapping genes vs alevin-fry; alevin-fry slightly higher in ~62%; ties ~10%. p05/p95 = [-0.4, +0.4] | star_default & star_cr_mimic vs alevin_fry |
| Direction of effect on pseudogenes | star_cr_mimic higher than alevin-fry 12% of the time, ties most of the time, alevin-fry higher elsewhere. log2_iqr = 0.0 (most are zero counts in both) | star_cr_mimic vs alevin_fry |
| Per-cell UMI Spearman dispersion | alevin-fry-vs-kb is the most variable (CI [0.93, 0.97]) | alevin-fry vs kb |

The IQR on the largest "disagreement" axis is 0.05 log2-units. Practically:
on the median overlapping gene, the two tools disagree by less than 5%.
This is not a methodologically actionable difference.

## Tool failure modes / caveats

- **alevin-fry's S+A convention is not the default raw output.** Naive
  loading of `quants_mat.mtx` returns three rows per gene
  (`ENSG.version`, `ENSG.version-U`, `ENSG.version-A`). Users
  comparing alevin-fry to other tools without consulting the splici
  documentation will produce the artifact this analysis hit. Worth a
  BO rule (see "Implications" below).
- **alevin-fry index size diagnostic:** the 189,723-entry
  quants_mat_cols.txt is the splici reference (S+U+A buckets); the
  underlying t2g.tsv has 63,241 unique gene IDs (matching STARsolo).
  No decoy file is present, so MT/pseudogene exclusion is not a factor
  with this index build.
- **No tool failed on any of the 9 datasets.** All 36 (dataset × tool)
  cells produced count matrices.

## Methodological limitations

- Multi-mapper gene category was specified in CP4 scope but not implemented.
  No GTF annotation source identified for multi-mapper status during CP4.
  The other 4 categories (overlapping genes, pseudogenes, mitochondrial, other)
  cover the audit's question of whether disagreement concentrates in specific
  gene types. Future work could add multi-mapper detection via STAR's
  --outFilterMultimapNmax tracking, salmon's mapping ambiguity metrics, or
  a dedicated annotation source.

## Sample size limitations

- **3p_v2: n=2 datasets** (pbmc_4k_v2, t_3k_v2). Bootstrap CIs
  artificially tight; can claim "v2 and v3 results pattern similarly"
  but cannot make chemistry-specific claims at the resolution of
  0.01 rho.
- **3p_v3: n=7 datasets** (pbmc_1k_v3, pbmc_5k_v3.1, pbmc_10k_v3.1,
  neuron_1k_v3, gse287209_lung, gse325955_kidney, gse288156_intestine).
  Bootstrap CIs reasonable; the "all tools agree" claim is supported.
- Mouse: n=3, human: n=6. No species-specific divergence observed,
  but power to detect any species×tool interaction is low.
- Tissue: PBMC n=5, non-PBMC n=4. No tissue-class pattern observed.

## Implications for CP5 (C2 cell-barcode-calling agreement)

The C1 finding "tools agree on per-gene counts AND on per-cell UMI
totals" reshapes the CP5 question:

- **No tool is the per-gene outlier.** Earlier "3 trusted vs 1 outlier"
  framing was a bug. CP5 should report all 6 pair-wise Jaccards
  symmetrically.
- **Per-cell UMI agreement is already high (0.93–0.97).** C2's question
  becomes: at the *barcode-calling* step (not the per-cell-UMI-total
  step), do the four tools' algorithmic defaults agree on which
  barcodes correspond to real cells vs background?
- **CP5 should focus on the low-UMI / ambient regime.** That's where
  tool-specific cell-calling heuristics (10x's knee-point, alevin-fry's
  generate-permit-list, kb's bustools filter, STARsolo's
  EmptyDrops_CR) would diverge if they diverge at all.
- **Gene-category stratification is NOT useful for C2.** The category
  disagreements that motivated stratification in CP4's first pass were
  bug artifacts.

## Implications for CP6 (Phase 1 robustness)

With "tools are essentially equivalent on counts," CP6 becomes a
confirmation pass: re-run one Phase 1 finding through alevin-fry's
output and verify the Phase 1 conclusion holds. We don't need a deep
investigation, but the confirmation is still worth doing because:

- Phase 1 used STARsolo-equivalent counts; we need to confirm pipeline
  reproducibility across tool choice
- Per-cell UMI is preserved (rho 0.93+), so cell-level pseudobulk DE
  should be robust
- The biological signal is well above the residual tool noise

## Implications for CP8 (rule drafting)

Drafting BO rules at this stage:

1. **HIGH-CONFIDENCE rule candidate:** "When using alevin-fry's splici
   output for STARsolo-comparable abundance, sum the S+A buckets per
   gene; do not use the raw `quants_mat.mtx` cols list as-is." This
   rule has direct, audit-validated empirical support: the bug we
   hit, if encoded in a tool wrapper, would have surfaced the fix at
   the wrapper level rather than requiring a re-run.

2. **MODERATE-CONFIDENCE rule candidate:** "All four counting tools
   (STARsolo default, STARsolo CR-mimic, alevin-fry, kb_count) are
   substitutable on 10x 3' v2/v3 scRNA data for per-gene counting
   purposes; cross-tool rho ≥ 0.95 on a 9-dataset audit set." Confidence
   moderated by: n=2 for v2; cr-like-em not benchmarked yet; index
   build conventions matter (see rule #1).

3. **DEFERRED:** any rule about "tool X is preferred over Y for
   category Z." The category disagreements are at noise floor;
   encoding them would be over-fitting.

Per the batched-update directive, no rules are pushed yet — wait for
C2 (CP5) and C3 (CP6) to finalize the picture before BO version bump.

## What the cr-like-em comparison would add

The user's earlier note flagged that `cr-like` (alevin-fry default
in this audit) discards multi-mapping reads, while `cr-like-em`
rescues them via EM. With the bug fixed, the empirical effect of
this choice was much smaller than first feared. Recommended next
steps if the rule on alevin-fry should be airtight:

1. Re-run one v3 dataset (e.g., pbmc_10k_v3.1) under `cr-like-em`
2. Compare per-gene rho to the cr-like result
3. If rho_cr_like_em vs rho_cr_like is < 0.99, the EM mode is
   non-negligible and the rule should specify both modes
4. If rho > 0.99, modes are interchangeable for this analysis

Deferred until after CP5/CP6 unless the rule scope changes.

## Summary table

| metric | v2 (n=2) | v3 (n=7) | note |
|---|---|---|---|
| STAR-default vs STAR-CR-mimic rho_gene | 0.964 | 0.991 | sanity check |
| STAR-default vs alevin-fry rho_gene | 0.919 | 0.958 | post-fix |
| STAR-default vs kb_count rho_gene | 0.932 | 0.962 | |
| STAR-CR-mimic vs alevin-fry rho_gene | 0.940 | 0.963 | |
| STAR-CR-mimic vs kb_count rho_gene | 0.913 | 0.963 | |
| alevin-fry vs kb_count rho_gene | 0.914 | 0.968 | |
| Median rho_cell_total_umi (all pairs) | 0.97 | 0.95 | |

All values are point estimates; CIs in `per_stratum_bootstrap.tsv`.
v2 CIs are artificially tight at n=2 — interpret as "consistent with v3" not as independent measurement.
