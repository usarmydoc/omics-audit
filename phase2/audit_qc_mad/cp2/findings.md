# QC-MAD Audit Findings

_Closed CP2: 2026-05-23. Lock: `phase2/repro.lock`. CP1 metrics: `audit_qc_mad/cp1/`._
_prior_audit_relationship: extends_prior (Phase 1 p1 mito-threshold QC)._

## Question

Does the choice between MAD-based outlier detection (Heumos 2023, citing
Germain 2020 / pipeComp) and quantile-based filtering produce materially
different cell sets in scRNA-seq QC?

## Working set

8 Census datasets (cellxgene_census, stable 2025-11-08; exact subset of the
Phase 1 p1 working set), spanning median-mito 0.07–10.83%, human + mouse:
heart, bone_marrow, liver, blood, lung, small_intestine, pancreas,
large_intestine. Each capped at 50,000 cells (seed 0). Blood + bone_marrow are
the immune/PBMC-like low-mito anchors; small/large_intestine + pancreas are the
high-mito anchors. Working set: `cp1/working_set.tsv`.

## Methods compared

- **C1 (quantile-data-driven):** n_genes < 5th pct OR total_counts < 5th pct
  OR pct_mt > 95th pct (all per-dataset).
- **C2 (quantile-fixed-floors):** n_genes < 200 OR total_counts < 500 OR
  pct_mt > 95th pct (typical scanpy defaults + per-dataset mito).
- **MAD k=3:** `scuttle::isOutlier(nmads=3)` on n_genes + total_counts
  (log=TRUE, lower tail) and pct_mt (raw, upper tail); outlier on any → filtered.
- **MAD k=5:** same, nmads=5.

6 pairwise comparisons. Bootstrap B=1000, dataset-level resampling (§2.4).

## Headline finding — CONDITIONAL tier (§5.3.2 equivalence)

**All 4 filtering methods produce largely equivalent cell sets across 8 Census
datasets.** Per-pair mean Jaccard 0.903–0.958 (medians 0.895–0.969); all
bootstrap CIs overlap; every per-pair CI lower bound ≥ 0.88.

| pair | mean Jaccard | 95% CI |
|---|---|---|
| C2 vs MAD3 | 0.938 | [0.880, 0.983] |
| C2 vs MAD5 | 0.938 | [0.887, 0.973] |
| MAD3 vs MAD5 | 0.958 | [0.939, 0.973] |
| C1 vs C2 | 0.931 | [0.911, 0.948] |
| C1 vs MAD3 | 0.910 | [0.882, 0.928] |
| C1 vs MAD5 | 0.903 | [0.893, 0.916] |

**Tier rationale (conditional, not hard_default):** the numerical bar for
§5.3.2 hard_default is met (≥8 datasets ✓, multiple tissues ✓, CIs overlap ✓,
mean ρ-equivalent ≥0.90 ✓, CI lower bounds ≥0.80 ✓). But hard_default also
requires equivalence to hold *uniformly, without exception*, and there is a
documented exception: small_intestine (C2 vs MAD5 = 0.774, C2 vs MAD3 = 0.786
per-dataset). The exception is feature-dependent (dataset gene-count
distribution, see Sub-finding 3). Strong equivalence + documented
feature-dependence is the §5.3.2 **conditional** case by definition. Tiering
this hard_default would require treating small_intestine as noise rather than a
real feature-exception, which the mechanism contradicts.

## Sub-finding 1 — C2 ≈ MAD3 (conditional, §5.3.2 equivalence)

C2 (typical scanpy defaults) is the *closest* method to MAD3 (the Heumos
recommendation): mean Jaccard 0.938, median 0.969, CI [0.880, 0.983]. **Common
fixed-floor defaults already approximate MAD-based outlier detection** more
closely than pure-quantile filtering does. Same conditional tier + same
small_intestine exception (0.786) as the headline.

## Sub-finding 2 — C1 is the methodologically distinct option (descriptive)

Pure-quantile C1 sits further from both C2 and MAD than they sit from each
other (C1 vs MAD3 0.910, C1 vs MAD5 0.903, C1 vs C2 0.931, vs C2 vs MAD3 0.938).
Comparative observation across the 3-way space, not a recommendation — no rule.

## Sub-finding 3 — disagreement driver is gene-count distribution, not tissue/mito (conditional, §5.3.2 feature-dependence)

Cross-method disagreement tracks **dataset gene-count distribution**, not
tissue or mitochondrial fraction. C1-vs-MAD3 Jaccard is flat across the
mito range (0.82–0.94); blood/PBMC (0.918) is indistinguishable from high-mito
tissues. The largest gaps occur where many cells are genuinely low-gene
(small_intestine: C2's fixed 200-gene/500-count floor removed 11,300 cells that
MAD5 kept). Disagreement concentrates in low_counts + low_genes cells, not
high_mito. Robust across datasets but generalizes via the mechanism; needs
replication on more low-gene-cell datasets for full generalization.

## Sub-finding 4 — MAD k=3 vs k=5 (descriptive)

MAD3 filters 2–3× more than MAD5 (liver 9,035 vs 4,768; heart 4,243 vs 1,595);
MAD5 filters ~0 on 3 datasets (bone_marrow, small/large_intestine). Confirms
Heumos's caveat that MAD is aggressive at low k. Parameter-sensitivity
observation; no rule unless paired with k-selection guidance this audit didn't test.

## Practical implication

- Common fixed floors (C2) and Heumos-recommended MAD3 produce substantially
  the same cell sets on most data (median Jaccard 0.969); switching is unlikely
  to change conclusions.
- Pure-quantile (C1) is the distinct option, differing from both.
- The choice matters most on datasets with many low-gene cells (small_intestine):
  fixed-floor harshness removes cells MAD retains.
- MAD k=3 aggressive, k=5 permissive — k choice materially changes aggressiveness.

## Heumos 2023 positioning — extends

Heumos recommends MAD-based filtering (Germain 2020) over fixed manual
thresholds as more robust to sample-level variation, without quantifying the
cell-set difference. This audit shows the difference from *typical defaults*
(C2) is small in most cases (median Jaccard 0.969) but larger on specific
dataset features (low-gene-cell populations). MAD is empirically distinct from
pure quantiles (C1), but the practical magnitude of switching from common
fixed-floor defaults (C2) to MAD is **smaller than the qualitative framing
suggests.**

## Sample size honesty (§3.1)

- n=8 datasets, 8 tissues (n=1 per tissue → tissue-specific claims not powered;
  tissue is used only as a mito-level proxy, and the finding is that mito level
  does *not* drive disagreement).
- One Census dataset per tissue; generalization assumes Census representativeness.
- The "low-gene-cell distribution drives disagreement" mechanism rests on
  small_intestine as the cleanest demonstration; needs replication on more
  low-gene tissues for a stronger generalization claim.
- §3.5 sanity check applied: filter fractions all in sane ranges; MAD5≈0 on 3
  datasets verified as genuine permissiveness (5 MADs), not degenerate failure.

## Limitations

- Census matrices are author-pre-filtered (already-deposited cells), not raw
  barcode×gene matrices — this audit compares QC *re-filtering* of deposited cells.
- 4 methods, not exhaustive (median+N*MAD hybrids untested).
- No downstream propagation tested (do QC-method differences cascade to
  clustering/DE? — future-audit candidate, like Audit 3 C3 did for cell-calling).

## Rules to draft (CP3 candidates — operator decision pending)

1. **scrna_qc_filtering_method_equivalence** (conditional, §5.3.2; severity
   info): all 4 methods produce largely equivalent cell sets (pair Jaccard
   0.90–0.97); QC-method choice is a modest variance source. Reassurance rule.
2. **scrna_qc_low_gene_dataset_caution** (conditional, §5.3.2 feature-dependence;
   severity warn): on datasets with many low-gene cells, fixed-floor methods (C2)
   can be harsh; prefer MAD-based or permissive-quantile filtering. Actionable.

Sub-findings 2 + 4 stay descriptive (no rule).
**Recommendation: draft both** — they mirror Audit 3's general-rule + exception-rule
pattern (the equivalence + the documented carve-out). Operator's call.
