# Audit QC-MAD — scRNA-seq Low-Quality Cell Filtering Methods

_Closed: 2026-05-23. Canonical drive `/mnt/nvme1/omics-audit/`. Lock: `phase2/repro.lock`._
_prior_audit_relationship: extends_prior (Phase 1 p1 mito-threshold QC)._

## Question

Does the choice between MAD-based outlier detection (Heumos 2023, citing
Germain 2020 / pipeComp) and quantile-based filtering produce materially
different cell sets in scRNA-seq QC?

## Working set

8 Census datasets (cellxgene_census, stable 2025-11-08; exact subset of Phase 1
p1), median-mito 0.07–10.83%, human + mouse: heart, bone_marrow, liver, blood,
lung, small_intestine, pancreas, large_intestine. Blood + bone_marrow are the
immune/PBMC-like anchors; small/large_intestine + pancreas the high-mito
anchors. Each capped at 50,000 cells (seed 0). Details: `cp1/working_set.tsv`.
Source matches Phase 1 p1 methodology.

## Finding 1 — Methods produce largely equivalent cell sets (Rule 1, conditional/info)

All 4 filtering methods (C1 pure-quantile, C2 fixed-floor, MAD k=3, MAD k=5)
produce per-pair mean Jaccard 0.90–0.96 (median 0.895–0.969) across 8 datasets;
all dataset-level bootstrap CIs overlap, all CI lower bounds ≥0.88. Method
choice is a modest variance source, not a major one. **Rule:**
`scrna_qc_filtering_method_equivalence` (conditional, info). Detail:
`cp2/findings.md#headline-finding`.

## Finding 2 — C2 ≈ MAD3 (the counter-intuitive result)

Typical scanpy defaults (C2: 200-gene floor, 500-count floor, 95th-pct mito)
approximate Heumos-recommended MAD k=3 closely — median Jaccard **0.969**, mean
0.938, CI [0.880, 0.983]. Conventional fixed-floor practice is already close to
MAD-based detection in cell-set terms. MAD is empirically distinct from pure
quantiles (C1), but the practical magnitude of switching from typical defaults
to MAD is small.

## Finding 3 — Disagreement driver is gene-count distribution, not tissue (Rule 2, conditional/warn)

Cross-method disagreement tracks dataset gene-count distribution, not tissue or
mitochondrial fraction. C1-vs-MAD3 Jaccard is flat across the mito range
(0.82–0.94); blood/PBMC is indistinguishable from high-mito tissues. Datasets
with many genuinely low-gene cells (small_intestine: C2-vs-MAD5 Jaccard 0.774,
C2 dropping 11,300 cells MAD kept) show the largest gaps. **Rule:**
`scrna_qc_low_gene_dataset_caution` (conditional, warn). Detail:
`cp2/findings.md#sub-finding-3`.

## Finding 4 — MAD parameter sensitivity (descriptive)

MAD k=3 filters 2–3× more cells than k=5 (liver 9,035 vs 4,768; heart 4,243 vs
1,595); k=5 filters ~0 on 3 datasets. Heumos's caveat that k=3 is aggressive is
empirically confirmed. Descriptive only — no rule.

## Synthesis

QC method choice matters less than discourse suggests. Conventional fixed-floor
defaults (C2) work well in most cases and approximate Heumos's MAD
recommendation. The real variance source is dataset gene-count distribution: on
tissues with many low-gene cells, prefer MAD or permissive-quantile over fixed
floors. This parallels Audit 3's arc — a convention (here, QC method) matters
less than assumed, while a dataset property (gene-count distribution, like
Audit 3's ambient burden) is the actual driver of divergence.

## Heumos 2023 positioning — extends

Heumos recommends MAD (Germain 2020) over fixed thresholds as more robust,
without quantifying the cell-set difference. This audit is the empirical
maintenance layer: the recommendation is directionally correct (MAD ≠ pure
quantile), but the practical magnitude of switching from common practice (C2)
to MAD is **smaller than the qualitative framing implies** (C2≈MAD3, 0.969).

## Methodological observations

- §3.5 candidate (tight CIs ≠ correctness; from Audit 3 CP4) applied as a
  sanity check — MAD5≈0 filtering verified as genuine, not degenerate.
- §3.4 companion-metrics requirement satisfied (Jaccard + directional agreement
  + per-metric distributions).
- §5.3.2 equivalence tiers (from Audit 3 CP7) applied cleanly — both rules
  conditional via the documented-feature-dependence clause.
- **Multiple audits now needed `best_practices_references` via workaround**
  (mechanism_notes + xref comment). A formal `references` field is a v1.0.4
  schema candidate (queued via the §3.5/standards pass, not adopted mid-audit).

## Limitations

- n=8 datasets, 8 tissues (n=1 per tissue; tissue used only as mito proxy).
- Census matrices are author-pre-filtered — this audit compares QC *re-filtering*
  of deposited cells, not raw barcode×gene matrices.
- 4 methods, not exhaustive (median+N*MAD hybrids untested).
- The low-gene-distribution mechanism rests primarily on small_intestine (n=1
  low-gene tissue); needs replication.
- No downstream propagation tested (future-audit candidate).

## Audits queued in DEFERRED.md

- QC method downstream propagation (does QC choice cascade to clustering/DE?).
- Additional low-gene-cell tissues to strengthen Rule 2 generalization.
- Hybrid filtering methods (median + N*MAD variants).

## Provenance

All outputs hash-registered in `phase2/repro.lock` (13 audit_qc_mad entries,
verify 13/13, 0 drift — `FINAL_LOCK_VERIFICATION.md`). Native runtimes
(scuttle::isOutlier via Rscript; scanpy + cellxgene_census in base). 2 rules
in `draft_rules/`, both validated `--strict-steps`.
