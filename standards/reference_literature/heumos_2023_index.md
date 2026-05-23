# Heumos et al. 2023 — recommendation → audit-coverage index

**Source:** Heumos, Schaar, Lance, et al. (Single-cell Best Practices
Consortium). "Best practices for single-cell analysis across modalities."
*Nature Reviews Genetics* 24, 550–572 (2023). doi:10.1038/s41576-023-00586-w.
PDF: `heumos_2023_best_practices.pdf`.

**Purpose:** position existing audit findings against this published
best-practices consensus. NOT an audit target list. **Temporal context:**
indexed 2026-05-23, paper is ~3 years old (synthesizes 2021–2022 benchmarks);
architectural recommendations hold, fast-moving subdomains have moved on.

## scRNA-seq main coverage

| Heumos section | Specific recommendation | Tool(s) named | Audit status | Relationship | Temporal note |
|---|---|---|---|---|---|
| QC — filtering low-quality cells | Permissive initial filtering; MAD-based outlier detection on QC covariates (counts, genes, mito%) rather than fixed cutoffs | scanpy, scater | covered_by_phase1_§1.1-1.3 | needs_review | durable (MAD/permissive-filtering is architectural) — needs_review: re-read Phase 1 §1.1-1.3 to confirm whether it used MAD or fixed/quantile thresholds |
| QC — ambient / empty droplets | Distinguish cells from empty droplets; correct ambient RNA | EmptyDrops, CellBender, SoupX | covered_by_audit3 | extends | durable — Audit 3 found permissiveness chain + ambient-burden-dependent biological propagation, beyond Heumos's qualitative guidance |
| QC — doublet detection | scDblFinder among top performers | scDblFinder, Scrublet, DoubletFinder | covered_by_phase1_§1.9c | needs_review | durable — needs_review: confirm Phase 1 §1.9c benchmarked scDblFinder as winner vs Heumos's claim |
| Normalization | Shifted-log for downstream DR; analytic Pearson residuals for HVG/feature selection; scran pooling where depth varies | scanpy, scran | partially_covered_by_phase1 | needs_review | durable — needs_review: Phase 1 normalization coverage incomplete vs Heumos's per-purpose split |
| Feature selection | Deviance / Pearson-residual HVG over dispersion-based | scry, scanpy | not_audited (used in Audit 3 CP6 pipeline, not itself audited) | confirms | durable — Audit 3 CP6 followed this (scry deviance) as method, did not benchmark it |
| Dimensionality reduction / clustering | Leiden over Louvain; evaluate at multiple resolutions | scanpy, leidenalg | covered_by_phase1_§1.5 | needs_review | durable — needs_review: confirm Phase 1 §1.5 Leiden-vs-Louvain + multi-resolution stance |
| Differential gene expression | Pseudobulk preferred over naive cell-level (Wilcoxon) tests to avoid pseudoreplication (cites Squair 2021) | muscat, DESeq2, edgeR, limma; (Wilcoxon as the cautioned baseline) | covered_by_phase1_§1.9a + phase2a_A6 | extends | durable — audit quantified the ~0.85x / ~0.78x cell-level inflation across 26 datasets that Heumos cites only generically via Squair 2021 |
| Gene set / pathway enrichment | Gene set (database) choice matters more than the statistical method (cites Holland 2020) | decoupleR, GSEA, fgsea, MSigDB | covered_by_audit1_E2 | extends | durable — Audit 1 quantified 24–124x variation across MSigDB collections that Holland's qualitative claim left unspecified |
| Read counting / quantification | Refers to Lafzi 2018 review; no specific cross-tool agreement recommendation | CellRanger, STARsolo, alevin, kallisto (review-level) | covered_by_audit3 | fills_gap | durable — Lafzi reviewed without benchmarking cross-tool agreement; Audit 3 produced the empirical per-gene + cell-calling agreement benchmark |
| From clusters to cell identities (annotation) | 3-step approach: automated annotation + manual marker-based + verification | CellTypist, scanpy, marker databases | not_audited (CellTypist used as Audit 3 CP6 method) | unknown | durable architecturally |
| Trajectory inference | Slingshot for simple topologies, PAGA / more complex methods otherwise (Saelens 2019 benchmark) | Slingshot, PAGA, RaceID | not_audited | unknown | moved on since 2023 (newer methods + topology-aware benchmarks) |
| RNA velocity | velocyto / scVelo; check phase portraits; dynamical model caveats | velocyto, scVelo | not_audited | unknown | moved on since 2023 (scVelo dynamical-model assumptions challenged; veloVI and successors) |
| Cell–cell communication | Use a framework that ranks across methods (Dimitrov 2022) | LIANA, CellPhoneDB, CellChat | not_audited | unknown | tool churn (LIANA+ and newer ensembles since 2023) |
| Perturbation modeling | Mixscape / Augur / MELD for perturbation effect analysis | Mixscape, Augur, MELD | not_audited | unknown | moved on since 2023 (foundation-model perturbation methods, e.g. GEARS/scGPT-class) |

## Modality-specific sections (one-line, not parsed in depth)

| Heumos section | Audit status | Relationship | Temporal note |
|---|---|---|---|
| Chromatin accessibility (scATAC-seq) | not_audited | unknown | moved on since 2023 |
| Surface protein (CITE-seq) | not_audited | unknown | moved on since 2023 |
| Adaptive immune receptor repertoires (AIRR / TCR-BCR) | not_audited | unknown | moved on since 2023 |
| Single-cell data resolved in space (spatial transcriptomics) | not_audited | unknown | moved on since 2023 (fastest-moving subdomain) |

## needs_review entries (4) — for a future Phase-1-re-read pass, not this one

1. **QC thresholds** — does Phase 1 §1.1-1.3 use MAD-based vs fixed/quantile? (mito_threshold rule used quantile-on-data.)
2. **Doublet detection** — does Phase 1 §1.9c name scDblFinder as the benchmarked winner?
3. **Clustering** — does Phase 1 §1.5 state Leiden>Louvain + multi-resolution?
4. **Normalization** — Phase 1 coverage of shifted-log vs Pearson-residual-by-purpose.

(Under the 5–7 cap; not resolved in this pass — each needs a Phase 1 findings re-read against Heumos's specific wording.)
