# Heumos et al. 2023 — recommendation → audit-coverage index (RESOLVED)

**Source:** Heumos, Schaar, Lance, et al. (Single-cell Best Practices
Consortium). "Best practices for single-cell analysis across modalities."
*Nature Reviews Genetics* 24, 550–572 (2023). doi:10.1038/s41576-023-00586-w.
PDF: `heumos_2023_best_practices.pdf`.

**Purpose:** position existing audit findings against this consensus. NOT an
audit target list. **Resolved 2026-05-23** against the corpus (Phase 1
scrnaseq_audit p1–p5, Phase 2a A2/A5/A6/A7, Audit 1 E2, Audit 3 C1–C3).
**Temporal context:** paper ~3 years old; architectural recommendations hold,
fast-moving subdomains (velocity, CCC, perturbation, spatial) have moved on.

## scRNA-seq main coverage

| Heumos section | Recommendation | Audit status | Relationship | Evidence | Temporal |
|---|---|---|---|---|---|
| QC — low-quality cell filtering | MAD-based outlier filtering on QC covariates; permissive initial filter | covered_by_audit_qc_mad | extends | Audit QC-MAD: QC method choice (MAD vs fixed-floor vs pure-quantile) produces largely equivalent cell sets (median pair Jaccard 0.90–0.97 across 8 Census datasets). Disagreement driven by dataset gene-count distribution, not tissue type. Heumos's MAD recommendation is empirically distinct from pure quantiles but very close to typical fixed-floor defaults (C2 vs MAD3 Jaccard 0.969). Phase 1 p1 (mito quantile-on-data) is the precursor. **Biological propagation also tested (Audit QC-MAD Propagation, 2026-05-25): QC method choice has a modest, tissue-INDEPENDENT downstream effect (ARI 0.80–0.91 vs C2 across blood/liver/small_intestine; annotation 88–98%) — methodological hygiene, not a substantive analytical decision, and notably does NOT show the edge-case amplification that cell-calling/ambient correction do.** | durable |
| QC — ambient RNA correction | Correct ambient RNA (SoupX, CellBender, DecontX) | covered_by_audit_ambient_correction | extends | Audit Ambient Correction (9 datasets + intestine/PBMC propagation): the 3 tools are NON-equivalent on per-gene contamination (cross-tool ρ 0.39–0.57; SoupX under-detects high ambient, CellBender tracks burden, DecontX most aggressive); only DecontX is ordering-sensitive; correct→QC is universally stricter; tool/ordering choices propagate to biology on high-ambient tissue (intestine ARI 0.50–0.70 vs baseline) and wash out on clean PBMC (0.85–0.90). 4 rules. Heumos names SoupX/CellBender without comparison, omits DecontX, doesn't address ordering | durable |
| QC — doublet detection | scDblFinder among top performers; consider multiple methods | confirms (+ extends) | confirms | Phase 1 p2_doublet_audit (scDblFinder vs Scrublet, 19 datasets) + Phase 2a A5 (Demuxlet ground truth: scDblFinder AUROC 0.834 > Scrublet 0.741). **Gap:** 2-tool, not the full multi-method ensemble | durable |
| Normalization & variance stabilization | Task-specific: shifted-log for DR, Pearson residuals for HVG, scran pooling for depth | not_audited | — | normalization-method choice is not a Phase 1/2a/3 project (used as method in Audit 3 CP6, not benchmarked) | durable → DEFERRED |
| Feature selection | Deviance / Pearson-residual HVG over dispersion | not_audited | — | scry deviance used as Audit 3 CP6 method, not benchmarked | durable → DEFERRED |
| Removing confounding variation (batch) | Integrate batches; scVI/Harmony/Scanorama; evaluate with ARI/LISI/kBET | partially_covers | confirms | Phase 1 p3_batch_correction (bbknn / harmony / scanorama / scvi via ARI/LISI/kBET). **Gap:** cell-cycle regression not covered; integration only on Census scRNA | durable |
| Clustering | Leiden over Louvain; sweep resolutions; select by metric | extends | extends | Phase 1 p9 (Leiden resolution sweep, 15 datasets) + Phase 2a A2 (ARI vs V-measure 5× divergence), rules `clustering_selection_metric_choice` + `..._homogeneity_or_completeness_alone_error`. **Gap:** Leiden-vs-Louvain itself not benchmarked | durable |
| Cell-type annotation | 3-step: automated + manual marker + verification | partially_covers | confirms | Phase 1 p5_annotation_tools (MarkerScore vs SingleR precision/recall/F1). **Gap:** CellTypist (used in CP6, not benchmarked); the 3-step workflow not audited as a workflow | durable |
| DGE | Pseudobulk preferred over cell-level Wilcoxon (Squair 2021) | extends | extends | Phase 1 p4_pseudobulk + Phase 2a A6 (muscat-dream inflation) + pseudobulk_vs_wilcoxon (26 datasets, quantified 0.78×/0.85× inflation), rule `scrna_muscat_dream_cell_level_de_inflation_warn` | durable |
| Gene set / pathway enrichment | Gene-set (database) choice matters more than method (Holland 2020); decoupleR ensemble | extends | extends | Audit 1 E2 (24–124× variation across MSigDB collections), rule `pathway_database_choice_warn`. **Gap:** decoupleR's multi-method ensemble not benchmarked | durable |
| Read counting / quantification | Lafzi 2018 review; no specific cross-tool agreement recommendation | covered_by_audit3 | fills_gap | Audit 3 C1 (per-gene Spearman ~0.96, 9 datasets), rule `scrna_counting_tool_per_gene_count_convergence` | durable |
| Cell calling (within QC) | Doesn't differentiate knee/EmptyDrops empirically | covered_by_audit3 | extends | Audit 3 C2/C3 (permissiveness chain + ambient-burden biological propagation), rules `scrna_cell_calling_permissiveness_chain`, `..._uniform_cell_caller_...`, `..._biological_propagation_high_ambient` | durable |
| Compositional / differential abundance | scCODA / MILO / DA-seq for composition shifts | not_audited | — | none | moving → DEFERRED |
| Trajectory inference | Slingshot (simple) / PAGA (complex), Saelens 2019 | not_audited | — | none | moved on since 2023 → DEFERRED |
| RNA velocity | velocyto / scVelo dynamical; phase-portrait checks | not_audited | — | none | moved on since 2023 → DEFERRED |
| Cell–cell communication | LIANA framework ranks across methods (Dimitrov 2022) | not_audited | — | none | tool churn → DEFERRED |
| Perturbation modeling | Mixscape / Augur / MELD | not_audited | — | none | moved on since 2023 → DEFERRED |

## Modality-specific sections (scanned, no corpus coverage)

| Section | Status | Temporal |
|---|---|---|
| Chromatin accessibility (scATAC-seq) | not_audited → DEFERRED | moved on since 2023 |
| Surface protein (CITE-seq) | not_audited | moved on since 2023 |
| Adaptive immune receptors (AIRR) | not_audited | moved on since 2023 |
| Spatial transcriptomics | not_audited | moved on since 2023 (fastest-moving) |

## Resolution note

All 4 prior `needs_review` entries are resolved: **QC** → partially_covers
(quantile-on-data, not MAD; ambient gap); **doublet** → confirms+extends
(p2 + A5 ground truth); **clustering** → extends (p9 + A2; Leiden-vs-Louvain
gap); **normalization** → not_audited (not a corpus project). The prior index
also under-counted coverage: **annotation** is partially_covers (p5 exists),
not not_audited; **batch integration** (p3) was not previously listed. Gaps
captured in `phase2/DEFERRED.md`; high-value candidates triaged in
`heumos_2023_corpus_coverage.md`.
