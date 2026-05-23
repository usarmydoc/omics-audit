# Phase 2 — corpus-level deferred audit candidates

Candidate future audits surfaced during audit work, NOT committed. Each runs
only if its re-surface trigger fires. Audit-specific scope-creep lives in the
per-audit DEFERRED.md files (audit1_main/, audit3_counting/); this file holds
corpus-level / new-modality candidates.

---

## New-modality audit candidates (from Heumos 2023 indexing, 2026-05-23)

Surfaced while indexing Heumos et al. 2023 (`standards/reference_literature/`).
These cover modalities/steps the corpus does not yet audit. **Re-surface
trigger (all): Audits 2 and 3 complete; clear bandwidth for a new-modality
audit; specific scientific motivation.** None is committed.

- **scRNA-seq trajectory inference tools** — Slingshot vs PAGA vs RaceID across
  topology types, per the Saelens 2019 benchmark framework.
- **RNA velocity inference** — velocyto vs scVelo (dynamical model) on real
  data, with phase-portrait diagnostics. (Note: scVelo dynamical-model
  assumptions have been challenged since 2023; scope against current methods.)
- **Cell–cell communication inference** — CellChat vs CellPhoneDB vs
  SingleCellSignalR, per the Dimitrov 2022 comparison framework (LIANA).
- **Perturbation modeling** — Mixscape vs Augur vs MELD on Perturb-seq data.
  (Note: foundation-model perturbation methods have emerged since 2023.)
- **scATAC-seq dimensionality reduction** — Signac vs cisTopic vs snapATAC,
  per the Chen 2019 benchmark.

These are bounded if scoped to a named benchmark + fixed tool set. Each would
follow the Audit 3 checkpoint pattern (CP0 inventory → … → rule drafting).
Temporal caveat: the fast-moving subdomains (velocity, CCC, perturbation,
spatial) will need a current-methods scan at scope time, since Heumos's 2023
tool lists have partially churned.

Sharpened scope (from Heumos 2023 corpus read, 2026-05-23):
- **Trajectory** — Slingshot vs PAGA vs Monocle3 across topology classes
  (linear / branching / cyclic) per Saelens 2019; metric = topology + ordering
  correlation. **Bounded** single audit. Temporal: moved on since 2023.
- **RNA velocity** — velocyto steady-state vs scVelo dynamical vs a current
  method (e.g. veloVI) on data with orthogonal lineage ground truth;
  phase-portrait + cross-boundary-direction diagnostic. **Bounded.** Temporal:
  scVelo dynamical assumptions challenged post-2023 — scope against current.
- **Cell–cell communication** — CellChat vs CellPhoneDB vs SingleCellSignalR
  via LIANA's consensus-rank framework (Dimitrov 2022); metric = rank
  agreement on known ligand-receptor pairs. **Bounded.** Temporal: LIANA+/
  newer ensembles since 2023.
- **Perturbation modeling** — Mixscape vs Augur vs MELD on Perturb-seq with
  known perturbation effects. **Medium / borderline open-ended** — overlaps a
  fast-moving foundation-model space; scope tightly or it sprawls.
- **scATAC-seq dim reduction** — Signac (LSI) vs cisTopic vs SnapATAC per
  Chen 2019. **Bounded** but requires scATAC hands-on context the operator
  may not have — deprioritize unless motivated.

---

## Partially-covered gaps (extensions of EXISTING audits — from Heumos read)

The corpus touches these areas but Heumos's recommendation has dimensions the
audits didn't address. Each is a natural EXTENSION of an existing audit, not a
new modality program. Re-surface trigger (all): the parent audit is revisited
OR a study needs the specific dimension. None committed.

- **Ambient RNA correction** — Phase 1 covered mito-QC; Heumos also recommends
  ambient correction (SoupX / CellBender / DecontX). Gap: do these correct the
  same contamination, and does the choice change DEGs/clusters? Natural
  extension of the QC audit + relevant to Audit 3's ambient finding.
  **Small-medium. Bounded.** Highest-value of the partials (ties to Audit 3
  C3's ambient story).
- **Normalization-by-purpose** — not a corpus project. Heumos: shifted-log for
  DR, Pearson residuals for HVG, scran for depth. Gap: does the
  normalization-method choice change downstream clustering/DE? **Medium.
  Bounded.** Architectural/durable.
- **Leiden vs Louvain** — Phase 1 p9 swept Leiden resolutions but never
  benchmarked Leiden-vs-Louvain. Gap: does the community-detection algorithm
  (not just resolution) change cluster assignments? **Small.** Natural P9
  extension.
- **Doublet multi-method ensemble** — Phase 1 p2 + A5 benchmarked scDblFinder
  vs Scrublet; Heumos suggests considering multiple methods. Gap: does an
  ensemble (e.g. scDblFinder ∪ DoubletFinder ∪ scds) beat the single best?
  **Small.** Extension of p2/A5.
- **Annotation 3-step workflow + CellTypist** — Phase 1 p5 benchmarked
  MarkerScore vs SingleR accuracy; Audit 3 CP6 used CellTypist without
  benchmarking it. Gap: CellTypist vs SingleR vs MarkerScore agreement, and
  does the automated→manual→verification workflow change calls? **Medium.**
- **decoupleR multi-method enrichment ensemble** — Audit 1 E2 covered database
  choice + paradigm; Heumos recommends decoupleR's run-many-methods-and-
  consensus approach. Gap: does the consensus beat any single enrichment
  method? **Small-medium.** Extension of Audit 1.
- **Compositional / differential-abundance analysis** — not_audited. scCODA /
  MILO / DA-seq for cell-composition shifts between conditions. **Medium.
  Bounded.** A distinct step (not DE, not clustering) the corpus skips entirely.
