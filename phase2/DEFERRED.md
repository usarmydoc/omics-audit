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

- **QC filtering method (MAD vs quantile)** — **RESOLVED 2026-05-23 by Audit
  QC-MAD** (`phase2/audit_qc_mad/AUDIT_QC_MAD_SYNTHESIS.md`). MAD vs fixed-floor
  vs pure-quantile produce largely equivalent cell sets (pair Jaccard 0.90–0.97,
  8 Census datasets); C2≈MAD3 (0.969); disagreement driven by gene-count
  distribution. 2 rules drafted. No longer a gap. (Ambient RNA correction,
  below, remains separate and open.)

### Audit QC-MAD follow-up candidates (surfaced at CP4 close, 2026-05-23)

Extensions of the closed QC-MAD audit. Re-surface trigger (all): bandwidth +
specific motivation. None committed.
- **QC method downstream propagation** — **RESOLVED 2026-05-25 (audit
  qc-mad-propagation CP0–CP2, CLOSED)** → `audit_qc_mad_propagation/`. 3 datasets
  (blood/liver/small_intestine, soma_joinid-matched to QC-MAD) × 4 QC methods ×
  downstream pipeline. **Finding: QC method choice is a modest, tissue-INDEPENDENT
  effect (ARI 0.80–0.91 vs C2 everywhere; annotation 88–98%) — the edge-case
  amplification of C3/Ambient-CP3 does NOT replicate; QC method is operational, not
  analytical.** 1 rule. Surfaced the cross-audit pattern (counts-reshaping choices
  propagate; cell-filtering choices don't). Synthesis:
  `audit_qc_mad_propagation/AUDIT_QC_MAD_PROPAGATION_SYNTHESIS.md`. Follow-ups below.

### QC-MAD-propagation follow-ups (surfaced at CP2 close, 2026-05-25)
- **Broader replication** — lift the tissue-independence claim from 3 to all 8
  QC-MAD datasets (flag_and_warn → stronger). **Small.**
- **Extreme high-ambient testing** — does the no-propagation finding break at
  ambient burden beyond small_intestine (e.g. tumor/organoid)? **Medium.**
- **Downstream-parameter sensitivity** — does QC-method propagation depend on
  clustering resolution / HVG choice (fixed at cp6 defaults here)? **Small-medium.**
- **Cross-audit pattern as a rule** — revisit "counts-reshaping propagates,
  cell-filtering doesn't" as a corpus-level rule after a 4th propagation test
  (3 audits = suggestive, not yet load-bearing).
- **Additional low-gene-cell tissue replication** — small_intestine demonstrated
  the gene-count-distribution mechanism (Rule 2); replicating on more low-gene
  tissues would lift Rule 2's generalization from conditional toward stronger.
  **Small, bounded.**
- **Hybrid filtering methods** — median + N*MAD variants and other approaches
  not tested in the 4-method comparison. **Small, bounded.**
- **Ambient RNA correction** — **RESOLVED 2026-05-25 (CP0–CP4 complete, audit CLOSED)** →
  `audit_ambient_correction/`. SoupX / CellBender / DecontX, 9 datasets +
  intestine/PBMC biological propagation. Findings: tools NON-equivalent
  (cross-tool ρ 0.39–0.57; SoupX under-detects high ambient, CellBender tracks
  burden, DecontX most aggressive); only DecontX is ordering-sensitive;
  correct→QC is universally stricter; tool/ordering propagate to biology on
  high-ambient tissue (intestine ARI 0.50–0.70 vs baseline) and wash out on
  clean PBMC (0.85–0.90). 4 rules contributed. Synthesis:
  `audit_ambient_correction/AUDIT_AMBIENT_CORRECTION_SYNTHESIS.md`.
  Follow-up candidates surfaced (below).

### Ambient correction follow-ups (surfaced by Audit Ambient Correction CP4)
- **Ambient correction parameter sensitivity** — this audit used defaults only.
  Gap: do tuned parameters (SoupX manual rho, CellBender --expected-cells /
  FPR, DecontX priors) close or widen the cross-tool gap? **Bounded.**
- **Tool–counting interaction** — ambient correction was run only on STARsolo
  upstream. Gap: does correction behave differently on alevin-fry / kb_count
  outputs (different ambient-gene profiles)? Cross with Audit 3. **Medium.**
- **High-ambient tissue replication** — biological propagation rests on n=1
  high-ambient (intestine). Gap: replicate on additional high-ambient tissues
  to move Finding 4 generalization from flag_and_warn toward hard_default.
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
