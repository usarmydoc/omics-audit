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
