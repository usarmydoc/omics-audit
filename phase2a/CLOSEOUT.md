# Phase 2a Closeout

**Status:** Phase 2a substantively complete. Schema v1.0.2 migration
delivered. Audit standards bedded in. Only A6 (NEBULA + muscat full
21 datasets) still running in background.

**Last updated:** 2026-05-16

---

## Deliverables — final state

### Rule corpus (post-categorization + migration)

`/mnt/nvme1/omics-audit/phase2/draft_rules/`:

| Path | Rules | Status | Validator |
|---|---|---|---|
| `bioc_version_sensitivity_v2.yaml` | 1 — `bulk_deg_bioc_version_warn` | active draft | ✓ PASS v1.0.2 |
| `clustering_metric_selection.yaml` | 2 — `clustering_selection_metric_choice`, `clustering_homogeneity_or_completeness_alone_error` | active draft | ✓ PASS v1.0.2 (both) |
| `mito_threshold_quantile.yaml` | 1 — `mito_threshold_compute_from_data` | active draft | ✓ PASS v1.0.2 |
| `sva_preprocessing_sensitivity.yaml` | 1 — `sva_preprocessing_choice_documented` | active draft | ✓ PASS v1.0.2 |
| `superseded/deg_tool_version_sensitivity.yaml` | 1 — `bulk_deg_tool_version_pin` | superseded | ✓ light_superseded |
| `pending_engine_support/bioc_release_transition_info.yaml` | 1 — `bulk_deg_bioc_release_transition_info` | pending engine feature | ✓ light_pending_engine_support |

**5 active rules + 1 superseded rule + 1 pending-engine-support rule.**
4 of the original 11 artifacts were evicted as non-rules (1 data table,
2 findings, 1 meta/standard) per `standards/RULE_CATEGORIZATION.md`.

### Schema and standards

`/mnt/nvme1/omics-audit/standards/`:

| File | Purpose |
|---|---|
| `AUDIT_STANDARDS.md` | v1.0.2 schema spec + standards (provenance, sample size, bootstrap units, confidence tiers, closeout amendments). Symlinked at `~/omics-audit/AUDIT_STANDARDS.md` |
| `SCHEMA_AMENDMENT_PROPOSAL_v1.0.1.md` | Friction amendments (status, created_date, out_of_scope list); folded into v1.0.2 |
| `SCHEMA_REVIEW_v1.0.0_against_phase2_rules.md` | Per-rule schema review; superseded by RULE_CATEGORIZATION.md |
| `RULE_CATEGORIZATION.md` | Reclassification of the 11 drafts into rules vs findings vs data vs meta |
| `CLOSEOUT_v1.0.2_migration.md` | Pre-migration closeout (Steps 1-7 + worked example) |

### Findings + reference docs

| File | Purpose |
|---|---|
| `phase2/STATUS.md` | Phase 2 main status |
| `phase2/PHASE2_FULL_SUMMARY.md` | Comprehensive Phase 2 + 2a writeup |
| `phase2/version_sensitivity/findings.md` | Step 5b deep dive |
| `phase2a/findings.md` | Phase 2a top-level findings (A1-A5, A7) |
| `phase2a/CLOSEOUT.md` | (this file) |
| `phase1_references.md` | Anchor sections for Phase 1 outputs cited by rules: `#p1-mito-threshold`, `#b1-deg-tool-agreement-published`, `#b4-batch-correction-published`, `#p9-clustering-resolution` |
| `phase2/MIGRATION_NOTES.md` | Running migration log; **4 Phase 2b candidates collected** |

### Validator + tooling

| File | Purpose |
|---|---|
| `phase2/scripts/validate_rules.py` | Schema v1.0.0-v1.0.2 validator. Strict pass for active rules; light pass for superseded/ and pending_engine_support/. Companion-metrics watchdog for tool-comparison rules |
| `phase2/draft_rules/VALIDATION_REPORT.txt` | Original v1.0.0 baseline report (all 11 fail); kept for provenance |

---

## Schema v1.0.2 in effect

Versioning history:

- **v1.0.0** — initial; required `last_reviewed`/`reviewer`/`revision_history` on every rule; `out_of_scope` string only; no `status` or `severity` fields
- **v1.0.1** — friction amendments: `status` enum (gates review fields), `created_date` split, `out_of_scope` list form
- **v1.0.2** — added `severity` enum (info/warn/error/reject), distinct from `confidence_tier`

NOT adopted (rejected per categorization pass):
`rule_type` enum, structured `recommendation.parameters`, expanded
`condition_type` (for meta-rules), `literature_note` / `data_table`
rule types. The schema stays narrow; the rule corpus stays focused
on actionable pipeline rules.

---

## Lock state

`/mnt/nvme1/omics-audit/phase2/repro.lock`:

- **411 verified entries, 0 missing, 0 mismatched**
- All Phase 2 + Phase 2a outputs hash-registered
- 5 active rule YAMLs registered with their latest hashes
- 4 Phase 1 outputs retroactively registered when needed for rule
  evidence (P1 mito, B1 DEG agreement, B4 batch correction, P9
  clustering resolution)
- Phase 1's own lock files (`bulk_audit/output/repro.lock` and
  `scrnaseq_audit/output/repro.lock`) remain immutable
  (`verified_outputs = {}`, preserving the documented Phase 1 repro gap)
- Stale entries from deleted/moved files archived under
  `stale_outputs` (2 entries: deleted `tool_concordance_reporting.yaml`
  and moved `deg_tool_version_sensitivity.yaml`)
- `environments` registry contains one env_id
  (`464cdbd07774154b`) covering R 4.5.3 / Bioc 3.22 / edgeR 4.8.2 /
  limma 3.66.0 / pyDESeq2 0.5.4 / scanpy / scipy / pandas

---

## Phase 2b candidates (4 collected, none acted on)

From `phase2/MIGRATION_NOTES.md`:

| # | Candidate | Source |
|---|---|---|
| B2b-1 | Pipeline step name registry | Rules 7, 3, 4 each chose a step name with no central registry; future engine work needs to lock this |
| B2b-2 | `action_type` enum extension for `reject_pipeline` | Rule 4 workaround: `severity: reject` + `action_type: flag_only` cosmetic placeholder |
| B2b-3 | Structured `recommendation.parameters` for `compute_from_data` rules | Currently 1 such rule (mito); revisit when 3+ accumulate |
| B2b-4 | Schema v1.0.3 amendment bundle | Bundle any future amendments into a single pass rather than incremental patches |

To be triaged after A6 completes and Phase 2a fully closes.

---

## A6 status — COMPLETE (2026-05-16)

A6 finished after 235.6 min (3h 56min). Final state: 26/26 datasets
(20 H. sapiens + 6 Mus musculus), all with full 3-tool output
(nebula + muscat_pb_DESeq2 + muscat_mm_dream).

Bug fixes applied during the run, documented in
`phase2/scripts/a6_nebula_muscat.R`:
1. NEBULA reported intercept coefficient instead of group coefficient
2. `group_cell()` returned NULL when cells already donor-sorted
3. `mmDS()` n_threads removed in muscat 1.24 → `BPPARAM=MulticoreParam(6)`
4. `mmDS()` output column names are `p_val/p_adj.loc`, not limma's
5. `pbDS()` requires `pb$group_id` colData column propagated from
   `metadata(pb)$experiment_info` after aggregation
6. `make.unique()` dedup of gene names (mmDS crashes on duplicates)

### A6 verdict vs §1.9a

`prior_audit_relationship: refines_prior`

§1.9a 5-dataset claim: muscat-dream inflated on lymph node + skin
fibroblast specifically (median 0.85× vs pseudobulk).

A6 26-dataset reality:
- Center holds: median inflation 0.78× (close to 0.85×)
- But range is enormous: 0.00× to 57.83×
- 7 of 26 (27%) show >2× inflation; 4 of 26 (15%) show >10×

HIGH inflation set (>50% genes called sig): 5 datasets
- c7775e88 (97%) — blood/COVID-19, naive CD4 T cells
- b617ee1b (89%) — breast/multi-cancer, T cells
- 16023185 (87%) — colon adenocarcinoma, stem cells
- a19d1667 (77%) — skin fibroblast (the §1.9a case ✓)
- a48343a2 (60%) — skeletal muscle

LOW activity set (<10% sig): 13 datasets including
b2dd6bc9 (MOUSE skin fibroblast) at 0.03% — same tissue, opposite
behavior. The inflation is dataset-driven, not tissue-driven.

Companion metrics:
- top-100 Jaccard muscat-dream vs pseudobulk: median 0.176, 21/26 < 0.30
- log2FC Spearman muscat-dream vs pseudobulk: median 0.841 (direction
  and magnitude agree on genes where both tools call a result; the
  disagreement is in WHICH genes pass FDR)

### A6 rule encoded

`phase2/draft_rules/muscat_dream_inflation_warn.yaml`:
- `rule_id: scrna_muscat_dream_cell_level_de_inflation_warn`
- `severity: warn`
- `confidence_tier: conditional`
- `prior_audit_relationship: refines_prior`
- Recommendation: report muscat-dream side-by-side with pseudobulk;
  flag dataset as INFLATED if mm/pb sig ratio > 2 AND top-100
  Jaccard < 0.30 on that specific dataset

Validates clean against schema v1.0.2.

A6 analysis output:
`phase2a/a6_analysis_full.tsv` — 26 rows × 15 cols; sha256
`5988c3c136a815f9...` registered in `phase2/repro.lock` under
kind `phase2a_a6_full_analysis`.

---

## What I did NOT do (per spec)

- No BioOrchestrator integration
- No Audit 2 (spatial) or Audit 3 (counting tools) work
- No new schema amendments (v1.1.0 deferred indefinitely;
  v1.0.3 bundle deferred to post-A6 Phase 2b triage)
- No new validator checks beyond companion-metrics + light pass
- No new pipeline step names beyond what was needed for the 4 migrated rules
- No retroactive change to Phase 1 lock files (kept empty per documented gap)

---

## Open questions awaiting supervised review

Carried from `standards/CLOSEOUT_v1.0.2_migration.md` open questions list, plus new ones surfaced during this migration:

1. **Retroactive Phase 1 output registration policy.** Pattern used: register
   Phase 1 outputs in Phase 2 lock when rules need to cite them; add a
   `Retroactive registration note` paragraph to `phase1_references.md`.
   4 Phase 1 outputs registered this way during migration. Confirm
   pattern.

2. **Pipeline step name registry.** Three step names used so far:
   `scrnaseq_qc_filtering` (mito), `scrnaseq_clustering_resolution_selection`
   (clustering). No central registry. Phase 2b candidate B2b-1
   recommends locking this.

3. **Rule 4 `severity: reject` + `action_type: flag_only` workaround.**
   Schema v1.0.2 lacks `reject_pipeline` action_type. Confirm severity
   drives engine behavior (rather than action_type) when these
   disagree.

4. **Companion-metrics watchdog.** Triggered on Rules 1, 8 (Jaccard
   mentions); both satisfied by adding direction agreement + log2FC
   references to description/recommendation. Working as intended.

5. **`bioc_release_transition_info.yaml` in pending_engine_support/.**
   Rule 2 was held because BioOrchestrator's L3 engine doesn't
   support multi-run history (`condition_type: analysis_context`).
   When the engine adds support, this rule needs final migration to
   strict v1.0.2 and promotion out of pending.

6. **Rule 4 confidence tier `hard_default`.** Used because the
   monotonicity finding is structural (constructional property of
   the metrics) rather than purely statistical. 14/15 dataset
   confirmation is at the lower bound of the hard_default
   "≥15 datasets" threshold. Confirm acceptable.

7. **A6 outcome (pending).** §1.9a re-check verdict will determine
   the next rule's `prior_audit_relationship` tag and whether any
   rule is encodable at all.

---

## Phase 2a — exit criteria

| Criterion | Status |
|---|---|
| Closeout amendments Issues 1-5, 7 addressed | ✓ done in earlier closeout |
| Issue 6 (A6 §1.9a re-check) | ✓ done 2026-05-16 — refines_prior, rule encoded |
| Schema v1.0.2 adopted as operating standard | ✓ done |
| All non-rules evicted from `draft_rules/` | ✓ done (4 evicted) |
| All active rules pass strict v1.0.2 validation | ✓ done (5/5) |
| Superseded rule passes light validation | ✓ done |
| Pending-engine-support rule passes light validation | ✓ done |
| Phase 2b candidates documented (not acted on) | ✓ done (4 in MIGRATION_NOTES.md) |
| All outputs hash-registered | ✓ done (411 entries) |
| Phase 1 lock immutability preserved | ✓ verified |
| BioOrchestrator integration untouched | ✓ verified |

Phase 2a closes when A6 lands and either confirms or contradicts the
§1.9a 5-dataset framing.

---

## What happens next (after A6)

Per the user's standing direction:

1. **A6 result** triggers either `as_original` confirmation, `refines_prior`,
   or `contradicts_prior` per the §1.9a 5-dataset findings.
2. **If contradicts:** stop, surface, await direction.
3. **If confirms or refines:** draft an A6 rule under schema v1.0.2,
   add `phase1_references.md#p4-pseudobulk` anchor, run validator.
4. **Phase 2b candidate triage:** review the 4 collected candidates
   in MIGRATION_NOTES.md; decide which become real work vs which get
   rejected.
5. **BioOrchestrator integration:** only after Phase 2a fully closes
   and supervised review of the 5 active rules.
6. **Resume original Phase 2:** pathway enrichment (Audit 1 main),
   spatial TX (Audit 2), counting tools (Audit 3).
