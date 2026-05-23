# §5.3.2 back-fill review — existing rules re-tiering candidates

_Generated: 2026-05-23, in support of the §5.3.2 equivalence-finding amendment._
_READ-ONLY review. No rule was re-tiered. Re-tiering, if any, is a separate decision._

## Method

Ran the new equivalence-tier watchdog across all rule YAMLs in
`phase2/draft_rules/` (10 files, 13 rules) and
`bioorchestrator/src/bioorchestrator/knowledge/rules/` (BO deployed copies +
v0 legacy). 27 rule instances scanned; the BO copies duplicate the phase2
drafts (pathway_*, clustering, mito, muscat, sva, bioc), so ~13 unique rules.

A rule was flagged if its evidence/title contained equivalence vocabulary
(converge, agreement, equivalent, nest, Jaccard, correlation, overlap,
Spearman, substitutable, interchangeable) AND it sits at an evidence-bearing
tier. Flag ≠ re-tier recommendation — most flags are disagreement/sensitivity
findings that §5.3.1 already handles.

## Classification

### A. Genuine equivalence-finding candidates (worth §5.3.2 review)

1. **`pathway_tool_paradigm_choice_warn`** (currently `hard_default`, §5.3.1)
   — **partial-equivalence finding.** Audit 1 E1 found ORA×ORA tools roughly
   equivalent but GSEA×ORA paradigms not. This is exactly the
   "equivalence holds on a subset of stratifications" case → §5.3.2
   `flag_and_warn (equivalence)` for the cross-paradigm part, with the
   within-ORA equivalence potentially `conditional`/`hard_default (equivalence)`.
   The current single `hard_default` tier conflates the two. **Candidate for
   split + §5.3.2 re-tier.** Flagged in the CP7 spec.

2. **`scrna_muscat_dream_cell_level_de_inflation_warn`** (currently
   `conditional`, §5.3.1) — **non-equivalence finding** (dream cell-level DE
   inflates relative to pseudobulk). This is the inverse of equivalence: it
   documents that two approaches do *not* agree. §5.3.2's `flag_and_warn
   (equivalence)` "equivalence is directional only / holds only on a subset"
   language may fit better than §5.3.1 `conditional`, since the finding is
   about a systematic discrepancy, not a feature-dependent winner. **Candidate
   for review.** Flagged in the CP7 spec.

### B. Watchdog WARN fired (review for evidence-language completeness)

3. **`clustering_selection_metric_choice`** (`conditional`) — watchdog WARN:
   evidence mentions "agreement" but not bootstrap CIs / correlation /
   stratification. Likely a **false positive** — this is a selection finding
   (which clustering metric to use); "agreement" refers to ARI-type metrics,
   not tool equivalence. No re-tier; the WARN just notes the evidence summary
   could be clearer. Stays §5.3.1.

### C. Disagreement / sensitivity findings — stay §5.3.1 (no re-tier)

These contain equivalence vocabulary because they *measure* (dis)agreement,
but the finding is "this choice matters / sources differ," which §5.3.1
`flag_and_warn` / `hard_default` already covers correctly:

- `pathway_database_choice_warn` (hard_default) — databases share few top genes; a *disagreement* finding.
- `pathway_mt_correction_method_info` (hard_default) — MT-correction method comparison.
- `sva_preprocessing_choice_documented` (flag_and_warn) — preprocessing changes results.
- `bulk_deg_bioc_version_warn` (flag_and_warn) — Bioconductor version sensitivity.

None of these assert tool *equivalence*; §5.3.2 does not apply.

### D. Pure selection / reject findings — §5.3.1 (no equivalence signal)

`clustering_homogeneity_or_completeness_alone_error`,
`pathway_no_uncorrected_pvalues_reject`, `pathway_ora_background_default_info`,
`mito_threshold_compute_from_data` — selection/threshold/reject rules, no
equivalence character.

## Recommendation

- **Two rules merit deliberate §5.3.2 review**: `pathway_tool_paradigm_choice_warn`
  (split partial-equivalence; the spec anticipated this) and
  `scrna_muscat_dream_cell_level_de_inflation_warn` (directional non-equivalence).
- **No automatic re-tiering.** Both are deployed BO rules; any change is a
  separate batched BioOrchestrator update decision, not part of this amendment.
- The watchdog is warn-only and did not block any rule.
- Audit 1's tier assignments are **not** back-filled here (out of scope unless
  a conflict surfaces; none blocks).

## Disposition

Queued for the next BioOrchestrator batched update review (not acted on now).
Recorded so the §5.3.2 criteria are applied consistently when those two rules
are next revised.
