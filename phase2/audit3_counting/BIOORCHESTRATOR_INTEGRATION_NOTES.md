# Audit 3 → BioOrchestrator integration notes

_For the next batched BioOrchestrator update. Documentation only — NO BO code
changes made in CP8. BioOrchestrator is personal/internal._

## 1. Four Audit 3 rules ready for the next batched update

Staged in `phase2/audit3_counting/draft_rules/`, all validated
(`--strict-steps`, schema v1.0.3), all `prior_audit_relationship: novel`,
all hash-registered:

| rule_id | tier | severity | pipeline_step |
|---|---|---|---|
| `scrna_counting_tool_per_gene_count_convergence` | hard_default | info | scrnaseq_counting |
| `scrna_cell_calling_permissiveness_chain` | hard_default | warn | scrnaseq_cell_calling |
| `scrna_uniform_cell_caller_eliminates_disagreement` | hard_default | warn | scrnaseq_cell_calling |
| `scrna_cell_calling_biological_propagation_high_ambient` | hard_default | warn | scrnaseq_cell_calling |

Promote with a single version bump (per the batched-update discipline — do not
push individually). The actionable rule for users is
`scrna_uniform_cell_caller_eliminates_disagreement` ("pin the caller").

## 2. Two new pipeline_step names needing engine support

Added to `standards/pipeline_step_registry.yaml` (schema 1.0.3):
- **`scrnaseq_counting`** — read counting / quantification (STARsolo,
  alevin-fry, kb-python, CellRanger).
- **`scrnaseq_cell_calling`** — cell-vs-empty-droplet calling (distinct from
  `scrnaseq_qc_filtering`).

The BO rule engine's trigger matching must recognize these step names when the
4 rules are loaded. Confirm the engine reads the registry (or its own copy)
and that `condition_type: pipeline_step` resolves these.

## 3. §5.3.2 equivalence-finding tier criteria

`AUDIT_STANDARDS.md` §5.3 was split into §5.3.1 (tool-selection) + §5.3.2
(equivalence-finding). All 4 Audit 3 rules are tiered under §5.3.2. If the BO
engine has any tier-aware logic (e.g., how strongly a `hard_default` rule is
applied), it should be aware that a §5.3.2 `hard_default` is an
*equivalence* claim (tools interchangeable) rather than a §5.3.1 *selection*
claim (one tool wins) — the recommended action differs (operational-choice
freedom vs pick-this-tool). The rule `severity` field already distinguishes
this (info/warn), so engine logic keyed on severity is fine as-is.

## 4. Back-fill candidates from the §5.3.2 review (not Audit 3 rules)

`standards_amendment_review.md` flagged two *existing deployed* BO rules for
possible §5.3.2 re-tiering when next revised (NOT auto-re-tiered):
- `pathway_tool_paradigm_choice` — partial-equivalence (ORA×ORA equivalent,
  GSEA×ORA not); candidate for split + §5.3.2.
- `scrna_muscat_dream_cell_level_de_inflation_warn` — directional
  non-equivalence.
Handle in the same batched review.

## 5. Future BO validation work — end-to-end validation candidate

Captured per CP8 spec: a **brain n=2 dataset end-to-end validation** as future
BioOrchestrator validation work — run a full pipeline through BO with the
Audit 3 rules active and confirm the engine surfaces them correctly on a
held-out tissue.

_Note: this item is referenced in the CP8 spec as arising from an
"idle-curiosity question"; the original framing is not in my current session
context (likely pre-compaction). Confirm the intended dataset/scope before
acting — recorded here so it is not lost._

## Status

None of the above is acted on now. This file is the hand-off for whenever the
next BioOrchestrator integration session happens.
