# Audit 3 (scRNA-seq counting tools) — Status

Canonical drive: `/mnt/nvme1/omics-audit/`. Lock: `phase2/repro.lock`.

## AUDIT STATUS: CLOSED (2026-05-23)

- **Checkpoints:** CP0–CP8 all complete (work window 2026-05-16 → 2026-05-23).
- **Datasets:** 9 (2 chemistries: 3' v2/v3; tissues: PBMC, T-cells, neuron,
  lung, kidney, intestine; human + mouse).
- **Tools:** 4 configs (STARsolo default, STARsolo CR-mimic, alevin-fry, kb-python).
- **Rules contributed:** 4 (all hard_default, §5.3.2 equivalence-tiered, novel).
- **Standards contributed:** §5.3.2 equivalence-finding tier criteria (adopted);
  §3.5 + clean-control candidates (queued).
- **Lock:** 92 entries, 0 drift, 0 missing (`FINAL_LOCK_VERIFICATION.md`).
- **Synthesis:** `AUDIT3_SYNTHESIS.md`. **BO hand-off:** `BIOORCHESTRATOR_INTEGRATION_NOTES.md`.
- **Outstanding decision:** §3.5 amendment adopt-now vs queue (CP8 Step 5) —
  awaiting explicit user input; currently queued.

## Checkpoints

- **CP0–CP2** — env, inventory, references, tool installs. COMPLETE.
- **CP3** — FASTQ acquisition + 4-tool counting (9 datasets × 4 configs). COMPLETE.
- **CP4 / C1** — per-gene count agreement. COMPLETE. Tools converge ρ~0.96.
  (USA-suffix bug found+fixed; gap fixes applied; verified.)
- **CP5 / C2** — cell-calling agreement (A native-3-tool, B common-caller,
  C native-4-tool). COMPLETE. Permissiveness chain STAR ⊂ alevin-fry ⊂ kb;
  uniform caller removes 94% of disagreement.
- **CP6 / C3** — biological propagation. COMPLETE. Ambient-burden-dependent:
  propagates on intestine (ARI 0.61, annotation 0.34), washes out on clean PBMC.
- **Standards amendment §5.3.2** — equivalence-finding tier criteria. COMPLETE
  2026-05-23 (see below).
- **CP7** — rule drafting. COMPLETE 2026-05-23. 4 rules drafted, all pass
  validator `--strict-steps` (schema v1.0.3); see below.
- **CP8** — close-out: synthesis + lock verify + index + DEFERRED triage.
  COMPLETE 2026-05-23. **AUDIT 3 COMPLETE.** See below.

## CP8 — close-out (COMPLETE 2026-05-23) — AUDIT 3 COMPLETE

- **Lock verification:** all 91 `audit3_counting/` entries verify, 0 drift,
  0 missing (inventory TSV re-registered after CP3-completion edits).
- **Audit synthesis:** `AUDIT3_SYNTHESIS.md` — C1+C2+C3 + 4 rules in one
  narrative (counts converge → cells nest → biology diverges on high-ambient).
- **AUDIT_INDEX.md** (repo root) created — Audit 3 listed COMPLETE, 4 rules.
- **DEFERRED.md** final triage: 3c/3d/3e/3f captured with triggers; 3g
  resolved in-audit (CP5 Deliverable C).
- **§3.5 decision:** QUEUED for next batched standards pass (non-blocking;
  travels with the clean-control closeout candidate). Not adopted at close.
- BioOrchestrator integration deferred to next batched BO update (4 Audit 3
  rules + the 2a/2b standards-review candidates).

## CP7 — rule drafting (COMPLETE 2026-05-23)

4 rules in `audit3_counting/draft_rules/`, all PASS `validate_rules.py
--strict-steps`, all hash-registered (verify 4/4, 0 drift), all
`prior_audit_relationship: novel`, tiered under §5.3.2 (equivalence):

| rule | tier | severity | claim |
|---|---|---|---|
| `scrna_counting_tool_per_gene_count_convergence` | hard_default | info | counting tools agree on per-gene counts (ρ~0.96) |
| `scrna_cell_calling_permissiveness_chain` | hard_default | warn | native callers nest STAR⊂alevin⊂kb (9/9); magnitude 3× on intestine = flag_and_warn (§5.3.1, n=1) |
| `scrna_uniform_cell_caller_eliminates_disagreement` | hard_default | warn | uniform EmptyDrops_CR → mean Jaccard 0.99 (94% divergence removed) — the actionable recommendation |
| `scrna_cell_calling_biological_propagation_high_ambient` | hard_default | warn | high-ambient: caller choice changes clustering/markers/annotation (existence hard_default by §5.3.2 contrast; generalization flag_and_warn, n=1) |

Original Rule 2 was split into 2a (permissiveness chain, descriptive) + 2b
(uniform caller, actionable) so the three sub-claims (nesting direction,
uniform-caller convergence, magnitude) carry distinct tiers cleanly.

2 new pipeline_step names (`scrnaseq_counting`, `scrnaseq_cell_calling`) added
to the registry via the §5.3.2 standards amendment. Companion metrics (§3.4):
per-gene log2 ratio + direction (Rule 1); per-cell UMI Spearman + set
differences (Rule 2a); 3 biological readouts (Rule 3) — surfaced in-rule.

Ready for CP8 (final synthesis + BO integration prep) on approval. No BO
changes made; rules staged in draft_rules/ for review.

## Standards amendment — §5.3.2 equivalence-finding tier criteria (2026-05-23)

Added in response to CP7 surfacing a structural gap: the original §5.3 tier
boundaries assume tool-selection audits ("a single tool dominates: wins
>60%"), so a robust equivalence/convergence finding (Audit 3: counting tools
agree, n=9, overlapping CIs, ρ~0.96) had no applicable tier — n=9 falls below
the §5.3.1 hard_default (≥15) and conditional (≥10) dataset floors, and the
dominance criterion is structurally inapplicable to an "all tools agree"
result.

Changes (no schema version bump; tier values unchanged):
- `AUDIT_STANDARDS.md` §5.3 → split into §5.3.1 (tool-selection, existing,
  retroactively numbered) + §5.3.2 (equivalence-finding, new). §5.3.2 sets
  equivalence tiers (hard_default: ≥8 datasets + overlapping CIs + ρ≥0.90;
  conditional: ≥6 + feature-dependence; flag_and_warn; insufficient_data <6),
  plus nested-effect and contrast-based (n=1-per-condition) tiering rules.
- `validate_rules.py` → added warn-only `_check_equivalence_tier_evidence`
  watchdog (flags equivalence-tier rules whose evidence omits CI/correlation/
  stratification language). Does not block. Validator passes: 13 rules,
  0 errors, 1 warn (the watchdog, as designed).
- `PENDING_AMENDMENTS.md` → §5.3.2 marked ADOPTED; registry-extension
  (scrnaseq_counting, scrnaseq_cell_calling) marked ADOPTED.
- Back-fill review: `audit3_counting/standards_amendment_review.md` —
  2 existing rules (pathway_tool_paradigm_choice, scrna_muscat_dream) flagged
  as §5.3.2 candidates for later review; **none re-tiered**.

CP7 now resumes: the 3 Audit 3 rules will be tiered under §5.3.2 (equivalence
criteria) rather than §5.3.1.
