# Audit 3 (scRNA-seq counting tools) — Status

Canonical drive: `/mnt/nvme1/omics-audit/`. Lock: `phase2/repro.lock`.

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
- **CP7** — rule drafting. PAUSED → resuming against amended §5.3.2.
- **CP8** — final synthesis + BO integration prep. NOT STARTED.

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
