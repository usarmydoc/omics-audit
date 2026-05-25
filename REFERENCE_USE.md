# Using the audit corpus as reference material

The audit corpus produces structured findings encoded as rule YAMLs in each
audit's `draft_rules/` directory. These are **reference documentation** for
scRNA / bulk analyses, not engine input.

## Workflow for using audit findings in an analysis

1. Identify the analysis step (counting, cell-calling, QC filtering,
   normalization, clustering, DE, pathway enrichment).
2. Open the corresponding audit's `findings.md` and `draft_rules/`.
3. Apply the recommendation; deviate only with a documented reason.
4. Cite the audit in methods if findings inform substantive choices.

## Current rule inventory

**11 scRNA audit rules across 4 audits** (Audit 3: 4 + QC-MAD: 2 + Ambient
Correction: 4 + QC-MAD Propagation: 1), plus the Audit 1 pathway rule set.

- **Audit 3 (counting tools):** 4 rules — `phase2/audit3_counting/draft_rules/`
  (per-gene convergence, permissiveness chain, uniform caller, biological
  propagation).
- **Audit QC-MAD (low-quality cell filtering):** 2 rules —
  `phase2/audit_qc_mad/draft_rules/` (method equivalence, low-gene caution).
- **Audit Ambient Correction (SoupX/CellBender/DecontX):** 4 rules —
  `phase2/audit_ambient_correction/draft_rules/` (tool non-equivalence,
  DecontX ordering sensitivity, correct→QC stricter, high-ambient biology
  propagation).
- **Audit QC-MAD Propagation (QC method → biology):** 1 rule —
  `phase2/audit_qc_mad_propagation/draft_rules/` (QC method choice has modest,
  tissue-independent downstream effect — operational, not analytical).
- **Audit 1 main (pathway enrichment):** rules — `phase2/draft_rules/`.
- **Phase 2a (Audit 1 rules + A6):** deployed copies —
  `bioorchestrator/src/bioorchestrator/knowledge/rules/`.

## How to read a rule

Rules use schema v1.0.3. Read `mechanism_notes` and `recommendation` for the
practical takeaway; `severity` (info / warn / error / reject) indicates whether
the rule is informational or actionable; `confidence_tier` (per
AUDIT_STANDARDS §5.3.1 selection / §5.3.2 equivalence) indicates evidence
strength. References to published consensus (e.g. Heumos 2023) appear in
`mechanism_notes` + a top-of-file `# xref:` comment (no dedicated field yet —
v1.0.4 candidate). Evidence entries cite hash-registered output paths.

## BioOrchestrator status

BioOrchestrator — the rule engine that would consume these structurally — is
**tabled**. The rules function as reference material for manual application.
The `BIOORCHESTRATOR_INTEGRATION_NOTES.md` files in audit dirs describe what a
future integration would entail, but no integration is planned at this time.

## See also

- `AUDIT_INDEX.md` — registry of audits, statuses, and findings.
- `standards/AUDIT_STANDARDS.md` — provenance, sample-size, tier, and rule
  schema standards.
- `standards/reference_literature/` — published best-practices sources
  (Heumos 2023) mapped to audit coverage.
