# Pending AUDIT_STANDARDS.md amendment candidates

Methodological observations surfaced during audit work that are candidates
for the next AUDIT_STANDARDS.md amendment pass. Queued here; **not**
implemented mid-audit (schema changes mid-audit are forbidden).

---

## §5.3.2 equivalence-finding tier criteria — ADOPTED 2026-05-23

**STATUS: ADOPTED** into AUDIT_STANDARDS.md §5.3.2 (2026-05-23). Tier values
unchanged (no schema bump). Surfaced when Audit 3 CP7 found that a robust
"counting tools converge" finding (n=9, overlapping CIs, ρ~0.96) had no tier
under the original §5.3 boundaries, which assume tool-selection audits
("a single tool dominates: wins >60%"). §5.3.2 adds equivalence/agreement
tier criteria (≥8 datasets + overlapping CIs + ρ≥0.90 for hard_default, etc.),
plus rules for nested-effect and contrast-based findings. Validator gained a
warn-only equivalence-tier watchdog. Back-fill review:
`audit3_counting/standards_amendment_review.md` (2 existing rules flagged for
later §5.3.2 review; none re-tiered). Existing §5.3 content retroactively
numbered §5.3.1 (tool-selection criteria).

---

## §3.5 candidate — "Bootstrap CIs reflect sampling variance, not correctness"

**DISPOSITION (Audit 3 CP8 close, 2026-05-23): QUEUED for next batched
standards pass — NOT adopted at Audit 3 close.** Rationale: unlike §5.3.2
(which blocked CP7 tiering and was adopted mid-process out of necessity), §3.5
blocks nothing. It is process guidance (a pre-publication sanity-check norm),
best written alongside the related `clean-control-for-stress-test` candidate
in a single deliberate pass rather than adopted one-off. Both candidates
below travel together to that pass.

**Surfaced from:** Audit 3 CP4 verification, 2026-05-22.

**Content:** Tight bootstrap CIs measure internal consistency of the sampling
distribution. They do not detect systematic errors in tool configuration,
gene ID handling, or other upstream issues that affect all bootstrap
replicates equally. Before treating tight CIs as strong evidence, verify
that the underlying data does not have systematic configuration errors.

**Audit 3 CP4 example:** bootstrap CIs were tight (B=1000, narrow intervals)
on data with an unhandled gene-ID suffix issue (alevin-fry USA-mode -U/-A
suffixes collapsed incorrectly, leaving only the ~10% ambiguous bucket as
alevin-fry's per-gene count). The CIs appeared to strongly support a
"tools disagree" headline (cross-tool rho ~0.58). After fixing the gene-ID
collapse (sum S+A per gene), the headline became "tools converge"
(rho ~0.96) — the CIs were tight in both cases, but only one was correct.

**Recommended amendment:** when an audit produces high-confidence findings
(tight CIs, low p-values), require a "sanity check" step that verifies
the underlying data is correctly configured before treating the findings
as load-bearing. This is a pre-publication check, not a methodological
change to bootstrap procedure.

**Status:** queue for next AUDIT_STANDARDS.md amendment pass; do not
implement schema changes mid-audit.

---

## Registry-extension proposal — two new pipeline_step names — ADOPTED 2026-05-23

**STATUS: ADOPTED.** User approved; `scrnaseq_counting` and
`scrnaseq_cell_calling` added to `pipeline_step_registry.yaml` (schema stays
1.0.3). CP7 unblocked.

**Surfaced from:** Audit 3 CP7 rule drafting, 2026-05-23.

**Need:** The 3 Audit 3 rules (C1/C2/C3) reference two pipeline steps not yet
in `standards/pipeline_step_registry.yaml`. Per the registry's own rule
("New step names must be added here BEFORE being used in rules") and
`--strict-steps` (errors on unknown names), CP7 cannot proceed until these
are registered. NOT added ad-hoc — surfaced here for deliberate approval.

This is a **registry content extension, not a schema-structure change**
(registry stays at schema_version 1.0.3; no v1.0.4 needed). Approving these
two entries unblocks CP7.

**Proposed entries (drop into `steps:` in pipeline_step_registry.yaml):**

```yaml
  - name: scrnaseq_counting
    description: >
      scRNA-seq read counting / quantification. Maps reads to features and
      produces a gene-by-barcode count matrix. Covers STARsolo, alevin-fry
      (salmon), kallisto|bustools (kb-python), CellRanger. Cell-calling is a
      SEPARATE downstream step (scrnaseq_cell_calling).
    modalities: ["scRNA-seq"]
    parameters_used:
      - counting_tool            # STARsolo | alevin-fry | kb-python | CellRanger
      - chemistry                # 10x 3' v2/v3/v3.1, 5' v2, etc.
      - reference                # genome + annotation build
      - multimapper_handling     # discard | EM-rescue (alevin cr-like vs cr-like-em)
    referenced_by_rules:
      - scrna_counting_tool_per_gene_count_convergence
    first_referenced_in: "Audit 3 / 2026-05-23"

  - name: scrnaseq_cell_calling
    description: >
      scRNA-seq cell-vs-empty-droplet calling on the raw count matrix.
      Distinguishes real cells from ambient/empty barcodes. Covers STARsolo
      --soloCellFilter (EmptyDrops_CR), alevin-fry generate-permit-list
      (knee / unfiltered), kb/bustools knee, and uniform downstream callers
      (DropletUtils emptyDropsCellRanger). DISTINCT from scrnaseq_qc_filtering,
      which filters already-called cells by gene/UMI/mito thresholds.
    modalities: ["scRNA-seq"]
    parameters_used:
      - cell_caller              # EmptyDrops_CR | bustools_knee | alevin_knee | unfiltered
      - caller_parameters        # FDR, expected-cells, knee params
      - ambient_burden           # tissue ambient-RNA level (modulates effect severity)
    referenced_by_rules:
      - scrna_cell_calling_permissiveness_chain
      - scrna_cell_calling_biological_propagation_high_ambient
    first_referenced_in: "Audit 3 / 2026-05-23"
```

**Status:** AWAITING APPROVAL. CP7 rule drafting is paused on this. Once
approved, add the two entries to the registry (schema_version stays 1.0.3),
then the 3 rules validate cleanly under `--strict-steps`.

---

## Closeout-amendment candidate — "Clean-control comparison for stress-test findings"

**Surfaced from:** Audit 3 CP6/C3, 2026-05-23.

**Content:** When an audit finding rests on a stress-test case (high-ambient
sample, unusually large dataset, edge condition), include a clean-control
comparison where feasible. The clean-control contrast distinguishes "real
effect on hard cases" from "effect specific to one weird dataset." Without it,
a dramatic single-dataset result is ambiguous between a genuine
condition-dependent effect and an idiosyncrasy of that dataset.

**Audit 3 C3 example:** C3 ran the downstream pipeline on the high-ambient
mouse intestine dataset (where C2 found a 3× cell-count spread across callers)
AND on a clean PBMC control. Intestine showed clustering ARI 0.61 / annotation
agreement 0.34 across tools; the PBMC control showed ARI 0.88 / agreement 0.90+.
The contrast established that the propagation is **ambient-burden-dependent**,
not intestine-specific — a much stronger and more generalizable claim than
"intestine clusters inconsistently." A single-dataset C3 (intestine only)
could not have made that distinction.

**Recommended amendment:** add to the closeout-amendment / audit-design
guidance: for any finding load-bearing on a stress-test condition, document
whether a clean control was run, and if not, why it was infeasible. The
control need not be elaborate — one well-chosen low-stress dataset processed
identically is enough to anchor the contrast.

**Status:** queue for next AUDIT_STANDARDS.md amendment pass; do not
implement schema changes mid-audit.
