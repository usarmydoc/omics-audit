# Audit Ambient Correction — Final Lock Verification

_2026-05-25, CP4 closeout. Lock: `phase2/repro.lock` (path-keyed SHA256 via dge_native)._

## Result: CLEAN

Full-audit re-hash (`dge_native.verify_outputs(prefix="audit_ambient_correction")`):

| metric | count |
|--------|-------|
| **ok (re-hash matches)** | **216** |
| drift (hash changed) | 0 |
| missing (registered, absent on disk) | 0 |

All registered outputs across CP0–CP4 re-hash to their stored SHA256. No drift, no missing entries.

## Coverage by checkpoint
- **CP1** — per_dataset_metrics, per_stratum_bootstrap, per_gene_contamination, findings, tool_failure_modes + 54 per-dataset (pergene + summary). Registered, kind `ambient_contamination_cp1`.
- **CP2** — per_dataset_ordering_metrics, per_dataset_cross_tool_metrics, cell_disposition, per_stratum_bootstrap, findings, tool_failure_modes + 108 per-(dataset×tool×ordering). Kind `ambient_ordering_cp2`.
- **CP3** — all_per_condition_metrics, cross_condition_comparison, contested_cell_disposition, findings, tool_failure_modes + 30 per-condition (obs/markers/summary). Kinds `ambient_propagation_cp3`, `cp3_*`.
- **CP4** — AUDIT_AMBIENT_CORRECTION_SYNTHESIS.md + 4 draft rule YAMLs. Kind `ambient_cp4_synthesis_rules`.

## Lock integrity (no clobbering of prior audits)
Registration used `dge_native.register_output` (path-keyed append), never the
`repro` CLI `verify` (whose directory-summary schema replaces `verified_outputs`
wholesale and is incompatible with the dge_native flat SHA256 format). Each
checkpoint's `git diff --stat phase2/repro.lock` showed insertions only, zero
deletions — the 1,678 pre-audit entries (Audit 1/3/QC-MAD/BioOrchestrator) and
all earlier-checkpoint entries remained intact throughout.

## Standards provenance
- Registry extension `scrnaseq_ambient_correction` added to
  `standards/pipeline_step_registry.yaml` (no schema change; logged ADOPTED in
  `standards/PENDING_AMENDMENTS.md` 2026-05-25). All 4 rules validate under
  `validate_rules.py --strict-steps` (0 errors).
