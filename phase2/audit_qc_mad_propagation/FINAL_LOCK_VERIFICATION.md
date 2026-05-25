# Audit QC-MAD Propagation — Final Lock Verification

_2026-05-25, CP2 closeout. Lock: `phase2/repro.lock` (path-keyed SHA256 via dge_native)._

## Result: CLEAN
Full-audit re-hash (`dge_native.verify_outputs(prefix="audit_qc_mad_propagation")`):

| metric | count |
|--------|-------|
| **ok (re-hash matches)** | **44** |
| drift | 0 |
| missing | 0 |

All registered outputs across CP0–CP2 re-hash to their stored SHA256.

## Coverage
- **CP1** — per-condition obs/markers/summary (36), per_comparison_metrics, per_stratum_bootstrap, method_specific_cell_disposition, findings, tool_failure_modes. Kinds `qcprop_cp1`, `qcprop_cp1_compare`.
- **CP2** — AUDIT_QC_MAD_PROPAGATION_SYNTHESIS.md + the 1 draft rule. Kind `qcprop_cp2_synthesis_rule`.

## Integrity
Registration used `dge_native.register_output` (path-keyed append), never the `repro` CLI `verify` (which clobbers). Per-checkpoint `git diff --stat phase2/repro.lock` was additions-only. The parallel pipeline batch's registration was done SERIALLY (post-batch) to avoid concurrent lock-write races.

## Standards
No AUDIT_STANDARDS.md or registry change (used existing §5.3.2 + the `scrnaseq_qc_filtering` step registered by Audit QC-MAD). The 1 rule validates under `validate_rules.py --strict-steps` (0 errors).
