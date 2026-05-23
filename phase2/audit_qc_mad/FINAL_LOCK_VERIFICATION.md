# Audit QC-MAD — Final lock state verification (CP4)

_Verified: 2026-05-23. Lock: `phase2/repro.lock`. Verifier: `dge_native.verify_outputs` (SHA256 re-hash)._

## Result: CLEAN

```
audit_qc_mad entries: 13
OK = 13    DRIFT = 0    MISSING = 0
```

Every registered Audit QC-MAD output re-hashes to its recorded SHA256.

## Entries by area

| area | entries | covers |
|---|---|---|
| `cp1/` | 9 | per-dataset metrics, jaccards, bootstrap, disagreement concentration, tool_failure_modes, working_set, 3 scripts |
| `cp2/` | 1 | findings.md |
| `draft_rules/` | 2 | the 2 CP3 rule YAMLs |
| root | 1 | STATUS.md (synthesis + final-lock registered at CP4 commit) |

## Coverage checks (all present)

- [OK] CP1 per_dataset_metrics + jaccards + per_stratum_bootstrap
- [OK] CP2 findings.md
- [OK] 2 draft rule YAMLs (validated --strict-steps)

## Notes

- Per-cell QC metric TSVs (`cp1/qc_metrics/`) and MAD flag files
  (`cp1/mad_flags/`) are intermediate, regenerable from `cp1/pull_qc.py` +
  `cp1/mad_filter.R`, and intentionally not committed/registered (kept local).
- Census version pinned (2025-11-08 stable) for reproducibility; working set is
  an exact subset of Phase 1 p1's dataset_ids.
- Native runtimes: scuttle::isOutlier via Rscript subprocess; scanpy +
  cellxgene_census in base. No rpy2.

## Conclusion

Audit QC-MAD is reproducibly closed: all metrics, findings, and rules are
hash-registered with zero drift.
