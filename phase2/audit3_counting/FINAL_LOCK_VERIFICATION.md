# Audit 3 — Final lock state verification (CP8)

_Verified: 2026-05-23. Lock: `phase2/repro.lock`. Verifier: `dge_native.verify_outputs` (SHA256 re-hash)._

## Result: CLEAN

```
TOTAL audit3_counting entries: 92
OK = 92    DRIFT = 0    MISSING = 0
```

Every registered Audit 3 output re-hashes to its recorded SHA256. Git tree
for `audit3_counting/` matches the lock (working files committed at close).

## Entries by area

| area | entries | covers |
|---|---|---|
| `c1/` | 12 | per-gene count agreement (metrics, bootstrap, findings, verification, superseded buggy set) |
| `c2/` | 12 | cell-calling agreement (A native-3, B common-caller, C native-4, findings) |
| `c3/` | 38 | biological propagation (per-tool pipeline outputs intestine + control, comparison TSVs, findings) |
| `draft_rules/` | 4 | the 4 CP7 rule YAMLs |
| `inventory/` | 8 | dataset inventory + read-structure verification |
| `scripts/` | 7 | pipeline + analysis scripts |
| `environment/` | 5 | tool versions / configs |
| root | 6 | STATUS, DEFERRED, synthesis, reports |

## Coverage checks (all present)

- [OK] C1 findings + per-dataset metrics + bootstrap
- [OK] C2 findings + A/B/C per-dataset metrics
- [OK] C3 findings + intestine cluster-agreement metrics
- [OK] 4 draft rule YAMLs
- [OK] CP5 regenerated alevin-fry (knee) outputs
- [OK] AUDIT3_SYNTHESIS.md

## Provenance notes

- Superseded buggy CP4 outputs retained at
  `c1/superseded_2026-05-18_buggy_usa_strip/` per AUDIT_STANDARDS §1.5
  (kept for reproducibility of the USA-mode-bug correction; not in the
  active verify set count for findings).
- Native runtimes throughout (scry deviance / scDblFinder / EmptyDrops_CR via
  Rscript subprocess; STARsolo / alevin-fry / kb-python; scanpy / CellTypist).
  No rpy2.
- `repro` CLI schema-validates the lock (repro_schema_version backfilled);
  the substantive hash check is `dge_native.verify_outputs` (SHA256), since
  the CLI's MD5/`files`-subkey format differs from the dge_native flat format.

## Conclusion

Audit 3 is reproducibly closed: all analytical outputs, findings, rules, and
the synthesis are hash-registered with zero drift.
