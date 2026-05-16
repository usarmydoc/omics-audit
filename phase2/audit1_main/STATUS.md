# Audit 1 main — STATUS

| Checkpoint | Status | Date |
|---|---|---|
| CP1 Input inventory + verification | ✓ complete | 2026-05-16 |
| CP2 Environment setup | pending | — |
| CP3 E1 Tool agreement | pending | — |
| CP4 E2 Database choice | pending | — |
| CP5 E3 Background gene set | pending | — |
| CP6 E4 Multiple testing | pending | — |
| CP7 Findings synthesis + rules | pending | — |

## CP1 summary
- 125 input DEG TSVs verified clean (48 TCGA + 47 Census + 30 GTEx)
- All in repro.lock with verified sha256
- 3 Census pseudobulk inputs below 50-sig sanity threshold (data-genuine,
  pseudobulk method just conservative on those datasets)
- Gene ID format heterogeneity flagged for CP2 normalization plan
- Sample sizes meet `conditional` confidence tier for all 3 categories

## CP1 deliverables
- `datasets/audit1_main_inputs.tsv` — inventory, 125 rows × 14 cols
- `datasets/CP1_report.md` — summary report
- Both registered in `~/omics-audit/phase2/repro.lock`
