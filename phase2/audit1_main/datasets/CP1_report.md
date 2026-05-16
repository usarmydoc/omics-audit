# Audit 1 main — Checkpoint 1 report

_Generated: 2026-05-16T17:16:51_

## Counts by input category

| Category | Total files | Verified | Hash mismatch | Not in lock | n with ≥50 sig | Median sig genes |
|---|---:|---:|---:|---:|---:|---:|
| tcga_cancer | 48 | 48 | 0 | 0 | 48 | 13754 |
| census_scrna | 47 | 47 | 0 | 0 | 44 | 4585 |
| gtex_tissue_pair | 30 | 30 | 0 | 0 | 30 | 18996 |

## Sample size adequacy per sub-audit (per AUDIT_STANDARDS.md §5.3)

Conditional confidence tier requires ≥10 per category. Hard default requires ≥15 across ≥3 tissues/categories.

| Category | n inputs | n comparisons (unique) | n tools | OK for conditional (≥10)? |
|---|---:|---:|---:|:---:|
| tcga_cancer | 48 | 16 | 3 | ✓ |
| census_scrna | 47 | 26 | 2 | ✓ |
| gtex_tissue_pair | 30 | 10 | 3 | ✓ |

## Verification: PASS — all files parse cleanly, no hash mismatches.

## Gene ID format coverage

| Category | symbol_human | ensembl_human | ensembl_mouse | symbol_mouse | other |
|---|---:|---:|---:|---:|---:|
| tcga_cancer | 31 | 0 | 0 | 0 | 17 |
| census_scrna | 7 | 14 | 0 | 5 | 21 |
| gtex_tissue_pair | 0 | 30 | 0 | 0 | 0 |
