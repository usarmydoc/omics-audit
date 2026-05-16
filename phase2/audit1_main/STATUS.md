# Audit 1 main — STATUS

| Checkpoint | Status | Date |
|---|---|---|
| CP1 Input inventory + verification | ✓ complete | 2026-05-16 |
| CP2 Environment setup | ✓ complete | 2026-05-16 |
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

## CP2 summary
- Python tools: gseapy 1.1.13 (existing), goatools 1.6.5 + mygene 3.2.2 (new)
- R tools: fgsea, GSVA, EGSEA, clusterProfiler, msigdbr 26.1.0, ReactomePA,
  rWikiPathways, org.{Hs,Mm}.eg.db installed via BiocManager (Bioc 3.22).
  limma + KEGGREST already present.
- 14 cached databases (~125 MB total): MSigDB Hallmark/C2-KEGG/C2-Reactome/
  C2-WikiPathways/C5-GO:BP for human + mouse, plus reactome.db full extract
  and fresh WikiPathways GMTs (2026-05-10).
- Both Python + R sessionInfo captured for env_id stamps.
- Gene ID heterogeneity from CP1 carried forward into CP2_report.md
  "Known issues" section — handled per-sub-audit, not pre-built infrastructure.

## CP2 deliverables
- `databases/` — 14 cached database TSVs/GMTs + 2 env-info files
- `CP2_report.md` — environment summary + Known issues
- All registered in `~/omics-audit/phase2/repro.lock` (16 cache files +
  CP2_report.md)
