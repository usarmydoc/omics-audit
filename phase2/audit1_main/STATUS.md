# Audit 1 main — STATUS

| Checkpoint | Status | Date |
|---|---|---|
| CP1 Input inventory + verification | ✓ complete | 2026-05-16 |
| CP2 Environment setup | ✓ complete | 2026-05-16 |
| CP3 E1 Tool agreement | ✓ complete | 2026-05-16 |
| CP4 E2 Database choice | ✓ complete | 2026-05-16 |
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

## CP3 summary
- 3 tools (fgsea, gseapy.enrichr, clusterProfiler ORA) × 125 inputs × 1 database (Hallmark) = 375 runs in 8.6 min
- 14 empty outputs from low-sig inputs (e.g., 5e57cd50 with 1 sig gene)
- 354 pairwise comparisons (after dropping 7 incomplete inputs)
- 9 strata = 3 input categories × 3 tool pairs; bootstrap B=1000
- Headline finding: ORA tools agree with each other (Spearman ~0.93, direction ~0.86); GSEA-vs-ORA disagrees fundamentally (Spearman ~0.30-0.44, direction ~0.50 = coin flip). Pattern holds across TCGA / Census / GTEx.

## CP3 deliverables
- `e1/runs/<input_id>__<tool>.tsv` (375 per-tool enrichment outputs)
- `e1/e1_run_summary.tsv` (timing + status)
- `e1/per_input_metrics.tsv` (354 input × tool-pair × metric rows)
- `e1/per_stratum_bootstrap.tsv` (9 strata with 95% bootstrap CIs)
- `e1/E1_findings.md`
- E1b queued in DEFERRED.md as a separate future audit
- All hash-registered in `~/omics-audit/phase2/repro.lock`

## CP4 summary
- fgsea held constant; 5 MSigDB databases × 125 inputs = 625 runs in 23.8 min
- 21 empty outputs (same low-sig inputs that produced empty in CP3)
- Within-database redundancy low (0.003-0.016 mean pairwise Jaccard) — not a confound
- Significant pathway count varies 24× across databases (Hallmark 20 → C5 GO:BP 476)
- Leading-edge gene Jaccard between database pairs: 0.08-0.43 across strata
- Pattern uniform across TCGA / Census / GTEx — database choice is a major effect
- Watchpoint: GSEA-only finding; ORA-companion run worth considering before CP7

## CP4 deliverables
- `e2/runs/<input_id>__<database>.tsv` (625 enrichment outputs)
- `e2/e2_run_summary.tsv`
- `e2/per_input_per_db_metrics.tsv` (604 rows)
- `e2/per_input_db_pair_metrics.tsv` (1202 rows)
- `e2/per_database_within_metrics.tsv` (5 rows)
- `e2/per_stratum_bootstrap.tsv` (30 strata)
- `e2/E2_findings.md`
- All hash-registered
