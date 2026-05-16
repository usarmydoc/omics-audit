# Audit 1 main — CP2 report

_Generated: 2026-05-16T18:21:32_

## Environments captured

- **Python:** see `databases/python_session_info.txt`
- **R:** see `databases/R_session_info.txt` (sessionInfo dump)

Both environments registered in `phase2/repro.lock` under kind
`audit1_main_database_cache`. environment_id stamps come from the
lock's `environments` registry (captured at write time).

## Tools installed

### Python (native, no rpy2)
- gseapy 1.1.13 (was already present)
- goatools 1.6.5 (installed CP2)
- mygene 3.2.2 (installed CP2)
- biothings-client 0.5.0 (dependency of mygene)

### R (via Rscript subprocess)
- fgsea 1.36.2
- GSVA 2.4.9
- EGSEA 1.38.0
- clusterProfiler 4.18.4
- limma 3.66.0 (was already present)
- msigdbr 26.1.0
- ReactomePA 1.54.0
- KEGGREST 1.50.0 (was already present)
- rWikiPathways 1.30.0
- org.Hs.eg.db 3.22.0
- org.Mm.eg.db 3.22.0

R version: 4.5.3; Bioconductor: 3.22; BiocManager: 1.30.27.

## Databases cached (14 files in `databases/`)

| Database | Human pathways | Mouse pathways |
|---|---:|---:|
| MSigDB Hallmark | 50 | 50 |
| MSigDB C2 KEGG_LEGACY | 186 | 186 |
| MSigDB C2 Reactome | 1839 | 1817 |
| MSigDB C2 WikiPathways | 925 | 922 |
| MSigDB C5 GO:BP | 7538 | 7402 |
| Reactome.db (full Bioc) | 2820 | 1815 |
| WikiPathways GMT (fresh, May 2026) | 985 | 214 |

## KEGG licensing note

KEGG REST API is free for academic/non-commercial use; rate limit
3 requests/sec. MSigDB C2 KEGG_LEGACY is the cached subset (186
pathways for both human and mouse) suitable for the E2 database-
choice audit. Direct KEGGREST calls deferred to sub-audit use only
if MSigDB coverage proves insufficient.

## Known issues (carried forward from CP1)

### Gene ID format heterogeneity in DEG inputs

Per CP1 (datasets/CP1_report.md), the 125 input DEG TSVs use a mix
of gene ID systems:
- TCGA: 31/48 symbol_human, 17/48 mixed_or_unknown
- Census: 7/47 symbol_human, 14/47 ensembl_human, 5/47 symbol_mouse,
  21/47 mixed_or_unknown
- GTEx: 30/30 ensembl_human (clean)

**Resolution strategy** (per CP1 follow-up guidance, 2026-05-16):
Gene ID mapping is handled at each sub-audit (E1-E4) using whatever
the chosen pathway tool natively expects. The cached MSigDB tables
include three gene ID columns (`gene_symbol`, `gene_ensembl`,
`gene_entrez`) to support cross-mapping per sub-audit need. If a
specific sub-audit hits a wall because of gene ID issues, it will
be surfaced then; no infrastructure pre-built.

### KEGG vs KEGG_LEGACY

MSigDB transitioned from `KEGG` to `KEGG_LEGACY` in 2024 when
KEGG's distribution license tightened. The cached collection is
the LEGACY snapshot, which is the most recent version Bioconductor
is permitted to redistribute. For pathway content this is
indistinguishable from current KEGG for the purposes of inter-tool
comparison (E1) and inter-database comparison (E2).

### Database staleness

MSigDB v2025.1 cached via msigdbr 26.1.0; WikiPathways GMT dated
2026-05-10 (fresh). Reactome.db is the Bioc 3.22 snapshot. These
are pinned via lock-file hashes — re-running CP2 with a newer
msigdbr release would produce different hashes and require a new
CP2 commit if the audit needs the updated content.
