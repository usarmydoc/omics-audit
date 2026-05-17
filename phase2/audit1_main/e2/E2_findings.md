# E2 — Database choice effect (fgsea held constant)

**Tool:** fgsea (GSEA-family)
**Databases (5):** MSigDB Hallmark, C2 KEGG_LEGACY, C2 Reactome, C2 WikiPathways, C5 GO:BP
**Bootstrap:** dataset-level, B=1000

**Caveat:** E2 is GSEA-paradigm-specific. ORA tools may show
different database sensitivity. If findings here are strong, an
ORA-companion run could be warranted (see DEFERRED.md, decision
deferred to CP7).

## Within-database structure (redundancy + size)

| Database | n pathways | mean pathway size | mean pairwise Jaccard (sampled) |
|---|---:|---:|---:|
| hallmark | 50 | 146.4 | 0.0138 |
| c2_kegg | 186 | 68.8 | 0.0160 |
| c2_reactome | 1839 | 55.7 | 0.0077 |
| c2_wikipathways | 925 | 44.6 | 0.0076 |
| c5_go_bp | 7538 | 82.4 | 0.0030 |

## Per-database significant-pathway burden (medians across all 125 inputs)

| Database | median n_sig (padj<0.05) | median frac_sig | median |NES| of sig |
|---|---:|---:|---:|
| c2_kegg | 36 | 0.191 | 1.755 |
| c2_reactome | 195 | 0.117 | 1.790 |
| c2_wikipathways | 76 | 0.088 | 1.767 |
| c5_go_bp | 476 | 0.071 | 1.742 |
| hallmark | 20 | 0.390 | 1.680 |

## Cross-database leading-edge gene Jaccard (per stratum, medians + 95% bootstrap CI)

| Category | DB-A | DB-B | n | sig-LE gene Jaccard | direction concordance (shared names) |
|---|---|---|---:|---:|---:|
| census_scrna | c2_kegg | c2_reactome | 42 | 0.229 [0.166, 0.281] | — |
| census_scrna | c2_kegg | c2_wikipathways | 44 | 0.434 [0.324, 0.518] | — |
| census_scrna | c2_kegg | c5_go_bp | 42 | 0.156 [0.118, 0.196] | — |
| census_scrna | c2_reactome | c2_wikipathways | 42 | 0.258 [0.164, 0.331] | — |
| census_scrna | c2_reactome | c5_go_bp | 42 | 0.321 [0.298, 0.368] | — |
| census_scrna | c2_wikipathways | c5_go_bp | 42 | 0.204 [0.163, 0.262] | — |
| census_scrna | hallmark | c2_kegg | 42 | 0.078 [0.053, 0.112] | — |
| census_scrna | hallmark | c2_reactome | 42 | 0.150 [0.101, 0.179] | — |
| census_scrna | hallmark | c2_wikipathways | 42 | 0.096 [0.058, 0.169] | — |
| census_scrna | hallmark | c5_go_bp | 42 | 0.126 [0.099, 0.182] | — |
| gtex_tissue_pair | c2_kegg | c2_reactome | 30 | 0.220 [0.200, 0.252] | — |
| gtex_tissue_pair | c2_kegg | c2_wikipathways | 30 | 0.320 [0.293, 0.350] | — |
| gtex_tissue_pair | c2_kegg | c5_go_bp | 30 | 0.163 [0.148, 0.178] | — |
| gtex_tissue_pair | c2_reactome | c2_wikipathways | 30 | 0.290 [0.257, 0.336] | — |
| gtex_tissue_pair | c2_reactome | c5_go_bp | 30 | 0.370 [0.316, 0.397] | — |
| gtex_tissue_pair | c2_wikipathways | c5_go_bp | 30 | 0.254 [0.229, 0.298] | — |
| gtex_tissue_pair | hallmark | c2_kegg | 30 | 0.196 [0.177, 0.213] | — |
| gtex_tissue_pair | hallmark | c2_reactome | 30 | 0.207 [0.194, 0.224] | — |
| gtex_tissue_pair | hallmark | c2_wikipathways | 30 | 0.217 [0.196, 0.230] | — |
| gtex_tissue_pair | hallmark | c5_go_bp | 30 | 0.184 [0.155, 0.191] | — |
| tcga_cancer | c2_kegg | c2_reactome | 48 | 0.277 [0.267, 0.288] | — |
| tcga_cancer | c2_kegg | c2_wikipathways | 48 | 0.361 [0.351, 0.378] | — |
| tcga_cancer | c2_kegg | c5_go_bp | 48 | 0.202 [0.195, 0.209] | — |
| tcga_cancer | c2_reactome | c2_wikipathways | 48 | 0.335 [0.316, 0.352] | — |
| tcga_cancer | c2_reactome | c5_go_bp | 48 | 0.398 [0.382, 0.413] | — |
| tcga_cancer | c2_wikipathways | c5_go_bp | 48 | 0.279 [0.267, 0.297] | — |
| tcga_cancer | hallmark | c2_kegg | 48 | 0.212 [0.203, 0.221] | — |
| tcga_cancer | hallmark | c2_reactome | 48 | 0.236 [0.227, 0.242] | — |
| tcga_cancer | hallmark | c2_wikipathways | 48 | 0.245 [0.234, 0.253] | — |
| tcga_cancer | hallmark | c5_go_bp | 48 | 0.207 [0.199, 0.221] | — |

## Interpretation

**Leading-edge gene Jaccard** is the most meaningful cross-database
metric here. Different databases test different pathway sets, so
pathway-name overlap is rare. But the *genes that drive*
significant pathways in each database can be compared — high
Jaccard means the databases capture the same biology through
different pathway lenses; low Jaccard means the databases
fundamentally see different signals.

**Direction concordance on shared pathway names** is included for
completeness but is data-thin: only Hallmark and select C2 pathways
have name overlap with other databases, so this metric has
small N at the stratum level.

**Within-database redundancy** (gene-set Jaccard between pathway
pairs in the same database) shows how much each database
double-counts the same biology under different pathway names.
High redundancy means a 'significant pathway count' from that
database is inflated by hierarchical / overlapping pathway
definitions (notably C5 GO:BP).

## Caveats

1. **GSEA paradigm only.** ORA tools may show different database
   sensitivity. ORA-companion deferred to CP7 decision per spec.
2. **Reactome.db full extract not included here.** This audit uses
   the MSigDB C2:CP:Reactome subset for cross-DB consistency
   (msigdbr-curated, gene-mapped uniformly). reactome.db full
   extract is available in CP2 cache for orthogonal validation.
3. **WikiPathways MSigDB subset vs fresh GMT not contrasted here.**
   The fresh May-2026 GMT is cached but not exercised by E2 since
   comparing MSigDB-snapshot vs fresh-GMT is database freshness,
   not database choice. Captured in DEFERRED.md if it becomes a
   question for E1b or a separate sub-audit.

