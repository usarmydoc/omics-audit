# Bioinformatics Audits

An empirical evidence base for default parameters and tool choices in bioinformatics pipelines. Anyone can use these audits, cite them, or contribute new ones in the same format. Every audit is fully reproducible via repro-lock.

## What this is

Bioinformatics pipelines rely on analytical defaults that are rarely empirically justified. Published best-practices papers (Heumos et al. 2023 for single-cell; equivalent papers for bulk RNA-seq, pathway enrichment, ATAC-seq, spatial) make recommendations based on benchmarks that are often run on a handful of datasets, sometimes years ago, sometimes by the tool authors themselves. The recommendations get treated as settled, propagate through tutorials and AI-generated pipelines, and the field rests on assumptions that were never tested at scale.

This repository is the ongoing empirical check on those recommendations. Each audit picks one default or tool choice, tests it across diverse public data with pinned versions, and reports whether the recommendation holds, where it holds conditionally, or where it needs to be flagged as genuinely uncertain. The goal is calibration, not contrarianism. Where the field is right, the audits confirm with broader evidence. Where the field is wrong, the audits surface specific failure modes.

The corpus is structured to grow. Older audits get refined or contradicted by newer ones, and the relationships are tracked explicitly. The intent is infrastructure the field can use at the moment someone is deciding what parameters to set, not a frozen snapshot.

## Current contents

**19 audits closed** across scRNA-seq and bulk RNA-seq. Phase 1 produced 13 audits as the initial corpus. Phase 2 and Phase 2a added 6 audits and refined or contradicted three Phase 1 claims.

The `relationship` column shows how each audit stands relative to earlier work in the corpus: **as_original** (stands as published), **refined** (later audit qualified the finding), **contradicted** (later audit overturned the specific claim), **extended** (later audit added new dimensions).

### Phase 1: scRNA-seq (P1–P9)

| Audit | Convention | Finding | Status | Relationship |
|---|---|---|---|---|
| P1 | Mito threshold (20%) | Tissue-dependent; removes 0% of cells in 55% of datasets, >10% in gut epithelia | Challenged | as_original |
| P2 | Doublet detection | scDblFinder vs Scrublet agree at chance (κ = 0.077); 20x rate difference | Challenged | as_original |
| P3 | Batch correction | scVI ≈ Harmony > BBKNN >> Scanorama (ARI 0.556 / 0.522 / 0.448 / 0.299) | Challenged | as_original |
| P4 | Cell-level vs pseudobulk DE | Wilcoxon inflates 3.8x; 9% overlap with pseudobulk | Challenged | **refined** (see Phase 2a) |
| P5 | Cell type annotation | Marker scoring (F1 = 0.56) outperforms SingleR (F1 = 0.25); CellTypist default fails on diverse tissue | Challenged | **extended** (see Annotation Methods) |
| P6 | Minimum cell filter (3) | 0% HVG loss across 30 dataset-tissue combinations | Validated | as_original |
| P7 | Normalization method | log-norm, scran, SCTransform produce identical AUC (0.739); log-norm 89x faster than scran | Validated | as_original |
| P8 | Integration method selection | scVI top in 19/30 datasets; decision tree predicts best method from n_donors and n_cell_types (80% accuracy) | Challenged | as_original |
| P9 | Clustering resolution (0.5) | Default over-clusters by ARI; optimal median = 0.15 | Challenged | **refined** (see A2) |

### Phase 1: bulk RNA-seq (B1–B4)

| Audit | Convention | Finding | Status | Relationship |
|---|---|---|---|---|
| B1 | DEG tool agreement | DESeq2-edgeR Jaccard = 0.52, 5.7x more concordant than scRNA-seq | Challenged | as_original |
| B2 | Normalization family | Length-based (TPM, FPKM) and count-based (TMM, VST) families share only 2% of top variable genes | Challenged | as_original |
| B3 | Low count filter (CPM > 1) | Robust gene retention (0.3% loss) but 2x DEG swing across thresholds | Challenged | as_original |
| B4 | Batch correction method | Method choice alters 80% of DEG lists | Challenged | **contradicted** (see A7) |

### Phase 2 / 2a: refinements to Phase 1

| Audit | Refines | Refined finding |
|---|---|---|
| A2 | P9 | Clustering optimal resolution depends on the selection metric. ARI gives median 0.15; V-measure gives median 0.50. The two metrics disagree by ~5x on 9/15 datasets. The choice of selection metric matters as much as the resolution itself. |
| A7 | B4 | The 80% DEG change result is preprocessing-dependent and SVA-specific. ComBat and limma findings hold under method-native preprocessing (DEG stability delta < 0.05). SVA on log-CPM produces 72% DEG change; svaseq on raw counts produces 43%. The cross-method "80% changes" headline overstates the result. |
| Phase 2a DE | P4 | The 3.8x Wilcoxon inflation is Wilcoxon-default specific. NEBULA's median inflation is ~0.85x (no inflation). muscat-dream is bimodal: conservative on most datasets, inflated on lymph node and skin fibroblast. The refined recommendation is to use NEBULA for cell-level DE; pseudoreplication is a Wilcoxon-default property, not a universal cell-level testing property. |

### Phase 2: additional audits

| Audit | Topic | Key finding | Rules contributed |
|---|---|---|---|
| Audit 1 main | Pathway enrichment | Tool and database choice both affect enrichment results; specific findings encoded as draft rules | Multiple |
| A6 | Ambient RNA correction | Method choice propagates to downstream biology on edge-case tissue; minimal effect on clean datasets | 4 |
| Audit 3 | Counting tools | Tool choice (CellRanger, STARsolo, Alevin-fry, kallisto) affects which cells are called and downstream biology in tissue-dependent ways | 4 |
| Audit QC-MAD | Low-quality cell filtering | MAD-based filtering (Heumos recommendation) and quantile-based filtering produce non-trivial disagreement on cell sets, worst in low-gene tissues (small intestine Jaccard 0.774) | 2 |
| Audit QC-MAD Propagation | Downstream consequences of QC method choice | QC filtering method choice produces a modest, tissue-INDEPENDENT effect on downstream biology (ARI 0.80–0.91). Cell-filtering choices don't propagate the way counts-reshaping choices do. | 1 |
| Audit Annotation Methods | Reference-based vs marker-based annotation | Reference-based tools (CellTypist, SingleR) agree more with each other than with marker-based tools, reflecting shared training-data assumptions rather than independent accuracy. For honest validation, pair one tool from each paradigm. | 1–2 |

### Cross-audit observation (documented, not yet rule-encoded)

Across three propagation audits (Audit 3 C3, Ambient Correction CP3, QC-MAD Propagation), a pattern emerged: choices that reshape counts (cell-calling, ambient correction) propagate to downstream biology on edge-case tissue; choices that only filter cells (QC method) do not. Three of three so far. Held back from rule encoding pending a fourth independent test, per the corpus's tier discipline.

### Companion: marker reliability atlas

264 aggregated marker scores across 80 datasets, 3.4M cells, 6 cell types (human and mouse). Used as input to P5 cell type annotation evaluation and as ongoing reference for marker-based annotation work.

## Methodology

Each audit is structured around the same template: question, working set, headline finding, specific findings, tier assignment (per §5.3.2 of the corpus standards), and a draft rule encoding the finding for downstream use. The §5.3.2 amendment provides numerical boundaries for confidence tiers (hard_default, conditional, flag_and_warn, literature_based) with explicit criteria for dataset count, tissue diversity, stability, and effect size. Tier assignments are checkable against the spec rather than left to judgment.

Every finding declares a `prior_audit_relationship` (as_original, refines_prior, contradicts_prior, extends_prior) and, where applicable, a `best_practices_relationship` mapping it to the relevant published best-practices paper (currently Heumos et al. 2023 for scRNA-seq; additional papers being added for bulk RNA-seq and other modalities).

See `METHODOLOGY.md` for the full specification and `BEST_PRACTICES_MAPPING.md` for the audit-to-paper relationship table.

## Data sources

- **scRNA-seq:** CELLxGENE Census v2025-11-08 (96.6M human + 18.4M mouse primary cells, queried via tiledbsoma lazy slices)
- **Bulk RNA-seq:** TCGA via GDC API (16 cancer types with ≥10 matched tumor-normal pairs); GTEx v10 raw counts (20 tissues)
- **Organisms:** Homo sapiens, Mus musculus

## Reproducibility

Every audit ships with a `repro.lock` file capturing the exact Python and R environment state at runtime, along with SHA-256 hashes of all output files. To verify reproduction:

```
repro restore
repro verify
```

`repro verify` exits non-zero on hash mismatch. Census queries are pinned to v2025-11-08. All result TSVs are committed alongside the code that produced them.

The repro tool (github.com/usarmydoc/repro) has been validated across the audits in this corpus by a single user. It has not yet been stress-tested by a wider contributor base. Outside contributions to this repo also exercise repro in the natural course of the audit work, which is one of the reasons contributions are welcome.

## Repository structure

```
bioinformatics-audits/
├── scrnaseq/
│   ├── shared/                 # Census access + metadata writer
│   ├── projects/p1–p9/         # Phase 1 projects: run_p*.py + emit_knowledge.py + repro.lock
│   ├── output/                 # Result TSVs
│   └── master_run.py           # Sequential orchestrator
├── bulk/
│   ├── shared/                 # Metadata writer
│   ├── projects/b1–b4/         # Phase 1 projects: run.py + data_access.py + analysis.py + repro.lock
│   └── output/                 # Result TSVs
├── phase2/
│   ├── audit1_main/            # Pathway enrichment
│   ├── audit3_counting/        # Counting tools (4 rules)
│   ├── audit_ambient_correction/  # Ambient RNA (4 rules)
│   ├── audit_qc_mad/           # QC method comparison (2 rules)
│   ├── audit_qc_mad_propagation/  # QC propagation (1 rule)
│   ├── audit_annotation_methods/  # Reference vs marker paradigms
│   └── DEFERRED.md             # Scoped follow-up candidates
├── marker_atlas/
│   ├── src/pipelines/          # Census access, scoring, reporting
│   └── output/                 # Atlas TSV + summary
├── BEST_PRACTICES_REGISTRY.md  # Canonical list of best-practices papers the corpus maps against
├── BEST_PRACTICES_MAPPING.md   # Audit-to-paper relationship table
├── METHODOLOGY.md              # Audit structure, tier criteria, schema
├── CONTRIBUTING.md             # How to submit an audit
├── requirements.txt
└── LICENSE
```

## Requirements

- Python 3.12+
- R 4.3+ with Seurat, scran, scDblFinder, SingleR, celldex, scuttle
- ~64 GB RAM for Census queries and scVI training
- GPU recommended for scVI (CUDA)

```
pip install -r requirements.txt
Rscript -e 'BiocManager::install(c("scDblFinder", "SingleR", "celldex", "scran", "scuttle"))'
```

## Running an audit

Each project runs independently:

```
cd scrnaseq
python projects/p1_mito_threshold/run_p1.py
python projects/p1_mito_threshold/emit_knowledge.py
```

Or all sequentially:

```
python scrnaseq/master_run.py
```

Runtimes range from minutes to several hours per audit; the full Phase 1 scRNA-seq corpus completes in ~12 hours on a workstation with GPU.

## Contributing

The corpus is designed to grow beyond a single author. New audits, refinements to existing audits, and extensions to new modalities are all welcome. See `CONTRIBUTING.md` for the audit template, quality bar, and review process.

Outside contributions also stress-test the repro tooling that ships with each audit, which has so far only been validated by a single user across this corpus. Failures, edge cases, and platform-specific issues surfaced by outside use are themselves valuable data.
