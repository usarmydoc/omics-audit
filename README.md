# An Empirical Cross-Atlas Audit of scRNA-seq Analytical Conventions

Code and data for the cross-modality audit of single-cell and bulk RNA-seq analytical conventions using CELLxGENE Census and public datasets.

## Overview

This repository contains 14 audit projects that empirically test commonly used defaults in RNA-seq analysis pipelines:

### scRNA-seq Audits (P1–P9)

| Project | Question | Key Finding |
|---------|----------|-------------|
| **P1** Mito Threshold | Is 20% mitochondrial threshold universal? | Tissue-dependent; removes 0% in 55% of datasets, >10% in gut |
| **P2** Doublet Tools | Do scDblFinder and Scrublet agree? | Cohen's kappa = 0.08; 20x rate difference |
| **P3** Batch Correction | Harmony vs Scanorama vs scVI? | scVI (ARI=0.56) >> Scanorama (ARI=0.30) |
| **P4** Pseudobulk | Pseudobulk vs Wilcoxon DE? | 9% DEG overlap; Wilcoxon inflates 3.8x |
| **P5** Annotation | Marker-score classifier accuracy? | F1=0.63; fibroblasts best, hepatocytes worst |
| **P6** Min Cell Filter | Is min_cells=3 justified? | **Validated**; threshold=10 loses only 2.1% HVGs |
| **P7** Normalization | log-norm vs SCTransform vs scran? | **All equivalent** (AUC ~0.82); log-norm 600x faster |
| **P8** Integration | Which method for which structure? | scVI preferred; decision tree inconclusive |
| **P9** Clustering Resolution | Is Leiden resolution=0.5 optimal? | Over-clusters; optimal median = 0.15 |

### Bulk RNA-seq Audits (B1–B4)

| Project | Question | Key Finding |
|---------|----------|-------------|
| **B1** DEG Tool Agreement | DESeq2 vs edgeR vs limma? | Jaccard = 0.52 (5.7x better than scRNA) |
| **B2** Normalization | TPM/FPKM vs TMM/VST? | Share only 2% of top genes |
| **B3** Low Count Filter | Is CPM>1 filter optimal? | Loses 0.3% tissue-specific genes |
| **B4** Batch Correction | ComBat vs limma vs SVA? | ComBat best (ARI 0.73→0.88) |

### Marker Reliability Atlas

Cross-dataset marker gene validation: 264 aggregated scores across 80 datasets, 3.4M cells, 6 cell types (human + mouse).

## Data Sources

- **scRNA-seq**: CELLxGENE Census v2025-11-08 (96.6M human + 18.4M mouse primary cells)
- **Bulk RNA-seq**: TCGA (16 cancer types) + GTEx (20 tissues)
- **Organisms**: Homo sapiens, Mus musculus

## Repository Structure

```
omics-audit/
├── scrnaseq/                  # scRNA-seq audits (P1–P9)
│   ├── shared/                # Census access + knowledge writer
│   ├── projects/p1–p9/        # Each project: run_p*.py + emit_knowledge.py
│   ├── output/                # Result TSVs + repro.lock files
│   └── master_run.py          # Sequential orchestrator
├── bulk/                      # Bulk RNA-seq audits (B1–B4)
│   ├── shared/                # Knowledge writer
│   ├── projects/b1–b4/        # Each project: run.py + data_access.py + analysis.py
│   └── output/                # Result TSVs
├── marker_atlas/              # Marker reliability atlas
│   ├── src/pipelines/         # Census access, scoring, reporting
│   └── output/                # Atlas TSV + summary
├── paper/                     # Manuscript drafts
├── requirements.txt
└── LICENSE
```

## Requirements

- Python 3.12+
- R 4.3+ with: Seurat, scran, scDblFinder, SingleR, celldex, scuttle
- ~64 GB RAM (for Census queries and scVI training)
- GPU recommended for scVI (CUDA)

```bash
pip install -r requirements.txt
Rscript -e 'BiocManager::install(c("scDblFinder", "SingleR", "celldex", "scran", "scuttle"))'
```

## Running

Each project runs independently:

```bash
cd scrnaseq
python projects/p1_mito_threshold/run_p1.py      # ~1 hour
python projects/p1_mito_threshold/emit_knowledge.py
```

Or run all sequentially:

```bash
python scrnaseq/master_run.py  # ~12 hours total
```

## Reproducibility

All outputs include `repro.lock` files (via [repro-lock](https://pypi.org/project/repro-lock/)) capturing the exact environment state at runtime. Census queries are pinned to version `2025-11-08`.

## BioOrchestrator Integration

Audit findings are encoded as machine-readable artifacts for [BioOrchestrator](https://github.com/usarmydoc/bioorchestrator):
- **L2**: 5 prompt enrichment paragraphs (empirical guidance for LLM pipeline design)
- **L3**: 18 validation rules appended to `scrna_seq.yaml`
- **L4**: 13 Snakefile templates with real shell commands

## Citation

If you use this code or data, please cite:

> [Manuscript in preparation]
