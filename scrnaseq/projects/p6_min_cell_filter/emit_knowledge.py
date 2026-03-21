#!/usr/bin/env python3
"""P6 — Emit BioOrchestrator knowledge artifacts for min cell filter audit."""

from __future__ import annotations

import logging
import os
import sys
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
BIOORCHESTRATOR_ROOT = Path(os.environ.get("BIOORCHESTRATOR_ROOT", str(Path.home() / "bioorchestrator" / "src" / "bioorchestrator")))
PROMPTS_DIR = BIOORCHESTRATOR_ROOT / "knowledge" / "prompts" / "scrna"
RULES_FILE = BIOORCHESTRATOR_ROOT / "knowledge" / "rules" / "scrna_seq.yaml"
SNAKEFILE_DIR = BIOORCHESTRATOR_ROOT / "knowledge" / "snakefile_templates"
OUTPUT_DIR = REPO_ROOT / "scrnaseq" / "output"


def emit_prompt_enrichment(result_df: pd.DataFrame) -> Path:
    path = PROMPTS_DIR / "v1.txt"
    n_datasets = result_df.groupby(["dataset_id", "tissue", "organism"]).ngroups
    t3 = result_df[result_df.threshold == 3]
    t10 = result_df[result_df.threshold == 10]
    t50 = result_df[result_df.threshold == 50]
    mean_hvg_lost_10 = t10.pct_hvg_lost.mean() if len(t10) > 0 else 0
    mean_hvg_lost_50 = t50.pct_hvg_lost.mean() if len(t50) > 0 else 0

    paragraph = f"""
## Empirical Guidance: Minimum Cell Filter Thresholds (P6 Audit)

Analysis of {n_datasets} datasets from CELLxGENE Census shows the default
min_cells=3 filter is generally conservative. Increasing to threshold=10
loses {mean_hvg_lost_10:.1f}% of HVGs on average, while threshold=50 loses
{mean_hvg_lost_50:.1f}%. For large datasets (>100K cells), a threshold of 10-20
removes noise genes without affecting HVG detection. For small datasets (<5K cells),
the default of 3 is appropriate. When designing pipelines, scale the min_cells
threshold with dataset size: use 3 for <5K cells, 10 for 5K-100K, and 20 for >100K.
Always explicitly set this parameter rather than relying on tool defaults.
"""

    if path.exists():
        existing = path.read_text()
        if "Minimum Cell Filter" not in existing:
            content = existing.rstrip() + "\n\n" + paragraph
        else:
            return path
    else:
        content = paragraph
    path.write_text(content)
    logger.info("Wrote P6 prompt enrichment to %s", path)
    return path


def emit_validation_rules() -> Path:
    new_rules = [
        {
            "id": "scrna_min_cells_explicit",
            "type": "conditional",
            "description": "min_cells gene filter must be explicitly set, not left to default",
            "condition": {
                "context": "experiment",
                "field": "constraints",
                "check": "key_value_gt",
                "key": "expected_cells",
                "threshold": 0,
                "fallback": True,
            },
            "required": {
                "scope": "parameters",
                "parameter_any": [
                    "qc_thresholds.min_cells_per_gene",
                    "qc_thresholds.min_cells",
                ],
            },
            "action": "warn",
            "message": (
                "No explicit min_cells gene filter parameter found. The default of 3 cells "
                "is appropriate for small datasets (<5K cells) but too permissive for large "
                "datasets (>100K cells), where it retains noise genes. Empirical audit across "
                "CELLxGENE Census shows increasing to 10-20 for large datasets removes noise "
                "without affecting HVG detection. Set min_cells_per_gene explicitly."
            ),
        },
        {
            "id": "scrna_min_cells_large_dataset",
            "type": "conditional",
            "description": "Flag large datasets (>500K cells) where min_cells=3 may be too permissive",
            "condition": {
                "context": "experiment",
                "field": "constraints",
                "check": "key_value_gt",
                "key": "expected_cells",
                "threshold": 500000,
                "fallback": False,
            },
            "required": {
                "scope": "parameters",
                "parameter_any": [
                    "qc_thresholds.min_cells_per_gene",
                ],
            },
            "action": "flag",
            "message": (
                "Dataset has >500K cells. The default min_cells=3 filter retains many noise "
                "genes in large datasets. Consider min_cells=20-50 to reduce dimensionality "
                "without losing biologically relevant highly variable genes."
            ),
        },
    ]

    path = RULES_FILE
    with open(path) as f:
        data = yaml.safe_load(f)

    existing_ids = {r["id"] for r in data.get("rules", [])}
    added = 0
    for rule in new_rules:
        if rule["id"] not in existing_ids:
            data["rules"].append(rule)
            added += 1
            logger.info("Added rule: %s", rule["id"])

    with open(path, "w") as f:
        f.write("# L3 validation rules for single-cell RNA-seq pipelines.\n")
        f.write("# Each rule is evaluated deterministically against a PipelineConfig + ExperimentContext.\n")
        f.write("# Rule types: range, conditional, compatibility\n")
        f.write("# Actions: reject (hard fail), warn (proceed with warning), flag (informational)\n")
        f.write("\n")
        yaml.dump(data, f, default_flow_style=False, sort_keys=False, width=120, allow_unicode=True)

    logger.info("Wrote %d new rules (total: %d)", added, len(data["rules"]))
    return path


def emit_snakefile_template() -> Path:
    path = SNAKEFILE_DIR / "min_cell_filter_audit.smk"
    path.parent.mkdir(parents=True, exist_ok=True)

    content = '''# Minimum Cell Filter Audit — Snakefile template for BioOrchestrator L4
# Tests multiple min_cells thresholds and reports HVG stability.
# Generated by P6 emit_knowledge.py

import os

WORKDIR = config.get("workdir", ".")
OUTPUT_DIR = os.path.join(WORKDIR, "output")
LOG_DIR = os.path.join(WORKDIR, "logs")
INPUT_H5AD = config.get("input_h5ad", os.path.join(WORKDIR, "input.h5ad"))
THRESHOLDS = config.get("thresholds", [1, 3, 5, 10, 20, 50])


rule all:
    input:
        os.path.join(OUTPUT_DIR, "repro_verified.flag"),


rule repro_snapshot:
    output:
        lockfile=os.path.join(OUTPUT_DIR, "repro.lock"),
    shell:
        """
        repro snapshot -o {output.lockfile} --quiet
        """


rule filter_sweep:
    """Apply each min_cells threshold and compute HVG overlap."""
    input:
        h5ad=INPUT_H5AD,
        lockfile=rules.repro_snapshot.output.lockfile,
    output:
        tsv=os.path.join(OUTPUT_DIR, "min_cell_filter_results.tsv"),
    log:
        os.path.join(LOG_DIR, "filter_sweep.log"),
    resources:
        mem_mb=32000,
    shell:
        """
        python -c "
import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_h5ad('{input.h5ad}')
thresholds = {THRESHOLDS}
results = []
# Baseline HVGs at threshold=3
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3')
hvg_baseline = set(adata.var_names[adata.var.highly_variable])

for t in thresholds:
    ad = sc.read_h5ad('{input.h5ad}')
    sc.pp.filter_genes(ad, min_cells=t)
    n_genes = ad.n_vars
    sc.pp.highly_variable_genes(ad, n_top_genes=min(2000, n_genes-1), flavor='seurat_v3')
    hvg_t = set(ad.var_names[ad.var.highly_variable])
    overlap = len(hvg_t & hvg_baseline) / len(hvg_t | hvg_baseline) if hvg_t | hvg_baseline else 1
    results.append(dict(threshold=t, genes_retained=n_genes, hvg_overlap=overlap))

pd.DataFrame(results).to_csv('{output.tsv}', sep='\\t', index=False)
        " 2>&1 | tee {log}
        """


rule repro_verify:
    input:
        tsv=rules.filter_sweep.output.tsv,
        lockfile=rules.repro_snapshot.output.lockfile,
    output:
        flag=os.path.join(OUTPUT_DIR, "repro_verified.flag"),
    shell:
        """
        repro verify {OUTPUT_DIR} -l {input.lockfile}
        touch {output.flag}
        """
'''

    path.write_text(content)
    logger.info("Wrote Snakefile template: %s", path)
    return path


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    tsv_path = OUTPUT_DIR / "p6_min_cell_filter.tsv"
    if not tsv_path.exists():
        logger.error("P6 output not found: %s", tsv_path)
        sys.exit(1)

    result_df = pd.read_csv(tsv_path, sep="\t")
    logger.info("Loaded %d rows from %s", len(result_df), tsv_path)

    emit_prompt_enrichment(result_df)
    emit_validation_rules()
    emit_snakefile_template()
    logger.info("P6 knowledge emission complete")


if __name__ == "__main__":
    main()
