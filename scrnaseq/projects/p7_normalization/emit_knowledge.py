#!/usr/bin/env python3
"""P7 — Emit BioOrchestrator knowledge artifacts for normalization comparison."""

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
    """L2: Append normalization guidance paragraph to v1.txt."""
    path = PROMPTS_DIR / "v1.txt"

    n_datasets = result_df.groupby(["dataset_id", "tissue", "organism"], observed=True).ngroups

    # Compute per-method summaries
    summaries = {}
    for method in ["log_normalize", "sctransform", "scran"]:
        sub = result_df[result_df.method == method]
        if len(sub) == 0:
            continue
        auc_vals = sub.mean_marker_auc.dropna()
        sil_vals = sub.silhouette_score.dropna()
        rt_vals = sub.runtime_seconds.dropna()
        summaries[method] = {
            "n": len(sub),
            "mean_auc": auc_vals.mean() if len(auc_vals) > 0 else 0.0,
            "mean_sil": sil_vals.mean() if len(sil_vals) > 0 else 0.0,
            "mean_rt": rt_vals.mean() if len(rt_vals) > 0 else 0.0,
        }

    ln = summaries.get("log_normalize", {})
    sct = summaries.get("sctransform", {})
    scr = summaries.get("scran", {})

    paragraph = f"""
## Empirical Guidance: Normalization Method Selection (P7 Audit)

Analysis of {n_datasets} datasets from CELLxGENE Census compared three normalization
methods — log-normalize (scanpy), SCTransform (Seurat), and scran (Bioconductor) —
on identical raw counts. Mean marker AUC scores were: log-normalize={ln.get('mean_auc', 0):.3f},
SCTransform={sct.get('mean_auc', 0):.3f}, scran={scr.get('mean_auc', 0):.3f}.
Mean silhouette scores were: log-normalize={ln.get('mean_sil', 0):.3f},
SCTransform={sct.get('mean_sil', 0):.3f}, scran={scr.get('mean_sil', 0):.3f}.
Runtime was fastest for log-normalize ({ln.get('mean_rt', 0):.0f}s), followed by
scran ({scr.get('mean_rt', 0):.0f}s) and SCTransform ({sct.get('mean_rt', 0):.0f}s).
When designing pipelines: (1) always explicitly specify the normalization method rather
than relying on tool defaults, (2) SCTransform should use sketch-based mode for datasets
>100K cells, (3) log-normalize is a safe default when runtime matters or when downstream
tools expect log-normalized input, (4) scran is preferred for heterogeneous datasets
where cell composition varies strongly across samples.
"""

    if path.exists():
        existing = path.read_text()
        if "Normalization Method Selection" not in existing:
            content = existing.rstrip() + "\n\n" + paragraph
        else:
            logger.info("P7 prompt enrichment already present, skipping")
            return path
    else:
        content = paragraph
    path.write_text(content)
    logger.info("Wrote P7 prompt enrichment to %s", path)
    return path


def emit_validation_rules() -> Path:
    """L3: Add 2 validation rules to scrna_seq.yaml."""
    new_rules = [
        {
            "id": "scrna_normalization_explicit",
            "type": "conditional",
            "description": "Normalization method must be explicitly specified, not left to tool defaults",
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
                    "normalization.method",
                    "normalization_method",
                    "normalize.method",
                ],
            },
            "action": "warn",
            "message": (
                "No explicit normalization method specified. Empirical comparison of "
                "log-normalize, SCTransform, and scran across CELLxGENE Census datasets "
                "shows all three produce similar marker AUC scores, but differ in runtime "
                "and HVG selection. Specify the normalization method explicitly to ensure "
                "reproducibility: use log-normalize for speed, SCTransform for Seurat pipelines "
                "with moderate cell counts, or scran for heterogeneous datasets."
            ),
        },
        {
            "id": "scrna_sctransform_large_dataset",
            "type": "conditional",
            "description": "SCTransform on >100K cells should use sketch-based approach",
            "condition": {
                "context": "experiment",
                "field": "constraints",
                "check": "key_value_gt",
                "key": "expected_cells",
                "threshold": 100000,
                "fallback": False,
            },
            "required": {
                "scope": "parameters",
                "parameter_any": [
                    "normalization.sketch",
                    "sctransform.sketch",
                    "SCTransform.ncells",
                ],
            },
            "action": "warn",
            "message": (
                "SCTransform is specified for a dataset with >100K cells but no sketch-based "
                "mode is configured. Standard SCTransform fits a regularized NB model on all "
                "cells, which is memory-intensive and slow at scale. Use "
                "SCTransform(vst.flavor='v2', ncells=5000) or the sketch-based workflow "
                "(Seurat v5) to subsample cells for model fitting while normalizing all cells. "
                "Alternatively, use log-normalize or scran which scale more linearly."
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
        else:
            logger.info("Rule %s already exists, skipping", rule["id"])

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
    """L4: Write normalization_comparison.smk Snakefile template."""
    path = SNAKEFILE_DIR / "normalization_comparison.smk"
    path.parent.mkdir(parents=True, exist_ok=True)

    content = '''# Normalization Method Comparison — Snakefile template for BioOrchestrator L4
# Compares log-normalize, SCTransform, and scran on identical raw counts.
# Generated by P7 emit_knowledge.py

import os

WORKDIR = config.get("workdir", ".")
OUTPUT_DIR = os.path.join(WORKDIR, "output")
LOG_DIR = os.path.join(WORKDIR, "logs")
SRC_DIR = os.path.join(WORKDIR, "src", "pipelines")
INPUT_H5AD = config.get("input_h5ad", os.path.join(WORKDIR, "input.h5ad"))
METHODS = config.get("methods", ["log_normalize", "sctransform", "scran"])


rule all:
    input:
        os.path.join(OUTPUT_DIR, "repro_verified.flag"),


rule repro_snapshot:
    """Capture environment state before any data access."""
    output:
        lockfile=os.path.join(OUTPUT_DIR, "repro.lock"),
    log:
        os.path.join(LOG_DIR, "repro_snapshot.log"),
    resources:
        mem_mb=512,
        runtime=5,
    shell:
        """
        repro snapshot -o {output.lockfile} --quiet 2>&1 | tee {log}
        """


rule normalize_lognorm:
    """Log-normalize: scanpy normalize_total + log1p."""
    input:
        h5ad=INPUT_H5AD,
        lockfile=rules.repro_snapshot.output.lockfile,
    output:
        h5ad=os.path.join(OUTPUT_DIR, "normalized_log_normalize.h5ad"),
    log:
        os.path.join(LOG_DIR, "normalize_lognorm.log"),
    resources:
        mem_mb=16000,
        runtime=30,
    shell:
        """
        python -c "
import scanpy as sc
adata = sc.read_h5ad('{input.h5ad}')
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.write('{output.h5ad}')
        " 2>&1 | tee {log}
        """


rule normalize_sctransform:
    """SCTransform via R/Seurat."""
    input:
        h5ad=INPUT_H5AD,
        lockfile=rules.repro_snapshot.output.lockfile,
    output:
        h5ad=os.path.join(OUTPUT_DIR, "normalized_sctransform.h5ad"),
    log:
        os.path.join(LOG_DIR, "normalize_sctransform.log"),
    resources:
        mem_mb=32000,
        runtime=60,
    shell:
        """
        python {SRC_DIR}/normalize_sctransform.py \\
            --input {input.h5ad} \\
            --output {output.h5ad} \\
            2>&1 | tee {log}
        """


rule normalize_scran:
    """scran deconvolution size factors + logNormCounts via R."""
    input:
        h5ad=INPUT_H5AD,
        lockfile=rules.repro_snapshot.output.lockfile,
    output:
        h5ad=os.path.join(OUTPUT_DIR, "normalized_scran.h5ad"),
    log:
        os.path.join(LOG_DIR, "normalize_scran.log"),
    resources:
        mem_mb=32000,
        runtime=60,
    shell:
        """
        python {SRC_DIR}/normalize_scran.py \\
            --input {input.h5ad} \\
            --output {output.h5ad} \\
            2>&1 | tee {log}
        """


rule compute_metrics:
    """Compute marker AUC, silhouette, HVG overlap for each normalization."""
    input:
        lognorm=rules.normalize_lognorm.output.h5ad,
        sctransform=rules.normalize_sctransform.output.h5ad,
        scran=rules.normalize_scran.output.h5ad,
    output:
        tsv=os.path.join(OUTPUT_DIR, "normalization_comparison.tsv"),
    log:
        os.path.join(LOG_DIR, "compute_metrics.log"),
    resources:
        mem_mb=32000,
        runtime=30,
    shell:
        """
        python {SRC_DIR}/normalization_metrics.py \\
            --lognorm {input.lognorm} \\
            --sctransform {input.sctransform} \\
            --scran {input.scran} \\
            --output {output.tsv} \\
            2>&1 | tee {log}
        """


rule repro_verify:
    """Hash all output files and verify reproducibility."""
    input:
        tsv=rules.compute_metrics.output.tsv,
        lockfile=rules.repro_snapshot.output.lockfile,
    output:
        flag=os.path.join(OUTPUT_DIR, "repro_verified.flag"),
    log:
        os.path.join(LOG_DIR, "repro_verify.log"),
    resources:
        mem_mb=512,
        runtime=5,
    shell:
        """
        repro verify {OUTPUT_DIR} -l {input.lockfile} 2>&1 | tee {log}
        VERIFY_EXIT=$?
        if [ $VERIFY_EXIT -ne 0 ]; then
            echo "ERROR: repro verify failed — output files may not be reproducible" | tee -a {log}
            exit 1
        fi
        touch {output.flag}
        echo "repro verify passed" | tee -a {log}
        """
'''

    path.write_text(content)
    logger.info("Wrote Snakefile template: %s", path)
    return path


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    tsv_path = OUTPUT_DIR / "p7_normalization.tsv"
    if not tsv_path.exists():
        logger.error("P7 output not found: %s", tsv_path)
        sys.exit(1)

    result_df = pd.read_csv(tsv_path, sep="\t")
    logger.info("Loaded %d rows from %s", len(result_df), tsv_path)

    emit_prompt_enrichment(result_df)
    emit_validation_rules()
    emit_snakefile_template()
    logger.info("P7 knowledge emission complete")


if __name__ == "__main__":
    main()
