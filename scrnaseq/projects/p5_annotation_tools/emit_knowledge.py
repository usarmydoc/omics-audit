#!/usr/bin/env python3
"""P5 — Emit BioOrchestrator knowledge artifacts for annotation tool agreement audit.

Generates:
- L2: Prompt enrichment paragraph for scrna/v1.txt
- L3: 2 new validation rules for scrna_seq.yaml
- L4: Snakefile template annotation_tools.smk
"""

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
    """Append P5 findings to scrna/v1.txt."""
    path = PROMPTS_DIR / "v1.txt"

    n_datasets = result_df.dataset_id.nunique()

    # Per-tool summary
    tool_summaries = []
    for tool, tdf in result_df.groupby("tool", observed=True):
        mean_f1 = tdf.f1.mean()
        mean_prec = tdf.precision.mean()
        mean_rec = tdf.recall.mean()
        tool_summaries.append(f"{tool}: mean F1={mean_f1:.3f}, P={mean_prec:.3f}, R={mean_rec:.3f}")

    tool_text = "; ".join(tool_summaries)

    # Best and worst cell types across all tools
    ct_f1 = result_df.groupby("cell_type", observed=True)["f1"].mean().sort_values()
    worst_ct = ct_f1.index[0] if len(ct_f1) > 0 else "unknown"
    worst_f1 = ct_f1.iloc[0] if len(ct_f1) > 0 else 0
    best_ct = ct_f1.index[-1] if len(ct_f1) > 0 else "unknown"
    best_f1 = ct_f1.iloc[-1] if len(ct_f1) > 0 else 0

    # MarkerScore-specific stats
    ms_df = result_df[result_df.tool == "MarkerScore"]
    ms_f1 = ms_df.f1.mean() if len(ms_df) > 0 else 0

    paragraph = f"""
## Empirical Guidance: Annotation Tool Accuracy (P5 Audit)

Analysis of {n_datasets} datasets from CELLxGENE Census comparing SingleR, CellTypist,
and a marker-score classifier against Census ground truth labels shows variable
accuracy across tools and cell types. Tool performance: {tool_text}. Accuracy
varies substantially by cell type: {best_ct} annotations are most reliable
(mean F1={best_f1:.3f}), while {worst_ct} is hardest to classify (mean
F1={worst_f1:.3f}). The marker-score classifier (mean F1={ms_f1:.3f}) uses
markers filtered by reliability_score > 0.5 from the Marker Reliability Atlas.
When designing annotation pipelines: (1) always specify the annotation tool
explicitly rather than relying on defaults, (2) if using a marker-score classifier,
apply a reliability threshold to exclude unreliable markers, and (3) consider
running multiple tools and taking consensus annotations for high-confidence results.
"""

    if path.exists():
        existing = path.read_text()
        if "Annotation Tool Accuracy" not in existing:
            content = existing.rstrip() + "\n\n" + paragraph
        else:
            logger.info("P5 guidance already present, skipping")
            return path
    else:
        content = paragraph

    path.write_text(content)
    logger.info("Wrote P5 prompt enrichment to %s", path)
    return path


def emit_validation_rules() -> Path:
    """Append 2 new P5 rules to scrna_seq.yaml."""
    new_rules = [
        {
            "id": "scrna_annotation_tool_explicit",
            "type": "conditional",
            "description": "Warn if no annotation tool is explicitly specified in the pipeline",
            "condition": {
                "context": "experiment",
                "field": "constraints",
                "check": "key_value_gt",
                "key": "expected_cells",
                "threshold": 0,
                "fallback": True,
            },
            "required": {
                "scope": "steps",
                "tool_name_any": [
                    "SingleR",
                    "CellTypist",
                    "celltypist",
                    "scType",
                    "marker_score",
                    "cell_type_annotation",
                ],
            },
            "action": "warn",
            "message": (
                "No cell type annotation tool is explicitly specified in the pipeline. "
                "Empirical audit of CELLxGENE Census data shows annotation accuracy varies "
                "substantially across tools and cell types. SingleR, CellTypist, and marker-score "
                "classifiers each have distinct strengths: reference-based methods (SingleR) "
                "excel on well-represented types, while marker-score classifiers depend on marker "
                "reliability. Always specify the annotation tool explicitly and consider running "
                "multiple tools for consensus annotations."
            ),
        },
        {
            "id": "scrna_marker_score_threshold",
            "type": "conditional",
            "description": "Flag if marker-score classifier is used without a reliability threshold",
            "condition": {
                "context": "experiment",
                "field": "steps",
                "check": "tool_name_present",
                "tool_names": [
                    "marker_score",
                    "marker_scoring",
                    "AddModuleScore",
                    "score_genes",
                ],
            },
            "required": {
                "scope": "parameters",
                "parameter_any": [
                    "annotation.reliability_threshold",
                    "annotation.min_reliability_score",
                    "marker_score.reliability_threshold",
                ],
            },
            "action": "flag",
            "message": (
                "A marker-score classifier is in the pipeline but no reliability threshold "
                "is specified. Empirical audit shows that markers with reliability_score <= 0.5 "
                "(from the Marker Reliability Atlas) degrade annotation accuracy significantly. "
                "Set annotation.reliability_threshold >= 0.5 to filter low-reliability markers "
                "before scoring. Markers like CD4 (T cells, AUC ~0.60) and MARCO (macrophages) "
                "are unreliable and should be excluded from score-based classifiers."
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
    """Write annotation_tools.smk Snakefile template."""
    path = SNAKEFILE_DIR / "annotation_tools.smk"
    path.parent.mkdir(parents=True, exist_ok=True)

    content = '''# Cell Type Annotation Tools — Snakefile template for BioOrchestrator L4
# Runs SingleR, CellTypist, and marker-score annotation, then compares to ground truth.
# Generated by P5 emit_knowledge.py

import os

WORKDIR = config.get("workdir", ".")
OUTPUT_DIR = os.path.join(WORKDIR, "output")
LOG_DIR = os.path.join(WORKDIR, "logs")
SRC_DIR = os.path.join(WORKDIR, "src")
INPUT_H5AD = config.get("input_h5ad", os.path.join(WORKDIR, "input.h5ad"))
ATLAS_TSV = config.get("marker_atlas", os.path.join(WORKDIR, "marker_reliability_atlas.tsv"))
ORGANISM = config.get("organism", "Homo sapiens")
RELIABILITY_THRESHOLD = config.get("reliability_threshold", 0.5)


rule all:
    input:
        os.path.join(OUTPUT_DIR, "repro_verified.flag"),


rule repro_snapshot:
    """Capture environment state before annotation."""
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


rule run_singler:
    """Annotate cell types using SingleR with celldex reference."""
    input:
        h5ad=INPUT_H5AD,
        lockfile=rules.repro_snapshot.output.lockfile,
    output:
        predictions=os.path.join(OUTPUT_DIR, "singler_predictions.tsv"),
    params:
        organism=ORGANISM,
    log:
        os.path.join(LOG_DIR, "singler.log"),
    resources:
        mem_mb=32000,
        runtime=60,
    shell:
        """
        Rscript -e "
library(SingleR)
library(SingleCellExperiment)
library(celldex)
library(zellkonverter)

adata <- readH5AD('{input.h5ad}')
sce <- as(adata, 'SingleCellExperiment')

if ('{params.organism}' == 'Homo sapiens') {{
    ref <- HumanPrimaryCellAtlasData()
}} else {{
    ref <- MouseRNAseqData()
}}

results <- SingleR(test=sce, ref=ref, labels=ref\\$label.main)
write.table(
    data.frame(barcode=rownames(results), singler_label=results\\$labels),
    '{output.predictions}', sep='\\t', row.names=FALSE, quote=FALSE
)
        " 2>&1 | tee {log}
        """


rule run_celltypist:
    """Annotate cell types using CellTypist default model."""
    input:
        h5ad=INPUT_H5AD,
        lockfile=rules.repro_snapshot.output.lockfile,
    output:
        predictions=os.path.join(OUTPUT_DIR, "celltypist_predictions.tsv"),
    log:
        os.path.join(LOG_DIR, "celltypist.log"),
    resources:
        mem_mb=16000,
        runtime=30,
    shell:
        """
        python -c "
import scanpy as sc
import celltypist
from celltypist import models
import pandas as pd

adata = sc.read_h5ad('{input.h5ad}')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

models.download_models(force_update=False)
model = models.Model.load()
predictions = celltypist.annotate(adata, model=model)

pred_df = predictions.predicted_labels.reset_index()
pred_df.columns = ['barcode', 'celltypist_label']
pred_df.to_csv('{output.predictions}', sep='\\t', index=False)
        " 2>&1 | tee {log}
        """


rule run_marker_score:
    """Annotate cell types using marker-score classifier with reliability filtering."""
    input:
        h5ad=INPUT_H5AD,
        atlas=ATLAS_TSV,
        lockfile=rules.repro_snapshot.output.lockfile,
    output:
        predictions=os.path.join(OUTPUT_DIR, "marker_score_predictions.tsv"),
    params:
        organism=ORGANISM,
        threshold=RELIABILITY_THRESHOLD,
    log:
        os.path.join(LOG_DIR, "marker_score.log"),
    resources:
        mem_mb=16000,
        runtime=20,
    shell:
        """
        python -c "
import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('{input.h5ad}')
atlas = pd.read_csv('{input.atlas}', sep='\\t')
atlas = atlas[(atlas.organism == '{params.organism}') & (atlas.reliability_score > {params.threshold})]

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

gene_names = list(adata.var_names)
scores = {{}}
for ct, ct_df in atlas.groupby('cell_type'):
    markers = [g for g in ct_df.gene.tolist() if g in gene_names]
    if markers:
        marker_idx = [gene_names.index(g) for g in markers]
        scores[ct] = np.asarray(adata.X[:, marker_idx].mean(axis=1)).ravel()

score_df = pd.DataFrame(scores, index=adata.obs_names)
predictions = score_df.idxmax(axis=1)

pred_df = pd.DataFrame({{'barcode': adata.obs_names, 'marker_score_label': predictions.values}})
pred_df.to_csv('{output.predictions}', sep='\\t', index=False)
        " 2>&1 | tee {log}
        """


rule compare_annotations:
    """Compare all tool predictions against ground truth and compute metrics."""
    input:
        singler=rules.run_singler.output.predictions,
        celltypist=rules.run_celltypist.output.predictions,
        marker_score=rules.run_marker_score.output.predictions,
        h5ad=INPUT_H5AD,
    output:
        metrics=os.path.join(OUTPUT_DIR, "annotation_metrics.tsv"),
        summary=os.path.join(OUTPUT_DIR, "annotation_summary.tsv"),
    log:
        os.path.join(LOG_DIR, "compare_annotations.log"),
    resources:
        mem_mb=8000,
        runtime=10,
    shell:
        """
        python -c "
import scanpy as sc
import pandas as pd
from sklearn.metrics import precision_recall_fscore_support

adata = sc.read_h5ad('{input.h5ad}')
gt = adata.obs['cell_type']

singler = pd.read_csv('{input.singler}', sep='\\t').set_index('barcode')
celltypist = pd.read_csv('{input.celltypist}', sep='\\t').set_index('barcode')
marker_score = pd.read_csv('{input.marker_score}', sep='\\t').set_index('barcode')

all_metrics = []
for tool_name, pred_col in [('SingleR', singler), ('CellTypist', celltypist), ('MarkerScore', marker_score)]:
    col = pred_col.columns[0]
    preds = pred_col.loc[gt.index, col]
    labels = sorted(set(gt.unique()))
    p, r, f, s = precision_recall_fscore_support(gt, preds, labels=labels, zero_division=0)
    for i, label in enumerate(labels):
        all_metrics.append(dict(tool=tool_name, cell_type=label, precision=p[i], recall=r[i], f1=f[i], n_cells=int(s[i])))

pd.DataFrame(all_metrics).to_csv('{output.metrics}', sep='\\t', index=False)

summary = pd.DataFrame(all_metrics).groupby('tool')[['precision','recall','f1']].mean().reset_index()
summary.to_csv('{output.summary}', sep='\\t', index=False)
        " 2>&1 | tee {log}
        """


rule repro_verify:
    """Hash outputs and verify reproducibility."""
    input:
        metrics=rules.compare_annotations.output.metrics,
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
            echo "ERROR: repro verify failed" | tee -a {log}
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
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    tsv_path = OUTPUT_DIR / "p5_annotation_tools.tsv"
    if not tsv_path.exists():
        logger.error("P5 output not found: %s — run run_p5.py first", tsv_path)
        sys.exit(1)

    result_df = pd.read_csv(tsv_path, sep="\t")
    logger.info("Loaded %d rows from %s", len(result_df), tsv_path)

    emit_prompt_enrichment(result_df)
    emit_validation_rules()
    emit_snakefile_template()

    logger.info("P5 knowledge emission complete")


if __name__ == "__main__":
    main()
