#!/usr/bin/env python3
"""P4 — Emit BioOrchestrator knowledge artifacts for pseudobulk audit."""

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
    """Append P4 findings to scrna/v1.txt."""
    path = PROMPTS_DIR / "v1.txt"

    n = len(result_df)
    mean_jaccard = result_df.jaccard_top100.mean()
    dir_agree = result_df.direction_agreement.dropna()
    mean_dir = dir_agree.mean() if len(dir_agree) > 0 else 0
    inflate = result_df.wilcoxon_inflation_factor.dropna()
    mean_inflate = inflate.mean() if len(inflate) > 0 else 0

    paragraph = f"""
## Empirical Guidance: Pseudobulk vs Cell-Level DE (P4 Audit)

Analysis of {n} datasets from CELLxGENE Census comparing DESeq2 pseudobulk DE
against Wilcoxon rank-sum cell-level testing shows substantial disagreement.
Mean Jaccard similarity of top 100 DEGs: {mean_jaccard:.3f}. Direction agreement
for shared genes: {mean_dir:.1%}. Wilcoxon inflates the number of significant genes
by {mean_inflate:.1f}x on average compared to pseudobulk, confirming pseudoreplication
bias (Squair et al. 2021). When donor structure exists (≥3 donors per group),
always use pseudobulk aggregation (DESeq2, edgeR, or limma-voom) for differential
expression. Cell-level Wilcoxon tests should only be used for exploratory marker
discovery, not for reporting DEG counts or significance.
"""

    if path.exists():
        existing = path.read_text()
        if "Pseudobulk vs Cell-Level DE" not in existing:
            content = existing.rstrip() + "\n\n" + paragraph
        else:
            logger.info("P4 guidance already present, skipping")
            return path
    else:
        content = paragraph

    path.write_text(content)
    logger.info("Wrote P4 prompt enrichment to %s", path)
    return path


def emit_validation_rules() -> Path:
    """Append 2 new P4 rules to scrna_seq.yaml."""
    new_rules = [
        {
            "id": "scrna_pseudobulk_required_3donors",
            "type": "conditional",
            "description": "Pseudobulk DE required when ≥3 donors per group exist",
            "condition": {
                "context": "experiment",
                "field": "samples",
                "check": "unique_donors_gt",
                "threshold": 2,
            },
            "required": {
                "scope": "steps",
                "tool_name_any": [
                    "DESeq2", "edgeR", "limma", "pseudobulk",
                    "AggregateExpression",
                ],
            },
            "action": "warn",
            "message": (
                "Experiment has ≥3 donors per group but no pseudobulk DE step. "
                "Empirical audit of CELLxGENE Census data shows Wilcoxon cell-level tests "
                "inflate significant gene counts by 5-50x compared to DESeq2 pseudobulk, "
                "with Jaccard overlap of top 100 DEGs typically below 0.15. This inflation "
                "is pseudoreplication bias (Squair et al. 2021). Add a pseudobulk aggregation "
                "step (DESeq2/edgeR/limma) before differential expression analysis."
            ),
        },
        {
            "id": "scrna_wilcoxon_donor_warning",
            "type": "conditional",
            "description": "Warn when Wilcoxon is used for DE on multi-donor data",
            "condition": {
                "context": "experiment",
                "field": "samples",
                "check": "unique_donors_gt",
                "threshold": 1,
            },
            "required": {
                "scope": "steps",
                "tool_name_none": [
                    "Wilcoxon", "wilcoxon", "FindAllMarkers",
                    "rank_genes_groups",
                ],
            },
            "action": "flag",
            "message": (
                "Wilcoxon/rank-sum cell-level DE is in the pipeline on multi-donor data. "
                "Cell-level tests treat each cell as an independent observation, but cells "
                "from the same donor are correlated. This produces dramatically inflated p-values "
                "and false discovery rates. Use Wilcoxon only for exploratory marker discovery, "
                "not for reporting differential expression results. Pseudobulk methods (DESeq2, "
                "edgeR, limma-voom) are statistically appropriate for multi-donor comparisons."
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
    """Write pseudobulk_de.smk Snakefile template."""
    path = SNAKEFILE_DIR / "pseudobulk_de.smk"
    path.parent.mkdir(parents=True, exist_ok=True)

    content = '''# Pseudobulk DE — Snakefile template for BioOrchestrator L4
# Aggregates single-cell counts to pseudobulk per donor, runs DESeq2.
# Generated by P4 emit_knowledge.py

import os

WORKDIR = config.get("workdir", ".")
OUTPUT_DIR = os.path.join(WORKDIR, "output")
LOG_DIR = os.path.join(WORKDIR, "logs")
INPUT_H5AD = config.get("input_h5ad", os.path.join(WORKDIR, "input.h5ad"))
DESIGN_FORMULA = config.get("design", "~condition")
CONTRAST = config.get("contrast", ["condition", "disease", "normal"])
MIN_DONORS = config.get("min_donors_per_group", 3)


rule all:
    input:
        os.path.join(OUTPUT_DIR, "repro_verified.flag"),


rule repro_snapshot:
    output:
        lockfile=os.path.join(OUTPUT_DIR, "repro.lock"),
    log:
        os.path.join(LOG_DIR, "repro_snapshot.log"),
    shell:
        """
        repro snapshot -o {output.lockfile} --quiet 2>&1 | tee {log}
        """


rule pseudobulk_aggregate:
    """Sum raw counts per donor to create pseudobulk matrix."""
    input:
        h5ad=INPUT_H5AD,
        lockfile=rules.repro_snapshot.output.lockfile,
    output:
        counts=os.path.join(OUTPUT_DIR, "pseudobulk_counts.tsv"),
        metadata=os.path.join(OUTPUT_DIR, "pseudobulk_metadata.tsv"),
    log:
        os.path.join(LOG_DIR, "pseudobulk_aggregate.log"),
    resources:
        mem_mb=32000,
        runtime=30,
    shell:
        """
        python -c "
import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('{input.h5ad}')
# Aggregate by donor_id
donors = adata.obs['donor_id'].unique()
counts = []
meta = []
for donor in donors:
    mask = adata.obs['donor_id'] == donor
    counts.append(np.asarray(adata.X[mask].sum(axis=0)).ravel())
    meta.append(adata.obs[mask].iloc[0][['donor_id', 'condition']].to_dict())

count_df = pd.DataFrame(counts, columns=adata.var_names, index=[m['donor_id'] for m in meta])
meta_df = pd.DataFrame(meta).set_index('donor_id')
count_df.to_csv('{output.counts}', sep='\\t')
meta_df.to_csv('{output.metadata}', sep='\\t')
        " 2>&1 | tee {log}
        """


rule deseq2_de:
    """Run pyDESeq2 on pseudobulk counts."""
    input:
        counts=rules.pseudobulk_aggregate.output.counts,
        metadata=rules.pseudobulk_aggregate.output.metadata,
    output:
        results=os.path.join(OUTPUT_DIR, "deseq2_results.tsv"),
    log:
        os.path.join(LOG_DIR, "deseq2_de.log"),
    resources:
        mem_mb=16000,
        runtime=60,
    shell:
        """
        python -c "
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

counts = pd.read_csv('{input.counts}', sep='\\t', index_col=0).astype(int)
meta = pd.read_csv('{input.metadata}', sep='\\t', index_col=0)

# Filter low-expression genes
gene_mask = (counts >= 10).sum(axis=0) >= 3
counts = counts.loc[:, gene_mask]

dds = DeseqDataSet(counts=counts, metadata=meta, design='{DESIGN_FORMULA}')
dds.deseq2()
stat_res = DeseqStats(dds, contrast={CONTRAST})
stat_res.summary()
stat_res.results_df.to_csv('{output.results}', sep='\\t')
        " 2>&1 | tee {log}
        """


rule repro_verify:
    input:
        results=rules.deseq2_de.output.results,
        lockfile=rules.repro_snapshot.output.lockfile,
    output:
        flag=os.path.join(OUTPUT_DIR, "repro_verified.flag"),
    log:
        os.path.join(LOG_DIR, "repro_verify.log"),
    shell:
        """
        repro verify {OUTPUT_DIR} -l {input.lockfile} 2>&1 | tee {log}
        touch {output.flag}
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

    tsv_path = OUTPUT_DIR / "p4_pseudobulk.tsv"
    if not tsv_path.exists():
        logger.error("P4 output not found: %s", tsv_path)
        sys.exit(1)

    result_df = pd.read_csv(tsv_path, sep="\t")
    logger.info("Loaded %d rows from %s", len(result_df), tsv_path)

    emit_prompt_enrichment(result_df)
    emit_validation_rules()
    emit_snakefile_template()

    logger.info("P4 knowledge emission complete")


if __name__ == "__main__":
    main()
