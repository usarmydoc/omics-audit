#!/usr/bin/env python3
"""P9 — Emit BioOrchestrator knowledge artifacts for clustering resolution audit."""

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
    """L2: Append paragraph about resolution selection to v1.txt."""
    path = PROMPTS_DIR / "v1.txt"

    n_datasets = result_df.groupby(
        ["dataset_id", "tissue", "organism"], observed=True,
    ).ngroups
    optimal = result_df[result_df.optimal_flag].copy()
    mean_opt = optimal.resolution.mean()
    median_opt = optimal.resolution.median()
    pct_05 = (optimal.resolution == 0.5).sum() / len(optimal) * 100 if len(optimal) > 0 else 0

    # ARI at res=0.5 vs best
    r05 = result_df[result_df.resolution == 0.5]
    mean_ari_05 = r05.ari.mean() if len(r05) > 0 else 0
    mean_ari_best = optimal.ari.mean() if len(optimal) > 0 else 0

    # Correlation between n_cell_types and optimal resolution
    from scipy import stats as scipy_stats
    if len(optimal) >= 5:
        rho, pval = scipy_stats.spearmanr(optimal.n_true_types, optimal.resolution)
        corr_text = (
            f"Spearman correlation between number of cell types and optimal "
            f"resolution is rho={rho:.2f} (p={pval:.2e}), "
        )
        if abs(rho) > 0.3:
            corr_text += (
                "confirming that datasets with more cell types benefit from "
                "higher resolution."
            )
        else:
            corr_text += (
                "suggesting resolution choice depends on factors beyond "
                "cell type count alone."
            )
    else:
        corr_text = ""

    paragraph = f"""
## Empirical Guidance: Clustering Resolution Selection (P9 Audit)

Analysis of {n_datasets} datasets from CELLxGENE Census shows that the default
Leiden resolution of 0.5 is optimal in only {pct_05:.0f}% of cases. The mean
optimal resolution is {mean_opt:.2f} (median {median_opt:.2f}). Mean ARI at
resolution 0.5 is {mean_ari_05:.3f} vs {mean_ari_best:.3f} at the dataset-specific
optimal. {corr_text} When designing scRNA-seq clustering pipelines, do NOT
blindly accept resolution=0.5. For datasets with few cell types (<5), lower
resolutions (0.1-0.3) often perform best. For complex tissues with many cell
types (>10), higher resolutions (0.8-1.5) recover more structure. Always
evaluate clustering quality against known markers or reference annotations,
and explicitly justify the chosen resolution.
"""

    if path.exists():
        existing = path.read_text()
        if "Clustering Resolution Selection" not in existing:
            content = existing.rstrip() + "\n\n" + paragraph
        else:
            logger.info("P9 prompt enrichment already present, skipping")
            return path
    else:
        content = paragraph
    path.write_text(content)
    logger.info("Wrote P9 prompt enrichment to %s", path)
    return path


def emit_validation_rules() -> Path:
    """L3: Add two validation rules for clustering resolution."""
    new_rules = [
        {
            "id": "scrna_resolution_explicit",
            "type": "conditional",
            "description": "Clustering resolution must be explicitly set, not left to tool default",
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
                    "qc_thresholds.clustering_resolution",
                    "clustering.resolution",
                    "leiden.resolution",
                ],
            },
            "action": "warn",
            "message": (
                "No explicit clustering resolution parameter found. The default "
                "resolution of 0.5 is optimal in only a minority of datasets. "
                "Empirical audit across CELLxGENE Census shows that optimal resolution "
                "varies from 0.1 to 1.5+ depending on tissue complexity. Datasets with "
                "few cell types (<5) need lower resolution (0.1-0.3); complex tissues "
                "(>10 types) need higher resolution (0.8-1.5). Always set resolution "
                "explicitly and justify the choice."
            ),
        },
        {
            "id": "scrna_resolution_default_05",
            "type": "range",
            "description": "Flag when clustering resolution is set to exactly 0.5, suggesting an unjustified default",
            "parameter": "qc_thresholds.clustering_resolution",
            "min": 0.49,
            "max": 0.51,
            "action": "flag",
            "invert": True,
            "message": (
                "clustering_resolution is set to exactly 0.5, which is the most common "
                "default but is rarely the optimal resolution for any specific dataset. "
                "Empirical audit shows the optimal resolution depends on the number of "
                "distinct cell types in the tissue: simple tissues need ~0.1-0.3, complex "
                "tissues need ~0.8-1.5. A resolution of 0.5 often over-clusters simple "
                "datasets and under-clusters complex ones. Consider evaluating multiple "
                "resolutions against known markers or reference annotations."
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
    """L4: Write clustering_resolution_audit.smk with repro rules."""
    path = SNAKEFILE_DIR / "clustering_resolution_audit.smk"
    path.parent.mkdir(parents=True, exist_ok=True)

    content = '''# Clustering Resolution Sensitivity Audit — Snakefile template for BioOrchestrator L4
# Sweeps Leiden resolutions and evaluates ARI vs ground truth cell type labels.
# Generated by P9 emit_knowledge.py

import os

WORKDIR = config.get("workdir", ".")
OUTPUT_DIR = os.path.join(WORKDIR, "output")
LOG_DIR = os.path.join(WORKDIR, "logs")
INPUT_H5AD = config.get("input_h5ad", os.path.join(WORKDIR, "input.h5ad"))
RESOLUTIONS = config.get("resolutions", [0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0])
CELL_TYPE_COL = config.get("cell_type_col", "cell_type")
N_PCS = config.get("n_pcs", 50)


rule all:
    input:
        os.path.join(OUTPUT_DIR, "repro_verified.flag"),


rule repro_snapshot:
    """Capture environment state before analysis."""
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


rule resolution_sweep:
    """Run Leiden at multiple resolutions and compute ARI vs ground truth."""
    input:
        h5ad=INPUT_H5AD,
        lockfile=rules.repro_snapshot.output.lockfile,
    output:
        tsv=os.path.join(OUTPUT_DIR, "clustering_resolution_results.tsv"),
    params:
        resolutions=RESOLUTIONS,
        cell_type_col=CELL_TYPE_COL,
        n_pcs=N_PCS,
    log:
        os.path.join(LOG_DIR, "resolution_sweep.log"),
    resources:
        mem_mb=32000,
        runtime=60,
    shell:
        """
        python -c "
import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

adata = sc.read_h5ad('{input.h5ad}')
cell_type_col = '{params.cell_type_col}'
n_pcs = {params.n_pcs}
resolutions = {params.resolutions}

# Preprocess
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
n_hvg = min(2000, adata.n_vars - 1)
sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)
adata = adata[:, adata.var.highly_variable].copy()
n_pcs = min(n_pcs, adata.n_vars - 1, adata.n_obs - 1)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=n_pcs, random_state=42)
sc.pp.neighbors(adata, n_pcs=n_pcs, random_state=42)

true_labels = adata.obs[cell_type_col].values
n_true = len(np.unique(true_labels))
results = []
for res in resolutions:
    sc.tl.leiden(adata, resolution=res, random_state=42, key_added='leiden')
    pred = adata.obs['leiden'].values
    n_clust = len(np.unique(pred))
    ari = adjusted_rand_score(true_labels, pred)
    over = n_clust / n_true if n_clust > n_true else 1.0
    under = n_true / n_clust if n_clust < n_true else 1.0
    results.append(dict(resolution=res, ari=round(ari, 4),
                        n_clusters=n_clust, n_true_types=n_true,
                        over_clustering=round(over, 3),
                        under_clustering=round(under, 3)))

df = pd.DataFrame(results)
best_ari = df.ari.max()
df['optimal_flag'] = df.ari == best_ari
df.to_csv('{output.tsv}', sep='\\t', index=False)
print(f'Best resolution: {{df.loc[df.ari.idxmax(), \"resolution\"]}} (ARI={{best_ari:.4f}})')
        " 2>&1 | tee {log}
        """


rule repro_verify:
    """Hash outputs and verify reproducibility."""
    input:
        tsv=rules.resolution_sweep.output.tsv,
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

    tsv_path = OUTPUT_DIR / "p9_clustering_resolution.tsv"
    if not tsv_path.exists():
        logger.error("P9 output not found: %s", tsv_path)
        sys.exit(1)

    result_df = pd.read_csv(tsv_path, sep="\t")
    logger.info("Loaded %d rows from %s", len(result_df), tsv_path)

    emit_prompt_enrichment(result_df)
    emit_validation_rules()
    emit_snakefile_template()
    logger.info("P9 knowledge emission complete")


if __name__ == "__main__":
    main()
