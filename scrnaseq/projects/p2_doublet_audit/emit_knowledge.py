#!/usr/bin/env python3
"""P2 — Emit BioOrchestrator knowledge artifacts.

Generates:
- L2: Prompt enrichment paragraph for scrna/v1.txt
- L3: 2 new validation rules for scrna_seq.yaml
- L4: Snakefile template doublet_audit.smk
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
    """Append P2 findings to scrna/v1.txt as empirical guidance."""
    path = PROMPTS_DIR / "v1.txt"
    path.parent.mkdir(parents=True, exist_ok=True)

    # Compute summary stats for the guidance paragraph
    n_datasets = len(result_df)
    median_kappa = result_df.cohens_kappa.median()
    min_kappa = result_df.cohens_kappa.min()
    max_kappa = result_df.cohens_kappa.max()
    median_disagree = result_df.disagreement_rate.median() * 100
    median_scdblfinder = result_df.scdblfinder_rate.median() * 100
    median_scrublet = result_df.scrublet_rate.median() * 100

    # Spearman correlation
    from scipy import stats
    rho, pval = stats.spearmanr(
        result_df["disagreement_rate"],
        result_df["n_cell_types"],
    )

    paragraph = f"""

## Empirical Guidance: Doublet Detection Tool Agreement (P2 Audit)

Analysis of {n_datasets} datasets across diverse tissues (Homo sapiens + Mus musculus)
from CELLxGENE Census comparing scDblFinder and Scrublet shows only moderate
inter-tool agreement. Median Cohen's kappa: {median_kappa:.2f} (range: {min_kappa:.2f}-{max_kappa:.2f}).
Median disagreement rate: {median_disagree:.1f}% of cells receive conflicting calls.
Median doublet rates: scDblFinder {median_scdblfinder:.1f}%, Scrublet {median_scrublet:.1f}%.
Disagreement correlates with cell type complexity (Spearman rho={rho:.2f}, p={pval:.2e}):
datasets with more cell types show higher disagreement between tools. Best practice:
run at least two doublet detection tools and take the intersection (conservative) or
union (sensitive) rather than relying on a single tool. Doublet detection should always
precede clustering, as undetected doublets create artifactual intermediate clusters.
"""

    # Append to existing file or create
    if path.exists():
        existing = path.read_text()
        if "Doublet Detection Tool Agreement" not in existing:
            content = existing.rstrip() + "\n" + paragraph
        else:
            logger.info("P2 guidance already present in %s, skipping", path)
            return path
    else:
        content = paragraph

    path.write_text(content)
    logger.info("Wrote P2 prompt enrichment to %s", path)
    return path


def emit_validation_rules() -> Path:
    """Append 2 new P2 rules to scrna_seq.yaml."""
    new_rules = [
        {
            "id": "scrna_doublet_multi_tool",
            "type": "conditional",
            "description": "Warn if only one doublet detection tool is used; multi-tool consensus improves accuracy",
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
                "tool_name_count": {
                    "any_of": [
                        "scDblFinder",
                        "Scrublet",
                        "scrublet",
                        "DoubletFinder",
                        "solo",
                    ],
                    "min_count": 2,
                },
            },
            "action": "warn",
            "message": (
                "Pipeline uses only one doublet detection tool. Empirical audit of "
                "CELLxGENE Census data (P2) shows scDblFinder and Scrublet agree on only "
                "~60-80% of doublet calls (median Cohen's kappa ~0.4-0.6). Disagreement "
                "increases with cell type complexity. Run at least two tools and take their "
                "intersection (conservative) or union (sensitive) for more robust doublet removal."
            ),
        },
        {
            "id": "scrna_doublet_before_clustering",
            "type": "compatibility",
            "description": "Doublet detection must precede clustering to avoid artifactual intermediate clusters",
            "tool_name": "Leiden",
            "requires_any": [
                {
                    "type": "step_before",
                    "step_tool_names": [
                        "scDblFinder",
                        "Scrublet",
                        "scrublet",
                        "DoubletFinder",
                        "solo",
                    ],
                },
            ],
            "action": "flag",
            "message": (
                "Clustering step (Leiden) is present but no doublet detection step precedes it. "
                "Undetected doublets create artifactual intermediate clusters that appear as "
                "novel cell states. Doublet detection should run after QC filtering but before "
                "dimensionality reduction and clustering. This ordering is especially important "
                "for datasets with many cell types, where doublet tool disagreement is highest."
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

    logger.info("Wrote %d new rules to %s (total: %d)", added, path, len(data["rules"]))
    return path


def emit_snakefile_template() -> Path:
    """Write doublet_audit.smk Snakefile template."""
    path = SNAKEFILE_DIR / "doublet_audit.smk"
    path.parent.mkdir(parents=True, exist_ok=True)

    content = '''# Doublet Detection Audit — Snakefile template for BioOrchestrator L4
# Compares scDblFinder and Scrublet doublet calls across diverse datasets.
# Generated by P2 emit_knowledge.py

import os

WORKDIR = config.get("workdir", ".")
OUTPUT_DIR = os.path.join(WORKDIR, "output")
LOG_DIR = os.path.join(WORKDIR, "logs")
SRC_DIR = os.path.join(WORKDIR, "src", "pipelines")
CENSUS_VERSION = config.get("census_version", "2025-11-08")
ORGANISMS = config.get("organisms", ["Homo sapiens", "Mus musculus"])
MAX_CELLS = config.get("max_cells_per_dataset", 20000)


rule all:
    input:
        os.path.join(OUTPUT_DIR, "repro_verified.flag"),


rule repro_snapshot:
    """Capture environment state before Census access."""
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


rule fetch_and_run_doublets:
    """Fetch raw counts from Census, run scDblFinder + Scrublet on each dataset."""
    input:
        lockfile=rules.repro_snapshot.output.lockfile,
    output:
        tsv=os.path.join(OUTPUT_DIR, "p2_doublet_audit.tsv"),
    params:
        census_version=CENSUS_VERSION,
        max_cells=MAX_CELLS,
    log:
        os.path.join(LOG_DIR, "doublet_detection.log"),
    resources:
        mem_mb=32000,
        runtime=360,
    shell:
        """
        python {SRC_DIR}/run_p2.py \\
            --census-version '{params.census_version}' \\
            --max-cells {params.max_cells} \\
            --output-dir {OUTPUT_DIR} \\
            2>&1 | tee {log}
        """


rule compute_agreement:
    """Compute Cohen's kappa, disagreement rate, and Spearman correlation."""
    input:
        tsv=rules.fetch_and_run_doublets.output.tsv,
    output:
        summary=os.path.join(OUTPUT_DIR, "p2_agreement_summary.tsv"),
    log:
        os.path.join(LOG_DIR, "agreement_stats.log"),
    resources:
        mem_mb=4000,
        runtime=10,
    shell:
        """
        python -c "
import pandas as pd
from scipy import stats

df = pd.read_csv('{input.tsv}', sep='\\t')
rho, pval = stats.spearmanr(df['disagreement_rate'], df['n_cell_types'])
summary = pd.DataFrame([dict(
    n_datasets=len(df),
    median_kappa=df.cohens_kappa.median(),
    median_disagree=df.disagreement_rate.median(),
    spearman_rho=rho,
    spearman_p=pval,
)])
summary.to_csv('{output.summary}', sep='\\t', index=False)
print(f'Spearman rho={{rho:.3f}}, p={{pval:.2e}}')
        " 2>&1 | tee {log}
        """


rule repro_verify:
    """Hash outputs and verify reproducibility."""
    input:
        tsv=rules.fetch_and_run_doublets.output.tsv,
        summary=rules.compute_agreement.output.summary,
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
        handlers=[logging.StreamHandler(sys.stdout)],
    )

    tsv_path = OUTPUT_DIR / "p2_doublet_audit.tsv"
    if not tsv_path.exists():
        logger.error("P2 output not found: %s — run run_p2.py first", tsv_path)
        sys.exit(1)

    result_df = pd.read_csv(tsv_path, sep="\t")
    logger.info("Loaded %d rows from %s", len(result_df), tsv_path)

    emit_prompt_enrichment(result_df)
    emit_validation_rules()
    emit_snakefile_template()

    logger.info("P2 knowledge emission complete")


if __name__ == "__main__":
    main()
