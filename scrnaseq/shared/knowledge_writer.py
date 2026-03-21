"""BioOrchestrator knowledge artifact writer.

Generates:
1. Prompt enrichment file (marker_validation_context.txt)
2. Validation rules (appended to scrna_seq.yaml)
3. Snakefile template (marker_atlas.smk)

Isolated from pipeline logic — takes aggregated scores as input.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

import yaml

from .scoring import AggregatedScore

logger = logging.getLogger(__name__)

BIOORCHESTRATOR_ROOT = Path(os.environ.get("BIOORCHESTRATOR_ROOT", str(Path.home() / "bioorchestrator" / "src" / "bioorchestrator")))
PROMPTS_DIR = BIOORCHESTRATOR_ROOT / "knowledge" / "prompts" / "scrna"
RULES_FILE = BIOORCHESTRATOR_ROOT / "knowledge" / "rules" / "scrna_seq.yaml"
SNAKEFILE_DIR = BIOORCHESTRATOR_ROOT / "knowledge" / "snakefile_templates"


def write_prompt_enrichment(
    aggregated: list[AggregatedScore],
    output_path: Path | None = None,
) -> Path:
    """Generate marker_validation_context.txt for BioOrchestrator LLM prompts.

    Structure:
    - Header
    - HIGH RELIABILITY markers per cell type
    - LOW RELIABILITY markers with caveats
    - TISSUE-SPECIFIC NOTES
    - GUIDANCE for LLM pipeline design
    """
    path = output_path or (PROMPTS_DIR / "marker_validation_context.txt")
    path.parent.mkdir(parents=True, exist_ok=True)

    high = [s for s in aggregated if s.mean_auc >= 0.80]
    low = [s for s in aggregated if s.mean_auc < 0.65]
    divergent = [s for s in aggregated if s.sd_auc > 0.10]

    lines = [
        "# Marker Validation Context — auto-generated from CELLxGENE Census",
        "# Append to scrna/v1.txt at prompt load time.",
        "",
        "## Marker Reliability Reference",
        "",
        "The following marker reliability scores were computed across multiple",
        "CELLxGENE Census datasets using ROC AUC, Mann-Whitney U tests, and",
        "bootstrap confidence intervals. Use these to inform marker selection",
        "in cell type annotation and QC steps.",
        "",
        "### High-Reliability Markers (mean AUC >= 0.80)",
        "",
        "These markers reliably distinguish their target cell type across datasets.",
        "Prefer these when designing marker-dependent pipeline steps:",
        "",
    ]

    # Group high by organism then cell type
    by_org_ct: dict[str, dict[str, list[AggregatedScore]]] = {}
    for s in high:
        by_org_ct.setdefault(s.organism, {}).setdefault(s.cell_type, []).append(s)

    for org in sorted(by_org_ct.keys()):
        lines.append(f"**{org}:**")
        for ct in sorted(by_org_ct[org].keys()):
            markers = sorted(by_org_ct[org][ct], key=lambda x: x.mean_auc, reverse=True)
            marker_strs = [f"{s.gene} (AUC={s.mean_auc:.2f}, N={s.n_datasets})" for s in markers]
            lines.append(f"- {ct}: {', '.join(marker_strs)}")
        lines.append("")

    lines.extend(["### Low-Reliability Markers (mean AUC < 0.65)", ""])
    lines.append("These markers show poor or inconsistent discrimination. Use with caution:")
    lines.append("")

    by_org_ct_low: dict[str, dict[str, list[AggregatedScore]]] = {}
    for s in low:
        by_org_ct_low.setdefault(s.organism, {}).setdefault(s.cell_type, []).append(s)

    for org in sorted(by_org_ct_low.keys()):
        lines.append(f"**{org}:**")
        for ct in sorted(by_org_ct_low[org].keys()):
            markers = sorted(by_org_ct_low[org][ct], key=lambda x: x.mean_auc)
            for s in markers:
                caveat = ""
                if s.sd_auc > 0.10:
                    caveat = " — highly variable across datasets"
                elif s.n_datasets < 3:
                    caveat = " — insufficient datasets for reliable estimate"
                lines.append(f"- {ct}: {s.gene} (AUC={s.mean_auc:.2f}{caveat})")
        lines.append("")

    lines.extend(["### Tissue-Specific Divergence Notes", ""])

    if divergent:
        lines.append("These markers behave inconsistently across tissues/datasets:")
        lines.append("")
        for s in sorted(divergent, key=lambda x: x.sd_auc, reverse=True):
            lines.append(
                f"- {s.gene} for {s.cell_type} ({s.organism}): sd_AUC={s.sd_auc:.2f} across "
                f"{s.n_datasets} datasets. Consider tissue-specific thresholds."
            )
    else:
        lines.append("No significant tissue-specific divergence detected.")

    lines.extend([
        "",
        "### Guidance for Pipeline Design",
        "",
        "When designing marker-dependent pipeline steps:",
        "1. Always prefer high-reliability markers for automated cell type annotation.",
        "2. If using low-reliability markers, add a QC warning step to flag uncertain annotations.",
        "3. Combine multiple markers per cell type — single markers are rarely sufficient.",
        "4. For tissue-specific analyses, verify marker performance in the target tissue.",
        "5. Require validation across at least 3 independent datasets before trusting a marker.",
        "6. Include a marker scoring/validation step before any cell type annotation step.",
        "7. When querying CELLxGENE Census, always include repro snapshot/verify for reproducibility.",
        "",
    ])

    content = "\n".join(lines)
    path.write_text(content)
    logger.info("Wrote prompt enrichment: %s", path)
    return path


def append_validation_rules(rules_path: Path | None = None) -> Path:
    """Append 5 new validation rules to scrna_seq.yaml.

    Uses existing BioOrchestrator rule schema (conditional, compatibility types).
    """
    path = rules_path or RULES_FILE

    new_rules = [
        {
            "id": "scrna_marker_min_datasets",
            "type": "conditional",
            "description": "Marker validation steps must use data from at least 3 independent datasets",
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
                "tool_name_any": ["marker_validation", "marker_scoring", "CELLxGENE_Census"],
            },
            "action": "warn",
            "message": (
                "Marker validation steps should use data from at least 3 independent datasets "
                "(n_datasets >= 3). Single-dataset marker scores are unreliable due to batch effects "
                "and annotation inconsistencies. Include a marker_validation or marker_scoring step "
                "that aggregates across multiple CELLxGENE Census datasets."
            ),
        },
        {
            "id": "scrna_marker_low_reliability_qc",
            "type": "conditional",
            "description": "Pipelines using low-reliability markers (AUC < 0.65) must include a QC warning step",
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
                "tool_name_any": ["marker_qc", "marker_validation", "QC_warning"],
            },
            "action": "warn",
            "message": (
                "Pipeline uses cell type markers but no marker QC or validation step is present. "
                "Low-reliability markers (AUC < 0.65) — such as CD4 for T cells or MARCO for macrophages — "
                "can produce unreliable annotations. Add a marker_validation or marker_qc step to flag "
                "uncertain cell type assignments."
            ),
        },
        {
            "id": "scrna_census_organism",
            "type": "compatibility",
            "description": "CELLxGENE Census data source requires organism to be human or mouse",
            "tool_name": "CELLxGENE_Census",
            "requires_any": [
                {
                    "type": "input_format",
                    "formats": ["human", "mouse", "h5ad", "cellxgene"],
                },
            ],
            "action": "reject",
            "message": (
                "CELLxGENE Census is specified as a data source but the organism must be human "
                "(Homo sapiens) or mouse (Mus musculus). Census does not contain data for other organisms. "
                "Use a different data source for non-human, non-mouse experiments."
            ),
        },
        {
            "id": "scrna_marker_before_annotation",
            "type": "compatibility",
            "description": "Marker scoring steps must precede any cell type annotation step",
            "tool_name": "cell_type_annotation",
            "requires_any": [
                {
                    "type": "step_before",
                    "step_tool_names": [
                        "marker_scoring", "marker_validation", "CELLxGENE_Census",
                        "SingleR", "CellTypist", "scType",
                    ],
                },
            ],
            "action": "flag",
            "message": (
                "Cell type annotation step found but no marker scoring or validation step precedes it. "
                "Running marker reliability assessment before annotation improves confidence in results "
                "and identifies potentially unreliable markers. Add a marker_scoring step before annotation."
            ),
        },
        {
            "id": "scrna_census_repro",
            "type": "conditional",
            "description": "Pipelines using CELLxGENE Census must include repro snapshot and verify steps",
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
                "tool_name_any": ["repro_snapshot", "repro_verify", "repro"],
            },
            "action": "warn",
            "message": (
                "Pipeline accesses external data (e.g., CELLxGENE Census) but no reproducibility "
                "tracking step is present. Add a repro_snapshot step before data access and a "
                "repro_verify step after all outputs are written to ensure reproducible results. "
                "Census data is versioned (LTS builds), but environment and output hashing via "
                "repro-lock provides end-to-end reproducibility proof."
            ),
        },
    ]

    # Read existing YAML
    with open(path) as f:
        data = yaml.safe_load(f)

    existing_ids = {r["id"] for r in data.get("rules", [])}

    # Append only new rules
    added = 0
    for rule in new_rules:
        if rule["id"] not in existing_ids:
            data["rules"].append(rule)
            added += 1
            logger.info("Added rule: %s", rule["id"])
        else:
            logger.info("Rule %s already exists, skipping", rule["id"])

    # Write back
    with open(path, "w") as f:
        # Write header comment
        f.write("# L3 validation rules for single-cell RNA-seq pipelines.\n")
        f.write("# Each rule is evaluated deterministically against a PipelineConfig + ExperimentContext.\n")
        f.write("# Rule types: range, conditional, compatibility\n")
        f.write("# Actions: reject (hard fail), warn (proceed with warning), flag (informational)\n")
        f.write("\n")
        yaml.dump(data, f, default_flow_style=False, sort_keys=False, width=120, allow_unicode=True)

    logger.info("Wrote %d new rules to %s (total: %d)", added, path, len(data["rules"]))
    return path


def write_snakefile_template(output_path: Path | None = None) -> Path:
    """Write marker_atlas.smk Snakefile template with real shell commands."""
    path = output_path or (SNAKEFILE_DIR / "marker_atlas.smk")
    path.parent.mkdir(parents=True, exist_ok=True)

    content = '''# Marker Reliability Atlas — Snakefile template for BioOrchestrator L4
# Real shell commands — no touch {output} placeholders.
# Generated by knowledge_writer.py

import os

WORKDIR = config.get("workdir", ".")
OUTPUT_DIR = os.path.join(WORKDIR, "output")
LOG_DIR = os.path.join(WORKDIR, "logs")
SRC_DIR = os.path.join(WORKDIR, "src", "pipelines")

ORGANISM = config.get("organism", "Homo sapiens")
MAX_WORKERS = config.get("max_workers", 8)
CENSUS_VERSION = config.get("census_version", "2025-11-08")


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


rule census_query:
    """Connect to CELLxGENE Census and fetch expression data for marker genes."""
    input:
        lockfile=rules.repro_snapshot.output.lockfile,
    output:
        expression_cache=os.path.join(OUTPUT_DIR, "expression_cache.npz"),
        dataset_metadata=os.path.join(OUTPUT_DIR, "dataset_metadata.tsv"),
    params:
        organism=ORGANISM,
        census_version=CENSUS_VERSION,
    log:
        os.path.join(LOG_DIR, "census_query.log"),
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        """
        python {SRC_DIR}/census_access.py \\
            --organism '{params.organism}' \\
            --census-version '{params.census_version}' \\
            --output-dir {OUTPUT_DIR} \\
            2>&1 | tee {log}
        """


rule marker_scoring:
    """Compute per-marker, per-cell-type, per-dataset reliability metrics."""
    input:
        expression_cache=rules.census_query.output.expression_cache,
        dataset_metadata=rules.census_query.output.dataset_metadata,
    output:
        raw_scores=os.path.join(OUTPUT_DIR, "raw_marker_scores.tsv"),
    params:
        max_workers=MAX_WORKERS,
    log:
        os.path.join(LOG_DIR, "marker_scoring.log"),
    resources:
        mem_mb=16000,
        runtime=60,
    threads: MAX_WORKERS
    shell:
        """
        python {SRC_DIR}/scoring.py \\
            --input {input.expression_cache} \\
            --metadata {input.dataset_metadata} \\
            --output {output.raw_scores} \\
            --max-workers {params.max_workers} \\
            2>&1 | tee {log}
        """


rule generate_report:
    """Generate TSV atlas and markdown summary report."""
    input:
        raw_scores=rules.marker_scoring.output.raw_scores,
    output:
        atlas_tsv=os.path.join(OUTPUT_DIR, "marker_reliability_atlas.tsv"),
        summary_md=os.path.join(OUTPUT_DIR, "summary_report.md"),
    log:
        os.path.join(LOG_DIR, "generate_report.log"),
    resources:
        mem_mb=4000,
        runtime=10,
    shell:
        """
        python {SRC_DIR}/reporting.py \\
            --input {input.raw_scores} \\
            --output-dir {OUTPUT_DIR} \\
            2>&1 | tee {log}
        """


rule write_knowledge:
    """Write BioOrchestrator knowledge artifacts (prompts, rules, Snakefile)."""
    input:
        atlas_tsv=rules.generate_report.output.atlas_tsv,
    output:
        prompt=os.path.join(OUTPUT_DIR, "marker_validation_context.txt"),
        rules_flag=os.path.join(OUTPUT_DIR, "rules_updated.flag"),
    log:
        os.path.join(LOG_DIR, "write_knowledge.log"),
    resources:
        mem_mb=2000,
        runtime=5,
    shell:
        """
        python {SRC_DIR}/knowledge_writer.py \\
            --atlas {input.atlas_tsv} \\
            --output-dir {OUTPUT_DIR} \\
            2>&1 | tee {log}
        """


rule repro_verify:
    """Hash all output files and verify reproducibility."""
    input:
        atlas_tsv=rules.generate_report.output.atlas_tsv,
        summary_md=rules.generate_report.output.summary_md,
        prompt=rules.write_knowledge.output.prompt,
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
