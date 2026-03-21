"""
BioOrchestrator knowledge artifact emitter for bulk RNA-seq audit.

Handles L2 (prompt enrichment), L3 (validation rules), and L4 (Snakefile
templates). All BioOrchestrator integration goes through this module —
analysis code must never write to BioOrchestrator paths directly.
"""

import os
from pathlib import Path
import yaml

BIOORCHESTRATOR_ROOT = Path(os.environ.get("BIOORCHESTRATOR_ROOT", str(Path.home() / "bioorchestrator" / "src" / "bioorchestrator")))
PROMPTS_DIR = BIOORCHESTRATOR_ROOT / "knowledge" / "prompts" / "bulk"
RULES_FILE = BIOORCHESTRATOR_ROOT / "knowledge" / "rules" / "bulk_rnaseq.yaml"
SNAKEFILE_DIR = BIOORCHESTRATOR_ROOT / "knowledge" / "snakefile_templates"


def _load_existing_rules() -> dict:
    """Load current bulk_rnaseq.yaml rules."""
    with open(RULES_FILE) as f:
        data = yaml.safe_load(f)
    return data if data else {"rules": []}


def _existing_rule_ids() -> set:
    data = _load_existing_rules()
    return {r["id"] for r in data.get("rules", [])}


def append_prompt_enrichment(text: str, source_project: str):
    """Append a paragraph to the bulk v1.txt prompt file (L2)."""
    prompt_file = PROMPTS_DIR / "v1.txt"
    content = prompt_file.read_text()
    marker = f"# --- {source_project} audit findings ---"
    if marker in content:
        return  # already appended
    block = f"\n\n{marker}\n{text}\n"
    prompt_file.write_text(content + block)
    print(f"  L2: appended {source_project} findings to {prompt_file}")


def append_validation_rules(rules: list[dict], source_project: str):
    """Append new validation rules to bulk_rnaseq.yaml (L3).

    Each rule dict must have: id, type, description, action, message,
    plus type-specific fields (parameter/min/max for range,
    condition/required for conditional).

    Skips rules whose id already exists.
    """
    data = _load_existing_rules()
    existing_ids = {r["id"] for r in data.get("rules", [])}
    added = 0
    for rule in rules:
        if rule["id"] not in existing_ids:
            data["rules"].append(rule)
            added += 1
    if added > 0:
        header = "# L3 validation rules for bulk RNA-seq pipelines.\n\n"
        with open(RULES_FILE, "w") as f:
            f.write(header)
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)
        print(f"  L3: added {added} rules from {source_project} to {RULES_FILE}")
    else:
        print(f"  L3: all {source_project} rules already present, skipping")


def write_snakefile_template(content: str, filename: str, source_project: str):
    """Write a Snakefile template to L4 (snakefile_templates/)."""
    outpath = SNAKEFILE_DIR / filename
    outpath.write_text(content)
    print(f"  L4: wrote {outpath} from {source_project}")
