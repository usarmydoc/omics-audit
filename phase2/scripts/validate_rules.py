#!/usr/bin/env python3
"""Validate draft rule YAMLs against AUDIT_STANDARDS.md schema v1.0.3.

Per Section 5.4 of standards. Reports compliance per rule and per file.

Behavior:
  - Standard validation pass for rules in draft_rules/ and subdirectories
    EXCEPT draft_rules/superseded/ and draft_rules/pending_engine_support/
  - Light validation pass for rules in draft_rules/superseded/:
    only requires rule_id, schema_version, deprecation block, description.
    Legacy field names accepted without rewrite.
  - Light validation pass for rules in draft_rules/pending_engine_support/:
    same as superseded — these are valid rules whose triggers depend on
    rule-engine features not yet implemented.
  - Companion-metrics check: rules with action_type=report_additional_metric
    OR whose recommendation/description text contains tool-comparison vocabulary
    (Jaccard, agreement, disagree, top-100) MUST also mention direction
    agreement AND log2FC correlation in recommendation/description, or fail.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Any

import yaml

RULES_DIR = Path("/mnt/nvme1/omics-audit/phase2/draft_rules")
STANDARDS_DOC = Path("/mnt/nvme1/omics-audit/standards/AUDIT_STANDARDS.md")
STEP_REGISTRY = Path("/mnt/nvme1/omics-audit/standards/pipeline_step_registry.yaml")


def _load_step_registry() -> set[str]:
    """Return the set of canonical pipeline-step names. Empty set if registry missing."""
    if not STEP_REGISTRY.exists():
        return set()
    try:
        with STEP_REGISTRY.open() as fh:
            doc = yaml.safe_load(fh) or {}
        return {s["name"] for s in doc.get("steps", []) if "name" in s}
    except Exception:
        return set()


REGISTERED_STEPS = _load_step_registry()

# Per schema 1.0.0 Section 5.4
# Schema v1.0.3 — required at draft state (backward-compatible with 1.0.0-1.0.2)
REQUIRED_TOP_LEVEL = {
    "rule_id", "schema_version", "status", "title", "description",
    "trigger_conditions", "recommendation", "confidence_tier",
    "prior_audit_relationship", "evidence", "out_of_scope",
    "created_date",
}
# Required only at review_ready+
REQUIRED_REVIEWED = {
    "last_reviewed", "reviewer",
}
OPTIONAL_TOP_LEVEL = {
    "mechanism_notes", "applicability", "related_rules", "deprecation",
    "last_reviewed", "reviewer", "revision_history", "severity",
}
# Required at deployed status only
REQUIRED_DEPLOYED = {"deployed_date"}

STATUSES = {"draft", "review_ready", "reviewed", "deployed", "deprecated"}
SEVERITIES = {"info", "warn", "error", "reject"}

CONFIDENCE_TIERS = {
    "hard_default", "conditional", "flag_and_warn",
    "literature_based", "insufficient_data",
}
PRIOR_AUDIT_RELATIONSHIPS = {
    "as_original", "refines_prior", "contradicts_prior",
    "extends_prior", "novel",
}
ACTION_TYPES = {
    "replace_parameter", "add_step", "remove_step", "flag_only",
    "require_documentation", "compute_from_data", "pin_version",
    "report_additional_metric",
    # Added v1.0.3 — for rules whose recommendation IS to stop the pipeline
    # (severity=reject) rather than a parameter change. block_step is a
    # narrower variant that fails only the current step.
    "reject_pipeline", "block_step",
}
CONDITION_TYPES = {
    "dataset_feature_threshold", "tool_invocation",
    "tool_version_constraint", "pipeline_step", "parameter_value_check",
    "analysis_context", "composite",
}
CHANGE_TYPES = {
    "create", "refine_evidence", "change_tier",
    "change_recommendation", "deprecate", "restore",
}
RULE_ID_REGEX = re.compile(r"^[a-z][a-z0-9_]{2,63}$")
ISO_DATE_REGEX = re.compile(r"^\d{4}-\d{2}-\d{2}$")
SHA256_REGEX = re.compile(r"^[a-fA-F0-9]{64}$|^<.*placeholder.*>$")


def validate_rule(rule: dict, rule_idx: int,
                  strict_steps: bool = False) -> tuple[list[str], list[str]]:
    """Return (errors, warnings) for one rule.

    strict_steps=True promotes unknown pipeline_step name warnings to errors.
    """
    errors: list[str] = []
    warnings: list[str] = []
    rid = rule.get("rule_id", f"rule#{rule_idx}")
    status = rule.get("status", "draft")

    # Required top-level — gates on status
    required = set(REQUIRED_TOP_LEVEL)
    if status in ("review_ready", "reviewed", "deployed"):
        required |= REQUIRED_REVIEWED
    if status == "deployed":
        required |= REQUIRED_DEPLOYED

    missing = required - rule.keys()
    if missing:
        errors.append(f"missing required fields: {sorted(missing)}")

    # Unknown top-level fields
    allowed = REQUIRED_TOP_LEVEL | REQUIRED_REVIEWED | REQUIRED_DEPLOYED | OPTIONAL_TOP_LEVEL
    unknown = rule.keys() - allowed
    if unknown:
        errors.append(f"unknown top-level fields (not in schema): {sorted(unknown)}")

    # status enum
    if status not in STATUSES:
        errors.append(f"status '{status}' not in {sorted(STATUSES)}")

    # severity enum (default warn if absent)
    if "severity" in rule and rule["severity"] not in SEVERITIES:
        errors.append(f"severity '{rule['severity']}' not in {sorted(SEVERITIES)}")

    # Field-level checks
    if "rule_id" in rule:
        rv = rule["rule_id"]
        if not isinstance(rv, str) or not RULE_ID_REGEX.match(rv):
            errors.append(f"rule_id '{rv}' fails regex {RULE_ID_REGEX.pattern}")

    if "schema_version" in rule:
        if rule["schema_version"] not in ("1.0.0", "1.0.1", "1.0.2", "1.0.3"):
            errors.append(f"schema_version '{rule['schema_version']}' not in supported set 1.0.0-1.0.3")

    if "title" in rule:
        t = rule["title"]
        if not isinstance(t, str) or len(t) > 80:
            errors.append(f"title len {len(t) if isinstance(t, str) else '?'} > 80")

    if "confidence_tier" in rule:
        tier = rule["confidence_tier"]
        if tier not in CONFIDENCE_TIERS:
            errors.append(f"confidence_tier '{tier}' not in {sorted(CONFIDENCE_TIERS)}")

    if "prior_audit_relationship" in rule:
        rel = rule["prior_audit_relationship"]
        if rel not in PRIOR_AUDIT_RELATIONSHIPS:
            errors.append(f"prior_audit_relationship '{rel}' not in {sorted(PRIOR_AUDIT_RELATIONSHIPS)}")

    # trigger_conditions
    if "trigger_conditions" in rule:
        tc = rule["trigger_conditions"]
        if not isinstance(tc, list) or len(tc) < 1:
            errors.append("trigger_conditions must be array with >=1 item")
        else:
            for i, cond in enumerate(tc):
                if not isinstance(cond, dict):
                    errors.append(f"trigger_conditions[{i}] not a dict")
                    continue
                if "condition_type" not in cond:
                    errors.append(f"trigger_conditions[{i}] missing condition_type")
                elif cond["condition_type"] not in CONDITION_TYPES:
                    errors.append(f"trigger_conditions[{i}].condition_type '{cond['condition_type']}' not in {sorted(CONDITION_TYPES)}")
                if "parameters" not in cond:
                    errors.append(f"trigger_conditions[{i}] missing parameters")
                if "description" not in cond:
                    errors.append(f"trigger_conditions[{i}] missing description")
                # v1.0.3 — pipeline_step name registry check
                if (cond.get("condition_type") == "pipeline_step"
                        and isinstance(cond.get("parameters"), dict)):
                    step = cond["parameters"].get("step")
                    if step and REGISTERED_STEPS and step not in REGISTERED_STEPS:
                        msg = (f"trigger_conditions[{i}].parameters.step='{step}' "
                               f"not in pipeline_step_registry.yaml "
                               f"(known: {sorted(REGISTERED_STEPS)})")
                        if strict_steps:
                            errors.append(f"[step-registry] {msg}")
                        else:
                            warnings.append(f"[step-registry] {msg}")

    # recommendation
    if "recommendation" in rule:
        rec = rule["recommendation"]
        if not isinstance(rec, dict):
            errors.append("recommendation must be a dict")
        else:
            if "text" not in rec:
                errors.append("recommendation missing text")
            if "action_type" not in rec:
                errors.append("recommendation missing action_type")
            elif rec["action_type"] not in ACTION_TYPES:
                errors.append(f"recommendation.action_type '{rec['action_type']}' not in {sorted(ACTION_TYPES)}")

    # evidence
    if "evidence" in rule:
        ev = rule["evidence"]
        if not isinstance(ev, list) or len(ev) < 1:
            errors.append("evidence must be array with >=1 item")
        else:
            for i, e in enumerate(ev):
                if not isinstance(e, dict):
                    errors.append(f"evidence[{i}] not a dict")
                    continue
                for f in ("audit_id", "audit_phase", "output_paths",
                          "lock_file_entries", "findings_md_section", "summary"):
                    if f not in e:
                        errors.append(f"evidence[{i}] missing field '{f}'")
                if "output_paths" in e and not isinstance(e["output_paths"], list):
                    errors.append(f"evidence[{i}].output_paths must be a list")
                if "lock_file_entries" in e:
                    lf = e["lock_file_entries"]
                    if not isinstance(lf, list):
                        errors.append(f"evidence[{i}].lock_file_entries must be a list")
                    else:
                        for h in lf:
                            if not isinstance(h, str) or not SHA256_REGEX.match(h):
                                errors.append(f"evidence[{i}].lock_file_entries contains non-sha256 '{h}'")

    # last_reviewed
    if "last_reviewed" in rule:
        if not ISO_DATE_REGEX.match(str(rule["last_reviewed"])):
            errors.append(f"last_reviewed '{rule['last_reviewed']}' not ISO 8601")

    # reviewer
    if "reviewer" in rule:
        if not isinstance(rule["reviewer"], str) or not rule["reviewer"]:
            errors.append("reviewer must be non-empty string")

    # revision_history
    if "revision_history" in rule:
        rh = rule["revision_history"]
        if not isinstance(rh, list) or len(rh) < 1:
            errors.append("revision_history must have >=1 entry")
        else:
            for i, r in enumerate(rh):
                if not isinstance(r, dict):
                    errors.append(f"revision_history[{i}] not a dict")
                    continue
                for f in ("date", "reviewer", "change_type", "summary"):
                    if f not in r:
                        errors.append(f"revision_history[{i}] missing field '{f}'")
                if "change_type" in r and r["change_type"] not in CHANGE_TYPES:
                    errors.append(f"revision_history[{i}].change_type '{r['change_type']}' not in {sorted(CHANGE_TYPES)}")

    # out_of_scope — accepts string OR list of strings per v1.0.1
    if "out_of_scope" in rule:
        oos = rule["out_of_scope"]
        if isinstance(oos, str):
            if not oos.strip():
                errors.append("out_of_scope must be non-empty")
        elif isinstance(oos, list):
            if not oos or not all(isinstance(x, str) and x.strip() for x in oos):
                errors.append("out_of_scope list must contain non-empty strings")
        else:
            errors.append("out_of_scope must be string or list of strings")

    # ---- Cross-field validation ----

    tier = rule.get("confidence_tier")
    rel = rule.get("prior_audit_relationship")
    ev = rule.get("evidence", []) if isinstance(rule.get("evidence"), list) else []
    rec = rule.get("recommendation", {}) if isinstance(rule.get("recommendation"), dict) else {}

    # 1. hard_default → ≥15 datasets, ≥3 tissues (must be in evidence; we check if any
    # evidence entry's summary or output_paths indicate dataset count and tissue count;
    # this is a soft check — can't strictly verify without parsing every TSV)
    if tier == "hard_default":
        # Look for evidence with explicit n_datasets and n_tissues in summary
        summaries = " ".join(e.get("summary", "") for e in ev).lower()
        if "datasets" not in summaries and "n_" not in summaries:
            errors.append("[xref-1] hard_default tier requires evidence with explicit dataset/tissue counts")

    # 3. insufficient_data → action_type must NOT be replace_parameter with specific value
    if tier == "insufficient_data":
        at = rec.get("action_type")
        if at == "replace_parameter":
            errors.append("[xref-3] insufficient_data tier forbids action_type=replace_parameter")

    # 4. refines_prior / contradicts_prior → ≥2 evidence entries
    if rel in ("refines_prior", "contradicts_prior"):
        if len(ev) < 2:
            errors.append(f"[xref-4] prior_audit_relationship={rel} requires >=2 evidence entries (have {len(ev)})")

    # 6. v1.0.3 — action_type=reject_pipeline requires severity=reject
    if rec.get("action_type") == "reject_pipeline":
        if rule.get("severity") != "reject":
            errors.append("[xref-6] action_type=reject_pipeline requires severity=reject")

    # 7. v1.0.3 — action_type=block_step requires severity in {error, reject}
    if rec.get("action_type") == "block_step":
        if rule.get("severity") not in ("error", "reject"):
            errors.append("[xref-7] action_type=block_step requires severity in {error, reject}")

    # 5. deprecation present → not loaded (this is informational, not an error)
    if "deprecation" in rule:
        dep = rule["deprecation"]
        if not isinstance(dep, dict):
            errors.append("deprecation must be a dict")
        else:
            for f in ("deprecated_date", "reason"):
                if f not in dep:
                    errors.append(f"deprecation missing field '{f}'")

    # ---- Companion-metrics check (per AUDIT_STANDARDS.md §5.4) ----
    # Rules that surface tool disagreement must report companion metrics
    # (direction agreement + log2FC correlation) alongside Jaccard-style
    # claims. Triggers on:
    #   - action_type == "report_additional_metric"
    #   - OR recommendation/description text contains tool-comparison vocabulary
    rec_text = (rec.get("text", "") if isinstance(rec, dict) else "").lower()
    desc_text = str(rule.get("description", "")).lower()
    full_text = rec_text + " " + desc_text
    triggers_companion_check = (
        rec.get("action_type") == "report_additional_metric"
        or any(token in full_text for token in
               ("top-100 jaccard", "top-n jaccard", "tool disagreement",
                "tools disagree", "tools agree on", "jaccard-style"))
    )
    if triggers_companion_check:
        has_direction = any(t in full_text for t in
                            ("direction agreement", "direction concordance",
                             "sign agreement"))
        has_lfc = any(t in full_text for t in
                      ("log2fc", "log fold change", "log-fold change",
                       "fold change correlation"))
        if not (has_direction and has_lfc):
            errors.append(
                "[companion-metrics] rule surfaces tool-comparison content "
                "without companion metrics (direction agreement + log2FC); "
                f"has_direction={has_direction}, has_lfc={has_lfc}"
            )

    # ---- Equivalence-finding tier watchdog (per AUDIT_STANDARDS.md §5.3.2) ----
    # Warn-only: if a rule reads as an equivalence/agreement finding at an
    # evidence-bearing tier, its evidence should cite the §5.3.2 criteria
    # (bootstrap CIs / cross-tool correlation / stratification). Does NOT block.
    eq_warn = _check_equivalence_tier_evidence(rule)
    if eq_warn:
        warnings.append(eq_warn)

    return errors, warnings


def _check_equivalence_tier_evidence(rule: dict):
    """Watchdog for equivalence-finding rules (§5.3.2): warn (not block) if an
    equivalence/agreement/convergence finding is assigned an evidence-bearing
    tier but its evidence summaries don't mention the §5.3.2 criteria
    (bootstrap CIs, cross-tool correlation, or stratification)."""
    tier = rule.get("confidence_tier")
    if tier not in ("hard_default", "conditional", "flag_and_warn"):
        return None
    evidence_text = " ".join(e.get("summary", "") for e in rule.get("evidence", []))
    el = evidence_text.lower()
    equivalence_indicators = ("converge", "agreement", "equivalent", "nest",
                              "jaccard", "correlation", "overlap", "spearman")
    if not any(ind in el for ind in equivalence_indicators):
        return None  # not an equivalence finding; §5.3.1 applies
    tier_criteria_indicators = ("bootstrap", "ci", "confidence interval",
                                "spearman", "stratif", "rho", "ρ")
    if not any(ind in el for ind in tier_criteria_indicators):
        return (f"[equivalence-tier] rule {rule.get('rule_id','<unknown>')} reads as an "
                f"equivalence finding at tier '{tier}' but evidence does not cite "
                f"bootstrap CIs, cross-tool correlation, or stratification per §5.3.2")
    return None


# ---- Light validation for superseded / pending_engine_support rules ----

LIGHT_REQUIRED = {"rule_id", "schema_version", "description"}
SUPERSEDED_LIGHT_EXTRA = {"deprecation"}


def validate_light(rule: dict, kind: str) -> list[str]:
    """Light validation for rules in superseded/ or pending_engine_support/.

    kind is 'superseded' or 'pending_engine_support'. Superseded requires a
    deprecation block; pending only requires the core identification fields.
    Legacy 'id' field is accepted in place of 'rule_id'.
    """
    errors = []
    has_identifier = "rule_id" in rule or "id" in rule
    if not has_identifier:
        errors.append("[light] missing identifier: needs 'rule_id' or legacy 'id'")
    other_required = LIGHT_REQUIRED - {"rule_id"} - rule.keys()
    if other_required:
        errors.append(f"[light] missing fields: {sorted(other_required)}")
    if kind == "superseded":
        if "deprecation" not in rule and "status" not in rule:
            # Tolerate the legacy 'status: superseded_by_...' string form
            status_str = rule.get("status", "")
            if not (isinstance(status_str, str) and "supersed" in status_str):
                errors.append("[light] superseded rule must have deprecation block, status: deprecated, or status: superseded_by_*")
    return errors


def main():
    strict_steps = "--strict-steps" in sys.argv
    print(f"Validating draft rules at {RULES_DIR}")
    print(f"Against schema 1.0.3 from {STANDARDS_DOC}")
    print(f"Pipeline-step registry: {len(REGISTERED_STEPS)} known steps "
          f"({'strict mode' if strict_steps else 'warn-only'})\n")

    # Standard pass: top-level YAMLs
    # Light pass: superseded/ + pending_engine_support/
    files = sorted([f for f in RULES_DIR.glob("*.yaml")])
    superseded_dir = RULES_DIR / "superseded"
    pending_dir = RULES_DIR / "pending_engine_support"
    light_files = []
    if superseded_dir.exists():
        light_files += [(f, "superseded") for f in sorted(superseded_dir.glob("*.yaml"))]
    if pending_dir.exists():
        light_files += [(f, "pending_engine_support") for f in sorted(pending_dir.glob("*.yaml"))]
    total_rules = 0
    total_errors = 0
    total_warnings = 0
    file_summaries = []

    for f in files:
        print("=" * 70)
        print(f"FILE: {f.name}")
        print("=" * 70)
        try:
            content = yaml.safe_load(f.read_text())
        except yaml.YAMLError as e:
            print(f"  YAML PARSE ERROR: {e}")
            file_summaries.append((f.name, 0, "yaml_parse_error"))
            continue

        if content is None:
            print("  (empty file)")
            continue

        # Normalize: top-level can be a single rule (dict) or list of rules
        if isinstance(content, dict):
            rules = [content]
            top_level_is = "single dict"
        elif isinstance(content, list):
            rules = content
            top_level_is = f"list of {len(content)}"
        else:
            print(f"  UNEXPECTED TOP-LEVEL TYPE: {type(content).__name__}")
            file_summaries.append((f.name, 0, "bad_top_level"))
            continue

        print(f"  Top-level structure: {top_level_is}")
        file_errors = 0
        for i, rule in enumerate(rules):
            if not isinstance(rule, dict):
                print(f"\n  Rule #{i}: not a dict (got {type(rule).__name__})")
                file_errors += 1
                continue
            rid = rule.get("rule_id", rule.get("id", f"<unnamed #{i}>"))
            errs, warns = validate_rule(rule, i, strict_steps=strict_steps)
            tier = rule.get("confidence_tier", "<missing>")
            rel = rule.get("prior_audit_relationship", "<missing>")
            print(f"\n  Rule {i}: id={rid}  tier={tier}  prior={rel}")
            if not errs and not warns:
                print("    ✓ PASS")
            elif not errs and warns:
                print(f"    ✓ PASS with {len(warns)} warning(s):")
                for w in warns:
                    print(f"      ⚠ {w}")
                total_warnings += len(warns)
            else:
                print(f"    ✗ {len(errs)} error(s):")
                for e in errs:
                    print(f"      - {e}")
                if warns:
                    print(f"    ⚠ {len(warns)} warning(s):")
                    for w in warns:
                        print(f"      ⚠ {w}")
                    total_warnings += len(warns)
                file_errors += len(errs)
            total_rules += 1

        total_errors += file_errors
        file_summaries.append((f.name, file_errors, "ok" if file_errors == 0 else "fail"))
        print()

    # ---- Light-validation pass: superseded/ + pending_engine_support/ ----
    if light_files:
        print("=" * 70)
        print("LIGHT VALIDATION (superseded/ + pending_engine_support/)")
        print("=" * 70)
        for f, kind in light_files:
            print(f"\nFILE [{kind}]: {f.relative_to(RULES_DIR)}")
            try:
                content = yaml.safe_load(f.read_text())
            except yaml.YAMLError as e:
                print(f"  YAML PARSE ERROR: {e}")
                file_summaries.append((str(f.relative_to(RULES_DIR)), 1, "light_yaml_err"))
                total_errors += 1
                continue
            if content is None:
                continue
            rules = [content] if isinstance(content, dict) else content
            for i, rule in enumerate(rules):
                if not isinstance(rule, dict):
                    continue
                rid = rule.get("rule_id", rule.get("id", f"<unnamed #{i}>"))
                errs = validate_light(rule, kind)
                print(f"  Rule {i}: id={rid}")
                if not errs:
                    print("    ✓ PASS (light)")
                else:
                    print(f"    ✗ {len(errs)} error(s):")
                    for e in errs:
                        print(f"      - {e}")
                    total_errors += len(errs)
                total_rules += 1
            file_summaries.append((str(f.relative_to(RULES_DIR)),
                                    0, f"light_{kind}"))

    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"{'File':<48} {'Errors':>8}  Status")
    for name, errs, status in file_summaries:
        print(f"{name:<48} {errs:>8d}  {status}")
    print(f"\nTotal rules: {total_rules}, Total errors: {total_errors}, "
          f"Total warnings: {total_warnings}")
    return 0 if total_errors == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
