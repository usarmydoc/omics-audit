#!/usr/bin/env python3
"""Unit tests for validate_rules.py — v1.0.3-specific checks.

Covers:
  - schema_version 1.0.3 accepted
  - pipeline_step name registry: warn by default, error under --strict-steps
  - action_type reject_pipeline + block_step accepted as valid enum values
  - cross-field validator [xref-6]: reject_pipeline requires severity=reject
  - cross-field validator [xref-7]: block_step requires severity in {error,reject}

Run: python3 -m pytest /mnt/nvme1/omics-audit/phase2/scripts/test_validate_rules.py -v
Or:  python3 /mnt/nvme1/omics-audit/phase2/scripts/test_validate_rules.py
"""
from __future__ import annotations

import sys
from pathlib import Path

# Import the validator module
SCRIPTS_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPTS_DIR))
import validate_rules as vr  # noqa: E402

# A 64-char hex placeholder that passes SHA256_REGEX
HASH = "a" * 64


def _baseline_rule(**overrides) -> dict:
    """Construct a minimal v1.0.3 rule dict that passes validation. Apply overrides on top."""
    rule = {
        "rule_id": "test_baseline_rule",
        "schema_version": "1.0.3",
        "status": "draft",
        "title": "Test baseline rule",
        "description": "A minimal valid rule used as a test baseline.",
        "created_date": "2026-05-17",
        "trigger_conditions": [
            {
                "condition_type": "pipeline_step",
                "parameters": {"step": "pathway_enrichment"},
                "description": "Fires at pathway enrichment step.",
            }
        ],
        "recommendation": {
            "text": "Test recommendation text.",
            "action_type": "flag_only",
        },
        "confidence_tier": "conditional",
        "prior_audit_relationship": "novel",
        "evidence": [
            {
                "audit_id": "TEST",
                "audit_phase": "Test",
                "output_paths": ["/tmp/test"],
                "lock_file_entries": [HASH],
                "findings_md_section": "test.md",
                "summary": "Test summary with n_datasets=10.",
            }
        ],
        "out_of_scope": ["test scope"],
    }
    rule.update(overrides)
    return rule


# ---- schema_version ----

def test_schema_version_103_accepted():
    errors, warnings = vr.validate_rule(_baseline_rule(), 0)
    assert errors == [], f"v1.0.3 baseline should pass: {errors}"


def test_schema_version_102_still_accepted():
    errors, _ = vr.validate_rule(_baseline_rule(schema_version="1.0.2"), 0)
    assert errors == [], f"v1.0.2 backward-compat should pass: {errors}"


def test_schema_version_unknown_rejected():
    errors, _ = vr.validate_rule(_baseline_rule(schema_version="9.9.9"), 0)
    assert any("schema_version" in e for e in errors), \
        f"Unknown schema_version should error: {errors}"


# ---- pipeline_step registry ----

def test_known_step_no_warning():
    """All 4 registry step names should pass without warning."""
    for step in ("pathway_enrichment", "scrnaseq_qc_filtering",
                 "scrnaseq_clustering_resolution_selection", "scrnaseq_de_test"):
        r = _baseline_rule()
        r["trigger_conditions"][0]["parameters"]["step"] = step
        errors, warnings = vr.validate_rule(r, 0)
        assert errors == [], f"Known step '{step}' should pass: {errors}"
        assert not any("step-registry" in w for w in warnings), \
            f"Known step '{step}' should not warn: {warnings}"


def test_unknown_step_warns_default():
    """Unknown step name yields a warning (not error) by default."""
    r = _baseline_rule()
    r["trigger_conditions"][0]["parameters"]["step"] = "totally_made_up_step"
    errors, warnings = vr.validate_rule(r, 0)
    assert errors == [], f"Unknown step should not error by default: {errors}"
    assert any("step-registry" in w for w in warnings), \
        f"Unknown step should warn: {warnings}"


def test_unknown_step_errors_under_strict():
    """Unknown step name becomes an error under --strict-steps."""
    r = _baseline_rule()
    r["trigger_conditions"][0]["parameters"]["step"] = "totally_made_up_step"
    errors, warnings = vr.validate_rule(r, 0, strict_steps=True)
    assert any("step-registry" in e for e in errors), \
        f"Strict mode should error: {errors}"
    assert not any("step-registry" in w for w in warnings), \
        f"Strict mode should not also warn: {warnings}"


def test_non_pipeline_step_condition_ignored():
    """Step-registry check fires ONLY on pipeline_step conditions."""
    r = _baseline_rule()
    r["trigger_conditions"] = [{
        "condition_type": "parameter_value_check",
        "parameters": {"step": "totally_made_up_step",  # not validated here
                       "parameter": "foo", "expected_value": "bar"},
        "description": "Not a pipeline_step condition.",
    }]
    errors, warnings = vr.validate_rule(r, 0)
    assert not any("step-registry" in w for w in warnings), \
        f"step check should be scoped to pipeline_step condition_type: {warnings}"


# ---- action_type enum extensions ----

def test_reject_pipeline_accepted():
    r = _baseline_rule()
    r["recommendation"]["action_type"] = "reject_pipeline"
    r["severity"] = "reject"
    errors, _ = vr.validate_rule(r, 0)
    assert errors == [], f"reject_pipeline + severity:reject should pass: {errors}"


def test_block_step_accepted():
    r = _baseline_rule()
    r["recommendation"]["action_type"] = "block_step"
    r["severity"] = "error"
    errors, _ = vr.validate_rule(r, 0)
    assert errors == [], f"block_step + severity:error should pass: {errors}"


# ---- Cross-field validators xref-6 and xref-7 ----

def test_xref6_reject_pipeline_requires_severity_reject():
    """[xref-6]: action_type=reject_pipeline + severity!=reject should error."""
    for sev in ("info", "warn", "error"):
        r = _baseline_rule()
        r["recommendation"]["action_type"] = "reject_pipeline"
        r["severity"] = sev
        errors, _ = vr.validate_rule(r, 0)
        assert any("xref-6" in e for e in errors), \
            f"reject_pipeline + severity:{sev} should trigger xref-6: {errors}"


def test_xref6_passes_with_severity_reject():
    r = _baseline_rule()
    r["recommendation"]["action_type"] = "reject_pipeline"
    r["severity"] = "reject"
    errors, _ = vr.validate_rule(r, 0)
    assert not any("xref-6" in e for e in errors), \
        f"reject_pipeline + severity:reject should pass xref-6: {errors}"


def test_xref7_block_step_requires_error_or_reject():
    """[xref-7]: action_type=block_step + severity in {info,warn} should error."""
    for sev in ("info", "warn"):
        r = _baseline_rule()
        r["recommendation"]["action_type"] = "block_step"
        r["severity"] = sev
        errors, _ = vr.validate_rule(r, 0)
        assert any("xref-7" in e for e in errors), \
            f"block_step + severity:{sev} should trigger xref-7: {errors}"


def test_xref7_passes_with_error_or_reject():
    for sev in ("error", "reject"):
        r = _baseline_rule()
        r["recommendation"]["action_type"] = "block_step"
        r["severity"] = sev
        errors, _ = vr.validate_rule(r, 0)
        assert not any("xref-7" in e for e in errors), \
            f"block_step + severity:{sev} should pass xref-7: {errors}"


# ---- Registry loaded sanity check ----

def test_registry_loaded():
    """The registry yaml should have loaded at module-import time."""
    assert len(vr.REGISTERED_STEPS) >= 4, \
        f"Expected ≥4 registered steps, got {sorted(vr.REGISTERED_STEPS)}"
    assert "pathway_enrichment" in vr.REGISTERED_STEPS
    assert "scrnaseq_qc_filtering" in vr.REGISTERED_STEPS


if __name__ == "__main__":
    # Allow direct execution without pytest
    import traceback
    tests = [(n, fn) for n, fn in globals().items()
             if n.startswith("test_") and callable(fn)]
    passed = failed = 0
    for name, fn in tests:
        try:
            fn()
            print(f"  PASS  {name}")
            passed += 1
        except AssertionError as e:
            print(f"  FAIL  {name}\n    {e}")
            failed += 1
        except Exception:
            print(f"  ERROR {name}")
            traceback.print_exc()
            failed += 1
    print(f"\n{passed} passed, {failed} failed (of {len(tests)} tests)")
    sys.exit(0 if failed == 0 else 1)
