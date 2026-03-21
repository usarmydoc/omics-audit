#!/usr/bin/env python3
"""Master sequential runner for all audit projects.

Order: shortest to longest estimated runtime.
Runs each project, then its knowledge emission, then BioOrchestrator tests.
Stops if tests fail.

Set SCRNASEQ_VENV to override the default venv path.
Set BIOORCHESTRATOR_ROOT to override the BioOrchestrator location.
"""

import os
import subprocess
import sys
import time
from pathlib import Path

MASTER_DIR = Path(__file__).resolve().parent
REPO_ROOT = MASTER_DIR.parent
PROJECTS_DIR = MASTER_DIR / "projects"

# Venv: override with SCRNASEQ_VENV env var, otherwise use .venv in repo root
VENV = os.environ.get("SCRNASEQ_VENV", str(REPO_ROOT / ".venv"))

# (label, run_script, emit_script, project_dir)
PROJECTS = [
    (
        "P6 Min Cell Filter",
        str(PROJECTS_DIR / "p6_min_cell_filter" / "run_p6.py"),
        str(PROJECTS_DIR / "p6_min_cell_filter" / "emit_knowledge.py"),
        str(PROJECTS_DIR / "p6_min_cell_filter"),
    ),
    (
        "P9 Clustering Resolution",
        str(PROJECTS_DIR / "p9_clustering_resolution" / "run_p9.py"),
        str(PROJECTS_DIR / "p9_clustering_resolution" / "emit_knowledge.py"),
        str(PROJECTS_DIR / "p9_clustering_resolution"),
    ),
    (
        "P4 Pseudobulk",
        str(PROJECTS_DIR / "p4_pseudobulk" / "run_p4.py"),
        str(PROJECTS_DIR / "p4_pseudobulk" / "emit_knowledge.py"),
        str(PROJECTS_DIR / "p4_pseudobulk"),
    ),
    (
        "P2 Doublet Audit",
        str(PROJECTS_DIR / "p2_doublet_audit" / "run_p2.py"),
        str(PROJECTS_DIR / "p2_doublet_audit" / "emit_knowledge.py"),
        str(PROJECTS_DIR / "p2_doublet_audit"),
    ),
    (
        "P7 Normalization",
        str(PROJECTS_DIR / "p7_normalization" / "run_p7.py"),
        str(PROJECTS_DIR / "p7_normalization" / "emit_knowledge.py"),
        str(PROJECTS_DIR / "p7_normalization"),
    ),
    (
        "P5 Annotation Tools",
        str(PROJECTS_DIR / "p5_annotation_tools" / "run_p5.py"),
        str(PROJECTS_DIR / "p5_annotation_tools" / "emit_knowledge.py"),
        str(PROJECTS_DIR / "p5_annotation_tools"),
    ),
    (
        "P3 Batch Correction",
        str(PROJECTS_DIR / "p3_batch_correction" / "run_p3.py"),
        str(PROJECTS_DIR / "p3_batch_correction" / "emit_knowledge.py"),
        str(PROJECTS_DIR / "p3_batch_correction"),
    ),
    (
        "P8 Integration Selection",
        str(PROJECTS_DIR / "p8_integration_selection" / "run_p8.py"),
        str(PROJECTS_DIR / "p8_integration_selection" / "emit_knowledge.py"),
        str(PROJECTS_DIR / "p8_integration_selection"),
    ),
]

BIOORCHESTRATOR_DIR = os.environ.get(
    "BIOORCHESTRATOR_ROOT",
    str(Path.home() / "bioorchestrator"),
).removesuffix("/src/bioorchestrator")


def run_cmd(label: str, cmd: list[str], cwd: str | None = None) -> int:
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}")
    t0 = time.time()
    result = subprocess.run(cmd, cwd=cwd)
    elapsed = time.time() - t0
    status = "PASS" if result.returncode == 0 else "FAIL"
    print(f"\n[{status}] {label} — {elapsed/60:.1f} min (exit {result.returncode})")
    return result.returncode


def run_tests() -> int:
    print("\n--- Running BioOrchestrator tests ---")
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "tests/", "-q"],
        cwd=BIOORCHESTRATOR_DIR,
        capture_output=True, text=True,
    )
    last_line = result.stdout.strip().split("\n")[-1] if result.stdout.strip() else ""
    print(f"  {last_line}")
    return result.returncode


def main():
    t_start = time.time()
    completed = []
    failed = []

    python = f"{VENV}/bin/python"

    for label, run_script, emit_script, project_dir in PROJECTS:
        # Run project
        rc = run_cmd(f"{label} — run", [python, run_script], cwd=project_dir)
        if rc != 0:
            print(f"\n*** {label} run FAILED (exit {rc}), skipping emit + continuing ***")
            failed.append(label)
            continue

        # Emit knowledge
        rc = run_cmd(f"{label} — emit knowledge", [python, emit_script], cwd=project_dir)
        if rc != 0:
            print(f"\n*** {label} emit FAILED (exit {rc}), continuing ***")
            failed.append(f"{label} (emit)")
            continue

        # Run BioOrchestrator tests
        rc = run_tests()
        if rc != 0:
            print(f"\n*** BioOrchestrator tests FAILED after {label} ***")
            print("*** STOPPING — fix rules before continuing ***")
            failed.append(f"{label} (tests)")
            break

        completed.append(label)
        print(f"\n✓ {label} complete ({len(completed)}/{len(PROJECTS)})")

    # Summary
    elapsed = time.time() - t_start
    print(f"\n{'='*60}")
    print(f"  MASTER RUN COMPLETE — {elapsed/60:.1f} minutes")
    print(f"{'='*60}")
    print(f"  Completed: {len(completed)}/{len(PROJECTS)}")
    for c in completed:
        print(f"    ✓ {c}")
    if failed:
        print(f"  Failed: {len(failed)}")
        for f in failed:
            print(f"    ✗ {f}")

    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
