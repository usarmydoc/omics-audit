#!/usr/bin/env python3
"""
B4 Batch Correction — main runner.

Compares ComBat, limma removeBatchEffect, and SVA on TCGA datasets
with plate-level batch structure.
"""

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from projects.b4_batch_correction.data_access import CANCER_TYPES, load_with_batch
from projects.b4_batch_correction.analysis import (
    _normalize_logcpm, _filter_low_counts,
    run_combat, run_limma_remove_batch, run_sva,
    evaluate_correction,
)
from shared.knowledge_writer import (
    append_prompt_enrichment,
    append_validation_rules,
    write_snakefile_template,
)

OUTPUT_DIR = PROJECT_ROOT / "output"


def _filter_low_counts(counts, min_cpm=1.0, min_samples=2):
    lib_sizes = counts.sum(axis=0)
    cpm = counts.div(lib_sizes, axis=1) * 1e6
    keep = (cpm >= min_cpm).sum(axis=1) >= min_samples
    return counts[keep]


def emit_knowledge(results_df: pd.DataFrame):
    """Emit L2, L3, L4 artifacts."""

    # --- L2 ---
    combat_results = results_df[results_df["method"] == "ComBat"]
    limma_results = results_df[results_df["method"] == "limma_removeBatch"]
    sva_results = results_df[results_df["method"] == "SVA"]

    mean_batch_before = results_df["batch_variance_before"].mean()
    combat_batch_after = combat_results["batch_variance_after"].mean()
    limma_batch_after = limma_results["batch_variance_after"].mean()
    sva_batch_after = sva_results["batch_variance_after"].mean()

    combat_bio = combat_results["biological_variance_retained"].mean()
    limma_bio = limma_results["biological_variance_retained"].mean()
    sva_bio = sva_results["biological_variance_retained"].mean()

    prompt_text = (
        "Batch correction comparison (B4 audit, 16 TCGA cancer types, plate-level batches):\n"
        f"Mean batch variance before correction: {mean_batch_before:.3f}.\n"
        f"Batch variance after: ComBat={combat_batch_after:.3f}, "
        f"limma removeBatchEffect={limma_batch_after:.3f}, SVA={sva_batch_after:.3f}.\n"
        f"Biological variance retained: ComBat={combat_bio:.3f}, "
        f"limma={limma_bio:.3f}, SVA={sva_bio:.3f}.\n"
        "Recommendation: Always check that batch correction does not remove "
        "biological signal. ComBat and limma removeBatchEffect require explicit "
        "batch labels; SVA discovers latent factors but may capture biological "
        "variation. Include the biological variable of interest in the model "
        "matrix when running ComBat or removeBatchEffect to protect it from "
        "removal. Verify correction by checking that batch variance decreases "
        "while tumor/normal separation is preserved in PCA."
    )
    append_prompt_enrichment(prompt_text, "B4_batch_correction")

    # --- L3 ---
    rules = [
        {
            "id": "bulk_batch_correction_requires_metadata",
            "type": "conditional",
            "description": "Batch correction requires batch metadata",
            "condition": {
                "context": "steps",
                "field": "tool_name",
                "check": "contains_any",
                "values": ["ComBat", "combat", "removeBatchEffect", "SVA", "sva"],
            },
            "required": {
                "scope": "experiment",
                "field": "batch_variable",
                "check": "exists",
            },
            "action": "warn",
            "message": (
                "Batch correction tool is used but no batch variable is specified "
                "in the experiment metadata. ComBat and removeBatchEffect require "
                "explicit batch labels. Specify the batch variable (e.g., plate, "
                "sequencing run, center)."
            ),
        },
        {
            "id": "bulk_batch_correction_bio_check",
            "type": "conditional",
            "description": "Require biological variance check after batch correction",
            "condition": {
                "context": "steps",
                "field": "tool_name",
                "check": "contains_any",
                "values": ["ComBat", "combat", "removeBatchEffect", "SVA", "sva"],
            },
            "required": {
                "scope": "steps",
                "step_type_any": ["pca_visualization", "variance_check"],
            },
            "action": "warn",
            "message": (
                "Batch correction is applied but no downstream variance check "
                "or PCA visualization step is included. Aggressive batch correction "
                "can remove biologically relevant variance. Add a PCA step after "
                "batch correction to verify that biological groups remain separated."
            ),
        },
        {
            "id": "bulk_sva_min_samples",
            "type": "conditional",
            "description": "SVA requires minimum samples per group",
            "condition": {
                "context": "steps",
                "field": "tool_name",
                "check": "contains_any",
                "values": ["SVA", "sva"],
            },
            "required": {
                "scope": "experiment",
                "field": "samples",
                "check": "min_group_size_gte",
                "threshold": 6,
            },
            "action": "warn",
            "message": (
                "SVA is used but at least one group has fewer than 6 samples. "
                "SVA requires sufficient samples to reliably estimate surrogate "
                "variables. Consider using ComBat with known batch labels instead."
            ),
        },
    ]
    append_validation_rules(rules, "B4_batch_correction")

    # --- L4 ---
    smk = '''# Batch correction comparison Snakefile template.
# Generated by B4_batch_correction audit.
# Real shell commands — no touch {output} placeholders.

configfile: "config.yaml"

WORKDIR = config.get("workdir", ".")
CANCER_TYPES = config.get("cancer_types", [])

rule all:
    input:
        expand("{workdir}/output/b4_batch_correction.tsv", workdir=WORKDIR),
        expand("{workdir}/output/repro_verified.flag", workdir=WORKDIR),

rule repro_snapshot:
    output:
        lockfile="{workdir}/output/repro.lock"
    shell:
        "repro snapshot -o {output.lockfile} --quiet"

rule batch_correction:
    input:
        counts="{workdir}/output/b3_cache/{cancer_type}_counts.parquet",
        meta="{workdir}/output/b3_cache/{cancer_type}_meta.parquet"
    output:
        results="{workdir}/output/b4_results/{cancer_type}_batch.tsv"
    params:
        cancer_type="{cancer_type}"
    resources:
        mem_mb=16000,
        runtime=30
    log:
        "{workdir}/logs/batch_{cancer_type}.log"
    shell:
        """
        python -c "
from projects.b4_batch_correction.run import process_cancer_type
process_cancer_type('{params.cancer_type}', '{output.results}')
" > {log} 2>&1
        """

rule aggregate_results:
    input:
        expand("{{workdir}}/output/b4_results/{ct}_batch.tsv",
               ct=CANCER_TYPES)
    output:
        "{workdir}/output/b4_batch_correction.tsv"
    shell:
        """
        head -1 {input[0]} > {output}
        for f in {input}; do tail -n +2 "$f" >> {output}; done
        """

rule repro_verify:
    input:
        results="{workdir}/output/b4_batch_correction.tsv",
        lockfile="{workdir}/output/repro.lock"
    output:
        "{workdir}/output/repro_verified.flag"
    shell:
        """
        repro verify -l {input.lockfile} --quiet
        if [ $? -ne 0 ]; then
            echo "REPRO VERIFY FAILED" >&2
            exit 1
        fi
        touch {output}
        """
'''
    write_snakefile_template(smk, "batch_correction_bulk.smk", "B4_batch_correction")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    all_results = []
    start_time = time.time()

    for i, project_id in enumerate(CANCER_TYPES):
        cancer_type = project_id.replace("TCGA-", "")
        print(f"\n{'='*60}")
        print(f"[{i+1}/16] {project_id}")
        print(f"{'='*60}")

        try:
            data = load_with_batch(project_id)
        except Exception as e:
            print(f"  FAILED to load: {e}")
            continue

        counts = data["counts"]
        sample_type = data["sample_type"]
        plate = data["plate"]
        n_plates = data["n_plates"]

        print(f"  Samples: {counts.shape[1]}, Plates: {n_plates}")

        if n_plates < 2:
            print(f"  SKIPPED: fewer than 2 plates")
            continue

        # Filter and normalize
        filtered = _filter_low_counts(counts)
        print(f"  Filtered genes: {filtered.shape[0]}")

        logcpm = _normalize_logcpm(filtered)
        print(f"  Log-CPM computed")

        # Run each correction method
        methods = {
            "ComBat": lambda: run_combat(logcpm, plate, sample_type),
            "limma_removeBatch": lambda: run_limma_remove_batch(logcpm, plate, sample_type),
            "SVA": lambda: run_sva(logcpm, sample_type),
        }

        for method_name, method_fn in methods.items():
            print(f"  Running {method_name}...")
            try:
                corrected = method_fn()
                result = evaluate_correction(
                    logcpm, corrected, filtered, plate, sample_type,
                    method_name, project_id, cancer_type,
                )
                all_results.append(result)
                print(f"    Batch var: {result['batch_variance_before']:.4f} → "
                      f"{result['batch_variance_after']:.4f}, "
                      f"Bio retained: {result['biological_variance_retained']:.4f}, "
                      f"ARI: {result['ari_before']:.4f} → {result['ari_after']:.4f}")
            except Exception as e:
                print(f"    {method_name} FAILED: {e}")
                import traceback
                traceback.print_exc()

    # Save results
    results_df = pd.DataFrame(all_results)
    output_file = OUTPUT_DIR / "b4_batch_correction.tsv"
    results_df.to_csv(output_file, sep="\t", index=False)
    print(f"\n{'='*60}")
    print(f"Results saved to {output_file}")
    print(f"Total rows: {len(results_df)}")

    # Emit BioOrchestrator artifacts
    print(f"\nEmitting BioOrchestrator artifacts...")
    emit_knowledge(results_df)

    elapsed = time.time() - start_time
    print(f"\nB4 completed in {elapsed/60:.1f} minutes")

    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    for method in results_df["method"].unique():
        subset = results_df[results_df["method"] == method]
        print(f"\n{method}:")
        print(f"  Batch variance: {subset['batch_variance_before'].mean():.4f} → "
              f"{subset['batch_variance_after'].mean():.4f}")
        print(f"  Bio variance retained: {subset['biological_variance_retained'].mean():.4f}")
        print(f"  ARI: {subset['ari_before'].mean():.4f} → {subset['ari_after'].mean():.4f}")
        if "deg_stability" in subset.columns:
            print(f"  DEG stability: {subset['deg_stability'].mean():.4f}")


if __name__ == "__main__":
    main()
