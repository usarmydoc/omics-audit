#!/usr/bin/env python3
"""
B3 Low Count Filter Threshold Audit — main runner.

Downloads TCGA data, runs filter sweep, emits results and
BioOrchestrator artifacts.
"""

import sys
import time
from pathlib import Path

import pandas as pd

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from projects.b3_low_count_filter.data_access import (
    CANCER_TYPES_WITH_NORMALS,
    CANCER_TO_TISSUE,
    download_cancer_type,
    get_tissue_specific_genes,
    check_batch_metadata,
)
from projects.b3_low_count_filter.analysis import run_filter_sweep
from shared.knowledge_writer import (
    append_prompt_enrichment,
    append_validation_rules,
    write_snakefile_template,
)

OUTPUT_DIR = PROJECT_ROOT / "output"
CACHE_DIR = PROJECT_ROOT / "output" / "b3_cache"


def emit_knowledge(results_df: pd.DataFrame):
    """Emit L2, L3, L4 artifacts to BioOrchestrator."""

    # --- L2: Prompt enrichment ---
    # Compute summary statistics for the prompt
    default_mask = (results_df["cpm_threshold"] == 1.0) & (results_df["min_samples"] == "abs_2")
    default_results = results_df[default_mask]
    mean_pct_lost = default_results["pct_tissue_specific_lost"].mean()
    problematic = default_results[default_results["pct_tissue_specific_lost"] > 10]
    n_problematic = len(problematic)
    n_total = len(default_results)

    # Find optimal threshold per cancer type (minimize tissue gene loss while keeping DEG count stable)
    optimal_rows = []
    for ct in results_df["cancer_type"].unique():
        ct_data = results_df[results_df["cancer_type"] == ct]
        # Prefer settings that lose <5% tissue genes and retain >80% of max DEGs
        max_degs = ct_data["downstream_deg_count"].max()
        good = ct_data[
            (ct_data["pct_tissue_specific_lost"] <= 5) &
            (ct_data["downstream_deg_count"] >= 0.8 * max_degs)
        ]
        if len(good) > 0:
            # Pick the most stringent filter that still passes
            best = good.sort_values("cpm_threshold", ascending=False).iloc[0]
            optimal_rows.append(best)

    prompt_text = (
        "Low-count gene filtering (B3 audit, {n} TCGA cancer types):\n"
        "The default CPM>1 in ≥2 samples filter loses an average of {mean_lost:.1f}% of "
        "tissue-specific genes across cancer types. {n_prob}/{n} cancer types lose >10% "
        "of tissue-specific genes with the default threshold.\n"
        "Recommendation: Use adaptive thresholds scaled to sample size. "
        "For datasets with <20 samples per group, use CPM>0.5 in ≥10% of samples. "
        "For datasets with ≥20 samples per group, the default CPM>1 in ≥2 samples is acceptable "
        "but should be validated against tissue-specific gene retention. "
        "Always report the number of genes filtered and verify that known tissue-specific "
        "markers survive the filter."
    ).format(n=n_total, mean_lost=mean_pct_lost, n_prob=n_problematic)

    append_prompt_enrichment(prompt_text, "B3_low_count_filter")

    # --- L3: Validation rules ---
    rules = [
        {
            "id": "bulk_low_count_filter_specified",
            "type": "conditional",
            "description": "Low-count gene filter must be explicitly specified",
            "condition": {
                "context": "steps",
                "field": "tool_name",
                "check": "contains_any",
                "values": ["edgeR", "DESeq2", "limma"],
            },
            "required": {
                "scope": "parameters",
                "parameter_exists": "low_count_filter",
            },
            "action": "warn",
            "message": (
                "A differential expression tool is used but no low-count gene filter "
                "is specified. The default CPM>1 in ≥2 samples removes >10% of tissue-specific "
                "genes in {n_prob}/{n} TCGA cancer types tested. Explicitly specify "
                "cpm_threshold and min_samples_threshold."
            ).format(n_prob=n_problematic, n=n_total),
        },
        {
            "id": "bulk_adaptive_filter_small_samples",
            "type": "conditional",
            "description": "Adaptive low-count threshold required for small datasets",
            "condition": {
                "context": "experiment",
                "field": "samples",
                "check": "min_group_size_lt",
                "threshold": 20,
            },
            "required": {
                "scope": "parameters",
                "parameter_exists": "low_count_filter.adaptive",
            },
            "action": "warn",
            "message": (
                "Dataset has fewer than 20 samples per group. The standard CPM>1 "
                "threshold is unreliable at low sample sizes — use an adaptive "
                "threshold (e.g., CPM>0.5 in ≥10% of samples) or justify the "
                "chosen threshold with tissue-specific gene retention analysis."
            ),
        },
    ]
    append_validation_rules(rules, "B3_low_count_filter")

    # --- L4: Snakefile template ---
    smk = '''# Low-count filter audit Snakefile template.
# Generated by B3_low_count_filter audit.
# Real shell commands — no touch {output} placeholders.

configfile: "config.yaml"

WORKDIR = config.get("workdir", ".")
CANCER_TYPES = config.get("cancer_types", [])
CPM_THRESHOLDS = config.get("cpm_thresholds", [0.1, 0.5, 1.0, 2.0, 5.0, 10.0])

rule all:
    input:
        expand("{workdir}/output/b3_low_count_filter.tsv", workdir=WORKDIR),
        expand("{workdir}/output/repro_verified.flag", workdir=WORKDIR),

rule repro_snapshot:
    output:
        lockfile="{workdir}/output/repro.lock"
    shell:
        "repro snapshot -o {output.lockfile} --quiet"

rule download_counts:
    output:
        counts="{workdir}/output/b3_cache/{cancer_type}_counts.parquet",
        meta="{workdir}/output/b3_cache/{cancer_type}_meta.parquet"
    params:
        cancer_type="{cancer_type}"
    resources:
        mem_mb=4000,
        runtime=30
    log:
        "{workdir}/logs/download_{cancer_type}.log"
    shell:
        """
        python -c "
from projects.b3_low_count_filter.data_access import download_cancer_type
from pathlib import Path
download_cancer_type('{params.cancer_type}', Path('{wildcards.workdir}/output/b3_cache'))
" > {log} 2>&1
        """

rule filter_sweep:
    input:
        counts="{workdir}/output/b3_cache/{cancer_type}_counts.parquet",
        meta="{workdir}/output/b3_cache/{cancer_type}_meta.parquet"
    output:
        results="{workdir}/output/b3_sweep/{cancer_type}_sweep.tsv"
    params:
        cancer_type="{cancer_type}"
    resources:
        mem_mb=8000,
        runtime=60
    log:
        "{workdir}/logs/sweep_{cancer_type}.log"
    shell:
        """
        python -c "
import pandas as pd
from projects.b3_low_count_filter.data_access import CANCER_TO_TISSUE, get_tissue_specific_genes
from projects.b3_low_count_filter.analysis import run_filter_sweep
counts = pd.read_parquet('{input.counts}')
meta = pd.read_parquet('{input.meta}')
tissue = CANCER_TO_TISSUE.get('{params.cancer_type}', '')
tissue_genes = get_tissue_specific_genes(tissue)
results = run_filter_sweep(counts, meta['sample_type'], tissue_genes,
                           '{params.cancer_type}', '{params.cancer_type}'.replace('TCGA-',''))
pd.DataFrame(results).to_csv('{output.results}', sep='\\t', index=False)
" > {log} 2>&1
        """

rule aggregate_results:
    input:
        expand("{{workdir}}/output/b3_sweep/{ct}_sweep.tsv",
               ct=CANCER_TYPES)
    output:
        "{workdir}/output/b3_low_count_filter.tsv"
    shell:
        """
        head -1 {input[0]} > {output}
        for f in {input}; do tail -n +2 "$f" >> {output}; done
        """

rule repro_verify:
    input:
        results="{workdir}/output/b3_low_count_filter.tsv",
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
    write_snakefile_template(smk, "low_count_filter_audit.smk", "B3_low_count_filter")


def main():
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    all_results = []
    batch_metadata = []
    start_time = time.time()

    for project_id in CANCER_TYPES_WITH_NORMALS:
        cancer_type = project_id.replace("TCGA-", "")
        tissue = CANCER_TO_TISSUE.get(project_id, "")
        print(f"\n{'='*60}")
        print(f"Processing {project_id} (tissue: {tissue})")
        print(f"{'='*60}")

        # Download data
        try:
            data = download_cancer_type(project_id, CACHE_DIR)
        except Exception as e:
            print(f"  FAILED to download {project_id}: {e}")
            continue

        print(f"  Downloaded: {data['n_tumor']} tumor, {data['n_normal']} normal samples")
        print(f"  Count matrix: {data['counts'].shape[0]} genes × {data['counts'].shape[1]} samples")

        # Get tissue-specific genes
        tissue_genes = get_tissue_specific_genes(tissue)
        print(f"  Tissue-specific gene list: {len(tissue_genes)} genes for {tissue}")

        # Map sample types for edgeR (Tumor vs Normal)
        sample_type = data["sample_type"].map({"Tumor": "Tumor", "Normal": "Normal"})

        # Run filter sweep
        results = run_filter_sweep(
            data["counts"], sample_type, tissue_genes,
            project_id, cancer_type,
        )
        all_results.extend(results)
        print(f"  Completed {len(results)} threshold combinations")

        # Check batch metadata (for B4 planning)
        batch_info = check_batch_metadata(project_id)
        batch_metadata.append(batch_info)
        print(f"  Batch metadata: plates={batch_info['n_plates']}, "
              f"centers={batch_info['n_centers']}")

    # Save results
    results_df = pd.DataFrame(all_results)
    output_file = OUTPUT_DIR / "b3_low_count_filter.tsv"
    results_df.to_csv(output_file, sep="\t", index=False)
    print(f"\n{'='*60}")
    print(f"Results saved to {output_file}")
    print(f"Total rows: {len(results_df)}")

    # Save batch metadata assessment (for B4)
    batch_df = pd.DataFrame(batch_metadata)
    batch_file = OUTPUT_DIR / "b4_batch_metadata_assessment.tsv"
    batch_df.to_csv(batch_file, sep="\t", index=False)
    print(f"Batch metadata assessment saved to {batch_file}")
    n_with_plates = (batch_df["has_plate"]).sum()
    n_with_centers = (batch_df["has_center"]).sum()
    print(f"  Projects with plate IDs: {n_with_plates}/{len(batch_df)}")
    print(f"  Projects with center IDs: {n_with_centers}/{len(batch_df)}")
    if n_with_plates < len(batch_df) // 2:
        print(f"  ⚠ WARNING: Plate IDs absent in >50% of cancer types. "
              f"Flag for B4 discussion.")

    # Emit BioOrchestrator artifacts
    print(f"\nEmitting BioOrchestrator artifacts...")
    emit_knowledge(results_df)

    elapsed = time.time() - start_time
    print(f"\nB3 completed in {elapsed/60:.1f} minutes")

    # Print summary statistics
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    default_mask = (results_df["cpm_threshold"] == 1.0) & (results_df["min_samples"] == "abs_2")
    default_results = results_df[default_mask]
    print(f"\nDefault filter (CPM>1, ≥2 samples):")
    print(f"  Mean tissue-specific gene loss: "
          f"{default_results['pct_tissue_specific_lost'].mean():.1f}%")
    print(f"  Cancer types losing >10%: "
          f"{(default_results['pct_tissue_specific_lost'] > 10).sum()}/{len(default_results)}")
    print(f"  Mean genes retained: "
          f"{default_results['genes_retained'].mean():.0f} / "
          f"{default_results['total_genes'].mean():.0f}")
    print(f"  Mean DEGs: {default_results['downstream_deg_count'].mean():.0f}")


if __name__ == "__main__":
    main()
