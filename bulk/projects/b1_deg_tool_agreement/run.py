#!/usr/bin/env python3
"""
B1 DEG Tool Agreement — main runner.

Runs DESeq2, edgeR, and limma-voom on TCGA tumor vs normal comparisons
and measures agreement.
"""

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from projects.b1_deg_tool_agreement.data_access import (
    CANCER_TYPES, CANCER_TO_TISSUE, load_cancer_type,
)
from projects.b1_deg_tool_agreement.analysis import (
    _filter_low_counts, run_deseq2, run_edger, run_limma_voom,
    compare_tools,
)
from shared.knowledge_writer import (
    append_prompt_enrichment,
    append_validation_rules,
    write_snakefile_template,
)

OUTPUT_DIR = PROJECT_ROOT / "output"


def emit_knowledge(results_df: pd.DataFrame):
    """Emit L2, L3, L4 artifacts."""

    # --- L2 ---
    mean_jaccard = results_df["jaccard_top100"].mean()
    mean_dir_agree = results_df["direction_agreement"].mean()
    mean_lfc_corr = results_df["log2fc_correlation"].mean()

    # Per-pair stats
    pair_stats = {}
    for pair in results_df["tool_pair"].unique():
        subset = results_df[results_df["tool_pair"] == pair]
        pair_stats[pair] = {
            "jaccard": subset["jaccard_top100"].mean(),
            "dir": subset["direction_agreement"].mean(),
            "lfc": subset["log2fc_correlation"].mean(),
        }

    # DEG counts
    mean_deseq2 = results_df["n_degs_deseq2"].mean()
    mean_edger = results_df["n_degs_edger"].mean()
    mean_limma = results_df["n_degs_limma"].mean()

    prompt_text = (
        "DEG tool agreement (B1 audit, 16 TCGA cancer types, tumor vs normal):\n"
        f"Mean Jaccard similarity of top 100 DEGs: {mean_jaccard:.3f} "
        f"(range across cancer types: {results_df['jaccard_top100'].min():.3f}-"
        f"{results_df['jaccard_top100'].max():.3f}).\n"
        f"Mean direction agreement: {mean_dir_agree:.3f}. "
        f"Mean log2FC correlation: {mean_lfc_corr:.3f}.\n"
        f"Mean DEGs per tool (FDR<0.05): DESeq2={mean_deseq2:.0f}, "
        f"edgeR={mean_edger:.0f}, limma-voom={mean_limma:.0f}.\n"
        f"By pair: DESeq2-edgeR Jaccard={pair_stats.get('DESeq2_vs_edgeR', {}).get('jaccard', 0):.3f}, "
        f"DESeq2-limma={pair_stats.get('DESeq2_vs_limma', {}).get('jaccard', 0):.3f}, "
        f"edgeR-limma={pair_stats.get('edgeR_vs_limma', {}).get('jaccard', 0):.3f}.\n"
        "Recommendation: DEG tool choice materially affects which genes are reported as "
        "significant. Always specify the DEG tool used and justify the choice. For "
        "high-confidence results, run at least two tools and report the intersection. "
        "edgeR and limma-voom typically agree most closely; DESeq2 tends to be more "
        "conservative with fewer DEGs."
    )
    append_prompt_enrichment(prompt_text, "B1_deg_tool_agreement")

    # --- L3 ---
    rules = [
        {
            "id": "bulk_deg_tool_specified",
            "type": "conditional",
            "description": "DE tool must be explicitly named",
            "condition": {
                "context": "steps",
                "field": "step_type",
                "check": "contains_any",
                "values": ["differential_expression"],
            },
            "required": {
                "scope": "parameters",
                "parameter_exists": "de_tool",
            },
            "action": "warn",
            "message": (
                "Differential expression step has no tool specified. "
                f"DESeq2, edgeR, and limma-voom agree on only {mean_jaccard:.0%} of "
                "top 100 DEGs (mean Jaccard across 16 TCGA cancer types). "
                "Explicitly specify the DE tool."
            ),
        },
        {
            "id": "bulk_deg_min_replicates_per_group",
            "type": "conditional",
            "description": "Require minimum 3 replicates per group for DE",
            "condition": {
                "context": "experiment",
                "field": "samples",
                "check": "min_group_size_lt",
                "threshold": 3,
            },
            "required": {
                "scope": "steps",
                "step_type_none": ["differential_expression"],
            },
            "action": "reject",
            "message": (
                "Differential expression is in the pipeline but at least one group "
                "has fewer than 3 biological replicates. This is insufficient for "
                "reliable dispersion estimation in any DE tool."
            ),
        },
        {
            "id": "bulk_single_deg_tool_warning",
            "type": "conditional",
            "description": "Warn if only one DE tool used without justification",
            "condition": {
                "context": "steps",
                "field": "step_type",
                "check": "count_eq",
                "value": "differential_expression",
                "count": 1,
            },
            "required": {
                "scope": "parameters",
                "parameter_exists": "de_tool_justification",
            },
            "action": "flag",
            "message": (
                "Only one DE tool is used. DEG tool agreement is low "
                f"(mean Jaccard={mean_jaccard:.3f} for top 100 DEGs). "
                "Consider running a second tool for validation or provide "
                "justification for the single-tool choice."
            ),
        },
    ]
    append_validation_rules(rules, "B1_deg_tool_agreement")

    # --- L4 ---
    smk = '''# DEG tool comparison Snakefile template.
# Generated by B1_deg_tool_agreement audit.
# Real shell commands — no touch {output} placeholders.

configfile: "config.yaml"

WORKDIR = config.get("workdir", ".")
CANCER_TYPES = config.get("cancer_types", [])

rule all:
    input:
        expand("{workdir}/output/b1_deg_tool_agreement.tsv", workdir=WORKDIR),
        expand("{workdir}/output/repro_verified.flag", workdir=WORKDIR),

rule repro_snapshot:
    output:
        lockfile="{workdir}/output/repro.lock"
    shell:
        "repro snapshot -o {output.lockfile} --quiet"

rule run_deg_comparison:
    input:
        counts="{workdir}/output/b3_cache/{cancer_type}_counts.parquet",
        meta="{workdir}/output/b3_cache/{cancer_type}_meta.parquet"
    output:
        results="{workdir}/output/b1_results/{cancer_type}_deg_comparison.tsv"
    params:
        cancer_type="{cancer_type}"
    resources:
        mem_mb=16000,
        runtime=30
    threads: 1
    log:
        "{workdir}/logs/deg_{cancer_type}.log"
    shell:
        """
        python -c "
import pandas as pd
from projects.b1_deg_tool_agreement.data_access import CANCER_TO_TISSUE
from projects.b1_deg_tool_agreement.analysis import (
    _filter_low_counts, run_deseq2, run_edger, run_limma_voom, compare_tools
)
counts = pd.read_parquet('{input.counts}')
meta = pd.read_parquet('{input.meta}')
filtered = _filter_low_counts(counts)
deseq2_res = run_deseq2(filtered, meta['sample_type'])
edger_res = run_edger(filtered, meta['sample_type'])
limma_res = run_limma_voom(filtered, meta['sample_type'])
tissue = CANCER_TO_TISSUE.get('{params.cancer_type}', '')
results = compare_tools(deseq2_res, edger_res, limma_res,
                        '{params.cancer_type}', '{params.cancer_type}'.replace('TCGA-',''),
                        tissue, len(meta))
pd.DataFrame(results).to_csv('{output.results}', sep='\\t', index=False)
" > {log} 2>&1
        """

rule aggregate_results:
    input:
        expand("{{workdir}}/output/b1_results/{ct}_deg_comparison.tsv",
               ct=CANCER_TYPES)
    output:
        "{workdir}/output/b1_deg_tool_agreement.tsv"
    shell:
        """
        head -1 {input[0]} > {output}
        for f in {input}; do tail -n +2 "$f" >> {output}; done
        """

rule repro_verify:
    input:
        results="{workdir}/output/b1_deg_tool_agreement.tsv",
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
    write_snakefile_template(smk, "deg_comparison.smk", "B1_deg_tool_agreement")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    all_results = []
    start_time = time.time()

    for i, project_id in enumerate(CANCER_TYPES):
        cancer_type = project_id.replace("TCGA-", "")
        tissue = CANCER_TO_TISSUE.get(project_id, "")
        print(f"\n{'='*60}")
        print(f"[{i+1}/16] {project_id} (tissue: {tissue})")
        print(f"{'='*60}")

        # Load cached data from B3
        try:
            data = load_cancer_type(project_id)
        except FileNotFoundError as e:
            print(f"  SKIPPED: {e}")
            continue

        counts = data["counts"]
        sample_type = data["sample_type"]
        n_total = counts.shape[1]

        print(f"  Samples: {data['n_tumor']} tumor, {data['n_normal']} normal")
        print(f"  Raw genes: {counts.shape[0]}")

        # Filter low-count genes
        filtered = _filter_low_counts(counts)
        print(f"  After CPM>1 filter: {filtered.shape[0]} genes")

        # Run DESeq2
        print(f"  Running DESeq2...")
        try:
            deseq2_res = run_deseq2(filtered, sample_type)
            n_deseq2 = (deseq2_res["padj"] < 0.05).sum()
            print(f"    DESeq2 DEGs: {n_deseq2}")
        except Exception as e:
            print(f"    DESeq2 FAILED: {e}")
            continue

        # Run edgeR
        print(f"  Running edgeR...")
        try:
            edger_res = run_edger(filtered, sample_type)
            n_edger = (edger_res["padj"] < 0.05).sum()
            print(f"    edgeR DEGs: {n_edger}")
        except Exception as e:
            print(f"    edgeR FAILED: {e}")
            continue

        # Run limma-voom
        print(f"  Running limma-voom...")
        try:
            limma_res = run_limma_voom(filtered, sample_type)
            n_limma = (limma_res["padj"] < 0.05).sum()
            print(f"    limma DEGs: {n_limma}")
        except Exception as e:
            print(f"    limma FAILED: {e}")
            continue

        # Compare
        print(f"  Comparing tools...")
        results = compare_tools(
            deseq2_res, edger_res, limma_res,
            project_id, cancer_type, tissue, n_total,
        )
        all_results.extend(results)

        for r in results:
            print(f"    {r['tool_pair']}: Jaccard={r['jaccard_top100']:.4f}, "
                  f"DirAgree={r['direction_agreement']}, LFC_r={r['log2fc_correlation']}")

    # Save results
    results_df = pd.DataFrame(all_results)
    output_file = OUTPUT_DIR / "b1_deg_tool_agreement.tsv"
    results_df.to_csv(output_file, sep="\t", index=False)
    print(f"\n{'='*60}")
    print(f"Results saved to {output_file}")
    print(f"Total rows: {len(results_df)}")

    # Emit BioOrchestrator artifacts
    print(f"\nEmitting BioOrchestrator artifacts...")
    emit_knowledge(results_df)

    elapsed = time.time() - start_time
    print(f"\nB1 completed in {elapsed/60:.1f} minutes")

    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    for pair in results_df["tool_pair"].unique():
        subset = results_df[results_df["tool_pair"] == pair]
        print(f"\n{pair}:")
        print(f"  Jaccard top100: {subset['jaccard_top100'].mean():.4f} "
              f"({subset['jaccard_top100'].min():.4f}-{subset['jaccard_top100'].max():.4f})")
        print(f"  Direction agreement: {subset['direction_agreement'].mean():.4f}")
        print(f"  Log2FC correlation: {subset['log2fc_correlation'].mean():.4f}")

    # DEG counts
    print(f"\nMean DEGs (FDR<0.05):")
    print(f"  DESeq2: {results_df['n_degs_deseq2'].mean():.0f}")
    print(f"  edgeR:  {results_df['n_degs_edger'].mean():.0f}")
    print(f"  limma:  {results_df['n_degs_limma'].mean():.0f}")


if __name__ == "__main__":
    main()
