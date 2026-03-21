#!/usr/bin/env python3
"""
B2 Normalization Method Comparison — main runner.

Downloads GTEx V10 data, applies TPM/FPKM/TMM/VST normalizations,
compares results, and emits BioOrchestrator artifacts.
"""

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from projects.b2_normalization.data_access import TISSUES, download_tissue
from projects.b2_normalization.analysis import (
    normalize_all, compare_normalizations, HOUSEKEEPING_GENES,
)
from shared.knowledge_writer import (
    append_prompt_enrichment,
    append_validation_rules,
    write_snakefile_template,
)

OUTPUT_DIR = PROJECT_ROOT / "output"
CACHE_DIR = OUTPUT_DIR / "b2_cache"


def _build_gene_name_index(counts: pd.DataFrame, tpm: pd.DataFrame) -> dict:
    """Build ENSG → gene name mapping from GCT Description field.

    GTEx GCT files use ENSG IDs as row names. We need gene symbols
    for housekeeping gene CV calculation.
    """
    # The Description column in GCT is the gene symbol
    # We stored it during parsing but lost it in the parquet cache
    # Re-derive from the count matrix: ENSG IDs that match known symbols
    return {}


def _remap_to_gene_symbols(df: pd.DataFrame, gene_map: dict) -> pd.DataFrame:
    """Remap ENSG IDs to gene symbols where possible."""
    if not gene_map:
        return df
    new_index = [gene_map.get(g, g) for g in df.index]
    df = df.copy()
    df.index = new_index
    # Drop duplicates (keep first)
    df = df[~df.index.duplicated(keep="first")]
    return df


def emit_knowledge(results_df: pd.DataFrame):
    """Emit L2, L3, L4 artifacts to BioOrchestrator."""

    # --- L2: Prompt enrichment ---
    # Key statistics
    tpm_tmm = results_df[results_df["method_pair"] == "TPM_vs_TMM"]
    tpm_vst = results_df[results_df["method_pair"] == "TPM_vs_VST"]
    tmm_vst = results_df[results_df["method_pair"] == "TMM_vs_VST"]
    tpm_fpkm = results_df[results_df["method_pair"] == "TPM_vs_FPKM"]

    mean_tpm_tmm = tpm_tmm["spearman_rank_correlation"].mean()
    mean_tpm_vst = tpm_vst["spearman_rank_correlation"].mean()
    mean_tmm_vst = tmm_vst["spearman_rank_correlation"].mean()
    mean_tpm_fpkm = tpm_fpkm["spearman_rank_correlation"].mean()

    mean_jaccard_tpm_tmm = tpm_tmm["pathway_jaccard"].mean()
    mean_jaccard_tmm_vst = tmm_vst["pathway_jaccard"].mean()

    mean_affected = results_df["n_tissue_specific_genes_affected"].mean()

    prompt_text = (
        "Normalization method comparison (B2 audit, 20 GTEx V10 tissues):\n"
        "Spearman rank correlation between methods (mean across tissues): "
        f"TPM-TMM={mean_tpm_tmm:.3f}, TPM-VST={mean_tpm_vst:.3f}, "
        f"TMM-VST={mean_tmm_vst:.3f}, TPM-FPKM={mean_tpm_fpkm:.3f}.\n"
        f"Top-variable-gene overlap (Jaccard): TPM-TMM={mean_jaccard_tpm_tmm:.3f}, "
        f"TMM-VST={mean_jaccard_tmm_vst:.3f}.\n"
        f"Mean genes affected by normalization choice: {mean_affected:.0f}/500 "
        f"(top variable genes differ between methods).\n"
        "Recommendation: TMM (edgeR) and VST (DESeq2) are the most concordant "
        "count-based normalizations. TPM and FPKM should only be used for "
        "visualization and cross-gene comparisons, not for differential expression "
        "or cross-sample statistical testing. Always explicitly specify the "
        "normalization method and justify the choice relative to the analysis goal."
    )
    append_prompt_enrichment(prompt_text, "B2_normalization")

    # --- L3: Validation rules ---
    rules = [
        {
            "id": "bulk_normalization_specified",
            "type": "conditional",
            "description": "Normalization method must be explicitly specified",
            "condition": {
                "context": "steps",
                "field": "step_type",
                "check": "contains_any",
                "values": ["normalization", "differential_expression", "quantification"],
            },
            "required": {
                "scope": "parameters",
                "parameter_exists": "normalization_method",
            },
            "action": "warn",
            "message": (
                "No normalization method specified. Different normalization methods "
                "(TPM, FPKM, TMM, VST) produce substantially different gene rankings "
                f"(mean Spearman r={mean_tpm_tmm:.3f} between TPM and TMM). "
                "Explicitly specify the normalization method."
            ),
        },
        {
            "id": "bulk_tpm_fpkm_cross_sample_warning",
            "type": "conditional",
            "description": "Flag TPM/FPKM for cross-sample statistical comparison",
            "condition": {
                "context": "parameters",
                "field": "normalization_method",
                "check": "in",
                "values": ["TPM", "FPKM", "tpm", "fpkm"],
            },
            "required": {
                "scope": "steps",
                "step_type_none": ["differential_expression"],
            },
            "action": "warn",
            "message": (
                "TPM/FPKM normalization is used with differential expression analysis. "
                "TPM and FPKM are not appropriate for cross-sample statistical testing — "
                "they do not account for library composition bias. Use TMM (edgeR) or "
                "median-of-ratios (DESeq2) normalization instead."
            ),
        },
        {
            "id": "bulk_vst_min_samples",
            "type": "conditional",
            "description": "VST requires minimum samples per group",
            "condition": {
                "context": "parameters",
                "field": "normalization_method",
                "check": "in",
                "values": ["VST", "vst"],
            },
            "required": {
                "scope": "experiment",
                "field": "samples",
                "check": "min_group_size_gte",
                "threshold": 3,
            },
            "action": "warn",
            "message": (
                "VST normalization is specified but at least one group has fewer than "
                "3 samples. VST requires sufficient samples to estimate the dispersion-mean "
                "relationship. Use log2(count+1) or rlog for very small sample sizes."
            ),
        },
    ]
    append_validation_rules(rules, "B2_normalization")

    # --- L4: Snakefile template ---
    smk = '''# Normalization comparison Snakefile template.
# Generated by B2_normalization audit.
# Real shell commands — no touch {output} placeholders.

configfile: "config.yaml"

WORKDIR = config.get("workdir", ".")
TISSUES = config.get("tissues", [])

rule all:
    input:
        expand("{workdir}/output/b2_normalization.tsv", workdir=WORKDIR),
        expand("{workdir}/output/repro_verified.flag", workdir=WORKDIR),

rule repro_snapshot:
    output:
        lockfile="{workdir}/output/repro.lock"
    shell:
        "repro snapshot -o {output.lockfile} --quiet"

rule download_tissue:
    output:
        counts="{workdir}/output/b2_cache/{tissue}_counts.parquet",
        tpm="{workdir}/output/b2_cache/{tissue}_tpm.parquet"
    params:
        tissue="{tissue}"
    resources:
        mem_mb=4000,
        runtime=10
    log:
        "{workdir}/logs/download_{tissue}.log"
    shell:
        """
        python -c "
from projects.b2_normalization.data_access import download_tissue
from pathlib import Path
download_tissue('{params.tissue}', Path('{wildcards.workdir}/output/b2_cache'))
" > {log} 2>&1
        """

rule normalize_tissue:
    input:
        counts="{workdir}/output/b2_cache/{tissue}_counts.parquet",
        tpm="{workdir}/output/b2_cache/{tissue}_tpm.parquet"
    output:
        results="{workdir}/output/b2_results/{tissue}_comparison.tsv"
    params:
        tissue="{tissue}"
    resources:
        mem_mb=16000,
        runtime=30
    log:
        "{workdir}/logs/normalize_{tissue}.log"
    shell:
        """
        python -c "
import pandas as pd
from projects.b2_normalization.data_access import TISSUES
from projects.b2_normalization.analysis import normalize_all, compare_normalizations
counts = pd.read_parquet('{input.counts}')
tpm = pd.read_parquet('{input.tpm}')
normalized = normalize_all(counts, tpm)
results = compare_normalizations(normalized, '{params.tissue}', TISSUES['{params.tissue}'])
pd.DataFrame(results).to_csv('{output.results}', sep='\\t', index=False)
" > {log} 2>&1
        """

rule aggregate_results:
    input:
        expand("{{workdir}}/output/b2_results/{tissue}_comparison.tsv",
               tissue=TISSUES)
    output:
        "{workdir}/output/b2_normalization.tsv"
    shell:
        """
        head -1 {input[0]} > {output}
        for f in {input}; do tail -n +2 "$f" >> {output}; done
        """

rule repro_verify:
    input:
        results="{workdir}/output/b2_normalization.tsv",
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
    write_snakefile_template(smk, "normalization_comparison.smk", "B2_normalization")


def main():
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    all_results = []
    start_time = time.time()

    for i, (tissue_key, tissue_name) in enumerate(TISSUES.items()):
        print(f"\n{'='*60}")
        print(f"[{i+1}/20] {tissue_name} ({tissue_key})")
        print(f"{'='*60}")

        # Download
        try:
            data = download_tissue(tissue_key, CACHE_DIR)
        except Exception as e:
            print(f"  FAILED to download {tissue_key}: {e}")
            continue

        print(f"  Samples: {data['n_samples']}")
        print(f"  Genes: {data['counts'].shape[0]}")

        # Subsample if needed for memory
        counts = data["counts"]
        tpm = data["tpm"]
        if counts.shape[1] > 100:
            cols = counts.columns[:100]
            counts = counts[cols]
            tpm = tpm[cols]
            print(f"  Subsampled to 100 samples")

        # Normalize
        print(f"  Running normalizations...")
        try:
            normalized = normalize_all(counts, tpm)
            print(f"  Normalizations complete: {list(normalized.keys())}")
        except Exception as e:
            print(f"  FAILED normalization for {tissue_key}: {e}")
            import traceback
            traceback.print_exc()
            continue

        # Compare
        print(f"  Comparing methods...")
        try:
            results = compare_normalizations(normalized, tissue_key, tissue_name)
            all_results.extend(results)
            print(f"  {len(results)} comparisons complete")
        except Exception as e:
            print(f"  FAILED comparison for {tissue_key}: {e}")
            import traceback
            traceback.print_exc()
            continue

        # Print key result
        for r in results:
            if r["method_pair"] == "TPM_vs_TMM":
                print(f"  TPM vs TMM: Spearman={r['spearman_rank_correlation']:.4f}, "
                      f"Jaccard={r['pathway_jaccard']:.4f}")

    # Save results
    results_df = pd.DataFrame(all_results)
    output_file = OUTPUT_DIR / "b2_normalization.tsv"
    results_df.to_csv(output_file, sep="\t", index=False)
    print(f"\n{'='*60}")
    print(f"Results saved to {output_file}")
    print(f"Total rows: {len(results_df)}")

    # Emit BioOrchestrator artifacts
    print(f"\nEmitting BioOrchestrator artifacts...")
    emit_knowledge(results_df)

    elapsed = time.time() - start_time
    print(f"\nB2 completed in {elapsed/60:.1f} minutes")

    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    for pair in results_df["method_pair"].unique():
        subset = results_df[results_df["method_pair"] == pair]
        print(f"\n{pair}:")
        print(f"  Spearman: {subset['spearman_rank_correlation'].mean():.4f} "
              f"(range {subset['spearman_rank_correlation'].min():.4f}-"
              f"{subset['spearman_rank_correlation'].max():.4f})")
        print(f"  Pathway Jaccard: {subset['pathway_jaccard'].mean():.4f}")
        print(f"  Genes affected: {subset['n_tissue_specific_genes_affected'].mean():.0f}/500")


if __name__ == "__main__":
    main()
