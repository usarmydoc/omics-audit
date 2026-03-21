"""Cell Type Marker Reliability Atlas — main entrypoint.

Orchestrates the full pipeline:
1. repro snapshot
2. Census connection + dataset discovery
3. Expression data fetching (sequential)
4. Marker scoring (parallel)
5. Score aggregation
6. TSV + markdown reporting
7. BioOrchestrator knowledge artifact generation
8. repro verify
"""

from __future__ import annotations

import logging
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import psutil

# Configure logging before imports that use it
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("logs/marker_atlas.log"),
    ],
)
logger = logging.getLogger(__name__)

from .census_access import (
    ALL_GENES,
    CENSUS_VERSION,
    MARKER_GENES,
    DatasetSlice,
    fetch_all_datasets,
    log_memory,
    open_census,
    discover_datasets,
    resolve_gene_ids,
)
from .scoring import AggregatedScore, MarkerScore, aggregate_scores, score_all_parallel
from .reporting import write_summary_markdown, write_tsv
from .knowledge_writer import (
    append_validation_rules,
    write_prompt_enrichment,
    write_snakefile_template,
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
LOGS_DIR = PROJECT_ROOT / "logs"


def run_repro_snapshot(output_dir: Path) -> bool:
    """Run repro snapshot to capture baseline environment."""
    lock_path = output_dir / "repro.lock"
    logger.info("Running repro snapshot -> %s", lock_path)
    result = subprocess.run(
        ["repro", "snapshot", "-o", str(lock_path), "--quiet"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        logger.warning("repro snapshot failed: %s", result.stderr)
        return False
    logger.info("repro snapshot complete: %s", lock_path)
    return True


def run_repro_verify(output_dir: Path) -> bool:
    """Run repro verify to hash outputs and verify reproducibility."""
    lock_path = output_dir / "repro.lock"
    logger.info("Running repro verify on %s", output_dir)
    result = subprocess.run(
        ["repro", "verify", str(output_dir), "-l", str(lock_path)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        logger.error("repro verify FAILED: %s", result.stderr)
        logger.error("stdout: %s", result.stdout)
        return False
    logger.info("repro verify passed")
    logger.info("repro.lock path: %s", lock_path)
    return True


def main() -> int:
    """Main pipeline orchestration."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    LOGS_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 70)
    logger.info("Cell Type Marker Reliability Atlas")
    logger.info("=" * 70)
    log_memory("startup")

    # --- Step 1: repro snapshot ---
    logger.info("Step 1/8: repro snapshot")
    run_repro_snapshot(OUTPUT_DIR)
    log_memory("after_snapshot")

    # --- Steps 2-5: Process each organism ---
    ORGANISMS = ["Homo sapiens", "Mus musculus"]

    all_dataset_slices: list[DatasetSlice] = []
    all_scores: list[MarkerScore] = []
    n_datasets_total = 0
    n_cells_total = 0

    census = open_census()
    try:
        for org_idx, organism in enumerate(ORGANISMS, 1):
            logger.info("=" * 50)
            logger.info("Organism %d/%d: %s", org_idx, len(ORGANISMS), organism)
            logger.info("=" * 50)

            # --- Resolve genes ---
            logger.info("Step 2: Resolving gene IDs for %s", organism)
            gene_df = resolve_gene_ids(census, organism)
            gene_joinids = gene_df.soma_joinid.tolist()
            gene_names = gene_df.feature_name.tolist()
            logger.info("Gene joinids resolved: %d genes", len(gene_names))
            log_memory(f"after_gene_resolve_{organism}")

            # --- Discover datasets ---
            logger.info("Step 3: Discovering datasets for %s", organism)
            dataset_info = discover_datasets(census, organism)
            n_dt = dataset_info[["dataset_id", "tissue"]].drop_duplicates().shape[0]
            n_datasets_total += n_dt
            logger.info("Found %d dataset-tissue combinations for %s", n_dt, organism)
            log_memory(f"after_discovery_{organism}")

            # --- Fetch expression data ---
            logger.info("Step 4: Fetching expression data for %s", organism)
            org_slices: list[DatasetSlice] = []
            for ds in fetch_all_datasets(census, dataset_info, gene_joinids, gene_names, organism):
                org_slices.append(ds)
                n_cells_total += ds.n_cells

            logger.info("Fetched %d datasets, %d cells for %s",
                        len(org_slices), sum(ds.n_cells for ds in org_slices), organism)
            all_dataset_slices.extend(org_slices)
            log_memory(f"after_fetch_{organism}")

            # --- Score markers (parallel) ---
            logger.info("Step 5: Scoring markers for %s", organism)
            target_types = list(MARKER_GENES.keys())
            work_items = [
                (ds.expr, ds.gene_names, ds.cell_types, target_types,
                 ds.dataset_id, ds.tissue, ds.organism)
                for ds in org_slices
            ]
            org_scores = score_all_parallel(work_items, max_workers=8)
            logger.info("Scored %d individual markers for %s", len(org_scores), organism)
            all_scores.extend(org_scores)
            log_memory(f"after_scoring_{organism}")

    finally:
        census.close()
        logger.info("Census connection closed")

    if not all_dataset_slices:
        logger.error("No datasets with usable data found. Exiting.")
        return 1

    # --- Step 6: Aggregate scores ---
    logger.info("Step 6/8: Aggregating scores across all organisms and datasets")
    aggregated: list[AggregatedScore] = aggregate_scores(all_scores)
    logger.info("Aggregated into %d gene-cell_type-organism scores", len(aggregated))

    # --- Step 7: Write outputs ---
    logger.info("Step 7/8: Writing outputs")

    tsv_path = write_tsv(aggregated, OUTPUT_DIR)

    run_metadata = {
        "date": datetime.now(timezone.utc).isoformat(),
        "census_version": CENSUS_VERSION,
        "organisms": ORGANISMS,
        "n_datasets": n_datasets_total,
        "n_datasets_with_data": len(all_dataset_slices),
        "n_cells_total": n_cells_total,
        "n_genes": len(ALL_GENES),
        "n_cell_types": len(set(s.cell_type for s in aggregated)),
    }

    md_path = write_summary_markdown(aggregated, all_scores, run_metadata, OUTPUT_DIR)
    log_memory("after_reporting")

    # --- BioOrchestrator knowledge artifacts ---
    logger.info("Writing BioOrchestrator knowledge artifacts")
    prompt_path = write_prompt_enrichment(aggregated)
    logger.info("Prompt enrichment: %s", prompt_path)

    rules_path = append_validation_rules()
    logger.info("Validation rules: %s", rules_path)

    smk_path = write_snakefile_template()
    logger.info("Snakefile template: %s", smk_path)
    log_memory("after_knowledge")

    # --- Step 8: repro verify ---
    logger.info("Step 8/8: repro verify")
    verified = run_repro_verify(OUTPUT_DIR)

    # --- Summary ---
    logger.info("=" * 70)
    logger.info("Pipeline complete!")
    logger.info("  Atlas TSV:      %s", tsv_path)
    logger.info("  Summary:        %s", md_path)
    logger.info("  Prompt:         %s", prompt_path)
    logger.info("  Rules:          %s", rules_path)
    logger.info("  Snakefile:      %s", smk_path)
    logger.info("  repro.lock:     %s", OUTPUT_DIR / "repro.lock")
    logger.info("  Verified:       %s", "YES" if verified else "NO")
    logger.info("  Organisms:      %s", ", ".join(ORGANISMS))
    logger.info("  Datasets:       %d queried, %d with data", n_datasets_total, len(all_dataset_slices))
    logger.info("  Cells scored:   %d", n_cells_total)
    logger.info("  Marker scores:  %d individual, %d aggregated", len(all_scores), len(aggregated))
    log_memory("final")
    logger.info("=" * 70)

    if not verified:
        logger.error("repro verify failed — exiting with non-zero code")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
