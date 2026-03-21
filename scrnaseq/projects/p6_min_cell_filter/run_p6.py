#!/usr/bin/env python3
"""P6 — Minimum Cell Filter Threshold Audit.

Question: Is the default "keep genes expressed in at least 3 cells" threshold
empirically justified, and how much does it affect downstream results?

Strategy:
- Query Census for 30+ datasets (diverse tissues, human + mouse)
- For each dataset, fetch sparse counts and compute gene-level stats
- Apply filter thresholds [1, 3, 5, 10, 20, 50]: genes retained, HVG overlap
- Spearman correlation: dataset size vs optimal threshold
"""

from __future__ import annotations

import gc
import logging
import sys
import time
import warnings
from pathlib import Path

import cellxgene_census
import numpy as np
import pandas as pd
import tiledbsoma
from scipy import stats as scipy_stats
from scipy.sparse import issparse
from tqdm import tqdm

warnings.filterwarnings("ignore", category=FutureWarning)

CENSUS_VERSION = "2025-11-08"
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = REPO_ROOT / "scrnaseq" / "output"
PROJECT_DIR = Path(__file__).resolve().parent

MAX_DATASETS_PER_ORGANISM = 15
MAX_CELLS_PER_DATASET = 30_000
MIN_CELLS_PER_DATASET = 500
THRESHOLDS = [1, 3, 5, 10, 20, 50]
N_HVG = 2000
RNG_SEED = 42

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(PROJECT_DIR / "p6_run.log"),
    ],
)
logger = logging.getLogger(__name__)


def _census_key(organism: str) -> str:
    return organism.lower().replace(" ", "_")


def discover_diverse_datasets(census, organism: str) -> pd.DataFrame:
    """Find diverse dataset-tissue combos."""
    exp = census["census_data"][_census_key(organism)]
    logger.info("Discovering datasets for %s...", organism)

    obs_df = exp.obs.read(
        column_names=["dataset_id", "tissue_general"],
        value_filter="is_primary_data == True",
    ).concat().to_pandas()

    grouped = obs_df.groupby(["dataset_id", "tissue_general"], observed=True).size().reset_index(name="n_cells")
    grouped = grouped[grouped.n_cells >= MIN_CELLS_PER_DATASET]
    grouped = grouped.rename(columns={"tissue_general": "tissue"})

    # Tissue diversity: 1 per tissue, largest first
    grouped = grouped.sort_values("n_cells", ascending=False)
    selected = grouped.drop_duplicates(subset="tissue", keep="first")
    selected = selected.head(MAX_DATASETS_PER_ORGANISM)

    logger.info(
        "%s: selected %d datasets across %d tissues",
        organism, len(selected), selected.tissue.nunique(),
    )
    del obs_df
    gc.collect()
    return selected.reset_index(drop=True)


def compute_hvgs(X_sparse, n_top: int = N_HVG) -> set[int]:
    """Compute top HVGs by normalized dispersion on sparse matrix.

    Returns set of gene indices (column indices).
    """
    # Mean and variance per gene
    if issparse(X_sparse):
        mean = np.asarray(X_sparse.mean(axis=0)).ravel()
        # Var = E[X^2] - (E[X])^2
        mean_sq = np.asarray(X_sparse.power(2).mean(axis=0)).ravel()
    else:
        mean = X_sparse.mean(axis=0)
        mean_sq = (X_sparse ** 2).mean(axis=0)

    var = mean_sq - mean ** 2
    var = np.maximum(var, 0)

    # Normalized dispersion (Seurat-style)
    # Avoid division by zero
    disp = np.where(mean > 0, var / mean, 0)
    # Log transform for ranking
    log_disp = np.log1p(disp)

    # Top n_top by dispersion (among expressed genes)
    expressed = mean > 0
    indices = np.where(expressed)[0]
    if len(indices) <= n_top:
        return set(indices)

    top_idx = indices[np.argsort(log_disp[indices])[-n_top:]]
    return set(top_idx.tolist())


def process_dataset(
    census, dataset_id: str, tissue: str, organism: str,
) -> list[dict] | None:
    """Fetch counts, apply filter thresholds, compute metrics."""
    exp = census["census_data"][_census_key(organism)]

    obs_filter = (
        f"dataset_id == '{dataset_id}' and tissue_general == '{tissue}' "
        f"and is_primary_data == True"
    )

    # Get obs joinids
    obs_df = exp.obs.read(
        column_names=["soma_joinid"],
        value_filter=obs_filter,
    ).concat().to_pandas()

    if len(obs_df) < MIN_CELLS_PER_DATASET:
        return None

    # Subsample if too large
    if len(obs_df) > MAX_CELLS_PER_DATASET:
        rng = np.random.default_rng(RNG_SEED)
        idx = rng.choice(len(obs_df), MAX_CELLS_PER_DATASET, replace=False)
        obs_df = obs_df.iloc[sorted(idx)].copy()

    obs_ids = obs_df.soma_joinid.tolist()

    # Fetch full expression (sparse)
    try:
        with exp.axis_query(
            measurement_name="RNA",
            obs_query=tiledbsoma.AxisQuery(coords=(obs_ids,)),
        ) as query:
            adata = query.to_anndata(X_name="raw")
    except Exception as e:
        logger.warning("Fetch failed for %s/%s: %s", dataset_id[:8], tissue, e)
        return None

    if adata.shape[0] == 0:
        return None

    X = adata.X
    n_cells_total = X.shape[0]
    n_genes_total = X.shape[1]

    # Compute cells-per-gene (how many cells express each gene)
    if issparse(X):
        cells_per_gene = np.asarray((X > 0).sum(axis=0)).ravel()
    else:
        cells_per_gene = (X > 0).sum(axis=0)

    # Compute genes-per-cell (before filtering)
    if issparse(X):
        genes_per_cell_raw = np.asarray((X > 0).sum(axis=1)).ravel()
    else:
        genes_per_cell_raw = (X > 0).sum(axis=1)

    # Baseline HVGs at threshold=3
    mask_3 = cells_per_gene >= 3
    if mask_3.sum() < 100:
        logger.info("  Too few genes pass threshold=3, skipping")
        return None

    if issparse(X):
        X_3 = X[:, mask_3]
    else:
        X_3 = X[:, mask_3]
    hvg_baseline = compute_hvgs(X_3)
    # Map back to original indices
    orig_indices_3 = np.where(mask_3)[0]
    hvg_baseline_orig = {orig_indices_3[i] for i in hvg_baseline}

    results = []
    for threshold in THRESHOLDS:
        mask = cells_per_gene >= threshold
        genes_retained = int(mask.sum())

        if genes_retained < 50:
            continue

        # Filter matrix
        if issparse(X):
            X_filt = X[:, mask]
            genes_per_cell = np.asarray((X_filt > 0).sum(axis=1)).ravel()
        else:
            X_filt = X[:, mask]
            genes_per_cell = (X_filt > 0).sum(axis=1)

        # Cells with 0 genes after filter
        cells_retained = int((genes_per_cell > 0).sum())

        # HVGs at this threshold
        hvg_this = compute_hvgs(X_filt)
        # Map to original indices for comparison
        orig_indices = np.where(mask)[0]
        hvg_this_orig = {orig_indices[i] for i in hvg_this}

        # Overlap with baseline (threshold=3)
        if hvg_baseline_orig:
            intersection = hvg_this_orig & hvg_baseline_orig
            union = hvg_this_orig | hvg_baseline_orig
            hvg_overlap = len(intersection) / len(union) if union else 1.0
            pct_hvg_lost = 1.0 - len(intersection) / len(hvg_baseline_orig) if hvg_baseline_orig else 0.0
        else:
            hvg_overlap = 1.0
            pct_hvg_lost = 0.0

        results.append({
            "dataset_id": dataset_id,
            "tissue": tissue,
            "organism": organism,
            "dataset_size_cells": n_cells_total,
            "threshold": threshold,
            "genes_total": n_genes_total,
            "genes_retained": genes_retained,
            "cells_retained": cells_retained,
            "hvg_overlap": round(hvg_overlap, 4),
            "pct_hvg_lost": round(pct_hvg_lost * 100, 2),
            "median_genes_per_cell": float(np.median(genes_per_cell)),
        })

    return results


def run_statistics(df: pd.DataFrame) -> None:
    """Spearman correlation: dataset size vs optimal threshold."""
    # For each dataset, find threshold that maximizes HVG overlap while retaining >95% genes
    optimal = []
    for (ds_id, tissue, org), grp in df.groupby(["dataset_id", "tissue", "organism"], observed=True):
        ds_size = grp.dataset_size_cells.iloc[0]
        # Optimal = highest threshold where pct_hvg_lost < 5%
        viable = grp[grp.pct_hvg_lost < 5.0]
        if len(viable) > 0:
            best = viable.loc[viable.threshold.idxmax()]
            optimal.append({"dataset_size": ds_size, "optimal_threshold": best.threshold})

    if len(optimal) < 5:
        logger.warning("Not enough datasets for correlation analysis")
        return

    opt_df = pd.DataFrame(optimal)
    rho, pval = scipy_stats.spearmanr(opt_df.dataset_size, opt_df.optimal_threshold)
    logger.info(
        "Spearman correlation (dataset_size vs optimal_threshold): rho=%.3f, p=%.3e (n=%d)",
        rho, pval, len(opt_df),
    )

    stats_path = PROJECT_DIR / "spearman_results.tsv"
    pd.DataFrame([{"test": "Spearman", "rho": rho, "p": pval, "n": len(opt_df)}]).to_csv(
        stats_path, sep="\t", index=False,
    )


def main():
    t0 = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PROJECT_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("=== P6: Minimum Cell Filter Threshold Audit ===")

    # Checkpoint
    checkpoint_path = PROJECT_DIR / "p6_checkpoint.tsv"
    if checkpoint_path.exists():
        existing = pd.read_csv(checkpoint_path, sep="\t")
        all_results = existing.to_dict("records")
        done_keys = {(r["dataset_id"], r["tissue"], r["organism"]) for r in all_results}
        logger.info("Resuming from checkpoint: %d rows", len(all_results))
    else:
        all_results = []
        done_keys = set()

    census = cellxgene_census.open_soma(census_version=CENSUS_VERSION)

    for organism in ["Homo sapiens", "Mus musculus"]:
        logger.info("\n--- Processing %s ---", organism)

        datasets = discover_diverse_datasets(census, organism)
        logger.info("Will process %d datasets", len(datasets))

        for _, row in tqdm(
            datasets.iterrows(), total=len(datasets),
            desc=f"P6 ({organism})",
        ):
            if pd.isna(row.dataset_id):
                continue

            key = (row.dataset_id, row.tissue, organism)
            if key in done_keys:
                logger.info("  Skipping %s/%s (done)", row.dataset_id[:8], row.tissue)
                continue

            try:
                results = process_dataset(census, row.dataset_id, row.tissue, organism)
            except Exception as e:
                logger.warning("  Error on %s/%s: %s", row.dataset_id[:8], row.tissue, e)
                continue

            if results:
                all_results.extend(results)
                done_keys.add(key)
                # Log summary for threshold=3
                t3 = [r for r in results if r["threshold"] == 3]
                if t3:
                    logger.info(
                        "  %s/%s: %d cells, %d→%d genes @3cells, median_gpd=%.0f",
                        row.dataset_id[:8], row.tissue,
                        t3[0]["dataset_size_cells"], t3[0]["genes_total"],
                        t3[0]["genes_retained"], t3[0]["median_genes_per_cell"],
                    )
                # Checkpoint
                pd.DataFrame(all_results).to_csv(checkpoint_path, sep="\t", index=False)
                gc.collect()

    census.close()

    result_df = pd.DataFrame(all_results)
    out_path = OUTPUT_DIR / "p6_min_cell_filter.tsv"
    result_df.to_csv(out_path, sep="\t", index=False)
    logger.info("\nSaved %d rows to %s", len(result_df), out_path)

    # Statistics
    if len(result_df) >= 10:
        run_statistics(result_df)

    # Summary
    elapsed = time.time() - t0
    n_datasets = result_df.groupby(["dataset_id", "tissue", "organism"]).ngroups
    logger.info("\n=== P6 COMPLETE ===")
    logger.info("Dataset-tissue combos: %d", n_datasets)
    logger.info("Total rows (datasets × thresholds): %d", len(result_df))

    # Key finding: how much do HVGs change at different thresholds?
    for t in THRESHOLDS:
        subset = result_df[result_df.threshold == t]
        if len(subset) > 0:
            logger.info(
                "  Threshold=%d: mean genes retained=%.0f, mean HVG lost=%.1f%%",
                t, subset.genes_retained.mean(), subset.pct_hvg_lost.mean(),
            )

    logger.info("Elapsed: %.1f minutes", elapsed / 60)


if __name__ == "__main__":
    main()
