#!/usr/bin/env python3
"""P9 — Clustering Resolution Sensitivity Audit.

Question: Is the default Leiden resolution of 0.5 empirically justified,
and which resolution best recovers known cell type structure?

Strategy:
- Query Census for 20+ datasets (diverse tissues, human + mouse, 10 per organism)
- For each dataset, fetch sparse counts (cap 30K cells), log-normalize, PCA, kNN, Leiden
- Sweep resolutions [0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
- Evaluate each resolution by ARI vs curated cell_type labels
- Spearman correlation: n_cell_types vs optimal_resolution
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
import scanpy as sc
import tiledbsoma
from scipy import stats as scipy_stats
from scipy.sparse import issparse
from sklearn.metrics import adjusted_rand_score
from tqdm import tqdm

warnings.filterwarnings("ignore", category=FutureWarning)

CENSUS_VERSION = "2025-11-08"
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = REPO_ROOT / "scrnaseq" / "output"
PROJECT_DIR = Path(__file__).resolve().parent

MAX_DATASETS_PER_ORGANISM = 10
MAX_CELLS_PER_DATASET = 30_000
MIN_CELLS_PER_DATASET = 500
MIN_CELL_TYPES = 3
RESOLUTIONS = [0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
N_PCS = 50
RNG_SEED = 42

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(PROJECT_DIR / "p9_run.log"),
    ],
)
logger = logging.getLogger(__name__)


def _census_key(organism: str) -> str:
    return organism.lower().replace(" ", "_")


def discover_diverse_datasets(census, organism: str) -> pd.DataFrame:
    """Find diverse dataset-tissue combos with >=MIN_CELL_TYPES curated cell types."""
    exp = census["census_data"][_census_key(organism)]
    logger.info("Discovering datasets for %s...", organism)

    obs_df = exp.obs.read(
        column_names=["dataset_id", "tissue_general", "cell_type"],
        value_filter="is_primary_data == True",
    ).concat().to_pandas()

    # Drop rows with NaN dataset_id
    obs_df = obs_df.dropna(subset=["dataset_id"])

    # Count cells and distinct cell types per dataset-tissue combo
    grouped = obs_df.groupby(["dataset_id", "tissue_general"], observed=True).agg(
        n_cells=("cell_type", "size"),
        n_cell_types=("cell_type", "nunique"),
    ).reset_index()

    grouped = grouped.rename(columns={"tissue_general": "tissue"})

    # Filter: enough cells and enough cell type diversity
    grouped = grouped[
        (grouped.n_cells >= MIN_CELLS_PER_DATASET)
        & (grouped.n_cell_types >= MIN_CELL_TYPES)
    ]

    # Tissue diversity: 1 per tissue, prefer datasets with more cell types
    grouped = grouped.sort_values(
        ["n_cell_types", "n_cells"], ascending=[False, False],
    )
    selected = grouped.drop_duplicates(subset="tissue", keep="first")
    selected = selected.head(MAX_DATASETS_PER_ORGANISM)

    logger.info(
        "%s: selected %d datasets across %d tissues (cell types range: %d-%d)",
        organism,
        len(selected),
        selected.tissue.nunique(),
        selected.n_cell_types.min() if len(selected) > 0 else 0,
        selected.n_cell_types.max() if len(selected) > 0 else 0,
    )
    del obs_df
    gc.collect()
    return selected.reset_index(drop=True)


def process_dataset(
    census, dataset_id: str, tissue: str, organism: str,
) -> list[dict] | None:
    """Fetch counts, log-normalize, PCA, kNN, Leiden sweep, ARI vs cell_type."""
    exp = census["census_data"][_census_key(organism)]

    obs_filter = (
        f"dataset_id == '{dataset_id}' and tissue_general == '{tissue}' "
        f"and is_primary_data == True"
    )

    # Get obs joinids with cell_type labels
    obs_df = exp.obs.read(
        column_names=["soma_joinid", "cell_type"],
        value_filter=obs_filter,
    ).concat().to_pandas()

    # Drop rows with missing cell_type
    obs_df = obs_df.dropna(subset=["cell_type"])

    if len(obs_df) < MIN_CELLS_PER_DATASET:
        return None

    n_true_types = obs_df["cell_type"].nunique()
    if n_true_types < MIN_CELL_TYPES:
        return None

    # Subsample if too large
    if len(obs_df) > MAX_CELLS_PER_DATASET:
        rng = np.random.default_rng(RNG_SEED)
        idx = rng.choice(len(obs_df), MAX_CELLS_PER_DATASET, replace=False)
        obs_df = obs_df.iloc[sorted(idx)].copy()

    obs_ids = obs_df.soma_joinid.tolist()
    cell_type_labels = obs_df["cell_type"].values

    # Fetch expression matrix
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

    # Attach ground truth labels (aligned by obs order)
    adata.obs["cell_type_true"] = cell_type_labels

    # Standard preprocessing: filter genes, log-normalize, HVG, PCA, neighbors
    sc.pp.filter_genes(adata, min_cells=3)
    if adata.n_vars < 100:
        logger.info("  Too few genes after filtering, skipping")
        return None

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # HVG selection
    n_hvg = min(2000, adata.n_vars - 1)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)
    adata = adata[:, adata.var.highly_variable].copy()

    # PCA
    n_pcs = min(N_PCS, adata.n_vars - 1, adata.n_obs - 1)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs, random_state=RNG_SEED)

    # kNN graph
    sc.pp.neighbors(adata, n_pcs=n_pcs, random_state=RNG_SEED)

    # Ground truth for ARI
    true_labels = adata.obs["cell_type_true"].values
    n_true = len(np.unique(true_labels))

    # Leiden sweep
    results = []
    best_ari = -1.0
    best_res = None

    for res in RESOLUTIONS:
        sc.tl.leiden(adata, resolution=res, random_state=RNG_SEED, key_added="leiden")
        pred_labels = adata.obs["leiden"].values
        n_clusters = len(np.unique(pred_labels))

        ari = adjusted_rand_score(true_labels, pred_labels)

        over_clustering = n_clusters / n_true if n_clusters > n_true else 1.0
        under_clustering = n_true / n_clusters if n_clusters < n_true else 1.0

        if ari > best_ari:
            best_ari = ari
            best_res = res

        results.append({
            "dataset_id": dataset_id,
            "tissue": tissue,
            "organism": organism,
            "resolution": res,
            "ari": round(ari, 4),
            "n_clusters": n_clusters,
            "n_true_types": n_true,
            "over_clustering": round(over_clustering, 3),
            "under_clustering": round(under_clustering, 3),
            "optimal_flag": False,  # will set below
        })

    # Mark optimal
    for r in results:
        if r["resolution"] == best_res:
            r["optimal_flag"] = True

    return results


def run_statistics(df: pd.DataFrame) -> None:
    """Spearman correlation: n_cell_types vs optimal_resolution."""
    optimal = df[df.optimal_flag].copy()

    if len(optimal) < 5:
        logger.warning("Not enough datasets for correlation analysis")
        return

    rho, pval = scipy_stats.spearmanr(optimal.n_true_types, optimal.resolution)
    logger.info(
        "Spearman correlation (n_cell_types vs optimal_resolution): "
        "rho=%.3f, p=%.3e (n=%d)",
        rho, pval, len(optimal),
    )

    stats_path = PROJECT_DIR / "spearman_results.tsv"
    pd.DataFrame([{
        "test": "Spearman",
        "variable_x": "n_true_types",
        "variable_y": "optimal_resolution",
        "rho": rho,
        "p": pval,
        "n": len(optimal),
    }]).to_csv(stats_path, sep="\t", index=False)
    logger.info("Spearman results saved to %s", stats_path)


def main():
    t0 = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PROJECT_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("=== P9: Clustering Resolution Sensitivity Audit ===")

    # Checkpoint
    checkpoint_path = PROJECT_DIR / "p9_checkpoint.tsv"
    if checkpoint_path.exists():
        existing = pd.read_csv(checkpoint_path, sep="\t")
        all_results = existing.to_dict("records")
        done_keys = {
            (r["dataset_id"], r["tissue"], r["organism"]) for r in all_results
        }
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
            desc=f"P9 ({organism})",
        ):
            if pd.isna(row.dataset_id):
                continue

            key = (row.dataset_id, row.tissue, organism)
            if key in done_keys:
                logger.info("  Skipping %s/%s (done)", row.dataset_id[:8], row.tissue)
                continue

            try:
                results = process_dataset(
                    census, row.dataset_id, row.tissue, organism,
                )
            except Exception as e:
                logger.warning(
                    "  Error on %s/%s: %s", row.dataset_id[:8], row.tissue, e,
                )
                continue

            if results:
                all_results.extend(results)
                done_keys.add(key)
                # Log summary for default resolution 0.5
                r05 = [r for r in results if r["resolution"] == 0.5]
                r_best = [r for r in results if r["optimal_flag"]]
                if r05 and r_best:
                    logger.info(
                        "  %s/%s: %d true types, "
                        "res=0.5 ARI=%.3f (%d clusters), "
                        "best res=%.1f ARI=%.3f (%d clusters)",
                        row.dataset_id[:8], row.tissue,
                        r05[0]["n_true_types"],
                        r05[0]["ari"], r05[0]["n_clusters"],
                        r_best[0]["resolution"],
                        r_best[0]["ari"], r_best[0]["n_clusters"],
                    )
                # Checkpoint
                pd.DataFrame(all_results).to_csv(
                    checkpoint_path, sep="\t", index=False,
                )
                gc.collect()

    census.close()

    result_df = pd.DataFrame(all_results)
    out_path = OUTPUT_DIR / "p9_clustering_resolution.tsv"
    result_df.to_csv(out_path, sep="\t", index=False)
    logger.info("\nSaved %d rows to %s", len(result_df), out_path)

    # Statistics
    if len(result_df) >= 10:
        run_statistics(result_df)

    # Summary
    elapsed = time.time() - t0
    n_datasets = result_df.groupby(
        ["dataset_id", "tissue", "organism"], observed=True,
    ).ngroups
    logger.info("\n=== P9 COMPLETE ===")
    logger.info("Dataset-tissue combos: %d", n_datasets)
    logger.info(
        "Total rows (datasets x resolutions): %d", len(result_df),
    )

    # Key finding: how does res=0.5 compare?
    for res in RESOLUTIONS:
        subset = result_df[result_df.resolution == res]
        if len(subset) > 0:
            logger.info(
                "  Resolution=%.1f: mean ARI=%.3f, mean clusters=%d, "
                "optimal in %d/%d datasets",
                res,
                subset.ari.mean(),
                int(subset.n_clusters.mean()),
                int(subset.optimal_flag.sum()),
                n_datasets,
            )

    # Overall: is 0.5 optimal?
    optimal = result_df[result_df.optimal_flag]
    if len(optimal) > 0:
        logger.info(
            "\nOptimal resolution distribution: mean=%.2f, median=%.2f, "
            "std=%.2f, range=[%.1f, %.1f]",
            optimal.resolution.mean(),
            optimal.resolution.median(),
            optimal.resolution.std(),
            optimal.resolution.min(),
            optimal.resolution.max(),
        )
        pct_05 = (optimal.resolution == 0.5).sum() / len(optimal) * 100
        logger.info(
            "Resolution 0.5 was optimal in %.1f%% of datasets", pct_05,
        )

    logger.info("Elapsed: %.1f minutes", elapsed / 60)


if __name__ == "__main__":
    main()
