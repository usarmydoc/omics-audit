#!/usr/bin/env python3
"""P1 — Mitochondrial Threshold Arbitrariness Audit.

Question: Is the ubiquitous 20% mitochondrial read threshold empirically
justified or tissue-dependent?

Strategy:
- Census obs has `raw_sum` (total UMIs per cell) — no matrix pull needed for denominator.
- Query only MT-gene expression via axis_query restricted to MT gene soma_joinids.
- Sum MT UMIs per cell, divide by raw_sum → mito fraction.
- Sample diverse datasets across tissues (human + mouse).
- Aggregate per dataset-tissue: median, P95, cells removed at 10/20/50% thresholds.
- Kruskal-Wallis + Dunn post-hoc across tissues.
"""

from __future__ import annotations

import gc
import logging
import sys
import time
from pathlib import Path

import cellxgene_census
import numpy as np
import pandas as pd
import tiledbsoma
from scipy import stats
from tqdm import tqdm

# ── Config ──────────────────────────────────────────────────────────
CENSUS_VERSION = "2025-11-08"
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = REPO_ROOT / "scrnaseq" / "output"
PROJECT_DIR = Path(__file__).resolve().parent

MAX_DATASETS_PER_ORGANISM = 30
MAX_CELLS_PER_DATASET = 50_000
MIN_CELLS_PER_DATASET = 200
RNG_SEED = 42
THRESHOLDS = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.50]

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(PROJECT_DIR / "p1_run.log"),
    ],
)
logger = logging.getLogger(__name__)


# ── Census helpers ──────────────────────────────────────────────────

def _census_key(organism: str) -> str:
    return organism.lower().replace(" ", "_")


def get_mt_gene_ids(census, organism: str) -> list[int]:
    """Get soma_joinids for mitochondrial genes."""
    prefix = "MT-" if "homo" in organism.lower() else "mt-"
    exp = census["census_data"][_census_key(organism)]
    var_df = exp.ms["RNA"].var.read(
        column_names=["soma_joinid", "feature_name"],
    ).concat().to_pandas()
    mt = var_df[var_df.feature_name.str.startswith(prefix)]
    logger.info("%s: found %d MT genes (prefix=%s)", organism, len(mt), prefix)
    return mt.soma_joinid.tolist()


def discover_diverse_datasets(
    census, organism: str, max_datasets: int = MAX_DATASETS_PER_ORGANISM,
) -> pd.DataFrame:
    """Find dataset-tissue combos with good tissue diversity.

    Returns DataFrame with columns: dataset_id, tissue, n_cells.
    Prioritizes tissue diversity — takes at most 3 datasets per tissue.
    """
    exp = census["census_data"][_census_key(organism)]

    logger.info("Discovering datasets for %s...", organism)
    obs_df = exp.obs.read(
        column_names=["dataset_id", "tissue_general"],
        value_filter="is_primary_data == True",
    ).concat().to_pandas()

    # Group by dataset-tissue, count cells
    grouped = obs_df.groupby(["dataset_id", "tissue_general"], observed=True).size().reset_index(name="n_cells")
    grouped = grouped[grouped.n_cells >= MIN_CELLS_PER_DATASET]
    grouped = grouped.rename(columns={"tissue_general": "tissue"})

    logger.info(
        "%s: %d dataset-tissue combos with >= %d cells across %d tissues",
        organism, len(grouped), MIN_CELLS_PER_DATASET, grouped.tissue.nunique(),
    )

    # Prioritize tissue diversity: up to 3 datasets per tissue, largest first
    selected = []
    for tissue, tdf in grouped.groupby("tissue", observed=True):
        top = tdf.nlargest(3, "n_cells")
        selected.append(top)

    result = pd.concat(selected).sort_values("n_cells", ascending=False)

    # Cap total
    if len(result) > max_datasets:
        # Ensure at least 1 per tissue, then fill with largest
        # Use drop_duplicates to get first per tissue (already sorted by n_cells desc)
        first_per_tissue = result.drop_duplicates(subset="tissue", keep="first")
        if len(first_per_tissue) >= max_datasets:
            result = first_per_tissue.head(max_datasets)
        else:
            remaining_idx = result.index.difference(first_per_tissue.index)
            remaining = result.loc[remaining_idx]
            budget = max_datasets - len(first_per_tissue)
            extras = remaining.nlargest(budget, "n_cells")
            result = pd.concat([first_per_tissue, extras])

    logger.info(
        "%s: selected %d dataset-tissue combos across %d tissues",
        organism, len(result), result.tissue.nunique(),
    )
    return result.reset_index(drop=True)


def compute_mito_fractions(
    census, dataset_id: str, tissue: str, organism: str, mt_joinids: list[int],
) -> pd.DataFrame | None:
    """Compute per-cell mito fraction for one dataset-tissue combo.

    Uses raw_sum from obs (total UMIs) and MT-gene sparse query for numerator.
    """
    exp = census["census_data"][_census_key(organism)]
    obs_filter = (
        f"dataset_id == '{dataset_id}' and tissue_general == '{tissue}' "
        f"and is_primary_data == True"
    )

    # Get obs with raw_sum
    obs_df = exp.obs.read(
        column_names=["soma_joinid", "raw_sum"],
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

    # Fetch MT-gene expression only (sparse)
    try:
        with exp.axis_query(
            measurement_name="RNA",
            obs_query=tiledbsoma.AxisQuery(coords=(obs_ids,)),
            var_query=tiledbsoma.AxisQuery(coords=(mt_joinids,)),
        ) as query:
            adata = query.to_anndata(X_name="raw")
    except Exception as e:
        logger.warning("Failed to fetch MT expression for %s/%s: %s", dataset_id[:8], tissue, e)
        return None

    if adata.shape[0] == 0:
        return None

    # Sum MT UMIs per cell
    X = adata.X
    if hasattr(X, "toarray"):
        mt_sum = np.asarray(X.sum(axis=1)).ravel()
    else:
        mt_sum = np.asarray(X).sum(axis=1)

    # Align with obs_df via soma_joinid
    adata_joinids = adata.obs["soma_joinid"].values
    obs_indexed = obs_df.set_index("soma_joinid")
    raw_sum = obs_indexed.loc[adata_joinids, "raw_sum"].values.astype(np.float64)

    # Compute mito fraction (guard against zero total)
    raw_sum = np.where(raw_sum == 0, 1.0, raw_sum)
    mito_frac = mt_sum / raw_sum

    return pd.DataFrame({
        "dataset_id": dataset_id,
        "tissue": tissue,
        "organism": organism,
        "mito_frac": mito_frac,
    })


def aggregate_dataset(cell_df: pd.DataFrame) -> dict:
    """Compute summary statistics for one dataset-tissue combo."""
    mf = cell_df.mito_frac.values
    n = len(mf)
    row = {
        "dataset_id": cell_df.dataset_id.iloc[0],
        "tissue": cell_df.tissue.iloc[0],
        "organism": cell_df.organism.iloc[0],
        "n_cells": n,
        "median_mito_pct": float(np.median(mf) * 100),
        "p95_mito_pct": float(np.percentile(mf, 95) * 100),
        "mean_mito_pct": float(np.mean(mf) * 100),
    }
    for t in THRESHOLDS:
        pct_name = f"pct_removed_at_{int(t*100):02d}pct"
        removed = np.sum(mf > t) / n * 100
        row[pct_name] = float(removed)
    return row


def run_statistics(df: pd.DataFrame) -> pd.DataFrame:
    """Kruskal-Wallis + Dunn post-hoc across tissues."""
    logger.info("Running Kruskal-Wallis test across tissues...")

    # Group by tissue, use median_mito_pct as the test variable
    tissues = df.tissue.unique()
    groups = [df[df.tissue == t].median_mito_pct.values for t in tissues]
    groups = [g for g in groups if len(g) >= 2]  # need at least 2 per group

    if len(groups) < 3:
        logger.warning("Not enough tissue groups (≥2 datasets each) for Kruskal-Wallis")
        return pd.DataFrame()

    stat, pval = stats.kruskal(*groups)
    logger.info("Kruskal-Wallis H=%.2f, p=%.2e", stat, pval)

    # Dunn post-hoc (manual pairwise Mann-Whitney with Bonferroni)
    tissue_names = [t for t, g in zip(tissues, [df[df.tissue == t].median_mito_pct.values for t in tissues]) if len(g) >= 2]
    tissue_groups = {t: df[df.tissue == t].median_mito_pct.values for t in tissue_names if len(df[df.tissue == t]) >= 2}

    comparisons = []
    tissue_list = list(tissue_groups.keys())
    n_comparisons = len(tissue_list) * (len(tissue_list) - 1) // 2

    for i in range(len(tissue_list)):
        for j in range(i + 1, len(tissue_list)):
            t1, t2 = tissue_list[i], tissue_list[j]
            g1, g2 = tissue_groups[t1], tissue_groups[t2]
            try:
                u_stat, u_pval = stats.mannwhitneyu(g1, g2, alternative="two-sided")
                bonf_pval = min(u_pval * n_comparisons, 1.0)
                comparisons.append({
                    "tissue_1": t1,
                    "tissue_2": t2,
                    "U_stat": float(u_stat),
                    "p_raw": float(u_pval),
                    "p_bonferroni": float(bonf_pval),
                    "median_1": float(np.median(g1)),
                    "median_2": float(np.median(g2)),
                })
            except Exception:
                pass

    stats_df = pd.DataFrame(comparisons)
    if not stats_df.empty:
        sig = stats_df[stats_df.p_bonferroni < 0.05]
        logger.info(
            "Dunn post-hoc: %d/%d pairwise comparisons significant (Bonferroni p<0.05)",
            len(sig), len(stats_df),
        )

    # Save stats
    stats_path = PROJECT_DIR / "kruskal_wallis_results.tsv"
    stats_meta = pd.DataFrame([{"test": "Kruskal-Wallis", "H": stat, "p": pval, "n_groups": len(groups)}])
    stats_meta.to_csv(stats_path, sep="\t", index=False)

    if not stats_df.empty:
        dunn_path = PROJECT_DIR / "dunn_posthoc_results.tsv"
        stats_df.to_csv(dunn_path, sep="\t", index=False)
        logger.info("Saved post-hoc results: %s", dunn_path)

    return stats_df


def identify_extreme_tissues(df: pd.DataFrame) -> None:
    """Log tissues where 20% removes <1% vs >30% of cells."""
    col = "pct_removed_at_20pct"
    if col not in df.columns:
        return

    low_removal = df[df[col] < 1.0]
    high_removal = df[df[col] > 30.0]

    logger.info("\n=== EXTREME TISSUES ===")
    if not low_removal.empty:
        logger.info("20%% threshold removes <1%% of cells in these tissues:")
        for _, row in low_removal.iterrows():
            logger.info("  %s (%s): %.1f%% removed", row.tissue, row.organism, row[col])
    else:
        logger.info("No tissues where 20%% threshold removes <1%% of cells.")

    if not high_removal.empty:
        logger.info("20%% threshold removes >30%% of cells in these tissues:")
        for _, row in high_removal.iterrows():
            logger.info("  %s (%s): %.1f%% removed", row.tissue, row.organism, row[col])
    else:
        logger.info("No tissues where 20%% threshold removes >30%% of cells.")


def main():
    t0 = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PROJECT_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("=== P1: Mitochondrial Threshold Audit ===")

    # Resume from checkpoint if available
    checkpoint_path = PROJECT_DIR / "p1_checkpoint.tsv"
    if checkpoint_path.exists():
        existing = pd.read_csv(checkpoint_path, sep="\t")
        all_results = existing.to_dict("records")
        done_keys = {(r["dataset_id"], r["tissue"], r["organism"]) for r in all_results}
        logger.info("Resuming from checkpoint: %d results loaded", len(all_results))
    else:
        all_results = []
        done_keys = set()

    census = cellxgene_census.open_soma(census_version=CENSUS_VERSION)

    for organism in ["Homo sapiens", "Mus musculus"]:
        logger.info("\n--- Processing %s ---", organism)

        # Get MT gene IDs
        mt_joinids = get_mt_gene_ids(census, organism)

        # Discover diverse datasets
        datasets = discover_diverse_datasets(census, organism)
        logger.info("Will process %d dataset-tissue combos", len(datasets))

        # Process each dataset
        for _, row in tqdm(
            datasets.iterrows(), total=len(datasets),
            desc=f"Mito audit ({organism})",
        ):
            # Skip NaN dataset_ids
            if pd.isna(row.dataset_id):
                continue

            # Skip if already done
            key = (row.dataset_id, row.tissue, organism)
            if key in done_keys:
                logger.info("  Skipping %s/%s (already done)", row.dataset_id[:8], row.tissue)
                continue

            try:
                cell_df = compute_mito_fractions(
                    census, row.dataset_id, row.tissue, organism, mt_joinids,
                )
            except Exception as e:
                logger.warning("  Error on %s/%s: %s", row.dataset_id[:8], row.tissue, e)
                continue

            if cell_df is not None and len(cell_df) > 0:
                agg = aggregate_dataset(cell_df)
                all_results.append(agg)
                done_keys.add(key)
                logger.info(
                    "  %s/%s: %d cells, median=%.1f%%, P95=%.1f%%, removed@20%%=%.1f%%",
                    row.dataset_id[:8], row.tissue, agg["n_cells"],
                    agg["median_mito_pct"], agg["p95_mito_pct"],
                    agg["pct_removed_at_20pct"],
                )
                del cell_df
                gc.collect()

                # Save checkpoint after each successful dataset
                pd.DataFrame(all_results).to_csv(checkpoint_path, sep="\t", index=False)
            else:
                logger.debug("  Skipped %s/%s", row.dataset_id[:8], row.tissue)

    census.close()

    # Build output TSV
    result_df = pd.DataFrame(all_results)
    out_path = OUTPUT_DIR / "p1_mito_threshold.tsv"
    result_df.to_csv(out_path, sep="\t", index=False)
    logger.info("\nSaved %d rows to %s", len(result_df), out_path)

    # Statistics
    if len(result_df) >= 3:
        run_statistics(result_df)
        identify_extreme_tissues(result_df)

    # Summary
    elapsed = time.time() - t0
    logger.info("\n=== P1 COMPLETE ===")
    logger.info("Datasets processed: %d", len(result_df))
    logger.info("Tissues: %d", result_df.tissue.nunique())
    logger.info("Organisms: %s", list(result_df.organism.unique()))
    logger.info(
        "Overall median mito %%: %.1f%% (range: %.1f%%–%.1f%%)",
        result_df.median_mito_pct.median(),
        result_df.median_mito_pct.min(),
        result_df.median_mito_pct.max(),
    )
    logger.info("Elapsed: %.1f minutes", elapsed / 60)

    return result_df


if __name__ == "__main__":
    main()
