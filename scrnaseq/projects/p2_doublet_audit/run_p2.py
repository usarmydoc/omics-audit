#!/usr/bin/env python3
"""P2 — Doublet Rate Audit.

Question: Do doublet detection tools agree on which cells are doublets,
and does disagreement correlate with cell type composition?

Strategy:
- Pull raw count matrices from CELLxGENE Census for 20+ datasets
  (10 human, 10 mouse, diverse tissues).
- For each dataset: run scDblFinder (via rpy2) and Scrublet (Python).
- Compute Cohen's kappa, doublet rates, disagreement rate.
- Correlate disagreement_rate with n_cell_types via Spearman.
"""

from __future__ import annotations

import gc
import logging
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import cellxgene_census
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scrublet
import tiledbsoma
from scipy import stats
from sklearn.metrics import cohen_kappa_score
from tqdm import tqdm

# ── Config ──────────────────────────────────────────────────────────
CENSUS_VERSION = "2025-11-08"
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = REPO_ROOT / "scrnaseq" / "output"
PROJECT_DIR = Path(__file__).resolve().parent

MAX_DATASETS_PER_ORGANISM = 12
MAX_CELLS_PER_DATASET = 20_000
MIN_CELLS_PER_DATASET = 500
RNG_SEED = 42

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(PROJECT_DIR / "p2_run.log"),
    ],
)
logger = logging.getLogger(__name__)


# ── Census helpers ──────────────────────────────────────────────────

def _census_key(organism: str) -> str:
    return organism.lower().replace(" ", "_")


def discover_diverse_datasets(
    census, organism: str, max_datasets: int = MAX_DATASETS_PER_ORGANISM,
) -> pd.DataFrame:
    """Find diverse dataset-tissue combos.

    Returns DataFrame with columns: dataset_id, tissue, n_cells.
    Prioritizes tissue diversity — takes at most 2 datasets per tissue.
    """
    exp = census["census_data"][_census_key(organism)]

    logger.info("Discovering datasets for %s...", organism)
    obs_df = exp.obs.read(
        column_names=["dataset_id", "tissue_general"],
        value_filter="is_primary_data == True",
    ).concat().to_pandas()

    # Skip NaN dataset_ids
    obs_df = obs_df.dropna(subset=["dataset_id"])

    # Group by dataset-tissue, count cells
    grouped = (
        obs_df
        .groupby(["dataset_id", "tissue_general"], observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    grouped = grouped[grouped.n_cells >= MIN_CELLS_PER_DATASET]
    grouped = grouped.rename(columns={"tissue_general": "tissue"})

    # Cap at MAX_CELLS to prefer moderately sized datasets (faster processing)
    grouped = grouped[grouped.n_cells <= MAX_CELLS_PER_DATASET * 2]

    logger.info(
        "%s: %d dataset-tissue combos with %d-%d cells across %d tissues",
        organism, len(grouped), MIN_CELLS_PER_DATASET,
        MAX_CELLS_PER_DATASET * 2, grouped.tissue.nunique(),
    )

    # Prioritize tissue diversity: up to 2 datasets per tissue, largest first
    selected = []
    for tissue, tdf in grouped.groupby("tissue", observed=True):
        top = tdf.nlargest(2, "n_cells")
        selected.append(top)

    if not selected:
        return pd.DataFrame(columns=["dataset_id", "tissue", "n_cells"])

    result = pd.concat(selected).sort_values("n_cells", ascending=False)

    # Cap total
    if len(result) > max_datasets:
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


def fetch_dataset(
    census, dataset_id: str, tissue: str, organism: str,
) -> tuple[sp.csr_matrix, np.ndarray] | None:
    """Fetch sparse raw counts and cell_type labels for one dataset-tissue.

    Returns (X_sparse, cell_types) or None on failure.
    Caps at MAX_CELLS_PER_DATASET cells.
    """
    exp = census["census_data"][_census_key(organism)]
    obs_filter = (
        f"dataset_id == '{dataset_id}' and tissue_general == '{tissue}' "
        f"and is_primary_data == True"
    )

    # Get obs with cell_type
    obs_df = exp.obs.read(
        column_names=["soma_joinid", "cell_type"],
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
    cell_types = obs_df.cell_type.values

    # Fetch full expression matrix (sparse)
    try:
        with exp.axis_query(
            measurement_name="RNA",
            obs_query=tiledbsoma.AxisQuery(coords=(obs_ids,)),
        ) as query:
            adata = query.to_anndata(X_name="raw")
    except Exception as e:
        logger.warning(
            "Failed to fetch expression for %s/%s: %s",
            dataset_id[:8], tissue, e,
        )
        return None

    if adata.shape[0] == 0:
        return None

    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    elif not sp.isspmatrix_csr(X):
        X = X.tocsr()

    # Align cell_types with adata obs order
    adata_joinids = adata.obs["soma_joinid"].values
    obs_indexed = obs_df.set_index("soma_joinid")
    cell_types = obs_indexed.loc[adata_joinids, "cell_type"].values

    return X, cell_types


# ── Doublet detection ───────────────────────────────────────────────

def run_scrublet(X: sp.csr_matrix, seed: int = RNG_SEED) -> np.ndarray:
    """Run Scrublet, return boolean array of predicted doublets."""
    scrub = scrublet.Scrublet(X, random_state=seed)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
    return predicted_doublets.astype(bool)


def run_scdblfinder_rpy2(X: sp.csr_matrix) -> np.ndarray | None:
    """Run scDblFinder via rpy2, return boolean array of predicted doublets."""
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import numpy2ri, pandas2ri
        from rpy2.robjects.packages import importr

        numpy2ri.activate()

        base = importr("base")
        sce_pkg = importr("SingleCellExperiment")
        scdblfinder_pkg = importr("scDblFinder")
        s4vectors = importr("S4Vectors")
        matrix_pkg = importr("Matrix")

        # Convert sparse matrix to R dgCMatrix
        csc = X.tocsc()

        r_sparse = matrix_pkg.sparseMatrix(
            i=ro.IntVector(csc.indices + 1),  # R is 1-indexed
            p=ro.IntVector(csc.indptr),
            x=ro.FloatVector(csc.data.astype(float)),
            dims=ro.IntVector(list(csc.shape)),
        )

        # Create SingleCellExperiment
        sce = sce_pkg.SingleCellExperiment(
            assays=base.list(counts=r_sparse),
        )

        # Run scDblFinder
        sce_result = scdblfinder_pkg.scDblFinder(sce)

        # Extract doublet class
        col_data = ro.r("colData")(sce_result)
        dbl_class = np.array(ro.r("as.character")(col_data.rx2("scDblFinder.class")))

        numpy2ri.deactivate()

        return dbl_class == "doublet"

    except Exception as e:
        logger.warning("rpy2 scDblFinder failed: %s", e)
        return None


def run_scdblfinder_subprocess(X: sp.csr_matrix) -> np.ndarray | None:
    """Fallback: run scDblFinder via Rscript subprocess."""
    try:
        import scipy.io

        with tempfile.TemporaryDirectory() as tmpdir:
            mtx_path = Path(tmpdir) / "counts.mtx"
            out_path = Path(tmpdir) / "doublet_calls.csv"

            # Write sparse matrix as MatrixMarket
            scipy.io.mmwrite(str(mtx_path), X.T)  # genes x cells for R

            r_script = f"""
library(Matrix)
library(scDblFinder)
library(SingleCellExperiment)

counts <- readMM("{mtx_path}")
sce <- SingleCellExperiment(assays=list(counts=counts))
sce <- scDblFinder(sce)
calls <- colData(sce)$scDblFinder.class
write.csv(data.frame(doublet=calls), "{out_path}", row.names=FALSE)
"""
            r_script_path = Path(tmpdir) / "run_scdblfinder.R"
            r_script_path.write_text(r_script)

            result = subprocess.run(
                ["Rscript", str(r_script_path)],
                capture_output=True, text=True, timeout=600,
            )

            if result.returncode != 0:
                logger.warning("Rscript scDblFinder failed: %s", result.stderr[:500])
                return None

            calls_df = pd.read_csv(out_path)
            return (calls_df["doublet"].values == "doublet").astype(bool)

    except Exception as e:
        logger.warning("Rscript subprocess scDblFinder failed: %s", e)
        return None


def run_scdblfinder(X: sp.csr_matrix) -> np.ndarray | None:
    """Run scDblFinder, trying rpy2 first, falling back to Rscript."""
    result = run_scdblfinder_rpy2(X)
    if result is not None:
        return result

    logger.info("  Falling back to Rscript subprocess for scDblFinder")
    return run_scdblfinder_subprocess(X)


# ── Per-dataset processing ──────────────────────────────────────────

def process_dataset(
    census, dataset_id: str, tissue: str, organism: str,
) -> dict | None:
    """Run both doublet tools on one dataset, compute agreement metrics."""
    fetch_result = fetch_dataset(census, dataset_id, tissue, organism)
    if fetch_result is None:
        return None

    X, cell_types = fetch_result
    n_cells = X.shape[0]

    # Number of unique cell types (excluding nan/unknown)
    ct_series = pd.Series(cell_types)
    ct_clean = ct_series.dropna()
    ct_clean = ct_clean[ct_clean != "unknown"]
    n_cell_types = ct_clean.nunique()

    # Run Scrublet
    try:
        scrublet_calls = run_scrublet(X)
    except Exception as e:
        logger.warning("  Scrublet failed for %s/%s: %s", dataset_id[:8], tissue, e)
        return None

    # Run scDblFinder
    scdblfinder_calls = run_scdblfinder(X)
    if scdblfinder_calls is None:
        logger.warning("  scDblFinder failed for %s/%s", dataset_id[:8], tissue)
        return None

    # Ensure same length
    if len(scrublet_calls) != len(scdblfinder_calls):
        logger.warning(
            "  Length mismatch: scrublet=%d, scdblfinder=%d",
            len(scrublet_calls), len(scdblfinder_calls),
        )
        return None

    # Compute metrics
    scrublet_rate = scrublet_calls.sum() / n_cells
    scdblfinder_rate = scdblfinder_calls.sum() / n_cells

    # Cohen's kappa
    kappa = cohen_kappa_score(
        scrublet_calls.astype(int),
        scdblfinder_calls.astype(int),
    )

    # Disagreement rate: fraction where tools disagree
    disagreement = (scrublet_calls != scdblfinder_calls).sum() / n_cells

    return {
        "dataset_id": dataset_id,
        "tissue": tissue,
        "organism": organism,
        "n_cells": n_cells,
        "scdblfinder_rate": round(float(scdblfinder_rate), 5),
        "scrublet_rate": round(float(scrublet_rate), 5),
        "cohens_kappa": round(float(kappa), 4),
        "n_cell_types": int(n_cell_types),
        "disagreement_rate": round(float(disagreement), 5),
    }


# ── Main ────────────────────────────────────────────────────────────

def main():
    t0 = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PROJECT_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("=== P2: Doublet Rate Audit ===")

    # Resume from checkpoint if available
    checkpoint_path = PROJECT_DIR / "p2_checkpoint.tsv"
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

        # Discover diverse datasets
        datasets = discover_diverse_datasets(census, organism)
        logger.info("Will process %d dataset-tissue combos", len(datasets))

        # Process each dataset
        for _, row in tqdm(
            datasets.iterrows(), total=len(datasets),
            desc=f"Doublet audit ({organism})",
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
                result = process_dataset(census, row.dataset_id, row.tissue, organism)
            except Exception as e:
                logger.warning("  Error on %s/%s: %s", row.dataset_id[:8], row.tissue, e)
                continue

            if result is not None:
                all_results.append(result)
                done_keys.add(key)
                logger.info(
                    "  %s/%s: %d cells, scDblFinder=%.1f%%, Scrublet=%.1f%%, "
                    "kappa=%.3f, n_types=%d, disagree=%.1f%%",
                    row.dataset_id[:8], row.tissue, result["n_cells"],
                    result["scdblfinder_rate"] * 100,
                    result["scrublet_rate"] * 100,
                    result["cohens_kappa"],
                    result["n_cell_types"],
                    result["disagreement_rate"] * 100,
                )

                # Save checkpoint after each successful dataset
                pd.DataFrame(all_results).to_csv(
                    checkpoint_path, sep="\t", index=False,
                )

                gc.collect()
            else:
                logger.debug("  Skipped %s/%s", row.dataset_id[:8], row.tissue)

    census.close()

    # Build output TSV
    result_df = pd.DataFrame(all_results)
    out_path = OUTPUT_DIR / "p2_doublet_audit.tsv"
    result_df.to_csv(out_path, sep="\t", index=False)
    logger.info("\nSaved %d rows to %s", len(result_df), out_path)

    # Spearman correlation: disagreement_rate vs n_cell_types
    if len(result_df) >= 5:
        rho, pval = stats.spearmanr(
            result_df["disagreement_rate"],
            result_df["n_cell_types"],
        )
        logger.info(
            "\n=== SPEARMAN CORRELATION ===\n"
            "disagreement_rate vs n_cell_types: rho=%.3f, p=%.2e",
            rho, pval,
        )
    else:
        logger.warning("Too few datasets (%d) for Spearman correlation", len(result_df))

    # Summary
    elapsed = time.time() - t0
    logger.info("\n=== P2 COMPLETE ===")
    logger.info("Datasets processed: %d", len(result_df))
    if not result_df.empty:
        logger.info("Tissues: %d", result_df.tissue.nunique())
        logger.info("Organisms: %s", list(result_df.organism.unique()))
        logger.info(
            "Median Cohen's kappa: %.3f (range: %.3f-%.3f)",
            result_df.cohens_kappa.median(),
            result_df.cohens_kappa.min(),
            result_df.cohens_kappa.max(),
        )
        logger.info(
            "Median disagreement rate: %.1f%% (range: %.1f%%-%.1f%%)",
            result_df.disagreement_rate.median() * 100,
            result_df.disagreement_rate.min() * 100,
            result_df.disagreement_rate.max() * 100,
        )
    logger.info("Elapsed: %.1f minutes", elapsed / 60)

    return result_df


if __name__ == "__main__":
    main()
