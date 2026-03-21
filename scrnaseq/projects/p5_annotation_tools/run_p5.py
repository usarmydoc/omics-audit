#!/usr/bin/env python3
"""P5 — Cell Type Annotation Tool Agreement Audit.

Question: How accurately do SingleR, CellTypist, and a marker-score classifier
annotate cell types vs Census ground truth?

Strategy:
- Select 15+ datasets (human + mouse) with >=3 distinct cell types matching
  the marker_reliability_atlas categories.
- For each dataset (cap 15K cells): fetch raw counts from Census.
- Run three annotation methods:
  1. SingleR (via rpy2) with celldex reference
  2. CellTypist with default model
  3. Marker-score classifier using marker_reliability_atlas.tsv
- Compute precision, recall, F1 per cell type vs Census ground truth.
- Output: p5_annotation_tools.tsv
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
from scipy import sparse
from sklearn.metrics import precision_recall_fscore_support
from tqdm import tqdm

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ── Config ──────────────────────────────────────────────────────────
CENSUS_VERSION = "2025-11-08"
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = REPO_ROOT / "scrnaseq" / "output"
PROJECT_DIR = Path(__file__).resolve().parent
ATLAS_PATH = REPO_ROOT / "marker_atlas" / "output" / "marker_reliability_atlas.tsv"

MAX_DATASETS = 20
MAX_CELLS_PER_DATASET = 15_000
MIN_CELLS_PER_DATASET = 200
MIN_CELL_TYPES = 3
RNG_SEED = 42

# Atlas cell types that we can map to
ATLAS_CELL_TYPES = {
    "Endothelial", "Fibroblasts", "Hepatocytes", "Macrophages", "NK cells", "T cells",
}

# Mapping from Census cell_type labels to atlas categories (many-to-one)
CENSUS_TO_ATLAS = {
    # Endothelial
    "endothelial cell": "Endothelial",
    "vascular endothelial cell": "Endothelial",
    "endothelial cell of lymphatic vessel": "Endothelial",
    "endothelial cell of vascular tree": "Endothelial",
    "capillary endothelial cell": "Endothelial",
    "blood vessel endothelial cell": "Endothelial",
    "endothelial cell of artery": "Endothelial",
    "endothelial cell of sinusoid": "Endothelial",
    "endothelial cell of hepatic sinusoid": "Endothelial",
    # Fibroblasts
    "fibroblast": "Fibroblasts",
    "myofibroblast cell": "Fibroblasts",
    "adventitial cell": "Fibroblasts",
    # Hepatocytes
    "hepatocyte": "Hepatocytes",
    "periportal region hepatocyte": "Hepatocytes",
    "pericentral region hepatocyte": "Hepatocytes",
    # Macrophages
    "macrophage": "Macrophages",
    "alveolar macrophage": "Macrophages",
    "kupffer cell": "Macrophages",
    "tissue-resident macrophage": "Macrophages",
    "inflammatory macrophage": "Macrophages",
    # NK cells
    "natural killer cell": "NK cells",
    "nk cell": "NK cells",
    "group 1 innate lymphoid cell": "NK cells",
    # T cells
    "t cell": "T cells",
    "cd4-positive, alpha-beta t cell": "T cells",
    "cd8-positive, alpha-beta t cell": "T cells",
    "regulatory t cell": "T cells",
    "effector memory cd4-positive, alpha-beta t cell": "T cells",
    "effector memory cd8-positive, alpha-beta t cell": "T cells",
    "naive thymus-derived cd4-positive, alpha-beta t cell": "T cells",
    "naive thymus-derived cd8-positive, alpha-beta t cell": "T cells",
    "gamma-delta t cell": "T cells",
    "mucosal invariant t cell": "T cells",
    "cd4-positive helper t cell": "T cells",
    "central memory cd4-positive, alpha-beta t cell": "T cells",
    "mature nk t cell": "T cells",
}

# Mapping from SingleR labels to atlas categories
SINGLER_TO_ATLAS = {
    # Human Primary Cell Atlas labels
    "endothelial_cells": "Endothelial",
    "endothelial cells": "Endothelial",
    "Endothelial_cells": "Endothelial",
    # Fibroblasts
    "fibroblasts": "Fibroblasts",
    "Fibroblasts": "Fibroblasts",
    # Hepatocytes
    "hepatocytes": "Hepatocytes",
    "Hepatocytes": "Hepatocytes",
    # Macrophages
    "macrophage": "Macrophages",
    "macrophages": "Macrophages",
    "Macrophage": "Macrophages",
    "Macrophages": "Macrophages",
    "monocyte": "Macrophages",
    "Monocyte": "Macrophages",
    "Monocytes": "Macrophages",
    # NK cells
    "nk_cell": "NK cells",
    "nk_cells": "NK cells",
    "NK_cell": "NK cells",
    "NK_cells": "NK cells",
    "NK cells": "NK cells",
    "NK cell": "NK cells",
    # T cells
    "t_cells": "T cells",
    "T_cells": "T cells",
    "T cells": "T cells",
    "T-cells": "T cells",
    "CD4+ T-cells": "T cells",
    "CD8+ T-cells": "T cells",
    "CD4+ T cells": "T cells",
    "CD8+ T cells": "T cells",
    "T_cell": "T cells",
}

# Mapping from CellTypist labels to atlas categories
CELLTYPIST_TO_ATLAS = {
    "endothelial cell": "Endothelial",
    "endothelial cells": "Endothelial",
    "Endothelial cells": "Endothelial",
    "fibroblast": "Fibroblasts",
    "fibroblasts": "Fibroblasts",
    "Fibroblasts": "Fibroblasts",
    "hepatocyte": "Hepatocytes",
    "hepatocytes": "Hepatocytes",
    "Hepatocytes": "Hepatocytes",
    "macrophage": "Macrophages",
    "macrophages": "Macrophages",
    "Macrophages": "Macrophages",
    "kupffer cell": "Macrophages",
    "Kupffer cells": "Macrophages",
    "alveolar macrophage": "Macrophages",
    "nk cell": "NK cells",
    "nk cells": "NK cells",
    "NK cells": "NK cells",
    "natural killer cell": "NK cells",
    "t cell": "T cells",
    "t cells": "T cells",
    "T cells": "T cells",
    "cd4+ t cell": "T cells",
    "cd8+ t cell": "T cells",
    "CD4+ T cells": "T cells",
    "CD8+ T cells": "T cells",
    "regulatory t cell": "T cells",
}

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(PROJECT_DIR / "p5_run.log"),
    ],
)
logger = logging.getLogger(__name__)


def _census_key(organism: str) -> str:
    return organism.lower().replace(" ", "_")


def verify_marker_atlas() -> pd.DataFrame:
    """Verify marker atlas exists with both organisms. Retry up to 10 times."""
    for attempt in range(10):
        if ATLAS_PATH.exists():
            df = pd.read_csv(ATLAS_PATH, sep="\t")
            if (
                "Homo sapiens" in df.organism.values
                and "Mus musculus" in df.organism.values
            ):
                logger.info(
                    "Marker atlas verified: %d rows, both organisms present",
                    len(df),
                )
                return df
            else:
                logger.warning(
                    "Marker atlas exists but missing organism(s). "
                    "Found: %s",
                    df.organism.unique().tolist(),
                )
        else:
            logger.warning("Marker atlas not found at %s", ATLAS_PATH)

        if attempt < 9:
            logger.info("Retrying in 5 minutes... (attempt %d/10)", attempt + 1)
            time.sleep(300)

    logger.error(
        "Marker atlas not ready after 10 attempts. "
        "Ensure marker_reliability_atlas.tsv exists with both organisms."
    )
    sys.exit(1)


def discover_eligible_datasets(census, organism: str) -> pd.DataFrame:
    """Find datasets with >=3 distinct cell types matching atlas categories."""
    exp = census["census_data"][_census_key(organism)]
    logger.info("Discovering eligible datasets for %s...", organism)

    obs_df = exp.obs.read(
        column_names=["dataset_id", "tissue_general", "cell_type"],
        value_filter="is_primary_data == True",
    ).concat().to_pandas()

    # Map Census cell types to atlas categories
    obs_df["atlas_type"] = obs_df.cell_type.str.lower().map(
        {k.lower(): v for k, v in CENSUS_TO_ATLAS.items()}
    )
    obs_df = obs_df.dropna(subset=["atlas_type"])

    results = []
    for ds_id, ds_df in obs_df.groupby("dataset_id", observed=True):
        if pd.isna(ds_id):
            continue
        distinct_types = ds_df.atlas_type.nunique()
        if distinct_types < MIN_CELL_TYPES:
            continue
        n_cells = len(ds_df)
        if n_cells < MIN_CELLS_PER_DATASET:
            continue
        # Skip very large datasets — obs query alone takes too long
        if n_cells > 100_000:
            continue
        tissue = (
            ds_df.tissue_general.mode().iloc[0]
            if len(ds_df.tissue_general.mode()) > 0
            else "unknown"
        )
        results.append({
            "dataset_id": ds_id,
            "tissue": tissue,
            "n_cells": n_cells,
            "n_atlas_types": distinct_types,
            "atlas_types": sorted(ds_df.atlas_type.unique().tolist()),
        })

    result_df = pd.DataFrame(results)
    if len(result_df) == 0:
        logger.warning("No eligible datasets found for %s", organism)
        return result_df

    # Sort by number of atlas types (desc) then by cell count (asc — prefer
    # smaller datasets that are faster to fetch from Census)
    result_df = result_df.sort_values(
        ["n_atlas_types", "n_cells"], ascending=[False, True]
    )
    # Tissue diversity: one per tissue
    result_df = result_df.drop_duplicates(subset="tissue", keep="first")
    result_df = result_df.head(MAX_DATASETS)

    logger.info(
        "%s: found %d eligible datasets across %d tissues",
        organism, len(result_df), result_df.tissue.nunique(),
    )

    del obs_df
    gc.collect()
    return result_df.reset_index(drop=True)


def fetch_dataset(
    census, dataset_id: str, organism: str,
) -> tuple[np.ndarray, pd.DataFrame, list[str]] | None:
    """Fetch raw counts for a dataset, capped at MAX_CELLS_PER_DATASET.

    Returns (count_matrix, obs_df, gene_names) or None.
    obs_df includes cell_type and atlas_type columns.
    """
    exp = census["census_data"][_census_key(organism)]

    obs_filter = (
        f"dataset_id == '{dataset_id}' and is_primary_data == True"
    )
    obs_df = exp.obs.read(
        column_names=["soma_joinid", "cell_type"],
        value_filter=obs_filter,
    ).concat().to_pandas()

    # Map to atlas types
    obs_df["atlas_type"] = obs_df.cell_type.str.lower().map(
        {k.lower(): v for k, v in CENSUS_TO_ATLAS.items()}
    )

    # Keep only cells with a valid atlas mapping
    obs_df = obs_df.dropna(subset=["atlas_type"]).copy()

    if len(obs_df) < MIN_CELLS_PER_DATASET:
        return None

    # Subsample if too large — proportional stratified sample
    if len(obs_df) > MAX_CELLS_PER_DATASET:
        frac = MAX_CELLS_PER_DATASET / len(obs_df)
        parts = []
        for _, grp in obs_df.groupby("atlas_type", observed=True):
            n = max(1, int(len(grp) * frac))
            parts.append(grp.sample(n=n, random_state=RNG_SEED, replace=False))
        obs_df = pd.concat(parts).reset_index(drop=True)

    obs_ids = obs_df.soma_joinid.tolist()

    try:
        with exp.axis_query(
            measurement_name="RNA",
            obs_query=tiledbsoma.AxisQuery(coords=(obs_ids,)),
        ) as query:
            adata = query.to_anndata(X_name="raw")
    except Exception as e:
        logger.warning("Failed to fetch expression for %s: %s", dataset_id[:8], e)
        return None

    if adata.shape[0] == 0:
        return None

    gene_names = (
        list(adata.var.feature_name)
        if "feature_name" in adata.var.columns
        else list(adata.var_names)
    )

    X = adata.X
    # Align obs_df with adata
    adata_joinids = adata.obs["soma_joinid"].values
    obs_indexed = obs_df.set_index("soma_joinid")
    aligned_obs = obs_indexed.loc[adata_joinids].reset_index()

    return X, aligned_obs, gene_names


def _map_labels(labels: pd.Series, mapping: dict[str, str]) -> pd.Series:
    """Map tool-specific labels to atlas categories using fuzzy matching."""
    # Build case-insensitive mapping
    lower_map = {k.lower(): v for k, v in mapping.items()}

    def _map_one(label):
        if pd.isna(label):
            return "unmapped"
        label_lower = str(label).lower()
        # Exact match
        if label_lower in lower_map:
            return lower_map[label_lower]
        # Substring match
        for key, val in lower_map.items():
            if key in label_lower or label_lower in key:
                return val
        # Keyword match
        for atlas_type in ATLAS_CELL_TYPES:
            keywords = atlas_type.lower().replace(" cells", "").replace("s", "")
            if keywords in label_lower:
                return atlas_type
        return "unmapped"

    return labels.apply(_map_one)


def run_singler(X, obs_df: pd.DataFrame, gene_names: list[str], organism: str) -> pd.Series | None:
    """Run SingleR via rpy2 and return mapped predictions."""
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri, numpy2ri
        from rpy2.robjects.packages import importr

        pandas2ri.activate()
        numpy2ri.activate()

        singler = importr("SingleR")
        sce_pkg = importr("SingleCellExperiment")
        s4vectors = importr("S4Vectors")

        if "homo" in organism.lower():
            celldex = importr("celldex")
            ref = celldex.HumanPrimaryCellAtlasData()
        else:
            celldex = importr("celldex")
            ref = celldex.MouseRNAseqData()

        # Convert to dense matrix for R
        if hasattr(X, "toarray"):
            X_dense = X.toarray().astype(np.float64)
        else:
            X_dense = np.asarray(X, dtype=np.float64)

        # Create R matrix (genes x cells — SingleR expects this orientation)
        nr, nc = X_dense.shape
        r_matrix = ro.r["matrix"](
            ro.vectors.FloatVector(X_dense.T.ravel()),
            nrow=len(gene_names),
            ncol=nr,
        )
        r_matrix.rownames = ro.vectors.StrVector(gene_names)

        # Run SingleR
        results = singler.SingleR(
            test=r_matrix,
            ref=ref,
            labels=ro.r("$")(ref, "label.main"),
        )

        # Extract predictions
        preds = list(ro.r("as.character")(ro.r("$")(results, "labels")))

        pandas2ri.deactivate()
        numpy2ri.deactivate()

        pred_series = pd.Series(preds, index=obs_df.index)
        return _map_labels(pred_series, SINGLER_TO_ATLAS)

    except Exception as e:
        logger.warning("SingleR failed: %s", e)
        try:
            pandas2ri.deactivate()
            numpy2ri.deactivate()
        except Exception:
            pass
        return None


def run_celltypist(X, obs_df: pd.DataFrame, gene_names: list[str]) -> pd.Series | None:
    """Run CellTypist and return mapped predictions."""
    try:
        import celltypist
        from celltypist import models
        import scanpy as sc
        import anndata as ad

        # Create AnnData object with log-normalized data (CellTypist expects this)
        if hasattr(X, "toarray"):
            X_dense = X.toarray().astype(np.float32)
        else:
            X_dense = np.asarray(X, dtype=np.float32)

        adata = ad.AnnData(
            X=X_dense,
            var=pd.DataFrame(index=gene_names),
            obs=obs_df[["cell_type", "atlas_type"]].reset_index(drop=True),
        )

        # Normalize for CellTypist
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        # Download default model if needed
        models.download_models(force_update=False)
        model = models.Model.load()

        # Annotate
        predictions = celltypist.annotate(adata, model=model)
        pred_labels = predictions.predicted_labels["predicted_labels"]

        return _map_labels(pred_labels, CELLTYPIST_TO_ATLAS)

    except Exception as e:
        logger.warning("CellTypist failed: %s", e)
        return None


def run_marker_score(
    X, obs_df: pd.DataFrame, gene_names: list[str],
    atlas_df: pd.DataFrame, organism: str,
) -> pd.Series | None:
    """Marker-score classifier: for each cell, compute mean expression of atlas
    markers per cell type, assign to highest-scoring type.

    Uses markers with reliability_score > 0.5.
    """
    try:
        # Filter atlas to this organism and high-reliability markers
        org_atlas = atlas_df[
            (atlas_df.organism == organism) & (atlas_df.reliability_score > 0.5)
        ].copy()

        if len(org_atlas) == 0:
            logger.warning("No high-reliability markers for %s", organism)
            return None

        # Build gene name to index mapping
        gene_to_idx = {g: i for i, g in enumerate(gene_names)}

        # For each atlas cell type, get marker gene indices
        ct_marker_indices = {}
        for ct, ct_df in org_atlas.groupby("cell_type", observed=True):
            markers = ct_df.gene.tolist()
            indices = [gene_to_idx[g] for g in markers if g in gene_to_idx]
            if indices:
                ct_marker_indices[ct] = indices

        if len(ct_marker_indices) < 2:
            logger.warning("Too few cell types with markers in gene list")
            return None

        # Convert to dense and log-normalize
        if hasattr(X, "toarray"):
            X_dense = X.toarray().astype(np.float32)
        else:
            X_dense = np.asarray(X, dtype=np.float32)

        lib_sizes = X_dense.sum(axis=1, keepdims=True)
        lib_sizes = np.where(lib_sizes == 0, 1, lib_sizes)
        X_norm = np.log1p(X_dense / lib_sizes * 1e4)

        # Score each cell for each cell type
        scores = {}
        for ct, indices in ct_marker_indices.items():
            scores[ct] = X_norm[:, indices].mean(axis=1)

        score_df = pd.DataFrame(scores, index=obs_df.index)
        predictions = score_df.idxmax(axis=1)

        return predictions

    except Exception as e:
        logger.warning("Marker-score classifier failed: %s", e)
        return None


def compute_metrics(
    ground_truth: pd.Series, predictions: pd.Series,
) -> list[dict]:
    """Compute per-cell-type precision, recall, F1."""
    # Only evaluate on cells where both GT and pred are in atlas types
    valid_mask = (
        ground_truth.isin(ATLAS_CELL_TYPES)
        & predictions.isin(ATLAS_CELL_TYPES)
    )
    gt = ground_truth[valid_mask]
    pred = predictions[valid_mask]

    if len(gt) == 0:
        return []

    labels = sorted(set(gt.unique()) | set(pred.unique()))
    prec, rec, f1, support = precision_recall_fscore_support(
        gt, pred, labels=labels, zero_division=0,
    )

    results = []
    for i, label in enumerate(labels):
        n = int(support[i]) if support[i] > 0 else int((gt == label).sum())
        results.append({
            "cell_type": label,
            "precision": float(prec[i]),
            "recall": float(rec[i]),
            "f1": float(f1[i]),
            "n_cells": n,
        })

    return results


def main():
    t0 = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PROJECT_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("=== P5: Cell Type Annotation Tool Agreement Audit ===")

    # Verify marker atlas
    atlas_df = verify_marker_atlas()

    # Resume from checkpoint
    checkpoint_path = PROJECT_DIR / "p5_checkpoint.tsv"
    if checkpoint_path.exists():
        existing = pd.read_csv(checkpoint_path, sep="\t")
        all_results = existing.to_dict("records")
        done_keys = {
            (r["dataset_id"], r["organism"], r["tool"])
            for r in all_results
        }
        logger.info("Resuming from checkpoint: %d rows", len(all_results))
    else:
        all_results = []
        done_keys = set()

    census = cellxgene_census.open_soma(census_version=CENSUS_VERSION)

    for organism in ["Homo sapiens", "Mus musculus"]:
        logger.info("\n--- Processing %s ---", organism)

        eligible = discover_eligible_datasets(census, organism)
        if len(eligible) == 0:
            continue

        logger.info("Will process %d datasets", len(eligible))

        for _, row in tqdm(
            eligible.iterrows(), total=len(eligible),
            desc=f"P5 ({organism})",
        ):
            ds_id = row.dataset_id
            if pd.isna(ds_id):
                continue

            # Check if all 3 tools are done for this dataset
            tools_done = {
                t for t in ["SingleR", "CellTypist", "MarkerScore"]
                if (ds_id, organism, t) in done_keys
            }
            if len(tools_done) == 3:
                logger.info("  Skipping %s (all tools done)", ds_id[:8])
                continue

            logger.info(
                "  Processing %s/%s, %d atlas types, %d cells",
                ds_id[:8], row.tissue, row.n_atlas_types, row.n_cells,
            )

            try:
                result = fetch_dataset(census, ds_id, organism)
            except Exception as e:
                logger.warning("  Fetch failed for %s: %s", ds_id[:8], e)
                continue

            if result is None:
                logger.info("  Skipped %s (insufficient data)", ds_id[:8])
                continue

            X, obs_df, gene_names = result
            ground_truth = obs_df["atlas_type"]
            n_cells = len(obs_df)

            logger.info(
                "  Fetched: %d cells, %d genes, types: %s",
                n_cells, len(gene_names),
                sorted(ground_truth.unique()),
            )

            # ── Tool 1: SingleR ──
            if "SingleR" not in tools_done:
                logger.info("  Running SingleR...")
                singler_preds = run_singler(X, obs_df, gene_names, organism)
                if singler_preds is not None:
                    metrics = compute_metrics(ground_truth, singler_preds)
                    for m in metrics:
                        record = {
                            "dataset_id": ds_id,
                            "tissue": row.tissue,
                            "organism": organism,
                            "tool": "SingleR",
                            **m,
                        }
                        all_results.append(record)
                    done_keys.add((ds_id, organism, "SingleR"))
                    logger.info(
                        "    SingleR: %d cell types scored",
                        len(metrics),
                    )
                else:
                    logger.info("    SingleR: skipped (failed)")

            # ── Tool 2: CellTypist ──
            if "CellTypist" not in tools_done:
                logger.info("  Running CellTypist...")
                ct_preds = run_celltypist(X, obs_df, gene_names)
                if ct_preds is not None:
                    metrics = compute_metrics(ground_truth, ct_preds)
                    for m in metrics:
                        record = {
                            "dataset_id": ds_id,
                            "tissue": row.tissue,
                            "organism": organism,
                            "tool": "CellTypist",
                            **m,
                        }
                        all_results.append(record)
                    done_keys.add((ds_id, organism, "CellTypist"))
                    logger.info(
                        "    CellTypist: %d cell types scored",
                        len(metrics),
                    )
                else:
                    logger.info("    CellTypist: skipped (failed)")

            # ── Tool 3: Marker-score classifier ──
            if "MarkerScore" not in tools_done:
                logger.info("  Running marker-score classifier...")
                ms_preds = run_marker_score(
                    X, obs_df, gene_names, atlas_df, organism,
                )
                if ms_preds is not None:
                    metrics = compute_metrics(ground_truth, ms_preds)
                    for m in metrics:
                        record = {
                            "dataset_id": ds_id,
                            "tissue": row.tissue,
                            "organism": organism,
                            "tool": "MarkerScore",
                            **m,
                        }
                        all_results.append(record)
                    done_keys.add((ds_id, organism, "MarkerScore"))
                    logger.info(
                        "    MarkerScore: %d cell types scored",
                        len(metrics),
                    )
                else:
                    logger.info("    MarkerScore: skipped (failed)")

            # Checkpoint after each dataset
            pd.DataFrame(all_results).to_csv(
                checkpoint_path, sep="\t", index=False,
            )

            del X, obs_df, gene_names
            gc.collect()

    census.close()

    # Save final output
    result_df = pd.DataFrame(all_results)
    out_path = OUTPUT_DIR / "p5_annotation_tools.tsv"
    result_df.to_csv(out_path, sep="\t", index=False)
    logger.info("\nSaved %d rows to %s", len(result_df), out_path)

    # Summary
    elapsed = time.time() - t0
    logger.info("\n=== P5 COMPLETE ===")
    logger.info("Total rows: %d", len(result_df))
    if len(result_df) > 0:
        for tool, tdf in result_df.groupby("tool", observed=True):
            mean_f1 = tdf.f1.mean()
            mean_prec = tdf.precision.mean()
            mean_rec = tdf.recall.mean()
            n_ds = tdf.dataset_id.nunique()
            logger.info(
                "  %s: mean F1=%.3f, P=%.3f, R=%.3f (%d datasets)",
                tool, mean_f1, mean_prec, mean_rec, n_ds,
            )
    logger.info("Elapsed: %.1f minutes", elapsed / 60)


if __name__ == "__main__":
    main()
