#!/usr/bin/env python3
"""P8 — Integration Method Selection Audit.

Question: Which batch correction method is best given specific dataset
characteristics, and can we build a decision rule?

Strategy:
- Wait for P3 batch correction results (dependency)
- If available, load P3 results as base dataset with ARI scores
- Augment with edge-case Census datasets (very few / many donors)
- Run Harmony + Scanorama on edge cases (skip scVI for speed)
- Record dataset features: n_donors, n_batches, n_cells, n_cell_types, batch_size_variance
- Identify best_method per dataset (highest ARI)
- Fit DecisionTreeClassifier (max_depth=4) predicting best_method from features
- Serialize tree rules to JSON
"""

from __future__ import annotations

import gc
import json
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
from scipy.sparse import issparse
from sklearn.metrics import adjusted_rand_score
from sklearn.tree import DecisionTreeClassifier, export_text
from tqdm import tqdm

warnings.filterwarnings("ignore", category=FutureWarning)

CENSUS_VERSION = "2025-11-08"
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = REPO_ROOT / "scrnaseq" / "output"
PROJECT_DIR = Path(__file__).resolve().parent
P3_RESULTS = REPO_ROOT / "scrnaseq" / "output" / "p3_batch_correction.tsv"

MAX_CELLS_PER_DATASET = 30_000
MIN_CELLS_PER_DATASET = 500
RNG_SEED = 42
P3_WAIT_INTERVAL = 300  # 5 minutes
P3_MAX_RETRIES = 20  # 100 minutes total

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(PROJECT_DIR / "p8_run.log"),
    ],
)
logger = logging.getLogger(__name__)


def _census_key(organism: str) -> str:
    return organism.lower().replace(" ", "_")


def wait_for_p3() -> pd.DataFrame | None:
    """Wait for P3 results to become available."""
    for attempt in range(1, P3_MAX_RETRIES + 1):
        if P3_RESULTS.exists():
            df = pd.read_csv(P3_RESULTS, sep="\t")
            if len(df) > 0:
                logger.info("P3 results loaded: %d rows from %s", len(df), P3_RESULTS)
                return df
        logger.info(
            "P3 results not ready (attempt %d/%d). Retrying in %d seconds...",
            attempt, P3_MAX_RETRIES, P3_WAIT_INTERVAL,
        )
        time.sleep(P3_WAIT_INTERVAL)

    logger.warning(
        "P3 results not available after %d minutes. Proceeding with edge-case datasets only.",
        P3_MAX_RETRIES * P3_WAIT_INTERVAL // 60,
    )
    return None


def discover_edge_case_datasets(census, category: str) -> pd.DataFrame:
    """Find datasets targeting specific edge cases.

    category: 'few_donors' (n=2) or 'many_donors' (n>10)
    """
    organism = "Homo sapiens"
    exp = census["census_data"][_census_key(organism)]
    logger.info("Discovering %s datasets for %s...", category, organism)

    obs_df = exp.obs.read(
        column_names=["dataset_id", "tissue_general", "donor_id"],
        value_filter="is_primary_data == True",
    ).concat().to_pandas()

    # Skip NaN dataset_ids
    obs_df = obs_df.dropna(subset=["dataset_id"])

    # Count donors per dataset
    donor_counts = (
        obs_df.groupby("dataset_id", observed=True)["donor_id"]
        .nunique()
        .reset_index(name="n_donors")
    )
    cell_counts = (
        obs_df.groupby("dataset_id", observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    dataset_info = donor_counts.merge(cell_counts, on="dataset_id")
    dataset_info = dataset_info[dataset_info.n_cells >= MIN_CELLS_PER_DATASET]

    # Get tissue for each dataset (most common)
    tissue_map = (
        obs_df.groupby("dataset_id", observed=True)["tissue_general"]
        .agg(lambda x: x.mode().iloc[0] if len(x.mode()) > 0 else x.iloc[0])
        .reset_index(name="tissue")
    )
    dataset_info = dataset_info.merge(tissue_map, on="dataset_id")

    if category == "few_donors":
        subset = dataset_info[dataset_info.n_donors == 2]
        # Sort by n_cells descending, take top 3
        subset = subset.sort_values("n_cells", ascending=False).head(3)
    elif category == "many_donors":
        subset = dataset_info[dataset_info.n_donors > 10]
        subset = subset.sort_values("n_donors", ascending=False).head(3)
    else:
        raise ValueError(f"Unknown category: {category}")

    logger.info(
        "%s: selected %d datasets (donors: %s)",
        category, len(subset),
        subset.n_donors.tolist(),
    )
    del obs_df
    gc.collect()
    return subset.reset_index(drop=True)


def run_integration(adata, method: str, batch_key: str = "donor_id") -> np.ndarray | None:
    """Run a single integration method. Returns corrected embedding or None."""
    try:
        ad = adata.copy()

        # Normalize + HVG + PCA
        sc.pp.normalize_total(ad, target_sum=1e4)
        sc.pp.log1p(ad)
        # Use batch_key for HVG if few batches, otherwise skip it
        n_batches = ad.obs[batch_key].nunique()
        try:
            if n_batches <= 50:
                sc.pp.highly_variable_genes(ad, n_top_genes=2000, batch_key=batch_key)
            else:
                sc.pp.highly_variable_genes(ad, n_top_genes=2000)
        except Exception:
            sc.pp.highly_variable_genes(ad, n_top_genes=2000)
        ad = ad[:, ad.var.highly_variable].copy()
        sc.pp.scale(ad, max_value=10)
        sc.tl.pca(ad, n_comps=30, random_state=RNG_SEED)

        if method == "harmony":
            import harmonypy
            ho = harmonypy.run_harmony(
                ad.obsm["X_pca"], ad.obs, batch_key,
                random_state=RNG_SEED, max_iter_harmony=20,
            )
            ad.obsm["X_corrected"] = ho.Z_corr.T
        elif method == "scanorama":
            import scanorama
            batches = ad.obs[batch_key].unique().tolist()
            adata_list = [ad[ad.obs[batch_key] == b].copy() for b in batches]
            corrected = scanorama.correct_scanpy(adata_list, return_dimred=True)
            # Reassemble
            import scipy.sparse as sp
            ad_merged = sc.concat(corrected)
            ad.obsm["X_corrected"] = ad_merged.obsm["X_scanorama"]
        elif method == "bbknn":
            import bbknn
            bbknn.bbknn(ad, batch_key=batch_key)
            # BBKNN modifies the neighbor graph directly; use PCA as embedding
            ad.obsm["X_corrected"] = ad.obsm["X_pca"]
        elif method == "scvi":
            import scvi as scvi_module
            ad_raw = adata.copy()
            sc.pp.normalize_total(ad_raw, target_sum=1e4)
            sc.pp.log1p(ad_raw)
            sc.pp.highly_variable_genes(ad_raw, n_top_genes=2000, batch_key=batch_key)
            ad_raw = ad_raw[:, ad_raw.var.highly_variable].copy()
            scvi_module.model.SCVI.setup_anndata(ad_raw, batch_key=batch_key)
            model = scvi_module.model.SCVI(ad_raw, n_latent=30)
            model.train(max_epochs=100, early_stopping=True, train_size=0.9)
            ad.obsm["X_corrected"] = model.get_latent_representation()
        else:
            raise ValueError(f"Unknown method: {method}")

        # Cluster on corrected embedding
        from sklearn.neighbors import NearestNeighbors

        rep_key = "X_corrected"
        nn = NearestNeighbors(n_neighbors=15, metric="euclidean")
        nn.fit(ad.obsm[rep_key])
        distances, indices = nn.kneighbors(ad.obsm[rep_key])

        # Build connectivity graph for Leiden
        import igraph as ig
        import leidenalg

        n_cells = ad.shape[0]
        sources = np.repeat(np.arange(n_cells), 15)
        targets = indices.ravel()
        weights = 1.0 / (1.0 + distances.ravel())
        g = ig.Graph(n=n_cells, edges=list(zip(sources.tolist(), targets.tolist())),
                     directed=False)
        g.es["weight"] = weights.tolist()
        g = g.simplify(combine_edges="max")

        partition = leidenalg.find_partition(
            g, leidenalg.RBConfigurationVertexPartition,
            weights="weight", seed=RNG_SEED, resolution_parameter=0.5,
        )
        ad.obs["leiden_corrected"] = pd.Categorical([str(x) for x in partition.membership])

        return ad.obs["leiden_corrected"].values

    except Exception as e:
        logger.warning("  %s failed: %s", method, e)
        return None


def compute_ari(labels_true, labels_pred) -> float:
    """Compute ARI between true cell type labels and integration-derived clusters."""
    mask = ~pd.isna(labels_true) & ~pd.isna(labels_pred)
    if mask.sum() < 10:
        return np.nan
    return adjusted_rand_score(labels_true[mask], labels_pred[mask])


def compute_dataset_features(
    census, dataset_id: str, organism: str = "Homo sapiens",
) -> dict | None:
    """Compute dataset-level features for the decision tree."""
    exp = census["census_data"][_census_key(organism)]

    obs_df = exp.obs.read(
        column_names=["soma_joinid", "dataset_id", "tissue_general",
                      "donor_id", "cell_type"],
        value_filter=f"dataset_id == '{dataset_id}' and is_primary_data == True",
    ).concat().to_pandas()

    if len(obs_df) < MIN_CELLS_PER_DATASET:
        return None

    n_cells = len(obs_df)
    n_donors = obs_df["donor_id"].nunique()
    n_batches = obs_df["donor_id"].nunique()  # Using donor as batch proxy
    n_cell_types = obs_df["cell_type"].nunique()
    tissue = obs_df["tissue_general"].mode().iloc[0] if len(obs_df) > 0 else "unknown"

    # Batch size variance (coefficient of variation of cells per donor)
    batch_sizes = obs_df.groupby("donor_id", observed=True).size()
    batch_size_variance = float(batch_sizes.std() / batch_sizes.mean()) if batch_sizes.mean() > 0 else 0.0

    return {
        "dataset_id": dataset_id,
        "tissue": tissue,
        "organism": organism,
        "n_donors": n_donors,
        "n_batches": n_batches,
        "n_cells": n_cells,
        "n_cell_types": n_cell_types,
        "batch_size_variance": round(batch_size_variance, 4),
    }


def process_edge_dataset(
    census, dataset_id: str, organism: str = "Homo sapiens",
    methods: list[str] | None = None,
) -> dict | None:
    """Fetch data and run integration methods on an edge-case dataset."""
    if methods is None:
        methods = ["scanorama", "bbknn"]

    exp = census["census_data"][_census_key(organism)]

    obs_df = exp.obs.read(
        column_names=["soma_joinid", "donor_id", "cell_type"],
        value_filter=f"dataset_id == '{dataset_id}' and is_primary_data == True",
    ).concat().to_pandas()

    if len(obs_df) < MIN_CELLS_PER_DATASET:
        return None

    # Subsample if too large
    if len(obs_df) > MAX_CELLS_PER_DATASET:
        rng = np.random.default_rng(RNG_SEED)
        idx = rng.choice(len(obs_df), MAX_CELLS_PER_DATASET, replace=False)
        obs_df = obs_df.iloc[sorted(idx)].copy()

    obs_ids = obs_df.soma_joinid.tolist()

    try:
        with exp.axis_query(
            measurement_name="RNA",
            obs_query=tiledbsoma.AxisQuery(coords=(obs_ids,)),
        ) as query:
            adata = query.to_anndata(X_name="raw")
    except Exception as e:
        logger.warning("Fetch failed for %s: %s", dataset_id[:8], e)
        return None

    if adata.shape[0] == 0:
        return None

    # Attach metadata
    adata.obs["donor_id"] = obs_df["donor_id"].values
    adata.obs["cell_type"] = obs_df["cell_type"].values

    # Need at least 2 donors for integration
    if adata.obs["donor_id"].nunique() < 2:
        logger.info("  %s: only 1 donor, skipping integration", dataset_id[:8])
        return None

    cell_types = adata.obs["cell_type"].values

    result = {}
    for method in methods:
        logger.info("  Running %s on %s (%d cells)...", method, dataset_id[:8], adata.shape[0])
        clusters = run_integration(adata, method, batch_key="donor_id")
        if clusters is not None:
            ari = compute_ari(cell_types, clusters)
            result[f"{method}_ari"] = round(ari, 4) if not np.isnan(ari) else np.nan
        else:
            result[f"{method}_ari"] = np.nan

    return result


def extract_tree_rules(tree: DecisionTreeClassifier, feature_names: list[str],
                       class_names: list[str]) -> dict:
    """Extract decision tree rules as JSON-serializable dict."""
    tree_text = export_text(tree, feature_names=feature_names)

    rules = {
        "model": "DecisionTreeClassifier",
        "max_depth": tree.max_depth,
        "n_features": len(feature_names),
        "feature_names": feature_names,
        "class_names": class_names,
        "tree_text": tree_text,
        "feature_importances": {
            name: round(float(imp), 4)
            for name, imp in zip(feature_names, tree.feature_importances_)
            if imp > 0
        },
        "rules": [],
    }

    # Extract leaf-level rules
    tree_ = tree.tree_
    feature = tree_.feature
    threshold = tree_.threshold
    children_left = tree_.children_left
    children_right = tree_.children_right
    value = tree_.value

    def _recurse(node_id, conditions):
        if children_left[node_id] == children_right[node_id]:
            # Leaf node
            class_idx = int(np.argmax(value[node_id]))
            n_samples = int(value[node_id].sum())
            rules["rules"].append({
                "conditions": list(conditions),
                "prediction": class_names[class_idx],
                "n_samples": n_samples,
                "confidence": round(float(value[node_id][0, class_idx] / n_samples), 3),
            })
            return

        feat_name = feature_names[feature[node_id]]
        thresh = round(float(threshold[node_id]), 4)

        # Left child: feature <= threshold
        _recurse(children_left[node_id],
                 conditions + [f"{feat_name} <= {thresh}"])
        # Right child: feature > threshold
        _recurse(children_right[node_id],
                 conditions + [f"{feat_name} > {thresh}"])

    _recurse(0, [])
    return rules


def normalize_p3_results(p3_df: pd.DataFrame) -> pd.DataFrame:
    """Normalize P3 batch correction results into our schema.

    P3 may have columns like: dataset_id, tissue, organism, method, ari, ...
    We need per-dataset rows with harmony_ari, scanorama_ari, scvi_ari columns.
    """
    # Determine actual column names in P3
    cols = p3_df.columns.tolist()
    logger.info("P3 columns: %s", cols)

    # Try to pivot if P3 has a 'method' column with per-method rows
    # Handle both 'ari' and 'ARI' column names
    ari_col = "ARI" if "ARI" in cols else "ari" if "ari" in cols else None
    if "method" in cols and ari_col:
        id_cols = [c for c in ["dataset_id", "tissue", "organism"] if c in cols]
        pivoted = p3_df.pivot_table(
            index=id_cols, columns="method", values=ari_col, aggfunc="first",
        ).reset_index()
        pivoted.columns = [
            f"{c.lower()}_ari" if c in ["harmony", "Harmony", "scanorama", "Scanorama",
                                         "scvi", "scVI"] else c
            for c in pivoted.columns
        ]
        return pivoted

    # If P3 already has harmony_ari, scanorama_ari, etc. columns, use as-is
    ari_cols = [c for c in cols if c.endswith("_ari")]
    if ari_cols:
        return p3_df

    logger.warning("Could not parse P3 results format. Columns: %s", cols)
    return pd.DataFrame()


def main():
    t0 = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PROJECT_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("=== P8: Integration Method Selection Audit ===")

    # Checkpoint
    checkpoint_path = PROJECT_DIR / "p8_checkpoint.tsv"
    if checkpoint_path.exists():
        existing = pd.read_csv(checkpoint_path, sep="\t")
        all_results = existing.to_dict("records")
        done_ids = {r["dataset_id"] for r in all_results}
        logger.info("Resuming from checkpoint: %d rows", len(all_results))
    else:
        all_results = []
        done_ids = set()

    # ── Step 1: Wait for P3 ──
    p3_df = wait_for_p3()
    p3_available = p3_df is not None and len(p3_df) > 0

    # ── Step 2: Open Census ──
    census = cellxgene_census.open_soma(census_version=CENSUS_VERSION)

    if p3_available:
        logger.info("── Loading P3 base results ──")
        p3_norm = normalize_p3_results(p3_df)

        if len(p3_norm) > 0:
            for _, row in p3_norm.iterrows():
                ds_id = row.get("dataset_id")
                if pd.isna(ds_id) or ds_id in done_ids:
                    continue

                # Compute dataset features from Census
                features = compute_dataset_features(census, ds_id)
                if features is None:
                    continue

                record = {**features}
                # Pull ARI values from P3
                for method in ["scanorama", "scvi"]:
                    col = f"{method}_ari"
                    record[col] = row.get(col, np.nan)

                # Run BBKNN on this dataset to add a third method
                logger.info("  Running BBKNN on P3 dataset %s...", ds_id[:8])
                organism = features.get("organism", "Homo sapiens")
                edge_result = process_edge_dataset(
                    census, ds_id, organism, methods=["bbknn"],
                )
                if edge_result and "bbknn_ari" in edge_result:
                    record["bbknn_ari"] = edge_result["bbknn_ari"]
                else:
                    record["bbknn_ari"] = np.nan

                all_results.append(record)
                done_ids.add(ds_id)
                logger.info(
                    "  P3 dataset %s: scanorama=%.3f scvi=%.3f bbknn=%.3f",
                    ds_id[:8],
                    record.get("scanorama_ari", np.nan),
                    record.get("scvi_ari", np.nan),
                    record.get("bbknn_ari", np.nan),
                )

            # Checkpoint
            pd.DataFrame(all_results).to_csv(checkpoint_path, sep="\t", index=False)
        else:
            logger.warning("P3 results could not be normalized")

    # ── Step 3: Edge-case datasets ──
    logger.info("── Discovering edge-case datasets ──")

    for category in ["few_donors", "many_donors"]:
        logger.info("\n--- Edge case: %s ---", category)
        try:
            edge_datasets = discover_edge_case_datasets(census, category)
        except Exception as e:
            logger.warning("Failed to discover %s datasets: %s", category, e)
            continue

        for _, row in tqdm(
            edge_datasets.iterrows(), total=len(edge_datasets),
            desc=f"P8 ({category})",
        ):
            ds_id = row["dataset_id"]
            if pd.isna(ds_id) or ds_id in done_ids:
                logger.info("  Skipping %s (done or NaN)", str(ds_id)[:8])
                continue

            # Compute features
            features = compute_dataset_features(census, ds_id)
            if features is None:
                continue

            # Run Scanorama + BBKNN (skip scVI for edge cases)
            integration_results = process_edge_dataset(
                census, ds_id, methods=["scanorama", "bbknn"],
            )

            record = {**features}
            if integration_results:
                record.update(integration_results)
            else:
                record["scanorama_ari"] = np.nan
                record["bbknn_ari"] = np.nan
            record["scvi_ari"] = np.nan  # Skipped for edge cases

            all_results.append(record)
            done_ids.add(ds_id)
            logger.info(
                "  Edge %s (%s): %d donors, scanorama=%.3f bbknn=%.3f",
                ds_id[:8], category,
                features["n_donors"],
                record.get("scanorama_ari", np.nan),
                record.get("bbknn_ari", np.nan),
            )

            # Checkpoint
            pd.DataFrame(all_results).to_csv(checkpoint_path, sep="\t", index=False)
            gc.collect()

    census.close()

    # ── Step 4: Determine best method per dataset ──
    result_df = pd.DataFrame(all_results)

    if len(result_df) == 0:
        logger.error("No results to process. Exiting.")
        sys.exit(1)

    ari_cols = [c for c in ["scanorama_ari", "scvi_ari", "bbknn_ari"]
                if c in result_df.columns]
    for col in ari_cols:
        if col not in result_df.columns:
            result_df[col] = np.nan

    # Best method = highest ARI
    ari_matrix = result_df[ari_cols].copy()
    method_names = [c.replace("_ari", "") for c in ari_cols]
    result_df["best_method"] = ari_matrix.idxmax(axis=1).str.replace("_ari", "")
    # Where all are NaN, best_method is NaN
    all_nan = ari_matrix.isna().all(axis=1)
    result_df.loc[all_nan, "best_method"] = np.nan

    # ── Step 5: Decision tree ──
    feature_cols = ["n_donors", "n_batches", "n_cells", "n_cell_types", "batch_size_variance"]
    # Ensure all feature columns exist
    for fc in feature_cols:
        if fc not in result_df.columns:
            result_df[fc] = np.nan

    # Filter to rows with valid best_method and features
    tree_df = result_df.dropna(subset=["best_method"] + feature_cols).copy()
    tree_rules = {}

    if len(tree_df) >= 5:
        X = tree_df[feature_cols].values
        y = tree_df["best_method"].values
        class_names = sorted(set(y))

        clf = DecisionTreeClassifier(max_depth=4, random_state=RNG_SEED, min_samples_leaf=2)
        clf.fit(X, y)

        # Predictions
        tree_df["decision_tree_prediction"] = clf.predict(X)
        tree_df["prediction_correct"] = tree_df["decision_tree_prediction"] == tree_df["best_method"]

        # Merge predictions back
        result_df = result_df.merge(
            tree_df[["dataset_id", "decision_tree_prediction", "prediction_correct"]],
            on="dataset_id", how="left",
        )

        # Extract rules
        tree_rules = extract_tree_rules(clf, feature_cols, class_names)
        accuracy = tree_df["prediction_correct"].mean()
        tree_rules["accuracy"] = round(float(accuracy), 4)
        tree_rules["n_training_samples"] = len(tree_df)

        logger.info("Decision tree accuracy: %.1f%% (%d samples)", accuracy * 100, len(tree_df))
        logger.info("Feature importances: %s", tree_rules["feature_importances"])
        logger.info("Tree rules:\n%s", tree_rules["tree_text"])
    else:
        logger.warning("Not enough valid rows (%d) for decision tree. Need >= 5.", len(tree_df))
        result_df["decision_tree_prediction"] = np.nan
        result_df["prediction_correct"] = np.nan
        tree_rules = {"error": "insufficient_data", "n_available": len(tree_df)}

    # ── Step 6: Output ──
    # Ensure output columns in correct order
    output_cols = [
        "dataset_id", "tissue", "organism", "n_donors", "n_batches", "n_cells",
        "n_cell_types", "best_method", "harmony_ari", "scanorama_ari", "scvi_ari",
        "decision_tree_prediction", "prediction_correct",
    ]
    for col in output_cols:
        if col not in result_df.columns:
            result_df[col] = np.nan

    out_path = OUTPUT_DIR / "p8_integration_selection.tsv"
    result_df[output_cols].to_csv(out_path, sep="\t", index=False)
    logger.info("Saved %d rows to %s", len(result_df), out_path)

    tree_path = OUTPUT_DIR / "p8_decision_tree.json"
    with open(tree_path, "w") as f:
        json.dump(tree_rules, f, indent=2)
    logger.info("Saved decision tree rules to %s", tree_path)

    # ── Summary ──
    elapsed = time.time() - t0
    logger.info("\n=== P8 COMPLETE ===")
    logger.info("Total datasets: %d", len(result_df))
    if "best_method" in result_df.columns:
        counts = result_df["best_method"].value_counts()
        for method, count in counts.items():
            logger.info("  %s: %d datasets (%.0f%%)", method, count, 100 * count / len(result_df))
    logger.info("P3 base results: %s", "loaded" if p3_available else "NOT AVAILABLE")
    logger.info("Edge-case datasets processed: %d",
                len([r for r in all_results if r.get("scvi_ari") is np.nan or pd.isna(r.get("scvi_ari", 0))]))
    logger.info("Elapsed: %.1f minutes", elapsed / 60)


if __name__ == "__main__":
    main()
