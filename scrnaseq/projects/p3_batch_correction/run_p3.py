#!/usr/bin/env python3
"""P3 — Batch Correction Consistency Audit.

Question: Do Harmony, Scanorama, and scVI produce consistent cell type
neighborhoods?

Strategy:
- Find Census datasets with >=2 donors (batch = donor_id), human + mouse
- For each dataset: log-normalize, PCA (50 PCs)
- Run 3 batch correction methods on the same PCA embedding
- Evaluate each with ARI, batch LISI, kBET acceptance rate
- Compare consistency across methods and datasets
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
from scipy.spatial import cKDTree
from sklearn.metrics import adjusted_rand_score
from tqdm import tqdm

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ── Config ──────────────────────────────────────────────────────────
CENSUS_VERSION = "2025-11-08"
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = REPO_ROOT / "scrnaseq" / "output"
PROJECT_DIR = Path(__file__).resolve().parent

MAX_DATASETS = 15
MAX_CELLS_PER_DATASET = 20_000
MAX_CELLS_SCVI = 50_000
MIN_CELLS_PER_DATASET = 500
MIN_DONORS = 2
N_PCS = 50
LEIDEN_RES = 0.5
LISI_K = 30
KBET_K = 25
KBET_N_SUBSETS = 50
RNG_SEED = 42

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(PROJECT_DIR / "p3_run.log"),
    ],
)
logger = logging.getLogger(__name__)


def _census_key(organism: str) -> str:
    return organism.lower().replace(" ", "_")


# ── Dataset Discovery ──────────────────────────────────────────────

def discover_eligible_datasets(census, organism: str) -> pd.DataFrame:
    """Find datasets with >=2 donors and known batch structure."""
    exp = census["census_data"][_census_key(organism)]
    logger.info("Discovering eligible datasets for %s...", organism)

    obs_df = exp.obs.read(
        column_names=["dataset_id", "tissue_general", "donor_id", "cell_type"],
        value_filter="is_primary_data == True",
    ).concat().to_pandas()

    # Group by dataset, count donors
    ds_stats = obs_df.groupby("dataset_id", observed=True).agg(
        tissue=("tissue_general", lambda x: x.mode().iloc[0] if len(x.mode()) > 0 else "unknown"),
        n_cells=("donor_id", "size"),
        n_donors=("donor_id", "nunique"),
        n_cell_types=("cell_type", "nunique"),
    ).reset_index()

    # Filter: >=2 donors, enough cells, multiple cell types for ARI
    ds_stats = ds_stats[
        (ds_stats.n_donors >= MIN_DONORS)
        & (ds_stats.n_cells >= MIN_CELLS_PER_DATASET)
        & (ds_stats.n_cell_types >= 3)
    ]

    # Skip NaN dataset_ids
    ds_stats = ds_stats.dropna(subset=["dataset_id"])

    # Prioritize tissue diversity: sort by n_donors (descending) then n_cells
    ds_stats = ds_stats.sort_values(
        ["n_donors", "n_cells"], ascending=[False, False],
    )

    # Cap at MAX_DATASETS, ensuring tissue diversity
    selected = []
    seen_tissues = set()
    # First pass: one per tissue
    for _, row in ds_stats.iterrows():
        if row.tissue not in seen_tissues:
            selected.append(row)
            seen_tissues.add(row.tissue)
        if len(selected) >= MAX_DATASETS:
            break
    # Second pass: fill remaining slots
    if len(selected) < MAX_DATASETS:
        for _, row in ds_stats.iterrows():
            if row.dataset_id not in {s.dataset_id for s in selected}:
                selected.append(row)
            if len(selected) >= MAX_DATASETS:
                break

    result = pd.DataFrame(selected).reset_index(drop=True)
    logger.info(
        "%s: found %d eligible datasets across %d tissues",
        organism, len(result), result.tissue.nunique(),
    )

    del obs_df
    gc.collect()
    return result


# ── Data Fetching ──────────────────────────────────────────────────

def fetch_dataset(
    census, dataset_id: str, organism: str,
) -> sc.AnnData | None:
    """Fetch raw counts for a dataset, subsample to MAX_CELLS_PER_DATASET."""
    exp = census["census_data"][_census_key(organism)]

    obs_filter = f"dataset_id == '{dataset_id}' and is_primary_data == True"

    obs_df = exp.obs.read(
        column_names=["soma_joinid", "donor_id", "cell_type"],
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

    try:
        with exp.axis_query(
            measurement_name="RNA",
            obs_query=tiledbsoma.AxisQuery(coords=(obs_ids,)),
        ) as query:
            adata = query.to_anndata(X_name="raw")
    except Exception as e:
        logger.warning("Failed to fetch expression for %s: %s", dataset_id[:8], e)
        return None

    if adata.shape[0] == 0 or adata.shape[1] < 500:
        return None

    # Align metadata
    adata_joinids = adata.obs["soma_joinid"].values
    obs_indexed = obs_df.set_index("soma_joinid")
    aligned = obs_indexed.loc[adata_joinids]
    adata.obs["donor_id"] = aligned["donor_id"].values
    adata.obs["cell_type"] = aligned["cell_type"].values

    return adata


# ── Preprocessing ──────────────────────────────────────────────────

def preprocess(adata: sc.AnnData) -> sc.AnnData:
    """Log-normalize, HVG selection, PCA (50 PCs)."""
    adata = adata.copy()

    # Store raw counts for scVI
    adata.layers["counts"] = adata.X.copy()

    # Log-normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # HVG selection
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3",
                                layer="counts")
    adata_hvg = adata[:, adata.var.highly_variable].copy()

    # Scale
    sc.pp.scale(adata_hvg, max_value=10)

    # PCA
    sc.tl.pca(adata_hvg, n_comps=N_PCS, random_state=RNG_SEED)

    # Copy PCA back to full adata
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.uns["pca"] = adata_hvg.uns.get("pca", {})

    del adata_hvg
    return adata


# ── Metrics ────────────────────────────────────────────────────────

def compute_lisi(embedding: np.ndarray, labels: np.ndarray, k: int = LISI_K) -> float:
    """Compute batch LISI (Local Inverse Simpson Index).

    For each cell, get k nearest neighbors, compute Simpson diversity
    of batch labels, invert. Higher LISI = better batch mixing.
    """
    tree = cKDTree(embedding)
    _, indices = tree.query(embedding, k=k + 1)  # +1 because query includes self
    indices = indices[:, 1:]  # remove self

    lisi_scores = []
    for i in range(len(labels)):
        neighbor_labels = labels[indices[i]]
        unique, counts = np.unique(neighbor_labels, return_counts=True)
        freqs = counts / counts.sum()
        simpson = np.sum(freqs ** 2)
        lisi_scores.append(1.0 / simpson)

    return float(np.mean(lisi_scores))


def compute_kbet(
    embedding: np.ndarray, batch_labels: np.ndarray,
    k: int = KBET_K, n_subsets: int = KBET_N_SUBSETS,
) -> float:
    """Approximate kBET acceptance rate.

    For random subsets of cells, test if batch proportions in k-neighborhood
    match global proportions (chi-square test). Report fraction accepted at p>0.05.
    """
    rng = np.random.default_rng(RNG_SEED)
    tree = cKDTree(embedding)

    # Global batch proportions
    unique_batches, global_counts = np.unique(batch_labels, return_counts=True)
    global_freqs = global_counts / global_counts.sum()

    n_cells = len(batch_labels)
    test_cells = rng.choice(n_cells, min(n_subsets, n_cells), replace=False)

    accepted = 0
    for cell_idx in test_cells:
        _, nn_idx = tree.query(embedding[cell_idx], k=k + 1)
        nn_idx = nn_idx[1:]  # remove self

        nn_labels = batch_labels[nn_idx]
        observed = np.array([np.sum(nn_labels == b) for b in unique_batches])
        expected = global_freqs * k

        # Skip if any expected < 1 (chi-square unreliable)
        if np.any(expected < 1):
            # Merge small categories for this test
            mask = expected >= 1
            if mask.sum() < 2:
                accepted += 1  # can't test, assume accepted
                continue
            observed = observed[mask]
            expected = expected[mask]

        try:
            _, pval = scipy_stats.chisquare(observed, f_exp=expected)
            if pval > 0.05:
                accepted += 1
        except Exception:
            accepted += 1  # degenerate case, assume accepted

    return float(accepted / len(test_cells))


def compute_ari(embedding: np.ndarray, cell_types: np.ndarray) -> float:
    """Leiden cluster on embedding, compute ARI vs ground truth cell_type."""
    # Build temporary AnnData for scanpy Leiden
    adata_tmp = sc.AnnData(X=np.zeros((len(cell_types), 1)))
    adata_tmp.obsm["X_emb"] = embedding
    sc.pp.neighbors(adata_tmp, use_rep="X_emb", n_neighbors=15, random_state=RNG_SEED)
    sc.tl.leiden(adata_tmp, resolution=LEIDEN_RES, random_state=RNG_SEED)

    leiden_labels = adata_tmp.obs["leiden"].values.astype(str)
    ari = adjusted_rand_score(cell_types, leiden_labels)

    del adata_tmp
    return float(ari)


# ── Batch Correction Methods ──────────────────────────────────────

def run_harmony(adata: sc.AnnData) -> np.ndarray | None:
    """Run Harmony on PCA embedding."""
    try:
        import harmonypy

        pca = adata.obsm["X_pca"].copy()
        # Ensure donor_id is string type to avoid tensor size mismatches
        meta = adata.obs[["donor_id"]].copy()
        meta["donor_id"] = meta["donor_id"].astype(str)

        # Cap nclust and force NumPy backend to avoid PyTorch tensor dimension bugs
        n_donors = meta["donor_id"].nunique()
        nclust = min(100, max(5, n_donors))
        ho = harmonypy.run_harmony(pca, meta, "donor_id",
                                   random_state=RNG_SEED, nclust=nclust,
                                   torch_device="cpu")
        return ho.Z_corr.T  # (n_cells, n_pcs)
    except Exception as e:
        logger.warning("Harmony failed: %s", e)
        return None


def run_scanorama(adata: sc.AnnData) -> np.ndarray | None:
    """Run Scanorama on per-batch expression matrices."""
    try:
        import scanorama

        donors = adata.obs["donor_id"].unique()
        # Get HVG indices for consistent feature space
        if "highly_variable" in adata.var.columns:
            hvg_mask = adata.var.highly_variable.values
        else:
            hvg_mask = np.ones(adata.shape[1], dtype=bool)

        gene_names = list(adata.var_names[hvg_mask])

        datasets = []
        genes_list = []
        batch_order = []
        for donor in donors:
            mask = adata.obs["donor_id"] == donor
            X_batch = adata[mask, hvg_mask].X
            if hasattr(X_batch, "toarray"):
                X_batch = X_batch.toarray()
            else:
                X_batch = np.asarray(X_batch)
            datasets.append(X_batch)
            genes_list.append(gene_names)
            batch_order.append(donor)

        corrected, _ = scanorama.correct(
            datasets, genes_list, return_dimred=False,
        )

        # Reconstruct in original order
        result = np.zeros((adata.shape[0], corrected[0].shape[1]), dtype=np.float32)
        idx = 0
        for i, donor in enumerate(donors):
            mask = adata.obs["donor_id"] == donor
            n = mask.sum()
            arr = corrected[i]
            if hasattr(arr, "toarray"):
                arr = arr.toarray()
            result[np.where(mask)[0]] = arr
            idx += n

        # PCA on corrected space for embedding
        from sklearn.decomposition import PCA
        pca = PCA(n_components=min(N_PCS, result.shape[1]), random_state=RNG_SEED)
        return pca.fit_transform(result)

    except Exception as e:
        logger.warning("Scanorama failed: %s", e)
        return None


def run_scvi(adata: sc.AnnData) -> np.ndarray | None:
    """Run scVI on raw counts."""
    try:
        import scvi as scvi_module

        adata_scvi = adata.copy()

        # Cap cells for scVI
        if adata_scvi.shape[0] > MAX_CELLS_SCVI:
            rng = np.random.default_rng(RNG_SEED)
            idx = rng.choice(adata_scvi.shape[0], MAX_CELLS_SCVI, replace=False)
            adata_scvi = adata_scvi[sorted(idx)].copy()

        # Use raw counts
        adata_scvi.X = adata_scvi.layers["counts"].copy()

        # Setup
        scvi_module.model.SCVI.setup_anndata(
            adata_scvi,
            batch_key="donor_id",
            layer=None,  # use .X which we set to counts
        )

        model = scvi_module.model.SCVI(
            adata_scvi,
            n_latent=30,
            gene_likelihood="nb",
        )

        # Train — try GPU first, fall back to CPU
        try:
            model.train(max_epochs=200, early_stopping=True)
        except Exception:
            logger.info("  scVI GPU unavailable, using CPU")
            try:
                model.train(max_epochs=200, early_stopping=True,
                            accelerator="cpu")
            except Exception as e2:
                logger.warning("  scVI CPU training failed: %s", e2)
                return None

        latent = model.get_latent_representation()

        # If we subsampled, we only have embeddings for the subset
        if adata_scvi.shape[0] < adata.shape[0]:
            # Subsampled: return embeddings for the subset only
            logger.info("  scVI ran on %d/%d cells (capped)", adata_scvi.shape[0], adata.shape[0])

        del adata_scvi, model
        gc.collect()
        return latent

    except ImportError:
        logger.warning("scVI not installed, skipping")
        return None
    except Exception as e:
        logger.warning("scVI failed: %s", e)
        return None


# ── Main Pipeline ──────────────────────────────────────────────────

METHODS = {
    "harmony": run_harmony,
    "scanorama": run_scanorama,
    "scvi": run_scvi,
}


def process_dataset(
    census, dataset_id: str, tissue: str, organism: str,
) -> list[dict]:
    """Process one dataset through all 3 methods, return metric rows."""
    adata = fetch_dataset(census, dataset_id, organism)
    if adata is None:
        logger.info("  Skipped %s (fetch returned None)", dataset_id[:8])
        return []

    n_cells = adata.shape[0]
    n_donors = adata.obs["donor_id"].nunique()
    logger.info(
        "  Fetched %s: %d cells, %d genes, %d donors",
        dataset_id[:8], n_cells, adata.shape[1], n_donors,
    )

    # Preprocess
    try:
        adata = preprocess(adata)
    except Exception as e:
        logger.warning("  Preprocessing failed for %s: %s", dataset_id[:8], e)
        return []

    cell_types = adata.obs["cell_type"].values.astype(str)
    batch_labels = adata.obs["donor_id"].values.astype(str)

    results = []
    for method_name, method_fn in METHODS.items():
        logger.info("  Running %s...", method_name)
        t_start = time.time()

        embedding = method_fn(adata)

        runtime = time.time() - t_start

        if embedding is None:
            logger.info("  %s returned None, skipping", method_name)
            continue

        # For scVI with subsampled data, adjust cell_types/batch_labels
        if method_name == "scvi" and embedding.shape[0] < len(cell_types):
            # scVI was capped — use the subset's labels
            n_sub = embedding.shape[0]
            rng = np.random.default_rng(RNG_SEED)
            sub_idx = sorted(rng.choice(len(cell_types), n_sub, replace=False))
            ct_sub = cell_types[sub_idx]
            bl_sub = batch_labels[sub_idx]
        else:
            ct_sub = cell_types
            bl_sub = batch_labels

        # Compute metrics
        try:
            ari = compute_ari(embedding, ct_sub)
        except Exception as e:
            logger.warning("  ARI failed for %s/%s: %s", dataset_id[:8], method_name, e)
            ari = np.nan

        try:
            lisi = compute_lisi(embedding, bl_sub)
        except Exception as e:
            logger.warning("  LISI failed for %s/%s: %s", dataset_id[:8], method_name, e)
            lisi = np.nan

        try:
            kbet = compute_kbet(embedding, bl_sub)
        except Exception as e:
            logger.warning("  kBET failed for %s/%s: %s", dataset_id[:8], method_name, e)
            kbet = np.nan

        record = {
            "dataset_id": dataset_id,
            "tissue": tissue,
            "organism": organism,
            "method": method_name,
            "ARI": ari,
            "LISI": lisi,
            "kBET": kbet,
            "runtime_seconds": round(runtime, 1),
        }
        results.append(record)

        logger.info(
            "  %s: ARI=%.3f, LISI=%.2f, kBET=%.3f, time=%.1fs",
            method_name, ari, lisi, kbet, runtime,
        )

        del embedding
        gc.collect()

    del adata
    gc.collect()
    return results


def main():
    t0 = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PROJECT_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("=== P3: Batch Correction Consistency Audit ===")

    # Resume from checkpoint
    checkpoint_path = PROJECT_DIR / "p3_checkpoint.tsv"
    if checkpoint_path.exists():
        existing = pd.read_csv(checkpoint_path, sep="\t")
        all_results = existing.to_dict("records")
        done_keys = {(r["dataset_id"], r["organism"]) for r in all_results}
        logger.info("Resuming from checkpoint: %d results", len(all_results))
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
            desc=f"P3 ({organism})",
        ):
            ds_id = row.dataset_id
            if pd.isna(ds_id):
                continue

            key = (ds_id, organism)
            if key in done_keys:
                logger.info("  Skipping %s (already done)", ds_id[:8])
                continue

            logger.info(
                "  Processing %s/%s, %d cells, %d donors",
                ds_id[:8], row.tissue, row.n_cells, row.n_donors,
            )

            try:
                records = process_dataset(census, ds_id, row.tissue, organism)
            except Exception as e:
                logger.warning("  Failed on %s: %s", ds_id[:8], e)
                continue

            if records:
                all_results.extend(records)
                done_keys.add(key)

                # Checkpoint after each dataset
                pd.DataFrame(all_results).to_csv(
                    checkpoint_path, sep="\t", index=False,
                )

    census.close()

    # Save final output
    result_df = pd.DataFrame(all_results)
    out_path = OUTPUT_DIR / "p3_batch_correction.tsv"
    result_df.to_csv(out_path, sep="\t", index=False)
    logger.info("\nSaved %d rows to %s", len(result_df), out_path)

    # Summary
    elapsed = time.time() - t0
    logger.info("\n=== P3 COMPLETE ===")
    logger.info("Total rows: %d", len(result_df))
    if len(result_df) > 0:
        for method in ["harmony", "scanorama", "scvi"]:
            mdf = result_df[result_df.method == method]
            if len(mdf) > 0:
                logger.info(
                    "%s: mean ARI=%.3f, mean LISI=%.2f, mean kBET=%.3f, mean time=%.1fs",
                    method,
                    mdf.ARI.mean(),
                    mdf.LISI.mean(),
                    mdf.kBET.mean(),
                    mdf.runtime_seconds.mean(),
                )
    logger.info("Elapsed: %.1f minutes", elapsed / 60)


if __name__ == "__main__":
    main()
