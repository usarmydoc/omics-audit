#!/usr/bin/env python3
"""P7 — Normalization Method Comparison.

Question: Do SCTransform, log-normalize, and scran produce meaningfully
different marker AUC scores?

Strategy:
- Load marker reliability atlas to identify cell types and marker genes
- Query Census for 15+ datasets with cell types matching the atlas
- For each dataset (capped at 15K cells), fetch raw counts
- Apply 3 normalization methods to identical raw counts:
  1. log-normalize (scanpy): normalize_total + log1p
  2. SCTransform (Seurat via rpy2): regularized NB regression
  3. scran (Bioconductor via rpy2): deconvolution size factors + logNormCounts
- Compute per-method: marker AUC, silhouette score, HVG overlap vs lognorm
- Record runtime per method
"""

from __future__ import annotations

import gc
import logging
import sys
import time
import warnings
from pathlib import Path

import anndata as ad
import cellxgene_census
import numpy as np
import pandas as pd
import scanpy as sc
import tiledbsoma
from scipy.sparse import issparse
from sklearn.metrics import roc_auc_score, silhouette_score
from tqdm import tqdm

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

CENSUS_VERSION = "2025-11-08"
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = REPO_ROOT / "scrnaseq" / "output"
PROJECT_DIR = Path(__file__).resolve().parent
ATLAS_PATH = REPO_ROOT / "marker_atlas" / "output" / "marker_reliability_atlas.tsv"

MAX_DATASETS_PER_ORGANISM = 8
MAX_CELLS_PER_DATASET = 15_000
MIN_CELLS_PER_DATASET = 500
MIN_CELLS_PER_TYPE = 30
N_HVG = 2000
N_PCS = 30
RNG_SEED = 42

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(PROJECT_DIR / "p7_run.log"),
    ],
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Marker atlas helpers
# ---------------------------------------------------------------------------

def load_atlas() -> pd.DataFrame:
    """Load marker atlas and keep only reliable markers (AUC > 0.65)."""
    atlas = pd.read_csv(ATLAS_PATH, sep="\t")
    reliable = atlas[atlas.mean_AUC > 0.65].copy()
    logger.info(
        "Marker atlas: %d reliable marker-celltype pairs (AUC > 0.65) from %d total",
        len(reliable), len(atlas),
    )
    return reliable


def atlas_cell_types(atlas: pd.DataFrame, organism: str) -> set[str]:
    """Unique cell types in atlas for an organism."""
    return set(atlas.loc[atlas.organism == organism, "cell_type"].unique())


def atlas_marker_genes(atlas: pd.DataFrame, organism: str) -> set[str]:
    """Unique marker genes in atlas for an organism."""
    return set(atlas.loc[atlas.organism == organism, "gene"].unique())


# ---------------------------------------------------------------------------
# Census dataset discovery
# ---------------------------------------------------------------------------

def _census_key(organism: str) -> str:
    return organism.lower().replace(" ", "_")


def _map_census_cell_type(ct: str) -> str | None:
    """Map Census fine-grained cell_type to atlas category via substring matching."""
    ct_lower = ct.lower()
    keywords = {
        "Hepatocytes": ["hepatocyte", "hepatoblast"],
        "Macrophages": ["macrophage", "kupffer"],
        "T cells": ["t cell", "cd4-positive", "cd8-positive", "regulatory t",
                     "memory t", "naive t", "effector t", "alpha-beta t",
                     "gamma-delta t"],
        "Endothelial": ["endothelial"],
        "Fibroblasts": ["fibroblast", "stellate", "myofibroblast"],
        "NK cells": ["natural killer", "nk cell"],
    }
    for category, kws in keywords.items():
        if any(kw in ct_lower for kw in kws):
            return category
    return None


def discover_datasets(census, organism: str, atlas: pd.DataFrame) -> pd.DataFrame:
    """Find datasets whose cell_type column has overlap with atlas cell types."""
    exp = census["census_data"][_census_key(organism)]
    target_types = atlas_cell_types(atlas, organism)
    logger.info("Discovering datasets for %s (target types: %s)...", organism, target_types)

    obs_df = exp.obs.read(
        column_names=["dataset_id", "tissue_general", "cell_type"],
        value_filter="is_primary_data == True",
    ).concat().to_pandas()

    # Map Census cell types to atlas categories via substring matching
    obs_df["mapped_type"] = obs_df.cell_type.apply(_map_census_cell_type)
    obs_df = obs_df.dropna(subset=["mapped_type"])
    # Keep only types present in atlas for this organism
    obs_df = obs_df[obs_df.mapped_type.isin(target_types)]
    if obs_df.empty:
        logger.warning("No matching cell types for %s", organism)
        return pd.DataFrame()

    # Group by dataset+tissue, count cells with matching types
    grouped = (
        obs_df.groupby(["dataset_id", "tissue_general"], observed=True)
        .agg(
            n_cells=("mapped_type", "size"),
            n_types=("mapped_type", "nunique"),
        )
        .reset_index()
    )
    grouped = grouped.rename(columns={"tissue_general": "tissue"})
    grouped = grouped[grouped.n_cells >= MIN_CELLS_PER_DATASET]

    # Prefer datasets with more cell types for better AUC evaluation
    grouped = grouped.sort_values(["n_types", "n_cells"], ascending=[False, False])

    # Tissue diversity: 1 per tissue
    selected = grouped.drop_duplicates(subset="tissue", keep="first")
    selected = selected.head(MAX_DATASETS_PER_ORGANISM)

    logger.info(
        "%s: selected %d datasets across %d tissues (cell types matched: %s)",
        organism, len(selected), selected.tissue.nunique(),
        sorted(target_types),
    )
    del obs_df
    gc.collect()
    return selected.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Normalization methods
# ---------------------------------------------------------------------------

def _lognorm(adata_raw: ad.AnnData) -> ad.AnnData:
    """Log-normalize: normalize_total + log1p."""
    adata = adata_raw.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def _sctransform(adata_raw: ad.AnnData) -> ad.AnnData | None:
    """SCTransform via rpy2 → Seurat → extract normalized data."""
    try:
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr

        seurat = importr("Seurat")
        seuratobj = importr("SeuratObject")
        base = importr("base")
        matrix_pkg = importr("Matrix")

        # Convert sparse to dense for transfer
        X = adata_raw.X
        if issparse(X):
            X_dense = np.asarray(X.todense())
        else:
            X_dense = np.asarray(X)

        # Filter out zero-UMI cells (SCTransform computes log(UMI) → NaN/Inf)
        cell_sums = np.asarray(X_dense.sum(axis=1)).ravel()
        nonzero_mask = cell_sums > 0
        if nonzero_mask.sum() < 100:
            logger.warning("SCTransform: too few non-zero cells")
            return None
        X_dense = X_dense[nonzero_mask]
        adata_filtered_obs = adata_raw.obs.iloc[nonzero_mask].copy()

        # Transpose: Seurat wants genes x cells
        X_t = X_dense.T

        # Create gene and cell names
        n_genes, n_cells = X_t.shape
        gene_names = list(adata_raw.var_names)
        cell_names = [f"Cell_{i}" for i in range(n_cells)]

        # Create R matrix with INTEGER counts (required by SCTransform)
        X_int = X_t.astype(int)
        r_matrix = ro.r["matrix"](
            ro.IntVector(X_int.flatten(order='F').tolist()),
            nrow=n_genes,
            ncol=n_cells,
        )
        ro.r.assign("r_mat", r_matrix)
        ro.r.assign("gene_names", ro.StrVector(gene_names))
        ro.r.assign("cell_names", ro.StrVector(cell_names))
        ro.r("rownames(r_mat) <- gene_names")
        ro.r("colnames(r_mat) <- cell_names")

        # Create Seurat object and run SCTransform
        # future.globals.maxSize must be increased for large datasets
        ro.r("""
            suppressWarnings(suppressMessages({
                options(future.globals.maxSize = 2 * 1024^3)
                sparse_mat <- Matrix::Matrix(r_mat, sparse = TRUE)
                sobj <- Seurat::CreateSeuratObject(counts = sparse_mat)
                sobj <- Seurat::SCTransform(sobj, verbose = FALSE)
                sct_data <- Seurat::GetAssayData(sobj, assay = "SCT", layer = "data")
                sct_dense <- as.matrix(sct_data)
            }))
        """)

        sct_dense = np.array(ro.r("sct_dense"))  # genes x cells
        sct_genes = list(ro.r("rownames(sct_dense)"))

        # Build anndata: cells x genes (transpose back)
        # Use filtered obs (zero-UMI cells removed)
        adata_sct = ad.AnnData(
            X=sct_dense.T,
            obs=adata_filtered_obs.reset_index(drop=True),
        )
        adata_sct.var_names = pd.Index(sct_genes)

        # Clean R memory
        ro.r("rm(r_mat, sparse_mat, sobj, sct_data, sct_dense); gc()")

        return adata_sct

    except Exception as e:
        logger.warning("SCTransform failed: %s", e)
        return None


def _scran(adata_raw: ad.AnnData) -> ad.AnnData | None:
    """scran via rpy2: computeSumFactors + logNormCounts."""
    try:
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr

        scran = importr("scran")
        scuttle = importr("scuttle")
        sce_pkg = importr("SingleCellExperiment")
        s4vec = importr("S4Vectors")
        base = importr("base")

        X = adata_raw.X
        if issparse(X):
            X_dense = np.asarray(X.todense())
        else:
            X_dense = np.asarray(X)

        # Filter out zero-UMI cells (scran's computeSumFactors fails on them)
        cell_sums = np.asarray(X_dense.sum(axis=1)).ravel()
        nonzero_mask = cell_sums > 0
        if nonzero_mask.sum() < 100:
            logger.warning("scran: too few non-zero cells")
            return None
        X_dense = X_dense[nonzero_mask]
        adata_filtered_obs = adata_raw.obs.iloc[nonzero_mask].copy()

        X_t = X_dense.T  # genes x cells
        n_genes, n_cells = X_t.shape
        gene_names = list(adata_raw.var_names)
        cell_names = [f"Cell_{i}" for i in range(n_cells)]

        # Integer counts required by scran's computeSumFactors
        X_int = X_t.astype(int)
        ro.r.assign("counts_mat", ro.r["matrix"](
            ro.IntVector(X_int.flatten(order='F').tolist()),
            nrow=n_genes,
            ncol=n_cells,
        ))
        ro.r.assign("gene_names", ro.StrVector(gene_names))
        ro.r.assign("cell_names", ro.StrVector(cell_names))
        ro.r("rownames(counts_mat) <- gene_names")
        ro.r("colnames(counts_mat) <- cell_names")

        ro.r("""
            suppressWarnings(suppressMessages({
                sce <- SingleCellExperiment::SingleCellExperiment(
                    assays = list(counts = counts_mat)
                )
                # Quick clustering for size factor estimation
                qclust <- scran::quickCluster(sce, min.size = 50)
                sce <- scran::computeSumFactors(sce, clusters = qclust)
                sce <- scuttle::logNormCounts(sce)
                norm_mat <- as.matrix(SummarizedExperiment::assay(sce, "logcounts"))
            }))
        """)

        norm_mat = np.array(ro.r("norm_mat"))  # genes x cells
        norm_genes = list(ro.r("rownames(norm_mat)"))

        adata_scran = ad.AnnData(
            X=norm_mat.T,
            obs=adata_filtered_obs.reset_index(drop=True),
        )
        adata_scran.var_names = pd.Index(norm_genes)

        ro.r("rm(counts_mat, sce, qclust, norm_mat); gc()")

        return adata_scran

    except Exception as e:
        logger.warning("scran failed: %s", e)
        return None


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def compute_marker_auc(
    adata_norm: ad.AnnData,
    atlas: pd.DataFrame,
    organism: str,
) -> float | None:
    """Mean one-vs-rest ROC AUC for atlas markers present in this dataset."""
    if "cell_type" not in adata_norm.obs.columns:
        return None

    gene_set = set(adata_norm.var_names)
    sub = atlas[(atlas.organism == organism) & (atlas.gene.isin(gene_set))]

    if sub.empty:
        return None

    aucs = []
    for _, row in sub.iterrows():
        gene = row["gene"]
        ct = row["cell_type"]

        if ct not in adata_norm.obs.cell_type.values:
            continue

        labels = (adata_norm.obs.cell_type == ct).astype(int).values
        if labels.sum() < MIN_CELLS_PER_TYPE or (1 - labels).sum() < MIN_CELLS_PER_TYPE:
            continue

        gene_idx = list(adata_norm.var_names).index(gene)
        X = adata_norm.X
        if issparse(X):
            expr = np.asarray(X[:, gene_idx].todense()).ravel()
        else:
            expr = X[:, gene_idx].ravel()

        try:
            auc = roc_auc_score(labels, expr)
            aucs.append(auc)
        except ValueError:
            continue

    return float(np.mean(aucs)) if aucs else None


def compute_silhouette(adata_norm: ad.AnnData) -> float | None:
    """Silhouette score of ground-truth cell_type labels on PCA of normalized data."""
    if "cell_type" not in adata_norm.obs.columns:
        return None

    labels = adata_norm.obs.cell_type.values
    unique_labels = np.unique(labels)
    if len(unique_labels) < 2:
        return None

    try:
        adata_tmp = adata_norm.copy()
        # HVG selection
        n_genes = adata_tmp.n_vars
        if n_genes > N_HVG:
            try:
                sc.pp.highly_variable_genes(adata_tmp, n_top_genes=N_HVG, flavor="seurat_v3")
            except Exception:
                sc.pp.highly_variable_genes(adata_tmp, n_top_genes=N_HVG)
            adata_tmp = adata_tmp[:, adata_tmp.var.highly_variable].copy()

        sc.pp.scale(adata_tmp, max_value=10)
        sc.tl.pca(adata_tmp, n_comps=min(N_PCS, adata_tmp.n_vars - 1))

        pca_coords = adata_tmp.obsm["X_pca"]

        # Subsample if too many cells for silhouette (>10K is slow)
        n = pca_coords.shape[0]
        if n > 10_000:
            rng = np.random.default_rng(RNG_SEED)
            idx = rng.choice(n, 10_000, replace=False)
            pca_coords = pca_coords[idx]
            labels = labels[idx]

        score = silhouette_score(pca_coords, labels, sample_size=min(n, 10_000), random_state=RNG_SEED)
        return float(score)
    except Exception as e:
        logger.warning("Silhouette computation failed: %s", e)
        return None


def compute_hvg_set(adata_norm: ad.AnnData, n_top: int = N_HVG) -> set[str]:
    """Top HVGs by seurat_v3 flavor."""
    try:
        adata_tmp = adata_norm.copy()
        n_genes = adata_tmp.n_vars
        if n_genes <= n_top:
            return set(adata_tmp.var_names)
        sc.pp.highly_variable_genes(adata_tmp, n_top_genes=n_top, flavor="seurat_v3")
        return set(adata_tmp.var_names[adata_tmp.var.highly_variable])
    except Exception:
        return set()


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 1.0
    union = a | b
    return len(a & b) / len(union) if union else 0.0


# ---------------------------------------------------------------------------
# Per-dataset processing
# ---------------------------------------------------------------------------

def process_dataset(
    census,
    dataset_id: str,
    tissue: str,
    organism: str,
    atlas: pd.DataFrame,
) -> list[dict] | None:
    """Fetch raw counts, apply 3 normalizations, compute metrics."""
    exp = census["census_data"][_census_key(organism)]

    obs_filter = (
        f"dataset_id == '{dataset_id}' and tissue_general == '{tissue}' "
        f"and is_primary_data == True"
    )

    # Fetch obs with cell_type for metrics
    obs_df = exp.obs.read(
        column_names=["soma_joinid", "cell_type"],
        value_filter=obs_filter,
    ).concat().to_pandas()

    # Map Census cell types to atlas categories via substring matching
    obs_df["cell_type"] = obs_df.cell_type.apply(_map_census_cell_type)
    obs_df = obs_df.dropna(subset=["cell_type"])
    target_types = atlas_cell_types(atlas, organism)
    obs_df = obs_df[obs_df.cell_type.isin(target_types)]

    if len(obs_df) < MIN_CELLS_PER_DATASET:
        logger.info("  %s/%s: too few atlas-matching cells (%d)", dataset_id[:8], tissue, len(obs_df))
        return None

    # Subsample if too large
    if len(obs_df) > MAX_CELLS_PER_DATASET:
        rng = np.random.default_rng(RNG_SEED)
        idx = rng.choice(len(obs_df), MAX_CELLS_PER_DATASET, replace=False)
        obs_df = obs_df.iloc[sorted(idx)].copy()

    obs_ids = obs_df.soma_joinid.tolist()

    # Fetch expression
    try:
        with exp.axis_query(
            measurement_name="RNA",
            obs_query=tiledbsoma.AxisQuery(coords=(obs_ids,)),
        ) as query:
            adata = query.to_anndata(X_name="raw")
    except Exception as e:
        logger.warning("  Fetch failed for %s/%s: %s", dataset_id[:8], tissue, e)
        return None

    if adata.shape[0] == 0:
        return None

    # Attach cell_type from obs_df (matched by soma_joinid)
    ct_map = obs_df.set_index("soma_joinid")["cell_type"]
    adata.obs["cell_type"] = adata.obs.index.map(
        lambda x: ct_map.get(int(x), None) if str(x).isdigit() else None
    )
    # If mapping by index didn't work, try via soma_joinid column
    if adata.obs.cell_type.isna().all() and "soma_joinid" in adata.obs.columns:
        adata.obs["cell_type"] = adata.obs["soma_joinid"].map(ct_map)

    # Set var_names to gene symbols for marker matching
    if "feature_name" in adata.var.columns:
        adata.var_names = adata.var.feature_name.astype(str)
        adata.var_names_make_unique()

    # Basic gene filter (min 3 cells)
    sc.pp.filter_genes(adata, min_cells=3)

    logger.info(
        "  %s/%s: %d cells, %d genes, %d cell types",
        dataset_id[:8], tissue, adata.n_obs, adata.n_vars,
        adata.obs.cell_type.nunique(),
    )

    # --- Apply normalizations ---
    methods = {}

    # 1. Log-normalize
    t0 = time.time()
    adata_ln = _lognorm(adata)
    rt_ln = time.time() - t0
    methods["log_normalize"] = (adata_ln, rt_ln)

    # 2. SCTransform
    t0 = time.time()
    adata_sct = _sctransform(adata)
    rt_sct = time.time() - t0
    if adata_sct is not None:
        # Attach cell_type metadata
        adata_sct.obs["cell_type"] = adata.obs["cell_type"].values
        methods["sctransform"] = (adata_sct, rt_sct)
    else:
        logger.info("  SCTransform failed for %s/%s, skipping", dataset_id[:8], tissue)

    # 3. scran
    t0 = time.time()
    adata_scr = _scran(adata)
    rt_scr = time.time() - t0
    if adata_scr is not None:
        adata_scr.obs["cell_type"] = adata.obs["cell_type"].values
        methods["scran"] = (adata_scr, rt_scr)
    else:
        logger.info("  scran failed for %s/%s, skipping", dataset_id[:8], tissue)

    if not methods:
        return None

    # HVG baseline from log-normalize
    hvg_lognorm = compute_hvg_set(adata_ln)

    results = []
    for method_name, (adata_norm, runtime) in methods.items():
        auc = compute_marker_auc(adata_norm, atlas, organism)
        sil = compute_silhouette(adata_norm)
        hvg_this = compute_hvg_set(adata_norm)
        hvg_j = jaccard(hvg_this, hvg_lognorm) if method_name != "log_normalize" else 1.0

        results.append({
            "dataset_id": dataset_id,
            "tissue": tissue,
            "organism": organism,
            "method": method_name,
            "mean_marker_auc": round(auc, 4) if auc is not None else None,
            "silhouette_score": round(sil, 4) if sil is not None else None,
            "hvg_jaccard_vs_lognorm": round(hvg_j, 4),
            "runtime_seconds": round(runtime, 2),
        })

        logger.info(
            "    %s: AUC=%.4f  sil=%.4f  HVG_J=%.3f  %.1fs",
            method_name,
            auc if auc is not None else 0.0,
            sil if sil is not None else 0.0,
            hvg_j,
            runtime,
        )

    # Free memory
    del adata, adata_ln
    if adata_sct is not None:
        del adata_sct
    if adata_scr is not None:
        del adata_scr
    gc.collect()

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    t0 = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PROJECT_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("=== P7: Normalization Method Comparison ===")

    atlas = load_atlas()

    # Checkpoint
    checkpoint_path = PROJECT_DIR / "p7_checkpoint.tsv"
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

        datasets = discover_datasets(census, organism, atlas)
        if datasets.empty:
            continue
        logger.info("Will process %d datasets", len(datasets))

        for _, row in tqdm(
            datasets.iterrows(), total=len(datasets),
            desc=f"P7 ({organism})",
        ):
            if pd.isna(row.dataset_id):
                continue

            key = (row.dataset_id, row.tissue, organism)
            if key in done_keys:
                logger.info("  Skipping %s/%s (done)", row.dataset_id[:8], row.tissue)
                continue

            try:
                results = process_dataset(census, row.dataset_id, row.tissue, organism, atlas)
            except Exception as e:
                logger.warning("  Error on %s/%s: %s", row.dataset_id[:8], row.tissue, e)
                continue

            if results:
                all_results.extend(results)
                done_keys.add(key)
                # Checkpoint
                pd.DataFrame(all_results).to_csv(checkpoint_path, sep="\t", index=False)
                gc.collect()

    census.close()

    result_df = pd.DataFrame(all_results)
    out_path = OUTPUT_DIR / "p7_normalization.tsv"
    result_df.to_csv(out_path, sep="\t", index=False)
    logger.info("\nSaved %d rows to %s", len(result_df), out_path)

    # --- Summary statistics ---
    elapsed = time.time() - t0
    n_datasets = result_df.groupby(["dataset_id", "tissue", "organism"], observed=True).ngroups
    logger.info("\n=== P7 COMPLETE ===")
    logger.info("Dataset-tissue combos: %d", n_datasets)
    logger.info("Total rows (datasets x methods): %d", len(result_df))

    for method in ["log_normalize", "sctransform", "scran"]:
        sub = result_df[result_df.method == method]
        if len(sub) > 0:
            auc_vals = sub.mean_marker_auc.dropna()
            sil_vals = sub.silhouette_score.dropna()
            hvg_vals = sub.hvg_jaccard_vs_lognorm.dropna()
            logger.info(
                "  %s: n=%d  mean_AUC=%.4f  mean_sil=%.4f  mean_HVG_J=%.3f  mean_runtime=%.1fs",
                method,
                len(sub),
                auc_vals.mean() if len(auc_vals) > 0 else 0.0,
                sil_vals.mean() if len(sil_vals) > 0 else 0.0,
                hvg_vals.mean() if len(hvg_vals) > 0 else 0.0,
                sub.runtime_seconds.mean(),
            )

    logger.info("Elapsed: %.1f minutes", elapsed / 60)


if __name__ == "__main__":
    main()
