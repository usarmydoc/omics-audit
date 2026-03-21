#!/usr/bin/env python3
"""P4 — Pseudobulk vs Pseudoreplication Audit.

Question: How much do DEG lists differ between DESeq2 pseudobulk and
Wilcoxon single-cell testing on the same data?

Strategy:
- Find Census datasets with disease vs normal AND ≥3 donors per group
- Pull count matrices for a cell type present in both groups
- Pseudobulk: aggregate counts per donor, run pyDESeq2
- Cell-level: Wilcoxon rank-sum on normalized single-cell data
- Compare top 100 DEGs: Jaccard, direction agreement, inflation factor
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
from tqdm import tqdm

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="pydeseq2")

# ── Config ──────────────────────────────────────────────────────────
CENSUS_VERSION = "2025-11-08"
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = REPO_ROOT / "scrnaseq" / "output"
PROJECT_DIR = Path(__file__).resolve().parent

MAX_DATASETS = 20
MAX_CELLS_PER_DATASET = 20_000
MIN_CELLS_PER_GROUP = 50
MIN_DONORS_PER_GROUP = 3
MIN_GENES = 500
TOP_N_DEGS = 100
RNG_SEED = 42

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(PROJECT_DIR / "p4_run.log"),
    ],
)
logger = logging.getLogger(__name__)


def _census_key(organism: str) -> str:
    return organism.lower().replace(" ", "_")


def discover_eligible_datasets(census, organism: str) -> pd.DataFrame:
    """Find datasets with disease vs normal, ≥3 donors per group."""
    exp = census["census_data"][_census_key(organism)]
    logger.info("Discovering eligible datasets for %s...", organism)

    # Get disease cells
    obs_disease = exp.obs.read(
        column_names=["dataset_id", "tissue_general", "donor_id", "disease", "cell_type"],
        value_filter="is_primary_data == True and disease != 'normal'",
    ).concat().to_pandas()

    # Get normal cells
    obs_normal = exp.obs.read(
        column_names=["dataset_id", "tissue_general", "donor_id", "disease", "cell_type"],
        value_filter="is_primary_data == True and disease == 'normal'",
    ).concat().to_pandas()

    # Find overlapping datasets
    disease_ds = set(obs_disease.dataset_id.unique())
    normal_ds = set(obs_normal.dataset_id.unique())
    overlap = disease_ds & normal_ds

    # Combine
    obs_d = obs_disease[obs_disease.dataset_id.isin(overlap)].copy()
    obs_n = obs_normal[obs_normal.dataset_id.isin(overlap)].copy()
    obs_d["group"] = "disease"
    obs_n["group"] = "normal"
    all_obs = pd.concat([obs_d, obs_n])

    # For each dataset, check donor counts per group
    results = []
    for ds_id, ds_df in all_obs.groupby("dataset_id", observed=True):
        donors_disease = ds_df[ds_df.group == "disease"].donor_id.nunique()
        donors_normal = ds_df[ds_df.group == "normal"].donor_id.nunique()
        if donors_disease < MIN_DONORS_PER_GROUP or donors_normal < MIN_DONORS_PER_GROUP:
            continue

        cells_disease = len(ds_df[ds_df.group == "disease"])
        cells_normal = len(ds_df[ds_df.group == "normal"])
        if cells_disease < MIN_CELLS_PER_GROUP or cells_normal < MIN_CELLS_PER_GROUP:
            continue

        tissue = ds_df.tissue_general.mode().iloc[0] if len(ds_df.tissue_general.mode()) > 0 else "unknown"
        diseases = ds_df[ds_df.group == "disease"].disease.unique()

        # Find a shared cell type with enough cells in both groups
        ct_disease = ds_df[ds_df.group == "disease"].cell_type.value_counts()
        ct_normal = ds_df[ds_df.group == "normal"].cell_type.value_counts()
        shared_cts = set(ct_disease.index) & set(ct_normal.index)
        best_ct = None
        best_n = 0
        for ct in shared_cts:
            n = min(ct_disease.get(ct, 0), ct_normal.get(ct, 0))
            if n > best_n:
                best_n = n
                best_ct = ct

        if best_ct is None or best_n < MIN_CELLS_PER_GROUP:
            continue

        results.append({
            "dataset_id": ds_id,
            "tissue": tissue,
            "disease": "; ".join(sorted(set(diseases))[:3]),
            "n_donors_disease": donors_disease,
            "n_donors_normal": donors_normal,
            "n_cells_disease": cells_disease,
            "n_cells_normal": cells_normal,
            "best_cell_type": best_ct,
            "best_ct_min_n": best_n,
            "total_donors": donors_disease + donors_normal,
        })

    result_df = pd.DataFrame(results)
    if len(result_df) == 0:
        logger.warning("No eligible datasets found for %s", organism)
        return result_df

    # Sort by total donors and take top MAX_DATASETS
    result_df = result_df.sort_values("total_donors", ascending=False).head(MAX_DATASETS)
    logger.info(
        "%s: found %d eligible datasets across %d tissues",
        organism, len(result_df), result_df.tissue.nunique(),
    )

    del obs_disease, obs_normal, all_obs
    gc.collect()
    return result_df.reset_index(drop=True)


def fetch_expression_for_de(
    census, dataset_id: str, organism: str, cell_type: str,
) -> tuple[np.ndarray, pd.DataFrame, list[str]] | None:
    """Fetch raw counts for a specific cell type in disease vs normal.

    Returns (count_matrix, obs_df, gene_names) or None.
    """
    exp = census["census_data"][_census_key(organism)]

    # Get obs for this dataset, cell type, disease/normal
    obs_filter = (
        f"dataset_id == '{dataset_id}' and is_primary_data == True "
        f"and cell_type == '{cell_type}'"
    )
    obs_df = exp.obs.read(
        column_names=["soma_joinid", "donor_id", "disease"],
        value_filter=obs_filter,
    ).concat().to_pandas()

    obs_df["group"] = obs_df.disease.apply(lambda x: "normal" if x == "normal" else "disease")

    # Check we have enough
    n_disease = (obs_df.group == "disease").sum()
    n_normal = (obs_df.group == "normal").sum()
    if n_disease < MIN_CELLS_PER_GROUP or n_normal < MIN_CELLS_PER_GROUP:
        return None

    # Check donor counts per group
    donors_d = obs_df[obs_df.group == "disease"].donor_id.nunique()
    donors_n = obs_df[obs_df.group == "normal"].donor_id.nunique()
    if donors_d < MIN_DONORS_PER_GROUP or donors_n < MIN_DONORS_PER_GROUP:
        return None

    # Subsample if too large — stratified by group
    if len(obs_df) > MAX_CELLS_PER_DATASET:
        rng = np.random.default_rng(RNG_SEED)
        per_group = MAX_CELLS_PER_DATASET // 2
        idx_d = obs_df[obs_df.group == "disease"].index
        idx_n = obs_df[obs_df.group == "normal"].index
        if len(idx_d) > per_group:
            idx_d = rng.choice(idx_d, per_group, replace=False)
        if len(idx_n) > per_group:
            idx_n = rng.choice(idx_n, per_group, replace=False)
        keep = np.concatenate([idx_d, idx_n])
        obs_df = obs_df.loc[sorted(keep)].copy()

    obs_ids = obs_df.soma_joinid.tolist()

    # Fetch ALL genes (need broad gene set for DE)
    try:
        with exp.axis_query(
            measurement_name="RNA",
            obs_query=tiledbsoma.AxisQuery(coords=(obs_ids,)),
        ) as query:
            adata = query.to_anndata(X_name="raw")
    except Exception as e:
        logger.warning("Failed to fetch expression for %s: %s", dataset_id[:8], e)
        return None

    if adata.shape[0] == 0 or adata.shape[1] < MIN_GENES:
        return None

    # Get gene names
    gene_names = list(adata.var.feature_name) if "feature_name" in adata.var.columns else list(adata.var_names)

    # Filter to top 10K genes by total expression to make DE tractable
    # This drops unexpressed/rare genes that would never appear in top DEGs
    X_raw = adata.X
    if hasattr(X_raw, "toarray"):
        gene_sums = np.asarray(X_raw.sum(axis=0)).ravel()
    else:
        gene_sums = np.asarray(X_raw).sum(axis=0)
    MAX_GENES_DE = 10_000
    if len(gene_sums) > MAX_GENES_DE:
        top_idx = np.argsort(gene_sums)[-MAX_GENES_DE:]
        top_idx = np.sort(top_idx)
        adata = adata[:, top_idx].copy()
        gene_names = [gene_names[i] for i in top_idx]
        logger.info("  Filtered to top %d genes by expression", MAX_GENES_DE)

    # Keep sparse for pseudobulk
    X = adata.X

    # Align obs_df with adata
    adata_joinids = adata.obs["soma_joinid"].values
    obs_indexed = obs_df.set_index("soma_joinid")
    aligned_obs = obs_indexed.loc[adata_joinids].reset_index()

    return X, aligned_obs, gene_names


def run_pseudobulk_deseq2(
    X, obs_df: pd.DataFrame, gene_names: list[str],
) -> pd.DataFrame | None:
    """Aggregate to pseudobulk per donor, run pyDESeq2."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    # Aggregate: sum counts per donor
    donors = obs_df.donor_id.values
    groups = obs_df.group.values
    unique_donors = sorted(set(donors))

    if len(unique_donors) < 6:
        return None

    # Build pseudobulk count matrix (donors x genes)
    donor_counts = []
    donor_meta = []
    for donor in unique_donors:
        mask = donors == donor
        group = groups[mask][0]  # all cells from same donor should have same group
        if hasattr(X, "toarray"):
            counts = np.asarray(X[mask].sum(axis=0)).ravel()
        else:
            counts = X[mask].sum(axis=0)
        donor_counts.append(counts)
        donor_meta.append({"donor_id": donor, "group": group})

    count_matrix = np.array(donor_counts, dtype=np.float64)
    meta_df = pd.DataFrame(donor_meta)

    # Filter low-expression genes (≥10 counts in ≥3 donors)
    gene_mask = (count_matrix >= 10).sum(axis=0) >= 3
    if gene_mask.sum() < 100:
        logger.warning("Too few genes pass filter for pseudobulk DE")
        return None

    count_filtered = count_matrix[:, gene_mask]
    genes_filtered = [g for g, m in zip(gene_names, gene_mask) if m]

    count_df = pd.DataFrame(count_filtered, columns=genes_filtered, index=meta_df.donor_id)
    count_df = count_df.astype(int)

    try:
        dds = DeseqDataSet(
            counts=count_df,
            metadata=meta_df.set_index("donor_id"),
            design="~group",
        )
        dds.deseq2()
        stat_res = DeseqStats(dds, contrast=["group", "disease", "normal"])
        stat_res.summary()
        results = stat_res.results_df.copy()
        results["gene"] = results.index
        results = results.rename(columns={"log2FoldChange": "log2fc", "padj": "padj", "pvalue": "pval"})
        return results[["gene", "log2fc", "pval", "padj"]].dropna(subset=["padj"])
    except Exception as e:
        logger.warning("pyDESeq2 failed: %s", e)
        return None


def run_wilcoxon_cellwise(
    X, obs_df: pd.DataFrame, gene_names: list[str],
) -> pd.DataFrame | None:
    """Wilcoxon rank-sum test on single-cell expression (cell-level DE)."""
    groups = obs_df.group.values
    mask_d = groups == "disease"
    mask_n = groups == "normal"

    if hasattr(X, "toarray"):
        X_dense = X.toarray().astype(np.float32)
    else:
        X_dense = np.asarray(X, dtype=np.float32)

    # Filter genes: expressed in ≥10 cells
    expressed = (X_dense > 0).sum(axis=0) >= 10
    if expressed.sum() < 100:
        return None

    X_filt = X_dense[:, expressed]
    genes_filt = [g for g, m in zip(gene_names, expressed) if m]

    # Simple log-normalize for fold change calculation
    lib_sizes = X_filt.sum(axis=1, keepdims=True)
    lib_sizes = np.where(lib_sizes == 0, 1, lib_sizes)
    X_norm = np.log1p(X_filt / lib_sizes * 10000)

    results = []
    X_d = X_norm[mask_d]
    X_n = X_norm[mask_n]

    for i, gene in enumerate(genes_filt):
        d_vals = X_d[:, i]
        n_vals = X_n[:, i]
        try:
            stat, pval = scipy_stats.mannwhitneyu(d_vals, n_vals, alternative="two-sided")
        except ValueError:
            continue
        log2fc = float(np.mean(d_vals) - np.mean(n_vals))  # log-space difference ≈ log2FC
        results.append({"gene": gene, "log2fc": log2fc, "pval": pval})

    if not results:
        return None

    result_df = pd.DataFrame(results)
    # BH correction
    from scipy.stats import false_discovery_control
    try:
        result_df["padj"] = false_discovery_control(result_df.pval.values, method="bh")
    except Exception:
        # Fallback manual BH
        n = len(result_df)
        ranked = result_df.pval.rank()
        result_df["padj"] = np.minimum(result_df.pval * n / ranked, 1.0)

    return result_df[["gene", "log2fc", "pval", "padj"]].dropna(subset=["padj"])


def compare_deg_lists(
    pb_results: pd.DataFrame, wc_results: pd.DataFrame, top_n: int = TOP_N_DEGS,
) -> dict:
    """Compare top DEGs from pseudobulk vs Wilcoxon."""
    # Top N by padj (ascending), then by abs(log2fc) descending
    pb_top = pb_results.nsmallest(top_n, "padj").head(top_n)
    wc_top = wc_results.nsmallest(top_n, "padj").head(top_n)

    pb_genes = set(pb_top.gene)
    wc_genes = set(wc_top.gene)

    # Jaccard
    intersection = pb_genes & wc_genes
    union = pb_genes | wc_genes
    jaccard = len(intersection) / len(union) if union else 0

    # Direction agreement (for shared genes)
    if intersection:
        pb_dir = pb_top.set_index("gene").loc[list(intersection), "log2fc"]
        wc_dir = wc_top.set_index("gene").loc[list(intersection), "log2fc"]
        common = pb_dir.index.intersection(wc_dir.index)
        if len(common) > 0:
            agree = sum(
                np.sign(pb_dir.loc[g]) == np.sign(wc_dir.loc[g])
                for g in common
            )
            direction_agreement = agree / len(common)
        else:
            direction_agreement = np.nan
    else:
        direction_agreement = np.nan

    # Inflation factor: ratio of significant genes (padj < 0.05)
    pb_sig = (pb_results.padj < 0.05).sum()
    wc_sig = (wc_results.padj < 0.05).sum()
    inflation = wc_sig / pb_sig if pb_sig > 0 else np.nan

    return {
        "jaccard_top100": float(jaccard),
        "direction_agreement": float(direction_agreement) if not np.isnan(direction_agreement) else None,
        "wilcoxon_sig_count": int(wc_sig),
        "pseudobulk_sig_count": int(pb_sig),
        "wilcoxon_inflation_factor": float(inflation) if not np.isnan(inflation) else None,
        "shared_in_top100": len(intersection),
    }


def main():
    t0 = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PROJECT_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("=== P4: Pseudobulk vs Pseudoreplication Audit ===")

    # Resume from checkpoint
    checkpoint_path = PROJECT_DIR / "p4_checkpoint.tsv"
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
            desc=f"P4 ({organism})",
        ):
            ds_id = row.dataset_id
            if pd.isna(ds_id):
                continue
            key = (ds_id, organism)
            if key in done_keys:
                logger.info("  Skipping %s (already done)", ds_id[:8])
                continue

            logger.info(
                "  Processing %s/%s (%s), cell_type=%s, donors=%d+%d",
                ds_id[:8], row.tissue, row.disease[:30],
                row.best_cell_type, row.n_donors_disease, row.n_donors_normal,
            )

            try:
                result = fetch_expression_for_de(
                    census, ds_id, organism, row.best_cell_type,
                )
            except Exception as e:
                logger.warning("  Fetch failed for %s: %s", ds_id[:8], e)
                continue

            if result is None:
                logger.info("  Skipped %s (insufficient data)", ds_id[:8])
                continue

            X, obs_df, gene_names = result
            n_donors = obs_df.donor_id.nunique()
            n_cells = len(obs_df)

            logger.info("  Fetched: %d cells, %d genes, %d donors", n_cells, len(gene_names), n_donors)

            # Run pseudobulk DESeq2
            pb_results = run_pseudobulk_deseq2(X, obs_df, gene_names)
            if pb_results is None or len(pb_results) < TOP_N_DEGS:
                logger.info("  Pseudobulk DE failed or too few results")
                del X, obs_df, gene_names
                gc.collect()
                continue

            # Run Wilcoxon cell-level
            wc_results = run_wilcoxon_cellwise(X, obs_df, gene_names)
            if wc_results is None or len(wc_results) < TOP_N_DEGS:
                logger.info("  Wilcoxon DE failed or too few results")
                del X, obs_df, gene_names, pb_results
                gc.collect()
                continue

            # Compare
            comparison = compare_deg_lists(pb_results, wc_results)

            record = {
                "dataset_id": ds_id,
                "tissue": row.tissue,
                "organism": organism,
                "disease": row.disease,
                "cell_type": row.best_cell_type,
                "n_donors": n_donors,
                "n_cells": n_cells,
                **comparison,
            }
            all_results.append(record)
            done_keys.add(key)

            logger.info(
                "  Result: Jaccard=%.3f, direction_agree=%s, inflation=%.1fx, shared=%d/100",
                comparison["jaccard_top100"],
                f"{comparison['direction_agreement']:.2f}" if comparison["direction_agreement"] else "N/A",
                comparison["wilcoxon_inflation_factor"] if comparison["wilcoxon_inflation_factor"] else float("nan"),
                comparison["shared_in_top100"],
            )

            # Checkpoint
            pd.DataFrame(all_results).to_csv(checkpoint_path, sep="\t", index=False)

            del X, obs_df, gene_names, pb_results, wc_results
            gc.collect()

    census.close()

    # Save final output
    result_df = pd.DataFrame(all_results)
    out_path = OUTPUT_DIR / "p4_pseudobulk.tsv"
    result_df.to_csv(out_path, sep="\t", index=False)
    logger.info("\nSaved %d rows to %s", len(result_df), out_path)

    # Summary
    elapsed = time.time() - t0
    logger.info("\n=== P4 COMPLETE ===")
    logger.info("Datasets processed: %d", len(result_df))
    if len(result_df) > 0:
        logger.info("Mean Jaccard (top 100): %.3f", result_df.jaccard_top100.mean())
        dir_agree = result_df.direction_agreement.dropna()
        if len(dir_agree) > 0:
            logger.info("Mean direction agreement: %.3f", dir_agree.mean())
        inflate = result_df.wilcoxon_inflation_factor.dropna()
        if len(inflate) > 0:
            logger.info("Mean Wilcoxon inflation: %.1fx", inflate.mean())
    logger.info("Elapsed: %.1f minutes", elapsed / 60)


if __name__ == "__main__":
    main()
