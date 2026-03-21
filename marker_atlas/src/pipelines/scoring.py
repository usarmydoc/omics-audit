"""Marker reliability scoring — sensitivity, specificity, AUC, Mann-Whitney, bootstrap CI.

All scoring runs on pre-fetched numpy arrays, so it's safe for multiprocessing.Pool.
"""

from __future__ import annotations

import logging
import multiprocessing as mp
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve

logger = logging.getLogger(__name__)

N_BOOTSTRAP = 200
RNG_SEED = 42
MIN_POSITIVE = 5   # minimum cells of target type to score
MIN_NEGATIVE = 5   # minimum cells of other types to score


@dataclass
class MarkerScore:
    """Per-gene, per-cell-type, per-dataset scoring result."""
    gene: str
    cell_type: str
    tissue: str
    organism: str
    dataset_id: str
    n_positive: int
    n_negative: int
    sensitivity: float
    specificity: float
    auc: float
    mw_statistic: float
    mw_pval: float
    ci_lower: float
    ci_upper: float


@dataclass
class AggregatedScore:
    """Cross-dataset aggregated score for one gene-cell_type pair."""
    gene: str
    cell_type: str
    tissue: str
    organism: str
    n_datasets: int
    mean_auc: float
    sd_auc: float
    ci_lower: float
    ci_upper: float
    mw_pval: float
    reliability_score: float


def _bootstrap_auc(
    y_true: np.ndarray,
    y_score: np.ndarray,
    n_boot: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    """Bootstrap 95% CI on AUC."""
    aucs = []
    n = len(y_true)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yt = y_true[idx]
        ys = y_score[idx]
        # Need both classes in bootstrap sample
        if yt.sum() == 0 or yt.sum() == n:
            continue
        try:
            aucs.append(roc_auc_score(yt, ys))
        except ValueError:
            continue

    if len(aucs) < 10:
        return (0.0, 1.0)  # uninformative CI

    return (float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5)))


def _optimal_threshold_metrics(
    y_true: np.ndarray, y_score: np.ndarray
) -> tuple[float, float]:
    """Find sensitivity and specificity at Youden's J optimal threshold."""
    fpr, tpr, _ = roc_curve(y_true, y_score)
    # Youden's J = sensitivity + specificity - 1 = tpr + (1 - fpr) - 1 = tpr - fpr
    j_scores = tpr - fpr
    best_idx = np.argmax(j_scores)
    sensitivity = float(tpr[best_idx])
    specificity = float(1.0 - fpr[best_idx])
    return sensitivity, specificity


def compute_marker_scores(
    expr: np.ndarray,
    gene_names: list[str],
    cell_types: np.ndarray,
    target_cell_type: str,
    dataset_id: str,
    tissue: str,
    organism: str = "Homo sapiens",
    n_bootstrap: int = N_BOOTSTRAP,
    rng_seed: int = RNG_SEED,
) -> list[MarkerScore]:
    """Compute all metrics for each gene against target_cell_type in one dataset.

    Binary label: 1 if cell_type == target_cell_type, else 0.
    """
    y_true = (cell_types == target_cell_type).astype(np.int32)
    n_pos = int(y_true.sum())
    n_neg = int(len(y_true) - n_pos)

    if n_pos < MIN_POSITIVE or n_neg < MIN_NEGATIVE:
        return []

    rng = np.random.default_rng(rng_seed)
    scores = []

    for gene_idx, gene in enumerate(gene_names):
        y_score = expr[:, gene_idx].astype(np.float64)

        # Skip if zero variance
        if y_score.std() == 0:
            continue

        # AUC
        try:
            auc = float(roc_auc_score(y_true, y_score))
        except ValueError:
            continue

        # Sensitivity / specificity at optimal threshold
        sensitivity, specificity = _optimal_threshold_metrics(y_true, y_score)

        # Mann-Whitney U
        pos_expr = y_score[y_true == 1]
        neg_expr = y_score[y_true == 0]
        try:
            mw_stat, mw_pval = stats.mannwhitneyu(pos_expr, neg_expr, alternative="two-sided")
        except ValueError:
            mw_stat, mw_pval = 0.0, 1.0

        # Bootstrap CI on AUC
        ci_lower, ci_upper = _bootstrap_auc(y_true, y_score, n_bootstrap, rng)

        scores.append(MarkerScore(
            gene=gene,
            cell_type=target_cell_type,
            tissue=tissue,
            organism=organism,
            dataset_id=dataset_id,
            n_positive=n_pos,
            n_negative=n_neg,
            sensitivity=sensitivity,
            specificity=specificity,
            auc=auc,
            mw_statistic=float(mw_stat),
            mw_pval=float(mw_pval),
            ci_lower=ci_lower,
            ci_upper=ci_upper,
        ))

    return scores


def _score_dataset_worker(args: tuple) -> list[MarkerScore]:
    """Worker function for multiprocessing.Pool.

    Scores all (gene, cell_type) combinations for one dataset.
    """
    expr, gene_names, cell_types, target_types, dataset_id, tissue, organism = args
    all_scores = []
    for target in target_types:
        scores = compute_marker_scores(
            expr, gene_names, cell_types, target,
            dataset_id, tissue, organism,
        )
        all_scores.extend(scores)
    return all_scores


def score_all_parallel(
    work_items: list[tuple],
    max_workers: int = 8,
) -> list[MarkerScore]:
    """Run scoring across datasets using multiprocessing.Pool.

    Each work_item is (expr, gene_names, cell_types, target_types, dataset_id, tissue).
    """
    if not work_items:
        return []

    n_workers = min(max_workers, len(work_items))
    logger.info("Scoring %d datasets with %d workers", len(work_items), n_workers)

    all_scores = []

    if n_workers <= 1:
        for item in tqdm(work_items, desc="Scoring datasets"):
            all_scores.extend(_score_dataset_worker(item))
    else:
        with mp.Pool(n_workers) as pool:
            results = list(tqdm(
                pool.imap_unordered(_score_dataset_worker, work_items),
                total=len(work_items),
                desc="Scoring datasets",
            ))
            for r in results:
                all_scores.extend(r)

    logger.info("Computed %d individual marker scores", len(all_scores))
    return all_scores


def aggregate_scores(scores: list[MarkerScore]) -> list[AggregatedScore]:
    """Aggregate per-dataset scores into per-(gene, cell_type) summaries.

    Groups by (gene, cell_type), computes mean/sd across datasets.
    reliability_score = mean_auc * (1 - sd_auc) * min(1.0, n_datasets / 5)
    """
    if not scores:
        return []

    # Convert to DataFrame for easy grouping
    records = [
        {
            "gene": s.gene,
            "cell_type": s.cell_type,
            "tissue": s.tissue,
            "organism": s.organism,
            "dataset_id": s.dataset_id,
            "auc": s.auc,
            "ci_lower": s.ci_lower,
            "ci_upper": s.ci_upper,
            "mw_pval": s.mw_pval,
            "sensitivity": s.sensitivity,
            "specificity": s.specificity,
        }
        for s in scores
    ]
    df = pd.DataFrame(records)

    aggregated = []

    for (gene, cell_type, organism), group in df.groupby(["gene", "cell_type", "organism"]):
        n_datasets = group.dataset_id.nunique()
        mean_auc = group.auc.mean()
        sd_auc = group.auc.std() if n_datasets > 1 else 0.0
        ci_lower = group.ci_lower.mean()
        ci_upper = group.ci_upper.mean()

        # Fisher's method for combining p-values
        pvals = group.mw_pval.values
        pvals = pvals[pvals > 0]  # avoid log(0)
        if len(pvals) > 0:
            # Fisher's method: -2 * sum(log(p)) ~ chi2(2k)
            chi2_stat = -2.0 * np.sum(np.log(pvals))
            combined_pval = float(stats.chi2.sf(chi2_stat, 2 * len(pvals)))
        else:
            combined_pval = 1.0

        reliability = mean_auc * (1.0 - sd_auc) * min(1.0, n_datasets / 5.0)

        # Get representative tissue (most common or "multiple")
        tissues = group.tissue.unique()
        tissue_str = tissues[0] if len(tissues) == 1 else "multiple"

        aggregated.append(AggregatedScore(
            gene=str(gene),
            cell_type=str(cell_type),
            tissue=tissue_str,
            organism=str(organism),
            n_datasets=n_datasets,
            mean_auc=float(mean_auc),
            sd_auc=float(sd_auc),
            ci_lower=float(ci_lower),
            ci_upper=float(ci_upper),
            mw_pval=float(combined_pval),
            reliability_score=float(reliability),
        ))

    # Sort by reliability score descending
    aggregated.sort(key=lambda x: x.reliability_score, reverse=True)

    logger.info("Aggregated into %d gene-cell_type scores", len(aggregated))
    return aggregated


# tqdm import needed for progress bars in parallel scoring
from tqdm import tqdm
