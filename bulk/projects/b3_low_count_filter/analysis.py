"""
B3 analysis: low-count filter threshold sweep.

For each cancer type, applies CPM thresholds × min-sample thresholds,
records genes retained, tissue-specific gene retention, and downstream
edgeR DEG count.
"""

import numpy as np
import pandas as pd

# Use rpy2 for edgeR
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri, pandas2ri, default_converter
from rpy2.robjects.packages import importr

# Build a converter that handles numpy + pandas ↔ R
_converter = default_converter + numpy2ri.converter + pandas2ri.converter

edgeR = importr("edgeR")
limma = importr("limma")
stats = importr("stats")
base = importr("base")

CPM_THRESHOLDS = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
MIN_SAMPLE_MODES = [
    ("abs_1", 1),
    ("abs_2", 2),
    ("abs_3", 3),
    ("pct_10", 0.10),
    ("pct_25", 0.25),
]


def _compute_min_samples(mode_value, n_samples: int) -> int:
    """Convert min_sample mode to absolute number."""
    if isinstance(mode_value, float) and mode_value < 1:
        return max(1, int(np.ceil(mode_value * n_samples)))
    return int(mode_value)


def _filter_by_cpm(counts: pd.DataFrame, cpm_threshold: float,
                   min_samples: int) -> pd.Index:
    """Return gene names passing the CPM filter."""
    lib_sizes = counts.sum(axis=0)
    cpm = counts.div(lib_sizes, axis=1) * 1e6
    keep = (cpm >= cpm_threshold).sum(axis=1) >= min_samples
    return counts.index[keep]


def _run_edger_deg(counts: pd.DataFrame, groups: pd.Series) -> int:
    """Run edgeR on filtered counts, return number of DEGs (FDR < 0.05).

    counts: genes × samples DataFrame (integer counts)
    groups: sample_id → group label Series
    """
    # Align samples
    common = counts.columns.intersection(groups.index)
    if len(common) < 4:
        return 0
    counts_aligned = counts[common]
    groups_aligned = groups[common]

    # Build R matrix outside converter context so it stays as R object
    np_counts = counts_aligned.values.astype(np.float64)
    nr, nc = np_counts.shape
    r_vec = ro.FloatVector(np_counts.T.flatten())  # R is column-major
    r_matrix = ro.r.matrix(r_vec, nrow=nr, ncol=nc)
    r_matrix.rownames = ro.StrVector(list(counts_aligned.index))
    r_matrix.colnames = ro.StrVector(list(counts_aligned.columns))
    r_groups = ro.StrVector(groups_aligned.values.tolist())

    # edgeR pipeline
    dge = edgeR.DGEList(counts=r_matrix, group=r_groups)
    dge = edgeR.calcNormFactors(dge)

    design = stats.model_matrix(ro.Formula("~group"),
                                data=ro.DataFrame({"group": r_groups}))
    dge = edgeR.estimateDisp(dge, design)
    fit = edgeR.glmQLFit(dge, design)
    qlf = edgeR.glmQLFTest(fit, coef=2)

    tt = edgeR.topTags(qlf, n=ro.IntVector([counts_aligned.shape[0]]))
    results = ro.r["as.data.frame"](tt)
    fdr_r = results.rx2("FDR")
    fdr = np.array(list(fdr_r))
    return int(np.sum(fdr < 0.05))


def run_filter_sweep(counts: pd.DataFrame, sample_type: pd.Series,
                     tissue_genes: list[str], project_id: str,
                     cancer_type: str) -> list[dict]:
    """Run the full CPM × min-samples sweep for one cancer type.

    Returns list of result dicts, one per threshold combination.
    """
    n_total = counts.shape[1]
    total_genes = counts.shape[0]

    # Identify which tissue-specific genes are in the count matrix
    tissue_genes_present = [g for g in tissue_genes if g in counts.index]
    n_tissue_genes = len(tissue_genes_present)

    results = []
    for cpm_thresh in CPM_THRESHOLDS:
        for mode_name, mode_value in MIN_SAMPLE_MODES:
            min_samp = _compute_min_samples(mode_value, n_total)
            passing_genes = _filter_by_cpm(counts, cpm_thresh, min_samp)
            n_retained = len(passing_genes)

            # Tissue-specific gene retention
            tissue_retained = [g for g in tissue_genes_present if g in passing_genes]
            n_tissue_retained = len(tissue_retained)
            pct_lost = 0.0
            if n_tissue_genes > 0:
                pct_lost = round(100.0 * (1 - n_tissue_retained / n_tissue_genes), 1)

            # Run downstream edgeR on filtered counts
            filtered_counts = counts.loc[passing_genes]
            try:
                n_degs = _run_edger_deg(filtered_counts, sample_type)
            except Exception as e:
                print(f"    edgeR failed for {project_id} cpm={cpm_thresh} "
                      f"min={mode_name}: {e}")
                n_degs = -1

            results.append({
                "dataset_id": project_id,
                "cancer_type": cancer_type,
                "cpm_threshold": cpm_thresh,
                "min_samples": mode_name,
                "min_samples_abs": min_samp,
                "genes_retained": n_retained,
                "total_genes": total_genes,
                "tissue_specific_genes_present": n_tissue_genes,
                "tissue_specific_genes_retained": n_tissue_retained,
                "pct_tissue_specific_lost": pct_lost,
                "downstream_deg_count": n_degs,
            })

    return results
