"""
B4 analysis: batch correction method comparison.

ComBat and SVA run via subprocess Rscript (rpy2 segfaults with sva).
limma removeBatchEffect runs via rpy2.
"""

import subprocess
import tempfile
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import adjusted_rand_score
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

import rpy2.robjects as ro
from rpy2.robjects.packages import importr

edgeR = importr("edgeR")
limma = importr("limma")
r_stats = importr("stats")


def _counts_to_r_matrix(counts: pd.DataFrame):
    np_counts = counts.values.astype(np.float64)
    nr, nc = np_counts.shape
    r_vec = ro.FloatVector(np_counts.T.flatten())
    r_matrix = ro.r.matrix(r_vec, nrow=nr, ncol=nc)
    r_matrix.rownames = ro.StrVector(list(counts.index))
    r_matrix.colnames = ro.StrVector(list(counts.columns))
    return r_matrix


def _filter_small_batches(batch: pd.Series, min_size: int = 2) -> pd.Index:
    counts = batch.value_counts()
    valid = counts[counts >= min_size].index
    return batch[batch.isin(valid)].index


def _normalize_logcpm(counts: pd.DataFrame) -> pd.DataFrame:
    r_matrix = _counts_to_r_matrix(counts)
    dge = edgeR.DGEList(counts=r_matrix)
    dge = edgeR.calcNormFactors(dge)
    log_cpm = ro.r["cpm"](dge, log=True)
    result = np.array(list(log_cpm))
    nr, nc = counts.shape
    result = result.reshape(nc, nr).T
    return pd.DataFrame(result, index=counts.index, columns=counts.columns)


def _filter_low_counts(counts, min_cpm=1.0, min_samples=2):
    lib_sizes = counts.sum(axis=0)
    cpm = counts.div(lib_sizes, axis=1) * 1e6
    keep = (cpm >= min_cpm).sum(axis=1) >= min_samples
    return counts[keep]


def _run_r_script(expr: pd.DataFrame, batch: pd.Series,
                  sample_type: pd.Series, method: str) -> pd.DataFrame:
    """Run ComBat or SVA via subprocess Rscript to avoid rpy2 segfault."""
    with tempfile.TemporaryDirectory() as tmpdir:
        expr_file = f"{tmpdir}/expr.tsv"
        meta_file = f"{tmpdir}/meta.tsv"
        out_file = f"{tmpdir}/corrected.tsv"

        expr.to_csv(expr_file, sep="\t")
        meta_df = pd.DataFrame({"batch": batch, "type": sample_type})
        meta_df.to_csv(meta_file, sep="\t")

        if method == "ComBat":
            r_code = f'''
library(sva)
expr <- as.matrix(read.table("{expr_file}", sep="\\t", header=TRUE, row.names=1, check.names=FALSE))
meta <- read.table("{meta_file}", sep="\\t", header=TRUE, row.names=1, check.names=FALSE)
# Align
common <- intersect(colnames(expr), rownames(meta))
expr <- expr[, common]
meta <- meta[common, ]
mod <- model.matrix(~type, data=meta)
corrected <- ComBat(dat=expr, batch=meta$batch, mod=mod)
write.table(corrected, "{out_file}", sep="\\t", quote=FALSE)
'''
        elif method == "SVA":
            r_code = f'''
library(sva)
expr <- as.matrix(read.table("{expr_file}", sep="\\t", header=TRUE, row.names=1, check.names=FALSE))
meta <- read.table("{meta_file}", sep="\\t", header=TRUE, row.names=1, check.names=FALSE)
common <- intersect(colnames(expr), rownames(meta))
expr <- expr[, common]
meta <- meta[common, ]
mod <- model.matrix(~type, data=meta)
mod0 <- model.matrix(~1, data=meta)
n.sv <- num.sv(expr, mod, method="be")
if (n.sv > 0) {{
  svobj <- sva(expr, mod, mod0, n.sv=n.sv)
  # Remove SV effects
  sv_matrix <- svobj$sv
  for (i in 1:nrow(expr)) {{
    fit <- lm(expr[i,] ~ sv_matrix)
    expr[i,] <- expr[i,] - sv_matrix %*% coef(fit)[-1]
  }}
}}
write.table(expr, "{out_file}", sep="\\t", quote=FALSE)
'''
        else:
            raise ValueError(f"Unknown method: {method}")

        result = subprocess.run(
            ["Rscript", "-e", r_code],
            capture_output=True, text=True, timeout=300,
        )
        if result.returncode != 0:
            raise RuntimeError(f"R script failed: {result.stderr[:500]}")

        corrected = pd.read_csv(out_file, sep="\t", index_col=0)
        return corrected


def run_combat(logcpm: pd.DataFrame, batch: pd.Series,
               sample_type: pd.Series) -> pd.DataFrame:
    common = logcpm.columns.intersection(batch.index).intersection(sample_type.index)
    valid = _filter_small_batches(batch[common])
    common = common.intersection(valid)
    if len(common) < 10:
        raise ValueError("Too few samples after filtering small batches")
    return _run_r_script(logcpm[common], batch[common], sample_type[common], "ComBat")


def run_limma_remove_batch(logcpm: pd.DataFrame, batch: pd.Series,
                            sample_type: pd.Series) -> pd.DataFrame:
    common = logcpm.columns.intersection(batch.index).intersection(sample_type.index)
    valid = _filter_small_batches(batch[common])
    common = common.intersection(valid)
    if len(common) < 10:
        raise ValueError("Too few samples after filtering small batches")
    expr = logcpm[common]
    batch_aligned = batch[common]
    type_aligned = sample_type[common]

    r_matrix = _counts_to_r_matrix(expr)
    r_batch = ro.StrVector(batch_aligned.values.tolist())
    r_type = ro.StrVector(type_aligned.values.tolist())
    design = r_stats.model_matrix(ro.Formula("~type"),
                                  data=ro.DataFrame({"type": r_type}))
    corrected = limma.removeBatchEffect(r_matrix, batch=r_batch, design=design)
    result = np.array(list(corrected)).reshape(expr.shape[1], expr.shape[0]).T
    return pd.DataFrame(result, index=expr.index, columns=expr.columns)


def run_sva(logcpm: pd.DataFrame, sample_type: pd.Series) -> pd.DataFrame:
    common = logcpm.columns.intersection(sample_type.index)
    if len(common) < 10:
        raise ValueError("Too few samples")
    return _run_r_script(logcpm[common], pd.Series("A", index=common),
                         sample_type[common], "SVA")


def _variance_explained_by_batch(expr: pd.DataFrame, batch: pd.Series) -> float:
    common = expr.columns.intersection(batch.index)
    expr_aligned = expr[common]
    batch_aligned = batch[common]
    batch_counts = batch_aligned.value_counts()
    valid_batches = batch_counts[batch_counts >= 2].index
    mask = batch_aligned.isin(valid_batches)
    if mask.sum() < 10:
        return np.nan
    expr_filtered = expr_aligned.loc[:, mask]
    batch_filtered = batch_aligned[mask]
    n_pcs = min(5, expr_filtered.shape[1] - 1)
    if n_pcs < 1:
        return np.nan
    pca = PCA(n_components=n_pcs)
    pcs = pca.fit_transform(expr_filtered.T)
    r2_total = 0
    for pc_idx in range(n_pcs):
        pc_values = pcs[:, pc_idx]
        groups = [pc_values[batch_filtered.values == b] for b in valid_batches]
        groups = [g for g in groups if len(g) >= 1]
        if len(groups) >= 2:
            grand_mean = np.mean(pc_values)
            ss_total = np.sum((pc_values - grand_mean) ** 2)
            ss_between = sum(len(g) * (np.mean(g) - grand_mean) ** 2 for g in groups)
            r2 = ss_between / ss_total if ss_total > 0 else 0
            r2_total += r2 * pca.explained_variance_ratio_[pc_idx]
    return r2_total


def _biological_variance_retained(expr_corrected: pd.DataFrame,
                                   sample_type: pd.Series) -> float:
    common = expr_corrected.columns.intersection(sample_type.index)
    expr_aligned = expr_corrected[common]
    type_aligned = sample_type[common]
    n_pcs = min(5, expr_aligned.shape[1] - 1)
    if n_pcs < 1:
        return np.nan
    pca = PCA(n_components=n_pcs)
    pcs = pca.fit_transform(expr_aligned.T)
    r2_total = 0
    for pc_idx in range(n_pcs):
        pc_values = pcs[:, pc_idx]
        groups_t = pc_values[type_aligned.values == "Tumor"]
        groups_n = pc_values[type_aligned.values == "Normal"]
        if len(groups_t) >= 2 and len(groups_n) >= 2:
            grand_mean = np.mean(pc_values)
            ss_total = np.sum((pc_values - grand_mean) ** 2)
            ss_between = (len(groups_t) * (np.mean(groups_t) - grand_mean) ** 2 +
                          len(groups_n) * (np.mean(groups_n) - grand_mean) ** 2)
            r2 = ss_between / ss_total if ss_total > 0 else 0
            r2_total += r2 * pca.explained_variance_ratio_[pc_idx]
    return r2_total


def _ari_clustering(expr: pd.DataFrame, sample_type: pd.Series) -> float:
    common = expr.columns.intersection(sample_type.index)
    expr_aligned = expr[common]
    type_aligned = sample_type[common]
    n_pcs = min(10, expr_aligned.shape[1] - 1)
    if n_pcs < 2:
        return np.nan
    pca = PCA(n_components=n_pcs)
    pcs = pca.fit_transform(expr_aligned.T)
    true_labels = (type_aligned == "Tumor").astype(int).values
    kmeans = KMeans(n_clusters=2, random_state=42, n_init=10)
    pred_labels = kmeans.fit_predict(pcs)
    return adjusted_rand_score(true_labels, pred_labels)


def _deg_stability(counts: pd.DataFrame, sample_type: pd.Series,
                   batch: pd.Series) -> float:
    common = counts.columns.intersection(sample_type.index).intersection(batch.index)
    counts_aligned = counts[common]
    type_aligned = sample_type[common]
    batch_aligned = batch[common]
    unique_batches = sorted(batch_aligned.unique())
    if len(unique_batches) < 2:
        return np.nan
    mid = len(unique_batches) // 2
    batch_set1 = set(unique_batches[:mid])
    batch_set2 = set(unique_batches[mid:])
    mask1 = batch_aligned.isin(batch_set1)
    mask2 = batch_aligned.isin(batch_set2)
    for mask in [mask1, mask2]:
        if len(type_aligned[mask].unique()) < 2:
            return np.nan
        if (type_aligned[mask] == "Tumor").sum() < 3 or (type_aligned[mask] == "Normal").sum() < 3:
            return np.nan

    def _quick_edger(c, g):
        r_matrix = _counts_to_r_matrix(c)
        r_groups = ro.StrVector(g.values.tolist())
        dge = edgeR.DGEList(counts=r_matrix, group=r_groups)
        dge = edgeR.calcNormFactors(dge)
        design = r_stats.model_matrix(ro.Formula("~group"),
                                      data=ro.DataFrame({"group": r_groups}))
        dge = edgeR.estimateDisp(dge, design)
        fit = edgeR.glmQLFit(dge, design)
        qlf = edgeR.glmQLFTest(fit, coef=2)
        tt = edgeR.topTags(qlf, n=ro.IntVector([500]))
        return set(ro.r["rownames"](tt))

    try:
        degs1 = _quick_edger(counts_aligned.loc[:, mask1], type_aligned[mask1])
        degs2 = _quick_edger(counts_aligned.loc[:, mask2], type_aligned[mask2])
        inter = len(degs1 & degs2)
        union = len(degs1 | degs2)
        return inter / union if union > 0 else 0
    except Exception:
        return np.nan


def evaluate_correction(expr_before: pd.DataFrame, expr_after: pd.DataFrame,
                        counts: pd.DataFrame, batch: pd.Series,
                        sample_type: pd.Series, method: str,
                        project_id: str, cancer_type: str) -> dict:
    batch_var_before = _variance_explained_by_batch(expr_before, batch)
    batch_var_after = _variance_explained_by_batch(expr_after, batch)
    bio_var = _biological_variance_retained(expr_after, sample_type)
    ari_before = _ari_clustering(expr_before, sample_type)
    ari_after = _ari_clustering(expr_after, sample_type)
    deg_stab = _deg_stability(counts, sample_type, batch)

    return {
        "dataset_id": project_id,
        "cancer_type": cancer_type,
        "method": method,
        "batch_variance_before": round(batch_var_before, 4) if not np.isnan(batch_var_before) else np.nan,
        "batch_variance_after": round(batch_var_after, 4) if not np.isnan(batch_var_after) else np.nan,
        "biological_variance_retained": round(bio_var, 4) if not np.isnan(bio_var) else np.nan,
        "ari_before": round(ari_before, 4) if not np.isnan(ari_before) else np.nan,
        "ari_after": round(ari_after, 4) if not np.isnan(ari_after) else np.nan,
        "deg_stability": round(deg_stab, 4) if not np.isnan(deg_stab) else np.nan,
    }
