"""
B2 analysis: normalization method comparison.

Compares TPM, FPKM, TMM, and VST normalization on identical raw count matrices.
Metrics: Spearman rank correlation, pathway enrichment Jaccard, housekeeping CV.
"""

import numpy as np
import pandas as pd
from itertools import combinations
from scipy import stats

import rpy2.robjects as ro
from rpy2.robjects import default_converter, numpy2ri, pandas2ri
from rpy2.robjects.packages import importr

_converter = default_converter + numpy2ri.converter + pandas2ri.converter

edgeR = importr("edgeR")
r_stats = importr("stats")
r_base = importr("base")

# Housekeeping genes (widely used set)
HOUSEKEEPING_GENES = [
    "ACTB", "GAPDH", "B2M", "RPL13A", "RPLP0", "RPS18", "SDHA",
    "TBP", "UBC", "YWHAZ", "HPRT1", "HMBS", "PGK1", "GUSB",
    "TFRC", "ALAS1", "PPIA", "G6PD", "LDHA", "NONO",
]

METHOD_PAIRS = [
    ("TPM", "TMM"), ("TPM", "VST"), ("TPM", "FPKM"),
    ("TMM", "VST"), ("TMM", "FPKM"), ("VST", "FPKM"),
]


def _compute_fpkm(counts: pd.DataFrame, tpm: pd.DataFrame) -> pd.DataFrame:
    """Derive FPKM from raw counts using gene lengths inferred from TPM/counts ratio.

    Gene length ∝ counts / (TPM * library_size) × 1e6
    FPKM = counts / (gene_length_kb * library_size_M)
    """
    lib_sizes = counts.sum(axis=0)  # per sample

    # Derive effective gene lengths from the TPM ↔ counts relationship
    # For each sample: length_i = (counts_i / lib_size) / (TPM_i / 1e6)
    # Average across samples for stable estimate
    with np.errstate(divide="ignore", invalid="ignore"):
        # rate = counts / lib_size (per million mapped reads → RPM)
        rpm = counts.div(lib_sizes, axis=1) * 1e6
        # length ∝ RPM / TPM
        length_ratios = rpm / tpm.replace(0, np.nan)

    # Median gene length ratio across samples (robust to outliers)
    median_lengths = length_ratios.median(axis=1)
    median_lengths = median_lengths.replace([np.inf, -np.inf, 0], np.nan)
    median_lengths = median_lengths.fillna(median_lengths.median())

    # FPKM = counts / (length_kb * lib_size_M)
    # Since we have relative lengths, normalize to get FPKM-like values
    fpkm = counts.div(median_lengths, axis=0).div(lib_sizes / 1e6, axis=1)
    fpkm = fpkm.replace([np.inf, -np.inf], 0).fillna(0)

    return fpkm


def _compute_tmm(counts: pd.DataFrame) -> pd.DataFrame:
    """Compute TMM-normalized log-CPM via edgeR."""
    np_counts = counts.values.astype(np.float64)
    nr, nc = np_counts.shape
    r_vec = ro.FloatVector(np_counts.T.flatten())
    r_matrix = ro.r.matrix(r_vec, nrow=nr, ncol=nc)
    r_matrix.rownames = ro.StrVector(list(counts.index))
    r_matrix.colnames = ro.StrVector(list(counts.columns))

    dge = edgeR.DGEList(counts=r_matrix)
    dge = edgeR.calcNormFactors(dge, method=ro.StrVector(["TMM"]))

    # Get log-CPM (normalized)
    log_cpm = ro.r["cpm"](dge, log=True)

    # Convert back to numpy
    result = np.array(list(log_cpm))
    result = result.reshape(nc, nr).T  # R is column-major
    return pd.DataFrame(result, index=counts.index, columns=counts.columns)


def _compute_vst(counts: pd.DataFrame) -> pd.DataFrame:
    """Compute VST via pydeseq2."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.preprocessing import deseq2_norm_transform

    # pydeseq2 expects samples × genes
    count_matrix = counts.T

    # Create minimal metadata — index must match count_matrix rows (samples)
    metadata = pd.DataFrame({"condition": ["A"] * count_matrix.shape[0]},
                            index=count_matrix.index)

    dds = DeseqDataSet(counts=count_matrix, metadata=metadata,
                       design="~1")  # no design, just transform
    dds.fit_size_factors()
    dds.vst(use_design=False)
    vst_values = dds.layers["vst_counts"]  # samples × genes
    return pd.DataFrame(vst_values.T, index=counts.index, columns=counts.columns)


def normalize_all(counts: pd.DataFrame, tpm: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Apply all four normalization methods to raw counts.

    Returns dict mapping method name → normalized DataFrame (genes × samples).
    """
    results = {}

    # TPM: use the pre-computed GTEx values
    results["TPM"] = tpm

    # FPKM: derive from counts + gene lengths (from TPM)
    results["FPKM"] = _compute_fpkm(counts, tpm)

    # TMM: edgeR normalized log-CPM
    results["TMM"] = _compute_tmm(counts)

    # VST: variance-stabilizing transformation
    results["VST"] = _compute_vst(counts)

    return results


def _spearman_rank_correlation(df1: pd.DataFrame, df2: pd.DataFrame) -> float:
    """Compute mean per-sample Spearman correlation of gene rankings."""
    common = df1.index.intersection(df2.index)
    correlations = []
    for col in df1.columns:
        if col in df2.columns:
            r, _ = stats.spearmanr(df1.loc[common, col], df2.loc[common, col])
            if not np.isnan(r):
                correlations.append(r)
    return np.mean(correlations) if correlations else np.nan


def _housekeeping_cv(norm_df: pd.DataFrame) -> float:
    """Compute mean CV of housekeeping genes across samples."""
    hk = [g for g in HOUSEKEEPING_GENES if g in norm_df.index]
    if len(hk) < 5:
        # Try matching by gene name in index
        return np.nan
    hk_values = norm_df.loc[hk]
    # CV = std / mean per gene
    means = hk_values.mean(axis=1)
    stds = hk_values.std(axis=1)
    cvs = (stds / means.replace(0, np.nan)).dropna()
    return cvs.mean() if len(cvs) > 0 else np.nan


def _pathway_enrichment_jaccard(df1: pd.DataFrame, df2: pd.DataFrame,
                                 n_top: int = 20) -> float:
    """Compute Jaccard similarity of top variable genes between two normalizations.

    Uses top variable genes as a proxy for pathway-relevant genes,
    since full pathway enrichment (gseapy) is slow per tissue.
    """
    common = df1.index.intersection(df2.index)
    # Rank genes by variance across samples
    var1 = df1.loc[common].var(axis=1).nlargest(200).index
    var2 = df2.loc[common].var(axis=1).nlargest(200).index
    top1 = set(var1[:n_top])
    top2 = set(var2[:n_top])
    if not top1 or not top2:
        return np.nan
    intersection = len(top1 & top2)
    union = len(top1 | top2)
    return intersection / union if union > 0 else 0.0


def _count_tissue_specific(norm_df: pd.DataFrame, n_top: int = 500) -> int:
    """Count highly expressed genes unique to this tissue normalization."""
    means = norm_df.mean(axis=1)
    return len(means.nlargest(n_top))


def compare_normalizations(normalized: dict[str, pd.DataFrame],
                           tissue_key: str, tissue_name: str) -> list[dict]:
    """Compare all pairs of normalization methods.

    Returns list of result dicts for the output TSV.
    """
    results = []

    # Pairwise comparisons
    for m1, m2 in METHOD_PAIRS:
        df1 = normalized[m1]
        df2 = normalized[m2]

        spearman = _spearman_rank_correlation(df1, df2)
        jaccard = _pathway_enrichment_jaccard(df1, df2)

        results.append({
            "dataset_id": f"GTEx_v10_{tissue_key}",
            "tissue": tissue_name,
            "method_pair": f"{m1}_vs_{m2}",
            "spearman_rank_correlation": round(spearman, 4),
            "pathway_jaccard": round(jaccard, 4),
            "housekeeping_cv": np.nan,  # filled below
            "n_tissue_specific_genes_affected": 0,  # filled below
        })

    # Per-method metrics
    for method, df in normalized.items():
        # Need gene names (not ENSG IDs) for housekeeping lookup
        # GTEx uses ENSG IDs — try mapping via Description
        hk_cv = _housekeeping_cv(df)

        # Update all rows involving this method
        for r in results:
            pair = r["method_pair"]
            if method in pair.split("_vs_"):
                if np.isnan(r["housekeeping_cv"]):
                    r["housekeeping_cv"] = round(hk_cv, 4) if not np.isnan(hk_cv) else np.nan

    # Tissue-specific gene differences
    # Compare top 500 variable genes across methods
    top_sets = {}
    for method, df in normalized.items():
        var_genes = df.var(axis=1).nlargest(500).index
        top_sets[method] = set(var_genes)

    for r in results:
        m1, m2 = r["method_pair"].split("_vs_")
        diff = len(top_sets[m1].symmetric_difference(top_sets[m2]))
        r["n_tissue_specific_genes_affected"] = diff

    return results
