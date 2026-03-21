"""
B1 analysis: DEG tool agreement between DESeq2, edgeR, and limma-voom.

Runs all three tools on identical count matrices and compares results.
"""

import numpy as np
import pandas as pd
from scipy import stats

import rpy2.robjects as ro
from rpy2.robjects import default_converter, numpy2ri, pandas2ri
from rpy2.robjects.packages import importr

_converter = default_converter + numpy2ri.converter + pandas2ri.converter

edgeR = importr("edgeR")
limma = importr("limma")
r_stats = importr("stats")
r_base = importr("base")

TOOL_PAIRS = [("DESeq2", "edgeR"), ("DESeq2", "limma"), ("edgeR", "limma")]


def _counts_to_r_matrix(counts: pd.DataFrame):
    """Convert pandas DataFrame to R matrix."""
    np_counts = counts.values.astype(np.float64)
    nr, nc = np_counts.shape
    r_vec = ro.FloatVector(np_counts.T.flatten())
    r_matrix = ro.r.matrix(r_vec, nrow=nr, ncol=nc)
    r_matrix.rownames = ro.StrVector(list(counts.index))
    r_matrix.colnames = ro.StrVector(list(counts.columns))
    return r_matrix


def _filter_low_counts(counts: pd.DataFrame, min_cpm: float = 1.0,
                       min_samples: int = 2) -> pd.DataFrame:
    """Apply standard low-count filter."""
    lib_sizes = counts.sum(axis=0)
    cpm = counts.div(lib_sizes, axis=1) * 1e6
    keep = (cpm >= min_cpm).sum(axis=1) >= min_samples
    return counts[keep]


def run_deseq2(counts: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    """Run DESeq2 via pydeseq2. Returns DataFrame with log2FC, pvalue, padj."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    common = counts.columns.intersection(groups.index)
    counts_aligned = counts[common]
    groups_aligned = groups[common]

    # pydeseq2 expects samples × genes
    count_matrix = counts_aligned.T
    metadata = pd.DataFrame({"condition": groups_aligned.values},
                            index=count_matrix.index)

    dds = DeseqDataSet(counts=count_matrix, metadata=metadata,
                       design="~condition")
    dds.deseq2()

    stat_res = DeseqStats(dds, contrast=["condition", "Tumor", "Normal"])
    stat_res.summary()

    results = stat_res.results_df
    return pd.DataFrame({
        "gene": results.index,
        "log2FoldChange": results["log2FoldChange"],
        "pvalue": results["pvalue"],
        "padj": results["padj"],
    }).set_index("gene")


def run_edger(counts: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    """Run edgeR glmQLFit. Returns DataFrame with log2FC, pvalue, padj."""
    common = counts.columns.intersection(groups.index)
    counts_aligned = counts[common]
    groups_aligned = groups[common]

    r_matrix = _counts_to_r_matrix(counts_aligned)
    r_groups = ro.StrVector(groups_aligned.values.tolist())

    dge = edgeR.DGEList(counts=r_matrix, group=r_groups)
    dge = edgeR.calcNormFactors(dge)

    design = r_stats.model_matrix(ro.Formula("~group"),
                                  data=ro.DataFrame({"group": r_groups}))
    dge = edgeR.estimateDisp(dge, design)
    fit = edgeR.glmQLFit(dge, design)
    qlf = edgeR.glmQLFTest(fit, coef=2)

    tt = edgeR.topTags(qlf, n=ro.IntVector([counts_aligned.shape[0]]))
    results = ro.r["as.data.frame"](tt)

    gene_names = list(ro.r["rownames"](results))
    logfc = np.array(list(results.rx2("logFC")))
    pval = np.array(list(results.rx2("PValue")))
    fdr = np.array(list(results.rx2("FDR")))

    return pd.DataFrame({
        "gene": gene_names,
        "log2FoldChange": logfc,
        "pvalue": pval,
        "padj": fdr,
    }).set_index("gene")


def run_limma_voom(counts: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    """Run limma-voom. Returns DataFrame with log2FC, pvalue, padj."""
    common = counts.columns.intersection(groups.index)
    counts_aligned = counts[common]
    groups_aligned = groups[common]

    r_matrix = _counts_to_r_matrix(counts_aligned)
    r_groups = ro.StrVector(groups_aligned.values.tolist())

    dge = edgeR.DGEList(counts=r_matrix, group=r_groups)
    dge = edgeR.calcNormFactors(dge)

    design = r_stats.model_matrix(ro.Formula("~group"),
                                  data=ro.DataFrame({"group": r_groups}))

    v = limma.voom(dge, design)
    fit = limma.lmFit(v, design)
    fit = limma.eBayes(fit)

    tt = limma.topTable(fit, coef=2, number=ro.IntVector([counts_aligned.shape[0]]),
                        sort_by=ro.StrVector(["none"]))

    gene_names = list(ro.r["rownames"](tt))
    logfc = np.array(list(tt.rx2("logFC")))
    pval = np.array(list(tt.rx2("P.Value")))
    adj_pval = np.array(list(tt.rx2("adj.P.Val")))

    return pd.DataFrame({
        "gene": gene_names,
        "log2FoldChange": logfc,
        "pvalue": pval,
        "padj": adj_pval,
    }).set_index("gene")


def compare_tools(deseq2_res: pd.DataFrame, edger_res: pd.DataFrame,
                  limma_res: pd.DataFrame, project_id: str,
                  cancer_type: str, tissue: str, n_samples: int) -> list[dict]:
    """Compare DEG results across all three tools.

    Returns list of result dicts (one per tool pair).
    """
    tool_results = {
        "DESeq2": deseq2_res,
        "edgeR": edger_res,
        "limma": limma_res,
    }

    # Get significant DEGs per tool (FDR < 0.05)
    sig_genes = {}
    n_degs = {}
    for tool, res in tool_results.items():
        sig = res[res["padj"] < 0.05].index
        sig_genes[tool] = set(sig)
        n_degs[tool] = len(sig)

    # Top 100 DEGs per tool (by adjusted p-value)
    top100 = {}
    for tool, res in tool_results.items():
        top = res.dropna(subset=["padj"]).nsmallest(100, "padj").index
        top100[tool] = set(top)

    results = []
    for t1, t2 in TOOL_PAIRS:
        res1 = tool_results[t1]
        res2 = tool_results[t2]
        common_genes = res1.index.intersection(res2.index)

        # Jaccard of top 100
        inter = len(top100[t1] & top100[t2])
        union = len(top100[t1] | top100[t2])
        jaccard = inter / union if union > 0 else 0

        # Direction agreement (among genes significant in both)
        both_sig = sig_genes[t1] & sig_genes[t2]
        both_sig_common = both_sig.intersection(common_genes)
        if len(both_sig_common) > 0:
            dir1 = np.sign(res1.loc[list(both_sig_common), "log2FoldChange"].values)
            dir2 = np.sign(res2.loc[list(both_sig_common), "log2FoldChange"].values)
            direction_agree = np.mean(dir1 == dir2)
        else:
            direction_agree = np.nan

        # Log2FC correlation (all common genes)
        lfc1 = res1.loc[common_genes, "log2FoldChange"].dropna()
        lfc2 = res2.loc[common_genes, "log2FoldChange"].dropna()
        lfc_common = lfc1.index.intersection(lfc2.index)
        if len(lfc_common) > 10:
            lfc_corr, _ = stats.pearsonr(lfc1[lfc_common], lfc2[lfc_common])
        else:
            lfc_corr = np.nan

        # Significance inflation (relative to DESeq2 as reference)
        if t1 == "DESeq2" and n_degs["DESeq2"] > 0:
            sig_inflation = n_degs[t2] / n_degs["DESeq2"]
        elif t2 == "DESeq2" and n_degs["DESeq2"] > 0:
            sig_inflation = n_degs[t1] / n_degs["DESeq2"]
        else:
            # For edgeR vs limma, report relative to DESeq2
            if n_degs["DESeq2"] > 0:
                sig_inflation = max(n_degs[t1], n_degs[t2]) / n_degs["DESeq2"]
            else:
                sig_inflation = np.nan

        results.append({
            "dataset_id": project_id,
            "cancer_type": cancer_type,
            "tissue": tissue,
            "n_samples": n_samples,
            "tool_pair": f"{t1}_vs_{t2}",
            "jaccard_top100": round(jaccard, 4),
            "direction_agreement": round(direction_agree, 4) if not np.isnan(direction_agree) else np.nan,
            "log2fc_correlation": round(lfc_corr, 4) if not np.isnan(lfc_corr) else np.nan,
            "significance_inflation": round(sig_inflation, 4) if not np.isnan(sig_inflation) else np.nan,
            "n_degs_deseq2": n_degs["DESeq2"],
            "n_degs_edger": n_degs["edgeR"],
            "n_degs_limma": n_degs["limma"],
        })

    return results
