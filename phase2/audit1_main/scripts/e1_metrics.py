#!/usr/bin/env python3
"""Audit 1 main CP3 — E1 metrics + bootstrap CIs.

Consumes e1/runs/<input_id>__<tool>.tsv (375 files from e1_run.py) and
computes pairwise tool-agreement metrics per input + bootstrap CIs per
stratum.

Metrics per (tool_pair × input):
  - top_10_jaccard, top_25_jaccard (Hallmark = 50 pathways; top-50 = 1.0)
  - fdr05_overlap_jaccard, fdr01_overlap_jaccard
  - spearman_full_rank — full 50-pathway ranking
  - pearson_neglog10padj — correlation of -log10(padj) values
  - direction_agreement — NES or log2OR sign concordance (where both
    tools emit a direction; NA when one or both are ORA-only without OR)

Bootstrap CIs at the stratum level (input_category × tool_pair):
  - dataset-level resampling, B=1000 per AUDIT_STANDARDS.md §2.4
  - report median + 2.5% / 97.5% percentiles

Outputs:
  e1/per_input_metrics.tsv (375 input_id × tool_pair × metric rows)
  e1/per_stratum_bootstrap.tsv (9 stratum × metric rows)
  e1/E1_findings.md
"""
from __future__ import annotations
from itertools import combinations
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

ROOT = Path("/mnt/nvme1/omics-audit/phase2/audit1_main")
RUNS = ROOT / "e1/runs"
INVENTORY = ROOT / "datasets/audit1_main_inputs.tsv"
PER_INPUT = ROOT / "e1/per_input_metrics.tsv"
PER_STRATUM = ROOT / "e1/per_stratum_bootstrap.tsv"
FINDINGS = ROOT / "e1/E1_findings.md"
TOOLS = ["fgsea", "gseapy_enrichr", "clusterProfiler_ORA"]
B = 1000
RNG = np.random.default_rng(42)


def load_run(input_id: str, tool: str) -> pd.DataFrame | None:
    p = RUNS / f"{input_id}__{tool}.tsv"
    if not p.exists() or p.stat().st_size < 100:
        return None
    df = pd.read_csv(p, sep="\t")
    return df if len(df) > 0 else None


def top_k_jaccard(a: pd.DataFrame, b: pd.DataFrame, k: int) -> float:
    ta = set(a.dropna(subset=["padj"]).sort_values("padj").head(k)["pathway_name"])
    tb = set(b.dropna(subset=["padj"]).sort_values("padj").head(k)["pathway_name"])
    u = ta | tb
    return len(ta & tb) / len(u) if u else float("nan")


def fdr_overlap_jaccard(a: pd.DataFrame, b: pd.DataFrame, thr: float) -> float:
    ta = set(a[a["padj"].fillna(1) < thr]["pathway_name"])
    tb = set(b[b["padj"].fillna(1) < thr]["pathway_name"])
    u = ta | tb
    return len(ta & tb) / len(u) if u else float("nan")


def spearman_full_rank(a: pd.DataFrame, b: pd.DataFrame) -> float:
    common = set(a["pathway_name"]) & set(b["pathway_name"])
    if len(common) < 10:
        return float("nan")
    ax = a.set_index("pathway_name").loc[list(common), "padj"]
    bx = b.set_index("pathway_name").loc[list(common), "padj"]
    if ax.isna().all() or bx.isna().all():
        return float("nan")
    rho, _ = scipy_stats.spearmanr(ax, bx, nan_policy="omit")
    return float(rho)


def pearson_neglog10padj(a: pd.DataFrame, b: pd.DataFrame) -> float:
    common = set(a["pathway_name"]) & set(b["pathway_name"])
    if len(common) < 10:
        return float("nan")
    ax = -np.log10(a.set_index("pathway_name").loc[list(common), "padj"].clip(1e-300))
    bx = -np.log10(b.set_index("pathway_name").loc[list(common), "padj"].clip(1e-300))
    r, _ = scipy_stats.pearsonr(ax, bx)
    return float(r)


def direction_agreement(a: pd.DataFrame, tool_a: str,
                        b: pd.DataFrame, tool_b: str) -> float:
    """Sign concordance of NES (fgsea) or log2_odds_ratio (ORA) on shared pathways."""
    def get_dir(df: pd.DataFrame, tool: str) -> pd.Series:
        if tool == "fgsea":
            col = "NES"
        else:
            col = "log2_odds_ratio"
        if col not in df.columns:
            return pd.Series(dtype=float)
        x = df.set_index("pathway_name")[col]
        return x.dropna()

    da = get_dir(a, tool_a)
    db = get_dir(b, tool_b)
    common = da.index.intersection(db.index)
    if len(common) < 5:
        return float("nan")
    sign_a = np.sign(da.loc[common])
    sign_b = np.sign(db.loc[common])
    nonzero = (sign_a != 0) & (sign_b != 0)
    if nonzero.sum() < 5:
        return float("nan")
    return float((sign_a[nonzero] == sign_b[nonzero]).mean())


METRIC_FNS = {
    "top_10_jaccard":          lambda a, b, ta, tb: top_k_jaccard(a, b, 10),
    "top_25_jaccard":          lambda a, b, ta, tb: top_k_jaccard(a, b, 25),
    "fdr05_overlap_jaccard":   lambda a, b, ta, tb: fdr_overlap_jaccard(a, b, 0.05),
    "fdr01_overlap_jaccard":   lambda a, b, ta, tb: fdr_overlap_jaccard(a, b, 0.01),
    "spearman_full_rank":      lambda a, b, ta, tb: spearman_full_rank(a, b),
    "pearson_neglog10padj":    lambda a, b, ta, tb: pearson_neglog10padj(a, b),
    "direction_agreement":     lambda a, b, ta, tb: direction_agreement(a, ta, b, tb),
}


def compute_per_input():
    inv = pd.read_csv(INVENTORY, sep="\t")
    rows = []
    skipped = 0
    for _, ir in inv.iterrows():
        input_id = ir["input_id"]
        loaded = {t: load_run(input_id, t) for t in TOOLS}
        if any(v is None for v in loaded.values()):
            skipped += 1
            continue
        for ta, tb in combinations(TOOLS, 2):
            row = {
                "input_id": input_id,
                "comparison_id": ir["comparison_id"],
                "input_category": ir["input_category"],
                "tool_a": ta,
                "tool_b": tb,
            }
            for name, fn in METRIC_FNS.items():
                row[name] = fn(loaded[ta], loaded[tb], ta, tb)
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(PER_INPUT, sep="\t", index=False)
    print(f"Wrote {len(df)} per-input metric rows to {PER_INPUT}; "
          f"skipped {skipped} inputs missing one or more tool outputs")
    return df


def bootstrap_ci(values: np.ndarray, b: int = B) -> tuple[float, float, float]:
    """Return (median, lo_2.5%, hi_97.5%) over b dataset-level bootstrap resamples."""
    vals = values[~np.isnan(values)]
    if len(vals) < 3:
        med = float(np.nanmedian(values)) if len(vals) > 0 else float("nan")
        return med, float("nan"), float("nan")
    n = len(vals)
    medians = np.empty(b)
    for i in range(b):
        idx = RNG.integers(0, n, n)
        medians[i] = np.nanmedian(vals[idx])
    return float(np.median(vals)), float(np.percentile(medians, 2.5)), float(np.percentile(medians, 97.5))


def compute_per_stratum(per_input: pd.DataFrame):
    rows = []
    metric_cols = list(METRIC_FNS.keys())
    for category, gcat in per_input.groupby("input_category"):
        for pair_key, gpair in gcat.groupby(["tool_a", "tool_b"]):
            ta, tb = pair_key
            row = {
                "input_category": category,
                "tool_a": ta,
                "tool_b": tb,
                "n_inputs": len(gpair),
            }
            for m in metric_cols:
                med, lo, hi = bootstrap_ci(gpair[m].values)
                row[f"{m}_median"] = med
                row[f"{m}_ci_lo"]  = lo
                row[f"{m}_ci_hi"]  = hi
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(PER_STRATUM, sep="\t", index=False)
    print(f"Wrote {len(df)} stratum rows (3 categories × 3 tool pairs) to {PER_STRATUM}")
    return df


def write_findings(per_stratum: pd.DataFrame):
    lines = ["# E1 — Tool agreement on identical DEG inputs\n",
             "**Database held constant:** MSigDB Hallmark (50 pathways, Hs + Mm)",
             "**Tools:** fgsea (GSEA-family), gseapy.enrichr (ORA), clusterProfiler ORA",
             "**Bootstrap:** dataset-level resampling, B=1000",
             ""]
    # Pair summary table
    lines.append("## Pairwise agreement summary (medians with 95% bootstrap CI)\n")
    lines.append("| Category | Tool pair | n | top_10 Jaccard | FDR<0.05 Jaccard | Spearman full | Direction agreement |")
    lines.append("|---|---|---:|---:|---:|---:|---:|")
    def fmt(med, lo, hi):
        if np.isnan(med): return "—"
        if np.isnan(lo) or np.isnan(hi): return f"{med:.3f}"
        return f"{med:.3f} [{lo:.3f}, {hi:.3f}]"
    for _, r in per_stratum.iterrows():
        pair = f"{r['tool_a']} × {r['tool_b']}"
        lines.append(
            f"| {r['input_category']} | {pair} | {r['n_inputs']} | "
            f"{fmt(r['top_10_jaccard_median'], r['top_10_jaccard_ci_lo'], r['top_10_jaccard_ci_hi'])} | "
            f"{fmt(r['fdr05_overlap_jaccard_median'], r['fdr05_overlap_jaccard_ci_lo'], r['fdr05_overlap_jaccard_ci_hi'])} | "
            f"{fmt(r['spearman_full_rank_median'], r['spearman_full_rank_ci_lo'], r['spearman_full_rank_ci_hi'])} | "
            f"{fmt(r['direction_agreement_median'], r['direction_agreement_ci_lo'], r['direction_agreement_ci_hi'])} |"
        )

    lines += [
        "",
        "## Interpretation (per AUDIT_STANDARDS.md §3.4 — companion metrics framing)",
        "",
        "Jaccard answers 'do the tools agree on which pathways pass FDR'.",
        "Spearman answers 'do the tools rank pathways the same way'.",
        "Direction agreement answers 'when both tools agree, do they assign",
        "the same biological direction (up/down)'.",
        "",
        "If Jaccard is low but Spearman + direction are high: tools largely",
        "agree on the underlying biology and rank order, but disagree on",
        "WHICH pathways pass the FDR cutoff — a multiple-testing-correction",
        "and threshold artifact, not a fundamental biological disagreement.",
        "",
        "If Jaccard AND Spearman are both low: tools fundamentally disagree;",
        "user must pick one and document the choice.",
        "",
        "## Caveats and gaps",
        "",
        "1. **Direction agreement is fgsea × clusterProfiler_ORA only.** ",
        "   gseapy.enrichr's log2_odds_ratio is currently uniformly positive ",
        "   (Enrichr-style OR > 1 for any enriched term), so its 'direction' ",
        "   carries no discrimination. Direction-agreement metrics with ",
        "   gseapy as one tool collapse to a constant.",
        "2. **top_50_jaccard omitted.** Hallmark has 50 pathways total, so ",
        "   top-50 == full pathway set; Jaccard would be 1.0 by construction. ",
        "   Replaced with top-10 / top-25.",
        "3. **GSVA, EGSEA, camera not included.** Those require expression ",
        "   matrices, not DEG TSVs; queued as E1b in DEFERRED.md.",
        "",
    ]
    FINDINGS.write_text("\n".join(lines) + "\n")
    print(f"Wrote {FINDINGS}")


def main():
    per_input = compute_per_input()
    per_stratum = compute_per_stratum(per_input)
    write_findings(per_stratum)
    sys.path.insert(0, "/mnt/nvme1/omics-audit/phase2/scripts")
    from dge_native import register_output
    for f in [PER_INPUT, PER_STRATUM, FINDINGS]:
        register_output(f, kind="audit1_main_e1_metric")
    print("Registered metrics + findings in lock.")


if __name__ == "__main__":
    main()
