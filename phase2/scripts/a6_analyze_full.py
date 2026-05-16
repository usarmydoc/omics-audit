#!/usr/bin/env python3
"""Phase 2a A6 — full 26-dataset verification of §1.9a (muscat-dream inflation).

Computes per-dataset:
  - n_sig per tool (nebula, muscat_pb_DESeq2, muscat_mm_dream)
  - inflation: each tool's n_sig / muscat_pb_DESeq2 n_sig (pseudobulk as reference)
  - top-100 Jaccard: each cell-level tool vs pseudobulk
  - log2FC Spearman: each cell-level tool vs pseudobulk

Aggregate stats:
  - Distribution of muscat-dream inflation across the 26-dataset corpus
  - n datasets with >2×, >5×, >10× inflation
  - Bimodality test (visual: bimodal distribution implies §1.9a refined)
  - Tissue/cell-type annotation of inflated datasets

Output: phase2a/a6_analysis_full.tsv + findings text for the rule.
"""
from __future__ import annotations
from pathlib import Path
import sys
import pandas as pd
import numpy as np
from scipy import stats as scipy_stats

A6_DIR = Path("/mnt/nvme1/omics-audit/phase2a/a6_mixed_model")
SUMMARY = A6_DIR / "a6_summary.tsv"
OUT = Path("/mnt/nvme1/omics-audit/phase2a/a6_analysis_full.tsv")


def jaccard_top_k(a: pd.DataFrame, b: pd.DataFrame, k: int = 100) -> float:
    ta = set(a.dropna(subset=["padj"]).sort_values("padj").head(k)["gene_id"])
    tb = set(b.dropna(subset=["padj"]).sort_values("padj").head(k)["gene_id"])
    u = ta | tb
    return len(ta & tb) / len(u) if u else float("nan")


def spearman_log2fc(a: pd.DataFrame, b: pd.DataFrame) -> float:
    common = set(a["gene_id"]) & set(b["gene_id"])
    if len(common) < 100:
        return float("nan")
    a_idx = a.set_index("gene_id")
    b_idx = b.set_index("gene_id")
    a_lfc = a_idx.loc[list(common), "log2FC"]
    b_lfc = b_idx.loc[list(common), "log2FC"]
    rho, _ = scipy_stats.spearmanr(a_lfc, b_lfc, nan_policy="omit")
    return float(rho)


def main():
    summary = pd.read_csv(SUMMARY, sep="\t")
    # Pull tissue + cell_type for each dataset from the orchestrator log
    log_text = Path("/mnt/nvme1/omics-audit/phase2a/logs/a6_rerun_dedup_2026-05-16.log").read_text()
    ds_meta = {}
    for line in log_text.split("\n"):
        # Format: "[short_id] tissue/disease, donors=X+Y, cell_type=Z"
        if "] " in line and "cell_type=" in line:
            import re
            m = re.search(r"\[([a-f0-9]{8})\] (.+?), donors=\d+\+\d+, cell_type=(.+)$", line)
            if m:
                ds_meta[m.group(1)] = (m.group(2), m.group(3))

    rows = []
    for ds_dir in sorted(A6_DIR.iterdir()):
        if not ds_dir.is_dir():
            continue
        short = ds_dir.name
        files = {
            "nebula":           ds_dir / f"{short}__nebula.tsv",
            "muscat_pb_DESeq2": ds_dir / f"{short}__muscat_pb_DESeq2.tsv",
            "muscat_mm_dream":  ds_dir / f"{short}__muscat_mm_dream.tsv",
        }
        if not all(p.exists() for p in files.values()):
            continue

        dfs = {tool: pd.read_csv(p, sep="\t") for tool, p in files.items()}
        sig = {tool: int((df["padj"].fillna(1) < 0.05).sum()) for tool, df in dfs.items()}
        n_tested = {tool: len(df) for tool, df in dfs.items()}

        pb_sig = sig["muscat_pb_DESeq2"]
        pb_df = dfs["muscat_pb_DESeq2"]
        tissue, cell_type = ds_meta.get(short, ("?", "?"))

        row = {
            "dataset_short": short,
            "tissue": tissue,
            "cell_type": cell_type,
            "nebula_sig": sig["nebula"],
            "muscat_pb_sig": pb_sig,
            "muscat_mm_dream_sig": sig["muscat_mm_dream"],
            "nebula_n": n_tested["nebula"],
            "muscat_pb_n": n_tested["muscat_pb_DESeq2"],
            "muscat_mm_dream_n": n_tested["muscat_mm_dream"],
            "muscat_mm_pct_sig": round(100 * sig["muscat_mm_dream"] / n_tested["muscat_mm_dream"], 2) if n_tested["muscat_mm_dream"] else 0,
            "nebula_inflation_vs_pb": round(sig["nebula"] / max(pb_sig, 1), 3),
            "muscat_mm_inflation_vs_pb": round(sig["muscat_mm_dream"] / max(pb_sig, 1), 3),
            "jaccard_top100_nebula_vs_pb": round(jaccard_top_k(dfs["nebula"], pb_df, 100), 4),
            "jaccard_top100_mm_vs_pb": round(jaccard_top_k(dfs["muscat_mm_dream"], pb_df, 100), 4),
            "spearman_log2fc_nebula_vs_pb": round(spearman_log2fc(dfs["nebula"], pb_df), 4),
            "spearman_log2fc_mm_vs_pb": round(spearman_log2fc(dfs["muscat_mm_dream"], pb_df), 4),
        }
        rows.append(row)

    df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, sep="\t", index=False)
    print(f"Wrote {len(df)} dataset rows to {OUT}\n")

    # Aggregate stats for muscat-dream inflation
    print(f"=== muscat-dream inflation vs pseudobulk (n={len(df)} datasets) ===")
    mm_infl = df["muscat_mm_inflation_vs_pb"].replace([np.inf], np.nan).dropna()
    print(f"  median: {mm_infl.median():.3f}×")
    print(f"  mean:   {mm_infl.mean():.3f}×")
    print(f"  range:  [{mm_infl.min():.3f}, {mm_infl.max():.3f}]")
    print(f"  n > 2×:  {(mm_infl > 2).sum()} of {len(mm_infl)}")
    print(f"  n > 5×:  {(mm_infl > 5).sum()} of {len(mm_infl)}")
    print(f"  n > 10×: {(mm_infl > 10).sum()} of {len(mm_infl)}")

    # Bimodal check: by % sig
    print(f"\n=== muscat-dream % significant genes (bimodal check) ===")
    pct = df["muscat_mm_pct_sig"]
    print(f"  median: {pct.median():.2f}%")
    print(f"  range:  [{pct.min():.2f}, {pct.max():.2f}]")
    high = df[df["muscat_mm_pct_sig"] > 50].sort_values("muscat_mm_pct_sig", ascending=False)
    low = df[df["muscat_mm_pct_sig"] < 10].sort_values("muscat_mm_pct_sig")
    mid = df[(df["muscat_mm_pct_sig"] >= 10) & (df["muscat_mm_pct_sig"] <= 50)]
    print(f"  HIGH (>50% sig): {len(high)} datasets")
    for _, r in high.iterrows():
        print(f"    {r['dataset_short']}  {r['tissue'][:35]:35s} {r['cell_type'][:25]:25s} pct={r['muscat_mm_pct_sig']:.1f}%")
    print(f"  MID (10-50% sig): {len(mid)} datasets")
    print(f"  LOW (<10% sig):  {len(low)} datasets")
    for _, r in low.iterrows():
        print(f"    {r['dataset_short']}  {r['tissue'][:35]:35s} {r['cell_type'][:25]:25s} pct={r['muscat_mm_pct_sig']:.2f}%")

    # Jaccard + Spearman distributions
    print(f"\n=== top-100 Jaccard: muscat-dream vs pseudobulk ===")
    jac = df["jaccard_top100_mm_vs_pb"].dropna()
    print(f"  median: {jac.median():.3f}  range: [{jac.min():.3f}, {jac.max():.3f}]")
    print(f"  n < 0.30: {(jac < 0.30).sum()} of {len(jac)}")

    print(f"\n=== log2FC Spearman: muscat-dream vs pseudobulk ===")
    rho = df["spearman_log2fc_mm_vs_pb"].dropna()
    print(f"  median: {rho.median():.3f}  range: [{rho.min():.3f}, {rho.max():.3f}]")

    print(f"\n=== verdict vs §1.9a (Phase 1 5-dataset subset claim) ===")
    print(f"  §1.9a claim: muscat-dream inflated on lymph node + skin fibroblast specifically")
    print(f"  §1.9a 5-dataset metric: muscat-dream median inflation 0.85× vs pseudobulk")
    print(f"  A6 full 26-dataset: median {mm_infl.median():.2f}×, range [{mm_infl.min():.2f}, {mm_infl.max():.2f}]")
    if abs(mm_infl.median() - 0.85) < 0.5 and len(high) <= 4:
        print(f"  → CONSISTENT with §1.9a framing (similar to as_original)")
    elif len(high) >= 3:
        print(f"  → REFINES_PRIOR: bimodal pattern confirmed broader than just lymph node + skin fibroblast")
    else:
        print(f"  → ambiguous — review per-dataset table")

    # Register output
    sys.path.insert(0, "/mnt/nvme1/omics-audit/phase2/scripts")
    from dge_native import register_output
    register_output(OUT, kind="phase2a_a6_full_analysis", n_rows=len(df))


if __name__ == "__main__":
    main()
