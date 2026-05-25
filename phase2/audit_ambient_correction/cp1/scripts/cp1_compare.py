#!/usr/bin/env python
"""CP1 Deliverable A — cross-tool per-gene contamination comparison + bootstrap.

Reads per_dataset/<ds>/{SoupX,DecontX,CellBender}_pergene.tsv (cols:
gene_id, gene_symbol, orig_total, corr_total), computes per-gene contamination
contam_g = (orig_total - corr_total)/orig_total on genes with orig_total>0,
aligns the 3 tools by Ensembl gene_id, and emits comparison + bootstrap tables.

Stop conditions (C1 USA-mode + CP4 §3.5 lessons):
  - gene-id overlap across tools < 0.99  -> abort (ID convention mismatch)
  - any tool global contamination ~0 or ~1 across all genes -> abort (config error)

Run with ambient_cb interpreter (numpy/scipy/pandas present).
"""
import sys, os, json
import numpy as np, pandas as pd
from scipy.stats import spearmanr

BASE = "/mnt/nvme1/omics-audit/phase2/audit_ambient_correction/cp1"
TOOLS = ["SoupX", "CellBender", "DecontX"]
PAIRS = [("SoupX", "CellBender"), ("SoupX", "DecontX"), ("CellBender", "DecontX")]
CHEM = {  # v2 vs v3 (collapsing v3.1 into v3)
    "10x_pbmc_1k_v3": "v3", "10x_pbmc_5k_v3.1": "v3", "10x_pbmc_10k_v3.1": "v3",
    "10x_neuron_1k_v3": "v3", "gse287209_human_lung_organoid": "v3",
    "gse325955_mouse_kidney_E18_5": "v3", "gse288156_mouse_intestine_scrna": "v3",
    "10x_t_3k_v2": "v2", "10x_pbmc_4k_v2": "v2",
}
B = 1000
RNG = np.random.default_rng(42)
EPS = 1e-6


def load_contam(ds, tool):
    """Return (full_gene_id_set, nonzero_df). full set = ID universe for the C1
    convention check; nonzero_df = genes with orig_total>0 (contam defined)."""
    p = f"{BASE}/per_dataset/{ds}/{tool}_pergene.tsv"
    if not os.path.exists(p):
        return None
    d = pd.read_csv(p, sep="\t")
    full_ids = set(d["gene_id"])
    nz = d[d["orig_total"] > 0].copy()
    nz["contam"] = ((nz["orig_total"] - nz["corr_total"]) / nz["orig_total"]).clip(lower=0)
    return full_ids, nz[["gene_id", "gene_symbol", "orig_total", "contam"]].set_index("gene_id")


def is_mito(sym):
    s = str(sym)
    return s.startswith("MT-") or s.startswith("mt-") or s.startswith("Mt-")


def boot_spearman(x, y, b=B):
    n = len(x)
    if n < 10:
        return (np.nan, np.nan)
    vals = np.empty(b)
    for i in range(b):
        idx = RNG.integers(0, n, n)
        r = spearmanr(x[idx], y[idx]).statistic
        vals[i] = r
    return (np.nanpercentile(vals, 2.5), np.nanpercentile(vals, 97.5))


def main(datasets):
    long_rows, metric_rows = [], []
    pooled = {}  # (chem, pair) -> list of (contamA, contamB) arrays

    for ds in datasets:
        loaded = {t: load_contam(ds, t) for t in TOOLS}
        present = [t for t in TOOLS if loaded[t] is not None]
        if len(present) < 2:
            print(f"[WARN] {ds}: <2 tools present ({present}); skipping")
            continue
        full = {t: loaded[t][0] for t in present}   # gene-id universe per tool
        df = {t: loaded[t][1] for t in present}      # nonzero contam frame per tool

        # C1-lesson check on ID CONVENTION: the full gene-id universes must match
        # (ENSG namespace, version suffixes, etc.). This is distinct from the
        # nonzero-expression set, which legitimately differs because CellBender
        # uses its own called cells vs STARsolo filtered.
        for a, b in PAIRS:
            if a in full and b in full:
                u = len(full[a] & full[b]) / max(len(full[a]), len(full[b]))
                if u < 0.999:
                    print(f"[ABORT] {ds} {a}vs{b}: gene-id UNIVERSE overlap {u:.4f} < 0.999 "
                          f"(real ID-convention mismatch — C1 USA-mode lesson)")
                    sys.exit(3)

        # canonical gene_id -> symbol map (always from SoupX/DecontX; CellBender's
        # symbol col is the Ensembl ID). Fixes mito strat + mouse/human prefixes.
        sym_src = df.get("SoupX") if df.get("SoupX") is not None else df.get("DecontX")
        canonical_sym = sym_src["gene_symbol"].to_dict() if sym_src is not None else {}

        # long format (canonical symbol for all tools)
        for t in present:
            for gid, row in df[t].iterrows():
                long_rows.append([ds, t, gid, canonical_sym.get(gid, row["gene_symbol"]),
                                  row["orig_total"], row["contam"]])

        # sanity: contamination not degenerate
        for t in present:
            mu = float(np.mean(df[t]["contam"].values))
            if mu < 1e-5 or mu > 0.999:
                print(f"[ABORT] {ds}/{t}: degenerate mean contamination {mu:.6f} (config error?)")
                sys.exit(2)

        for a, b in PAIRS:
            if a not in df or b not in df:
                continue
            da, db = df[a], df[b]
            common = da.index.intersection(db.index)  # genes nonzero in BOTH tools
            xa = da.loc[common, "contam"].values
            xb = db.loc[common, "contam"].values
            expr = da.loc[common, "orig_total"].values  # canonical expression (STARsolo filtered for SoupX/DecontX)
            sym = np.array([canonical_sym.get(g, "") for g in common])  # canonical symbols (not CellBender's Ensembl-as-symbol)

            rho = spearmanr(xa, xb).statistic
            lo, hi = boot_spearman(xa, xb)
            l2 = np.log2((xa + EPS) / (xb + EPS))
            mito = np.array([is_mito(s) for s in sym])
            rho_mt = spearmanr(xa[mito], xb[mito]).statistic if mito.sum() >= 5 else np.nan
            rho_nonmt = spearmanr(xa[~mito], xb[~mito]).statistic if (~mito).sum() >= 10 else np.nan
            # expression deciles
            dec = pd.qcut(pd.Series(expr).rank(method="first"), 10, labels=False)
            rho_lowexpr = spearmanr(xa[dec <= 1], xb[dec <= 1]).statistic if (dec <= 1).sum() >= 10 else np.nan
            rho_highexpr = spearmanr(xa[dec >= 8], xb[dec >= 8]).statistic if (dec >= 8).sum() >= 10 else np.nan

            metric_rows.append({
                "dataset": ds, "chemistry": CHEM.get(ds, "?"), "pair": f"{a}_vs_{b}",
                "n_genes": len(common), "spearman": rho, "boot_lo": lo, "boot_hi": hi,
                "log2ratio_median": float(np.median(l2)), "log2ratio_iqr": float(np.subtract(*np.percentile(l2, [75, 25]))),
                "mean_contam_A": float(np.mean(xa)), "mean_contam_B": float(np.mean(xb)),
                "directional": a if np.mean(xa) > np.mean(xb) else b,  # which estimates higher contamination
                "spearman_mito": rho_mt, "spearman_nonmito": rho_nonmt,
                "spearman_lowexpr": rho_lowexpr, "spearman_highexpr": rho_highexpr,
            })
            pooled.setdefault((CHEM.get(ds, "?"), f"{a}_vs_{b}"), []).append((xa, xb))

    # write long + per-dataset metrics
    pd.DataFrame(long_rows, columns=["dataset", "tool", "gene_id", "gene_symbol", "orig_total", "contam"]
                 ).to_csv(f"{BASE}/per_gene_contamination.tsv", sep="\t", index=False)
    mdf = pd.DataFrame(metric_rows)
    mdf.to_csv(f"{BASE}/per_dataset_metrics.tsv", sep="\t", index=False)

    # per-stratum bootstrap (chemistry x pair, genes pooled across datasets in stratum)
    strat_rows = []
    for (chem, pair), arrs in sorted(pooled.items()):
        xa = np.concatenate([a for a, _ in arrs]); xb = np.concatenate([b for _, b in arrs])
        rho = spearmanr(xa, xb).statistic
        lo, hi = boot_spearman(xa, xb)
        strat_rows.append({"chemistry": chem, "pair": pair, "n_datasets": len(arrs),
                           "n_genes_pooled": len(xa), "spearman": rho, "boot_lo": lo, "boot_hi": hi})
    pd.DataFrame(strat_rows).to_csv(f"{BASE}/per_stratum_bootstrap.tsv", sep="\t", index=False)
    print("[compare] wrote per_gene_contamination.tsv, per_dataset_metrics.tsv, per_stratum_bootstrap.tsv")
    print(mdf.to_string(index=False) if not mdf.empty else "(no metrics)")


if __name__ == "__main__":
    main(sys.argv[1:])
