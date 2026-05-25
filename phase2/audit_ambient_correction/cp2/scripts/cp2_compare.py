#!/usr/bin/env python
"""CP2 Deliverable B — ordering analysis comparison + bootstrap.

Reads per_dataset/<ds>/<tool>_<O1|O2>_{pergene.tsv,percell.tsv,survivors.txt,summary.json}.
Computes (a) per-tool O1-vs-O2 comparison and (b) cross-tool comparison within each
ordering, plus bootstrap CIs. Run with ambient_cb interpreter.

Sanity check (per spec): CellBender's corrected matrix is identical across orderings
by design, so O1 vs O2 differ ONLY via the QC basis -> cell sets should differ small
but NON-zero. If CellBender O1 and O2 survivor sets are LITERALLY identical, the QC
basis switch is broken -> flagged loudly.
"""
import sys, os, json, glob
import numpy as np, pandas as pd
from scipy.stats import spearmanr

BASE = "/mnt/nvme1/omics-audit/phase2/audit_ambient_correction/cp2"
TOOLS = ["SoupX", "CellBender", "DecontX"]
PAIRS = [("SoupX", "CellBender"), ("SoupX", "DecontX"), ("CellBender", "DecontX")]
CHEM = {"10x_pbmc_1k_v3":"v3","10x_pbmc_5k_v3.1":"v3","10x_pbmc_10k_v3.1":"v3","10x_neuron_1k_v3":"v3",
        "gse287209_human_lung_organoid":"v3","gse325955_mouse_kidney_E18_5":"v3",
        "gse288156_mouse_intestine_scrna":"v3","10x_t_3k_v2":"v2","10x_pbmc_4k_v2":"v2"}
B = 1000; RNG = np.random.default_rng(42); EPS = 1e-6


def contam(ds, tool, ordering):
    p = f"{BASE}/per_dataset/{ds}/{tool}_{ordering}_pergene.tsv"
    if not os.path.exists(p): return None
    d = pd.read_csv(p, sep="\t"); d = d[d.orig_total > 0].copy()
    d["contam"] = ((d.orig_total - d.corr_total) / d.orig_total).clip(lower=0)
    return d.set_index("gene_id")["contam"]


def survivors(ds, tool, ordering):
    p = f"{BASE}/per_dataset/{ds}/{tool}_{ordering}_survivors.txt"
    return set(l.strip() for l in open(p)) if os.path.exists(p) else None


def percell(ds, tool, ordering):
    p = f"{BASE}/per_dataset/{ds}/{tool}_{ordering}_percell.tsv"
    return pd.read_csv(p, sep="\t").set_index("barcode")["corr_total"] if os.path.exists(p) else None


def jacc(a, b): return len(a & b) / len(a | b) if (a and b) else np.nan


def boot_sp(x, y):
    n = len(x)
    if n < 10: return (np.nan, np.nan)
    v = np.array([spearmanr(x[i], y[i]).statistic for i in (RNG.integers(0, n, n) for _ in range(B))])
    return (np.nanpercentile(v, 2.5), np.nanpercentile(v, 97.5))


def sp_common(a, b):
    c = a.index.intersection(b.index)
    return (spearmanr(a[c], b[c]).statistic, a[c].values, b[c].values, len(c)) if len(c) >= 10 else (np.nan, np.array([]), np.array([]), len(c))


def main(datasets):
    ord_rows, cross_rows, disp_rows, strat = [], [], [], []
    flags = []
    for ds in datasets:
        # ---- per-tool O1 vs O2 ----
        for t in TOOLS:
            sO1, sO2 = survivors(ds, t, "O1"), survivors(ds, t, "O2")
            cO1, cO2 = contam(ds, t, "O1"), contam(ds, t, "O2")
            if sO1 is None or sO2 is None: continue
            j = jacc(sO1, sO2)
            o1only, o2only = len(sO1 - sO2), len(sO2 - sO1)
            disp_rows.append({"dataset": ds, "tool": t, "n_O1": len(sO1), "n_O2": len(sO2),
                              "shared": len(sO1 & sO2), "O1_only": o1only, "O2_only": o2only, "jaccard": j})
            # CellBender identical-cell-set sanity check
            if t == "CellBender" and sO1 == sO2:
                flags.append(f"[CHECK FAIL] {ds} CellBender O1==O2 survivor sets identical — QC basis switch may be broken")
            rho_g, xa, xb, ng = sp_common(cO1, cO2) if (cO1 is not None and cO2 is not None) else (np.nan, np.array([]), np.array([]), 0)
            lo, hi = boot_sp(xa, xb) if ng >= 10 else (np.nan, np.nan)
            pcO1, pcO2 = percell(ds, t, "O1"), percell(ds, t, "O2")
            rho_cell = np.nan
            if pcO1 is not None and pcO2 is not None:
                shared = pcO1.index.intersection(pcO2.index)
                if len(shared) >= 10: rho_cell = spearmanr(pcO1[shared], pcO2[shared]).statistic
            ord_rows.append({"dataset": ds, "chemistry": CHEM.get(ds, "?"), "tool": t,
                             "cell_jaccard": j, "O1_only": o1only, "O2_only": o2only,
                             "pergene_contam_spearman": rho_g, "boot_lo": lo, "boot_hi": hi,
                             "percell_total_spearman": rho_cell, "n_genes": ng,
                             "corrected_identical_by_design": (t == "CellBender")})
            if t != "CellBender" and ng >= 10:
                strat.append((CHEM.get(ds, "?"), f"ordering_{t}", xa, xb))
        # ---- cross-tool within each ordering ----
        for ordering in ("O1", "O2"):
            cs = {t: contam(ds, t, ordering) for t in TOOLS}
            ss = {t: survivors(ds, t, ordering) for t in TOOLS}
            for a, b in PAIRS:
                if cs[a] is None or cs[b] is None: continue
                rho, xa, xb, ng = sp_common(cs[a], cs[b])
                cross_rows.append({"dataset": ds, "chemistry": CHEM.get(ds, "?"), "ordering": ordering,
                                   "pair": f"{a}_vs_{b}", "pergene_contam_spearman": rho, "n_genes": ng,
                                   "cell_jaccard": jacc(ss[a], ss[b])})
                if ng >= 10: strat.append((CHEM.get(ds, "?"), f"crosstool_{a}_vs_{b}_{ordering}", xa, xb))

    pd.DataFrame(ord_rows).to_csv(f"{BASE}/per_dataset_ordering_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(cross_rows).to_csv(f"{BASE}/per_dataset_cross_tool_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(disp_rows).to_csv(f"{BASE}/cell_disposition.tsv", sep="\t", index=False)

    # pooled bootstrap by chemistry x comparison-type
    srows = {}
    for chem, ctype, xa, xb in strat:
        srows.setdefault((chem, ctype), [[], []])
        srows[(chem, ctype)][0].append(xa); srows[(chem, ctype)][1].append(xb)
    out = []
    for (chem, ctype), (XA, XB) in sorted(srows.items()):
        xa, xb = np.concatenate(XA), np.concatenate(XB)
        rho = spearmanr(xa, xb).statistic; lo, hi = boot_sp(xa, xb)
        out.append({"chemistry": chem, "comparison": ctype, "n_genes_pooled": len(xa),
                    "spearman": rho, "boot_lo": lo, "boot_hi": hi})
    pd.DataFrame(out).to_csv(f"{BASE}/per_stratum_bootstrap.tsv", sep="\t", index=False)

    print("[cp2_compare] wrote ordering/cross-tool/disposition/bootstrap tables")
    if flags:
        print("\n".join(flags))
    else:
        print("[sanity] CellBender O1/O2 cell sets differ (non-identical) on all datasets — QC basis switch OK")
    print("\n=== per-tool O1-vs-O2 (cell_jaccard, pergene contam spearman) ===")
    od = pd.DataFrame(ord_rows)
    if not od.empty:
        print(od[["dataset","tool","cell_jaccard","O1_only","O2_only","pergene_contam_spearman","percell_total_spearman"]].to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1:])
