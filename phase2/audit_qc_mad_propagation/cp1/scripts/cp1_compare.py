#!/usr/bin/env python
"""CP1 Step 4-5 — compare QC methods' downstream biology (C2 reference) + bootstrap.

Per dataset: each method vs C2 (C1/MAD3/MAD5 vs C2) + cross pairs. ARI/NMI on
leiden (r1.0/0.5/1.5) over shared cells, cluster-matched top-50 marker Jaccard,
annotation agreement. Method-specific cell disposition (cells C2 drops but MADk
keeps → annotated as what). Bootstrap CIs (B=1000). audit3_counting env.
Usage: cp1_compare.py <short> [...]
"""
import sys, os, itertools
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr  # noqa
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

CP1 = Path("/mnt/nvme1/omics-audit/phase2/audit_qc_mad_propagation/cp1")
METHODS = ["C1", "C2", "MAD3", "MAD5"]
REF = "C2"
B = 1000; RNG = np.random.default_rng(42)


def load(short, m):
    p = CP1 / "per_condition" / short / f"{m}.obs.tsv"
    o = pd.read_csv(p, sep="\t", index_col=0) if p.exists() else None
    mk = CP1 / "per_condition" / short / f"{m}.markers.tsv"
    return o, (pd.read_csv(mk, sep="\t") if mk.exists() else None)


def marker_jacc(mka, mkb, clua, club, shared):
    if mka is None or mkb is None: return np.nan
    la, lb = clua.loc[shared], club.loc[shared]; js = []
    for ca in la.unique():
        cells_a = set(shared[la == ca]); best, bo = None, 0
        for cb in lb.unique():
            ov = len(cells_a & set(shared[lb == cb]))
            if ov > bo: bo, best = ov, cb
        if best is None: continue
        ga = set(mka[mka.cluster.astype(str) == str(ca)].gene.head(50))
        gb = set(mkb[mkb.cluster.astype(str) == str(best)].gene.head(50))
        if ga or gb: js.append(len(ga & gb) / len(ga | gb))
    return float(np.mean(js)) if js else np.nan


def boot(a, b, fn):
    n = len(a)
    if n < 20: return (np.nan, np.nan)
    v = np.array([fn(a[i], b[i]) for i in (RNG.integers(0, n, n) for _ in range(B))])
    return (np.percentile(v, 2.5), np.percentile(v, 97.5))


def main(shorts):
    rows, disp = [], []
    pairs = [(m, REF) for m in METHODS if m != REF] + [("C1", "MAD3"), ("C1", "MAD5"), ("MAD3", "MAD5")]
    for short in shorts:
        D = {m: load(short, m) for m in METHODS}
        D = {m: v for m, v in D.items() if v[0] is not None}
        for a, b in pairs:
            if a not in D or b not in D: continue
            oa, mka = D[a]; ob, mkb = D[b]
            sh = oa.index.intersection(ob.index)
            if len(sh) < 20: continue
            res = {"dataset": short, "comparison": f"{a}_vs_{b}", "method_a": a, "method_b": b,
                   "n_intersection": len(sh)}
            for rk in ("r1.0", "r0.5", "r1.5"):
                ca = oa.loc[sh, f"leiden_{rk}"].astype(str).values
                cb = ob.loc[sh, f"leiden_{rk}"].astype(str).values
                res[f"ari_{rk}"] = adjusted_rand_score(ca, cb)
                if rk == "r1.0":
                    res["nmi_r1.0"] = normalized_mutual_info_score(ca, cb)
                    lo, hi = boot(ca, cb, adjusted_rand_score); res["ari_r1.0_lo"], res["ari_r1.0_hi"] = lo, hi
            res["marker_jaccard_median"] = marker_jacc(mka, mkb, oa.loc[sh, "leiden_r1.0"].astype(str),
                                                        ob.loc[sh, "leiden_r1.0"].astype(str), sh)
            la = oa.loc[sh, "celltypist_label"].astype(str).values
            lb = ob.loc[sh, "celltypist_label"].astype(str).values
            res["annotation_agreement_pct"] = float((la == lb).mean() * 100)
            rows.append(res)
        # method-specific cell disposition: cells C2 drops but MADk/C1 keeps -> annotated as what
        if REF in D:
            ref_cells = set(D[REF][0].index)
            for m in [x for x in METHODS if x != REF and x in D]:
                om = D[m][0]; extra = [c for c in om.index if c not in ref_cells]
                if not extra: continue
                vc = om.loc[extra, "celltypist_label"].astype(str).value_counts().head(5)
                for lbl, n in vc.items():
                    disp.append({"dataset": short, "method": m, "vs": REF,
                                 "cells_kept_not_in_C2": len(extra), "annotated_label": lbl, "n": int(n)})

    rdf = pd.DataFrame(rows)
    rdf.to_csv(CP1 / "per_comparison_metrics.tsv", sep="\t", index=False)
    # per_stratum_bootstrap: derived from the SAME inline CIs (no redundant re-bootstrap)
    rdf[["dataset", "comparison", "n_intersection", "ari_r1.0", "ari_r1.0_lo", "ari_r1.0_hi"]].rename(
        columns={"n_intersection": "n", "ari_r1.0_lo": "boot_lo", "ari_r1.0_hi": "boot_hi"}
    ).to_csv(CP1 / "per_stratum_bootstrap.tsv", sep="\t", index=False)
    pd.DataFrame(disp).to_csv(CP1 / "method_specific_cell_disposition.tsv", sep="\t", index=False)
    print("[cp1_compare] wrote per_comparison / bootstrap / disposition tables")
    r = pd.DataFrame(rows)
    if not r.empty:
        print(r[["dataset", "comparison", "n_intersection", "ari_r1.0", "nmi_r1.0",
                 "marker_jaccard_median", "annotation_agreement_pct"]].to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1:] or ["liver", "small_intestine", "blood"])
