#!/usr/bin/env python
"""CP3 Deliverable C — biological-propagation comparison.

Per dataset, compares each corrected condition to the no-correction baseline
(ARI/NMI on leiden_r1.0, annotation agreement, marker Jaccard on cluster-matched
top-50 sets) and runs cross-condition pairwise comparisons. Bootstrap CIs on
ARI/NMI (cell-level, B=1000). Run with audit3_counting interpreter (sklearn).
Usage: cp3_compare.py <intestine|pbmc> [...]
"""
import sys, os, json, itertools
import numpy as np, pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

BASE = "/mnt/nvme1/omics-audit/phase2/audit_ambient_correction/cp3"
CONDS = ["nocorr", "SoupX_O1", "SoupX_O2", "CellBender", "DecontX_O1", "DecontX_O2"]
B = 1000; RNG = np.random.default_rng(42)


def load(sub, cond):
    o = f"{BASE}/{sub}/{cond}.obs.tsv"; m = f"{BASE}/{sub}/{cond}.markers.tsv"
    if not os.path.exists(o): return None
    obs = pd.read_csv(o, sep="\t", index_col=0)
    mk = pd.read_csv(m, sep="\t") if os.path.exists(m) else None
    return {"obs": obs, "markers": mk}


def marker_jaccard(mk_a, mk_b, clu_a, clu_b, shared):
    """Match a-clusters to b-clusters by cell overlap on shared cells, then top-50 gene Jaccard."""
    if mk_a is None or mk_b is None: return np.nan
    la = clu_a.loc[shared]; lb = clu_b.loc[shared]
    ja = []
    for ca in la.unique():
        cells_a = set(shared[la == ca])
        # best-matching b cluster by cell overlap
        best, bestov = None, 0
        for cb in lb.unique():
            ov = len(cells_a & set(shared[lb == cb]))
            if ov > bestov: bestov, best = ov, cb
        if best is None: continue
        ga = set(mk_a[mk_a.cluster.astype(str) == str(ca)].gene.head(50))
        gb = set(mk_b[mk_b.cluster.astype(str) == str(best)].gene.head(50))
        if ga or gb: ja.append(len(ga & gb) / len(ga | gb))
    return float(np.mean(ja)) if ja else np.nan


def boot_ari(a, b):
    n = len(a)
    if n < 20: return (np.nan, np.nan)
    v = np.array([adjusted_rand_score(a[i], b[i]) for i in (RNG.integers(0, n, n) for _ in range(B))])
    return (np.percentile(v, 2.5), np.percentile(v, 97.5))


def main(subs):
    per_cond, cross, disp = [], [], []
    for sub in subs:
        D = {c: load(sub, c) for c in CONDS}
        D = {c: v for c, v in D.items() if v is not None}
        if "nocorr" not in D:
            print(f"[WARN] {sub}: no baseline"); continue
        base = D["nocorr"]; bclu = base["obs"]["leiden_r1.0"].astype(str); blab = base["obs"]["celltypist_label"].astype(str)
        for c, v in D.items():
            obs = v["obs"]; clu = obs["leiden_r1.0"].astype(str); lab = obs["celltypist_label"].astype(str)
            shared = base["obs"].index.intersection(obs.index)
            a, b = bclu.loc[shared].values, clu.loc[shared].values
            ari = adjusted_rand_score(a, b); nmi = normalized_mutual_info_score(a, b)
            lo, hi = boot_ari(a, b) if c != "nocorr" else (np.nan, np.nan)
            ann_agree = float((blab.loc[shared].values == lab.loc[shared].values).mean())
            mj = marker_jaccard(base["markers"], v["markers"], bclu, clu, shared) if c != "nocorr" else 1.0
            per_cond.append({"dataset": sub, "condition": c, "n_cells": len(obs),
                             "n_shared_with_base": len(shared),
                             "n_clusters_r1.0": clu.nunique(), "n_labels": lab.nunique(),
                             "ARI_vs_base": ari if c != "nocorr" else np.nan,
                             "NMI_vs_base": nmi if c != "nocorr" else np.nan,
                             "ari_boot_lo": lo, "ari_boot_hi": hi,
                             "annotation_agreement_vs_base": ann_agree if c != "nocorr" else 1.0,
                             "marker_jaccard_vs_base": mj})
        # cross-condition pairwise (ARI + annotation agreement on shared cells)
        for ca, cb in itertools.combinations([c for c in CONDS if c in D], 2):
            oa, ob = D[ca]["obs"], D[cb]["obs"]
            sh = oa.index.intersection(ob.index)
            if len(sh) < 20: continue
            a, b = oa["leiden_r1.0"].astype(str).loc[sh].values, ob["leiden_r1.0"].astype(str).loc[sh].values
            la, lb = oa["celltypist_label"].astype(str).loc[sh].values, ob["celltypist_label"].astype(str).loc[sh].values
            cross.append({"dataset": sub, "pair": f"{ca}_vs_{cb}", "n_shared": len(sh),
                          "ARI": adjusted_rand_score(a, b), "NMI": normalized_mutual_info_score(a, b),
                          "annotation_agreement": float((la == lb).mean())})
        # contested cells: where SoupX_O1/DecontX_O1/CellBender disagree with baseline annotation
        for c in [x for x in ("SoupX_O1", "DecontX_O1", "CellBender") if x in D]:
            obs = D[c]["obs"]; sh = base["obs"].index.intersection(obs.index)
            bl = blab.loc[sh].values; cl = obs["celltypist_label"].astype(str).loc[sh].values
            mism = sh[bl != cl]
            for bcat in pd.Series(blab.loc[mism].values).value_counts().head(5).items():
                disp.append({"dataset": sub, "condition": c, "baseline_label": bcat[0], "n_reannotated": int(bcat[1])})

    pd.DataFrame(per_cond).to_csv(f"{BASE}/{'intestine' if 'intestine' in subs else 'pbmc'}_per_condition_metrics.tsv", sep="\t", index=False) if len(subs) == 1 else None
    pd.DataFrame(per_cond).to_csv(f"{BASE}/all_per_condition_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(cross).to_csv(f"{BASE}/cross_condition_comparison.tsv", sep="\t", index=False)
    pd.DataFrame(disp).to_csv(f"{BASE}/contested_cell_disposition.tsv", sep="\t", index=False)
    print("[cp3_compare] wrote per_condition / cross_condition / contested tables")
    print(pd.DataFrame(per_cond)[["dataset","condition","n_cells","n_clusters_r1.0","ARI_vs_base","annotation_agreement_vs_base","marker_jaccard_vs_base"]].to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1:] or ["pbmc", "intestine"])
