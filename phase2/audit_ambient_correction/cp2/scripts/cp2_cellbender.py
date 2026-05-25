#!/usr/bin/env python
"""CP2 Deliverable B — CellBender under one ordering (reuses CP1 corrected output).

CellBender's correction is raw-based and ordering-invariant (it cannot be refit on
QC-passed cells — it needs the raw droplet distribution). So both orderings reuse
CP1's cb_filtered.h5; they differ ONLY in the QC-metric basis:
  O1 (correct->QC): C2 QC computed on CORRECTED counts.
  O2 (QC->correct): C2 QC computed on ORIGINAL (raw) counts for CellBender's cells.
This asymmetry vs SoupX/DecontX (which refit on survivors) is documented in findings.

Usage: cp2_cellbender.py <cb_filtered_h5> <starsolo_raw_dir> <starsolo_filt_features> <out_prefix> <O1|O2>
"""
import sys, json, csv
import numpy as np, scipy.io, scipy.sparse as sp
from cellbender.remove_background.downstream import anndata_from_h5

cb_h5, raw_dir, filt_features, out_prefix, ordering = sys.argv[1:6]

# canonical gene_id -> symbol (CellBender var has only Ensembl IDs; use STARsolo features)
sym = {}
for line in open(filt_features):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 2: sym[p[0]] = p[1]

a = anndata_from_h5(cb_h5)                              # cells x genes (corrected)
id_col = next((c for c in ("gene_id","gene_ids","id","ID") if c in a.var.columns), None)
gene_id = (a.var[id_col].astype(str).values if id_col else a.var_names.astype(str).values)
gene_sym = np.array([sym.get(g, g) for g in gene_id])
cb_bc = a.obs_names.astype(str).values
Xcorr = (a.X if sp.issparse(a.X) else sp.csr_matrix(a.X)).tocsr()   # cells x genes

# raw counts for CB cells
raw = scipy.io.mmread(f"{raw_dir}/matrix.mtx").tocsc()             # genes x cells
raw_gene = [l.split("\t")[0] for l in open(f"{raw_dir}/features.tsv")]
raw_bc = {b.strip(): i for i, b in enumerate(open(f"{raw_dir}/barcodes.tsv"))}
def norm(b): return b[:-2] if b.endswith("-1") else b
raw_bc_n = {norm(b): i for b, i in raw_bc.items()}
cols = [raw_bc.get(b, raw_bc_n.get(norm(b))) for b in cb_bc]
assert all(c is not None for c in cols), "barcode match failure"
rid = {g: i for i, g in enumerate(raw_gene)}
grows = [rid[g] for g in gene_id]                                  # align raw genes to CB order
Xorig = raw[np.ix_(grows, cols)].T.tocsr()                         # cells x genes (raw)

mito = np.array([str(s).startswith(("MT-", "mt-", "Mt-")) for s in gene_sym])
def c2_qc(X):                                                      # X: cells x genes
    tot = np.asarray(X.sum(axis=1)).ravel()
    ng = np.asarray((X > 0).sum(axis=1)).ravel()
    mt = np.asarray(X[:, mito].sum(axis=1)).ravel()
    mpct = np.where(tot > 0, 100 * mt / tot, 0)
    p95 = np.percentile(mpct, 95)
    return (ng >= 200) & (tot >= 500) & (mpct <= p95), p95

basis = Xcorr if ordering == "O1" else Xorig
keep, p95 = c2_qc(basis)
surv = cb_bc[keep]
orig_sub, corr_sub = Xorig[keep], Xcorr[keep]
n_in = len(cb_bc)
print(f"[CB {ordering}] QC on {'corrected' if ordering=='O1' else 'original'}: {n_in} -> {keep.sum()} cells (mito_p95={p95:.2f})")

og = np.asarray(orig_sub.sum(axis=0)).ravel(); cg = np.asarray(corr_sub.sum(axis=0)).ravel()
with open(f"{out_prefix}_pergene.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["gene_id","gene_symbol","orig_total","corr_total"])
    for g, s, o, c in zip(gene_id, gene_sym, og, cg): w.writerow([g, s, f"{o:.6g}", f"{c:.6g}"])
oc = np.asarray(orig_sub.sum(axis=1)).ravel(); cc = np.asarray(corr_sub.sum(axis=1)).ravel()
with open(f"{out_prefix}_percell.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["barcode","orig_total","corr_total"])
    for b, o, c in zip(surv, oc, cc): w.writerow([b, f"{o:.6g}", f"{c:.6g}"])
open(f"{out_prefix}_survivors.txt", "w").write("\n".join(surv) + "\n")
json.dump({"tool":"CellBender","ordering":ordering,"n_in":int(n_in),"n_survivors":int(keep.sum()),
           "retention":float(keep.sum()/n_in),"mito_p95":float(p95),
           "global_frac_removed":float(1 - cg.sum()/og.sum())},
          open(f"{out_prefix}_summary.json","w"), indent=2)
print(f"[CB {ordering}] DONE survivors={keep.sum()} retention={keep.sum()/n_in:.3f}")
