#!/usr/bin/env python
"""CP3 — convert CellBender cb_filtered.h5 to a condition matrix dir.
CellBender correction is ordering-invariant -> one matrix per dataset (the CellBender
condition). Rounds posterior-mean corrected counts to integers for downstream
(scDblFinder/scry expect counts), matching SoupX/DecontX roundToInt.
Usage: cp3_cellbender_matrix.py <cb_filtered_h5> <out_dir>
"""
import sys, os
import numpy as np, scipy.io as sio, scipy.sparse as sp
from cellbender.remove_background.downstream import anndata_from_h5

cb_h5, out_dir = sys.argv[1], sys.argv[2]
os.makedirs(out_dir, exist_ok=True)
a = anndata_from_h5(cb_h5)
id_col = next((c for c in ("gene_id", "gene_ids", "id", "ID") if c in a.var.columns), None)
ensembl = (a.var[id_col].astype(str).values if id_col else a.var_names.astype(str).values)
X = (a.X if sp.issparse(a.X) else sp.csr_matrix(a.X)).tocsr()      # cells x genes
X.data = np.rint(X.data)                                          # integer counts
X.eliminate_zeros()
sio.mmwrite(f"{out_dir}/matrix.mtx", X.tocoo())
open(f"{out_dir}/barcodes.txt", "w").write("\n".join(a.obs_names.astype(str)) + "\n")
open(f"{out_dir}/ensembl.txt", "w").write("\n".join(ensembl) + "\n")
print(f"[CB matrix] wrote {X.shape[0]} cells x {X.shape[1]} genes -> {out_dir}")
