#!/usr/bin/env python
"""CP1 Deliverable A — CellBender per-gene contamination.

Run with the ambient_cb interpreter. Computes the SAME matrix-difference metric
used for SoupX/DecontX:
    contam_g = (orig_total_g - corrected_total_g) / orig_total_g
on CellBender's own called-cell set: orig = raw STARsolo counts for CB's cell
barcodes; corrected = CB _filtered.h5.

Usage: cp1_extract_cellbender.py <cb_filtered_h5> <starsolo_raw_dir> <out_prefix>
"""
import sys, json
import numpy as np
import scipy.io, scipy.sparse as sp
from cellbender.remove_background.downstream import anndata_from_h5

cb_h5, raw_dir, out_prefix = sys.argv[1], sys.argv[2], sys.argv[3]

# --- CellBender corrected (cells x genes) ---
a = anndata_from_h5(cb_h5)
# locate the Ensembl gene-id column in var
id_col = next((c for c in ("gene_id", "gene_ids", "id", "ID", "ensembl_id")
               if c in a.var.columns), None)
cb_gene_id = (a.var[id_col].astype(str).values if id_col is not None
              else a.var_names.astype(str).values)
cb_sym = (a.var["gene_name"].astype(str).values if "gene_name" in a.var.columns
          else (a.var["feature_name"].astype(str).values if "feature_name" in a.var.columns
                else cb_gene_id))
cb_bc = a.obs_names.astype(str).values
Xcorr = a.X if sp.issparse(a.X) else sp.csr_matrix(a.X)        # cells x genes
corr_total_by_gene = np.asarray(Xcorr.sum(axis=0)).ravel()     # per gene

# --- raw STARsolo (genes x cells, CellRanger MTX) ---
raw = scipy.io.mmread(f"{raw_dir}/matrix.mtx").tocsc()         # genes x cells
raw_gene_id = [l.split("\t")[0] for l in open(f"{raw_dir}/features.tsv")]
raw_bc = [l.strip() for l in open(f"{raw_dir}/barcodes.tsv")]
raw_bc_idx = {b: i for i, b in enumerate(raw_bc)}

# --- match CB cell barcodes back to raw ---
def norm(b):  # tolerate a trailing -1 difference
    return b[:-2] if b.endswith("-1") else b
raw_bc_norm = {norm(b): i for i, b in enumerate(raw_bc)}
cols, matched = [], 0
for b in cb_bc:
    j = raw_bc_idx.get(b, raw_bc_norm.get(norm(b)))
    if j is not None:
        cols.append(j); matched += 1
frac_bc = matched / len(cb_bc)
print(f"[CB] barcode match: {matched}/{len(cb_bc)} ({frac_bc:.3f})")
if frac_bc < 0.95:
    sys.exit(f"FATAL: only {frac_bc:.3f} of CB barcodes matched raw — barcode convention mismatch.")

raw_cells = raw[:, cols]                                       # genes x CB-cells
orig_total_by_gene_raw = np.asarray(raw_cells.sum(axis=1)).ravel()

# --- align genes: CB gene order vs raw features order, on Ensembl ID ---
raw_id_idx = {g: i for i, g in enumerate(raw_gene_id)}
n_match = sum(1 for g in cb_gene_id if g in raw_id_idx)
frac_gene = n_match / len(cb_gene_id)
print(f"[CB] gene-id match to raw: {n_match}/{len(cb_gene_id)} ({frac_gene:.3f})")
if frac_gene < 0.99:
    sys.exit(f"FATAL: only {frac_gene:.3f} of CB gene IDs found in raw features — ID convention mismatch (C1 USA-mode lesson).")

orig_total = np.array([orig_total_by_gene_raw[raw_id_idx[g]] if g in raw_id_idx else np.nan
                       for g in cb_gene_id])

import csv
with open(f"{out_prefix}_pergene.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["gene_id", "gene_symbol", "orig_total", "corr_total"])
    for gid, gsym, o, c in zip(cb_gene_id, cb_sym, orig_total, corr_total_by_gene):
        w.writerow([gid, gsym, f"{o:.6g}", f"{c:.6g}"])

gfr = 1 - (float(np.nansum(corr_total_by_gene)) / float(np.nansum(orig_total)))
json.dump({"tool": "CellBender", "n_cells": int(len(cb_bc)), "n_genes": int(len(cb_gene_id)),
           "barcode_match_frac": frac_bc, "gene_id_match_frac": frac_gene,
           "global_frac_removed": gfr},
          open(f"{out_prefix}_summary.json", "w"), indent=2)
print(f"[CB] DONE global_frac_removed={gfr:.4f}")
