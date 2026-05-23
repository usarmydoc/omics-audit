#!/usr/bin/env python3
"""CP1 step 1: pull each working-set dataset from Census, compute the 3 QC
metrics per cell via scanpy, save per-cell metrics TSV. Only the QC metrics
are retained (the 4 filtering methods operate on them) — full matrices are
not kept. Census version pinned to 2025-11-08 for reproducibility.
Run in base env (cellxgene_census + scanpy)."""
import sys, time
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc
import cellxgene_census as cc

BASE = Path("/mnt/nvme1/omics-audit/phase2/audit_qc_mad/cp1")
OUT = BASE / "qc_metrics"; OUT.mkdir(parents=True, exist_ok=True)
CENSUS_VERSION = "2025-11-08"
CELL_CAP = 50000
SEED = 0

ws = pd.read_csv(BASE / "working_set.tsv", sep="\t")
census = cc.open_soma(census_version=CENSUS_VERSION)
print(f"census {CENSUS_VERSION} open", flush=True)

for _, row in ws.iterrows():
    short, did, tissue = row["short"], row["dataset_id"], row["tissue"]
    outf = OUT / f"{short}.tsv"
    if outf.exists():
        print(f"  [{short}] exists, skip", flush=True); continue
    t0 = time.time()
    # organism: census splits by experiment; try human first then mouse
    adata = None
    for org in ("homo_sapiens", "mus_musculus"):
        try:
            a = cc.get_anndata(census, organism=org,
                               obs_value_filter=f"dataset_id == '{did}'",
                               column_names={"var": ["feature_id", "feature_name"]})
            if a.n_obs > 0:
                adata = a; organism = org; break
        except Exception as e:
            print(f"    [{short}] {org} pull err: {e}", flush=True)
    if adata is None or adata.n_obs == 0:
        print(f"  [{short}] NO CELLS for {did} — FAIL", flush=True)
        (OUT / f"{short}.FAILED").write_text(f"no cells for {did}\n"); continue
    # cap cells
    if adata.n_obs > CELL_CAP:
        rng = np.random.default_rng(SEED)
        keep = rng.choice(adata.n_obs, CELL_CAP, replace=False)
        adata = adata[keep].copy()
    # mito genes by symbol prefix (human MT-, mouse mt-)
    sym = adata.var["feature_name"].astype(str)
    adata.var["mt"] = sym.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
    df = adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].copy()
    df.insert(0, "tissue", tissue); df.insert(0, "short", short)
    df.index.name = "cell"
    df.to_csv(outf, sep="\t")
    print(f"  [{short}] {organism} {adata.n_obs} cells, {int(adata.var['mt'].sum())} mt genes, "
          f"{time.time()-t0:.0f}s -> {outf.name}", flush=True)

census.close()
print("pull_qc done", flush=True)
