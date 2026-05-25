#!/usr/bin/env python
"""CP1 Step 1 — filtered Census pull, soma_joinid-matched to QC-MAD's exact cells.

Run in BASE env (cellxgene_census 1.17.0 + scanpy). Census version PINNED.
QC-MAD's qc_metrics `cell` column is the 0-based POSITION within the
dataset_id-filtered, soma_joinid-ordered full pull (NOT the global soma_joinid).
So: ordered obs read -> map QC-MAD positions to global soma_joinids -> obs_coords
pull exactly those 50k cells (efficient; never materializes blood's 2.6M X).
Sanity-check recomputed median QC metrics vs QC-MAD's saved values.
"""
import sys
import numpy as np, pandas as pd, scanpy as sc
import cellxgene_census as cc

CENSUS = "2025-11-08"   # PINNED
QCM = "/mnt/nvme1/omics-audit/phase2/audit_qc_mad/cp1/qc_metrics"
OUT = "/mnt/nvme1/omics-audit/phase2/audit_qc_mad_propagation/cp1/matrices"
DS = {"liver": "34f5307e-7b4d-4a48-b68f-2ba844c6414b",
      "small_intestine": "a37f857c-779f-464e-9310-3db43a1811e7",
      "blood": "46104f0b-9af5-466a-ae0f-56b8dc1969a2"}

census = cc.open_soma(census_version=CENSUS)
print(f"census {CENSUS} open", flush=True)
ok = True
for short, did in DS.items():
    qm = pd.read_csv(f"{QCM}/{short}.tsv", sep="\t")
    pos = qm["cell"].astype(int).values                    # 0-based positions, not soma_joinid
    # ordered obs read (soma_joinid ascending) -> map position -> global soma_joinid
    obs = census["census_data"]["homo_sapiens"].obs.read(
        value_filter=f"dataset_id == '{did}'", column_names=["soma_joinid"]).concat().to_pandas()
    obs = obs.sort_values("soma_joinid").reset_index(drop=True)
    gids = [int(x) for x in obs["soma_joinid"].values[pos]]
    a = cc.get_anndata(census, organism="Homo sapiens", obs_coords=gids)
    a.var_names = a.var["feature_id"].astype(str).values   # Ensembl
    a.obs_names = a.obs["soma_joinid"].astype(str).values
    # sanity: recompute QC metrics, compare medians to QC-MAD saved
    a.var["mt"] = a.var["feature_name"].astype(str).str.upper().str.startswith(("MT-", "MT."))
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None)
    got = (a.n_obs, np.median(a.obs.n_genes_by_counts), np.median(a.obs.pct_counts_mt))
    exp = (len(qm), qm.n_genes_by_counts.median(), qm.pct_counts_mt.median())
    match = (got[0] == exp[0]) and abs(got[1]-exp[1]) <= max(2, 0.01*exp[1]) and abs(got[2]-exp[2]) <= 0.05
    print(f"[{short}] pulled {a.n_obs} cells x {a.n_vars} genes | "
          f"median n_genes got/exp {got[1]:.0f}/{exp[1]:.0f}, mito {got[2]:.2f}/{exp[2]:.2f} | "
          f"{'MATCH' if match else 'MISMATCH!!'}", flush=True)
    if not match:
        ok = False; print(f"[{short}] SANITY FAIL — Census data may have changed", flush=True)
    a.write_h5ad(f"{OUT}/{short}_raw.h5ad")
census.close()
print("PULL_DONE" if ok else "PULL_SANITY_FAILED", flush=True)
sys.exit(0 if ok else 1)
