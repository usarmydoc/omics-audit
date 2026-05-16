#!/usr/bin/env python3
"""Run a6_nebula_muscat.R against ONE real dataset for debugging.

Fetches the smallest P4-eligible dataset from Census, writes inputs,
runs the R script with full stdout/stderr capture printed to terminal.
"""
from __future__ import annotations
import gzip, subprocess, sys, time, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy import io as sio
from scipy.sparse import csr_matrix

warnings.filterwarnings("ignore")
P4_ROOT = Path("/mnt/nvme1/omics-audit/scrnaseq_audit")
PHASE2_ROOT = Path("/mnt/nvme1/omics-audit/phase2")
sys.path.insert(0, str(P4_ROOT))
from projects.p4_pseudobulk.run_p4 import (  # noqa: E402
    discover_eligible_datasets, fetch_expression_for_de,
)

OUT = Path("/mnt/nvme1/omics-audit/phase2a/a6_debug_one")
OUT.mkdir(parents=True, exist_ok=True)
R_SCRIPT = PHASE2_ROOT / "scripts" / "a6_nebula_muscat.R"

import cellxgene_census
print("Opening Census...")
census = cellxgene_census.open_soma(census_version="2025-11-08")

print("Discovering eligible (human)...")
eligible = discover_eligible_datasets(census, "Homo sapiens")
eligible = eligible.assign(
    _cost=eligible["n_cells_disease"].fillna(0) + eligible["n_cells_normal"].fillna(0)
).sort_values("_cost")
row = eligible.iloc[0]
short = row.dataset_id[:8]
print(f"Smallest dataset: {short}, {row.tissue}, {row.best_cell_type}, "
      f"donors={row.n_donors_disease}+{row.n_donors_normal}")

print("Fetching expression...")
fetched = fetch_expression_for_de(census, row.dataset_id, "Homo sapiens", row.best_cell_type)
X, obs_df, gene_names = fetched
print(f"  {len(obs_df)} cells, {len(gene_names)} genes")

# Write inputs
ds_out = OUT / short
ds_out.mkdir(parents=True, exist_ok=True)
counts_p = ds_out / "counts.mtx.gz"
meta_p = ds_out / "meta.tsv"
genes_p = ds_out / "genes.tsv"
X_T = X.T.tocoo() if hasattr(X, "tocsc") else csr_matrix(np.asarray(X).T).tocoo()
with gzip.open(counts_p, "wb") as gz:
    sio.mmwrite(gz, X_T, field="integer")
pd.DataFrame({
    "cell_id": [f"cell_{i}" for i in range(len(obs_df))],
    "donor_id": obs_df["donor_id"].astype(str).values,
    "group": obs_df["group"].astype(str).values,
}).to_csv(meta_p, sep="\t", index=False)
pd.Series(gene_names).to_csv(genes_p, sep="\t", index=False, header=False)
print(f"Inputs written to {ds_out}")

print(f"\n=== Running R script (stdout + stderr live) ===\n")
t0 = time.time()
proc = subprocess.run(
    ["Rscript", "--vanilla", str(R_SCRIPT),
     str(counts_p), str(meta_p), str(genes_p), str(ds_out), short],
    timeout=14400,
)
print(f"\n=== Done in {time.time()-t0:.1f}s, rc={proc.returncode} ===")
print(f"\nOutputs:")
for p in sorted(ds_out.glob(f"{short}__*.tsv")):
    df = pd.read_csv(p, sep="\t")
    n_sig = int((df["padj"].fillna(1) < 0.05).sum())
    print(f"  {p.name}: {len(df)} genes, {n_sig} sig")
