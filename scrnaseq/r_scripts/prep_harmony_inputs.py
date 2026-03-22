#!/usr/bin/env python3
"""Prep PCA embeddings + donor labels for Harmony R script.

Exports the EXACT same datasets used in P3, with log-normalized + PCA
embeddings and donor_id labels. Harmony in R will correct the PCA space.

Usage:
    python prep_harmony_inputs.py
"""

import gc
from pathlib import Path

import cellxgene_census
import numpy as np
import pandas as pd
import scanpy as sc
import tiledbsoma

CENSUS_VERSION = "2025-11-08"
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = REPO_ROOT / "scrnaseq" / "projects" / "p3_batch_correction" / "harmony_inputs"
MAX_CELLS = 20_000
RNG_SEED = 42
N_PCS = 50
N_HVG = 2000

P3_TSV = REPO_ROOT / "scrnaseq" / "output" / "p3_batch_correction.tsv"


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    p3 = pd.read_csv(P3_TSV, sep="\t")
    datasets = p3.drop_duplicates(["dataset_id", "tissue", "organism"])[
        ["dataset_id", "tissue", "organism"]
    ].reset_index(drop=True)

    print(f"Exporting PCA embeddings for {len(datasets)} P3 datasets\n")

    census = cellxgene_census.open_soma(census_version=CENSUS_VERSION)

    from tqdm import tqdm
    for _, row in tqdm(datasets.iterrows(), total=len(datasets), desc="Prep Harmony"):
        ds_id = row.dataset_id
        organism = row.organism
        tissue = row.tissue
        org_key = organism.lower().replace(" ", "_")
        org_tag = "human" if "homo" in organism.lower() else "mouse"
        tag = f"{org_tag}_{ds_id[:8]}_{tissue.replace(' ', '_')}"

        pca_path = OUTPUT_DIR / f"{tag}_pca.tsv"
        meta_path = OUTPUT_DIR / f"{tag}_meta.tsv"

        if pca_path.exists() and meta_path.exists():
            print(f"  Skipping {tag} (exists)")
            continue

        print(f"  Processing {tag}...", flush=True)

        exp = census["census_data"][org_key]

        print(f"    Fetching obs metadata...", flush=True)
        obs_df = exp.obs.read(
            column_names=["soma_joinid", "donor_id", "cell_type"],
            value_filter=f"dataset_id == '{ds_id}' and is_primary_data == True",
        ).concat().to_pandas()

        if len(obs_df) < 500:
            print(f"    Too few cells ({len(obs_df)}), skipping")
            continue

        # Need ≥2 donors
        if obs_df.donor_id.nunique() < 2:
            print(f"    Only 1 donor, skipping")
            continue

        # Subsample
        if len(obs_df) > MAX_CELLS:
            rng = np.random.default_rng(RNG_SEED)
            idx = rng.choice(len(obs_df), MAX_CELLS, replace=False)
            obs_df = obs_df.iloc[sorted(idx)].copy()

        # Fetch expression
        print(f"    Fetching expression ({len(obs_df)} cells)...", flush=True)
        try:
            with exp.axis_query(
                measurement_name="RNA",
                obs_query=tiledbsoma.AxisQuery(coords=(obs_df.soma_joinid.tolist(),)),
            ) as q:
                adata = q.to_anndata(X_name="raw")
        except Exception as e:
            print(f"    Fetch failed: {e}", flush=True)
            continue

        print(f"    Fetched: {adata.shape[0]} cells x {adata.shape[1]} genes", flush=True)

        # Attach metadata
        ct_map = obs_df.set_index("soma_joinid")
        adata.obs["donor_id"] = adata.obs["soma_joinid"].map(ct_map["donor_id"])
        adata.obs["cell_type"] = adata.obs["soma_joinid"].map(ct_map["cell_type"])

        # Log-normalize + HVG + PCA (same pipeline as P3)
        print(f"    Normalizing + HVG + PCA...", flush=True)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        n_donors = adata.obs.donor_id.nunique()
        try:
            if n_donors <= 50:
                sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG,
                                            batch_key="donor_id")
            else:
                sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG)
        except Exception:
            sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG)

        adata_hvg = adata[:, adata.var.highly_variable].copy()
        sc.pp.scale(adata_hvg, max_value=10)
        sc.tl.pca(adata_hvg, n_comps=N_PCS, random_state=RNG_SEED)
        print(f"    PCA done: {adata_hvg.obsm['X_pca'].shape}", flush=True)

        pca = adata_hvg.obsm["X_pca"]

        # Save PCA embedding (cells x PCs)
        pca_df = pd.DataFrame(pca, columns=[f"PC{i+1}" for i in range(pca.shape[1])])
        pca_df.to_csv(pca_path, sep="\t", index=False)

        # Save metadata
        meta_df = pd.DataFrame({
            "donor_id": adata.obs["donor_id"].values,
            "cell_type": adata.obs["cell_type"].values,
        })
        meta_df.to_csv(meta_path, sep="\t", index=False)

        print(f"    Saved: {pca.shape[0]} cells, {pca.shape[1]} PCs, {n_donors} donors")

        del adata, adata_hvg
        gc.collect()

    census.close()
    print(f"\nDone. Inputs saved to {OUTPUT_DIR}")
    print("Now run run_harmony.R in RStudio.")


if __name__ == "__main__":
    main()
