#!/usr/bin/env python3
"""Prep input matrices for SingleR R script.

Exports each P5 dataset as a gzipped TSV (genes x cells) that
run_singler.R can load directly. Run this first, then run_singler.R
in RStudio.

Usage:
    python prep_singler_inputs.py
"""

import gzip
import sys
from pathlib import Path

import cellxgene_census
import numpy as np
import pandas as pd
import tiledbsoma

CENSUS_VERSION = "2025-11-08"
OUTPUT_DIR = Path(__file__).resolve().parent.parent / "projects" / "p5_annotation_tools" / "singler_inputs"
MAX_CELLS = 5000
MIN_CELLS = 500
RNG_SEED = 42

# Same cell type mapping as P5
def map_ct(ct):
    ct_l = ct.lower()
    if "endothelial" in ct_l: return "Endothelial"
    if "fibroblast" in ct_l or "stellate" in ct_l or "myofibroblast" in ct_l: return "Fibroblasts"
    if "macrophage" in ct_l or "kupffer" in ct_l: return "Macrophages"
    if "hepatocyte" in ct_l: return "Hepatocytes"
    if any(k in ct_l for k in ["t cell", "cd4-positive", "cd8-positive", "regulatory t",
                                 "memory t", "naive t", "effector t", "gamma-delta"]): return "T cells"
    if any(k in ct_l for k in ["natural killer", "nk cell"]): return "NK cells"
    return None


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    census = cellxgene_census.open_soma(census_version=CENSUS_VERSION)

    for organism in ["Homo sapiens", "Mus musculus"]:
        org_key = organism.lower().replace(" ", "_")
        org_tag = "human" if "homo" in organism.lower() else "mouse"
        exp = census["census_data"][org_key]

        print(f"\n=== {organism} ===")

        # Discover datasets (same logic as P5)
        obs_df = exp.obs.read(
            column_names=["dataset_id", "tissue_general", "cell_type"],
            value_filter="is_primary_data == True",
        ).concat().to_pandas()

        obs_df["mapped"] = obs_df.cell_type.apply(map_ct)
        obs_df = obs_df.dropna(subset=["mapped"])

        datasets = []
        for ds_id, grp in obs_df.groupby("dataset_id", observed=True):
            if pd.isna(ds_id): continue
            n_types = grp.mapped.nunique()
            n_cells = len(grp)
            if n_types < 3 or n_cells < MIN_CELLS or n_cells > 100_000: continue
            tissue = grp.tissue_general.mode().iloc[0]
            datasets.append({"dataset_id": ds_id, "tissue": tissue,
                           "n_cells": n_cells, "n_types": n_types})

        datasets = sorted(datasets, key=lambda x: (-x["n_types"], x["n_cells"]))
        # 1 per tissue, max 20
        seen_tissues = set()
        selected = []
        for d in datasets:
            if d["tissue"] not in seen_tissues:
                selected.append(d)
                seen_tissues.add(d["tissue"])
            if len(selected) >= 20:
                break

        print(f"Selected {len(selected)} datasets")

        for d in selected:
            ds_id = d["dataset_id"]
            tag = f"{org_tag}_{ds_id[:8]}_{d['tissue'].replace(' ', '_')}"
            out_path = OUTPUT_DIR / f"{tag}.tsv.gz"

            if out_path.exists():
                print(f"  Skipping {tag} (exists)")
                continue

            print(f"  Exporting {tag}...")

            # Fetch obs
            ds_obs = exp.obs.read(
                column_names=["soma_joinid", "cell_type"],
                value_filter=f"dataset_id == '{ds_id}' and is_primary_data == True",
            ).concat().to_pandas()

            ds_obs["mapped"] = ds_obs.cell_type.apply(map_ct)
            ds_obs = ds_obs.dropna(subset=["mapped"])

            # Subsample
            if len(ds_obs) > MAX_CELLS:
                rng = np.random.default_rng(RNG_SEED)
                frac = MAX_CELLS / len(ds_obs)
                parts = []
                for _, grp in ds_obs.groupby("mapped", observed=True):
                    n = max(1, int(len(grp) * frac))
                    parts.append(grp.sample(n=n, random_state=RNG_SEED, replace=False))
                ds_obs = pd.concat(parts)

            # Fetch expression
            try:
                with exp.axis_query(
                    measurement_name="RNA",
                    obs_query=tiledbsoma.AxisQuery(coords=(ds_obs.soma_joinid.tolist(),)),
                ) as q:
                    adata = q.to_anndata(X_name="raw")
            except Exception as e:
                print(f"    Fetch failed: {e}")
                continue

            # Set gene names
            if "feature_name" in adata.var.columns:
                adata.var_names = adata.var.feature_name.astype(str)
                adata.var_names_make_unique()

            # Filter genes expressed in ≥3 cells
            if hasattr(adata.X, "toarray"):
                X_dense = adata.X.toarray()
            else:
                X_dense = np.asarray(adata.X)

            gene_mask = (X_dense > 0).sum(axis=0) >= 3
            X_filt = X_dense[:, gene_mask].astype(int)
            genes_filt = [g for g, m in zip(adata.var_names, gene_mask) if m]

            # Save as genes x cells TSV (gzipped)
            # Cell IDs encode the mapped type for later use
            cell_ids = [f"{ds_obs.iloc[i]['mapped']}_{i}" for i in range(len(ds_obs))]
            if len(cell_ids) != X_filt.shape[0]:
                cell_ids = [f"cell_{i}" for i in range(X_filt.shape[0])]

            df = pd.DataFrame(X_filt.T, index=genes_filt,
                            columns=cell_ids[:X_filt.shape[0]])

            with gzip.open(out_path, "wt") as f:
                df.to_csv(f, sep="\t")

            print(f"    Saved: {X_filt.shape[0]} cells, {len(genes_filt)} genes → {out_path.name}")

            # Also save ground truth labels
            gt_path = OUTPUT_DIR / f"{tag}_labels.tsv"
            labels_df = pd.DataFrame({
                "cell_id": cell_ids[:X_filt.shape[0]],
                "ground_truth": ds_obs["mapped"].values[:X_filt.shape[0]],
            })
            labels_df.to_csv(gt_path, sep="\t", index=False)

    census.close()
    print(f"\nDone. Inputs saved to {OUTPUT_DIR}")
    print("Now run run_singler.R in RStudio.")


if __name__ == "__main__":
    main()
