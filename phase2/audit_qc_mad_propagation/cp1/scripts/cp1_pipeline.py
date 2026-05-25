#!/usr/bin/env python
"""CP1 Step 3 — downstream pipeline on one (dataset x QC method) condition.

Mirrors Audit 3 CP6 / Ambient CP3 EXACTLY, except: the QC method IS the cell
filter, so the pipeline SKIPS its own cell QC (min_genes/max_mito) — it runs from
gene-filter onward on the method's survivor cells. Census var supplies symbols.
Reuses cp6_scdblfinder.R + cp6_deviance_hvg.R. Run in audit3_counting env.
Usage: cp1_pipeline.py <short> <method>
"""
import sys, time, json, tempfile, subprocess, csv
from pathlib import Path
import numpy as np, scipy.io as sio, scipy.sparse as sp
import anndata as ad, scanpy as sc, pandas as pd

CP6 = Path("/mnt/nvme1/omics-audit/phase2/audit3_counting/scripts")
CP1 = Path("/mnt/nvme1/omics-audit/phase2/audit_qc_mad_propagation/cp1")
sc.settings.verbosity = 1
QC_MIN_CELLS, N_HVG, N_PCS, K = 3, 2000, 50, 15
LEIDEN = {"r1.0": 1.0, "r0.5": 0.5, "r1.5": 1.5}
MODEL = {"liver": "Healthy_Human_Liver.pkl", "small_intestine": "Cells_Intestinal_Tract.pkl",
         "blood": "Immune_All_Low.pkl"}


def run_r(script, *a):
    r = subprocess.run(["Rscript", str(CP6 / script), *map(str, a)], capture_output=True, text=True)
    if r.returncode != 0: raise RuntimeError(f"{script} failed:\n{r.stdout}\n{r.stderr}")
    return r.stdout


def main():
    short, method = sys.argv[1], sys.argv[2]
    t0 = time.time(); model = MODEL[short]
    outdir = CP1 / "per_condition" / short; outdir.mkdir(parents=True, exist_ok=True)
    A = sc.read_h5ad(CP1 / "matrices" / f"{short}_raw.h5ad")
    surv = set(int(x) for x in np.loadtxt(CP1 / "survivors" / f"{short}_{method}.txt", dtype=int).ravel())
    A = A[[b for b in A.obs_names if int(b) in surv]].copy()
    A.var_names = A.var["feature_name"].astype(str).values
    A.var["ensembl"] = A.var["feature_id"].astype(str).values
    A.var_names_make_unique()
    n0 = A.n_obs
    print(f"[{short}/{method}] {n0} survivor cells x {A.n_vars} genes; model={model}", flush=True)

    # NO cell QC here (the QC method already filtered). Gene filter only.
    sc.pp.filter_genes(A, min_cells=QC_MIN_CELLS)
    if not sp.issparse(A.X): A.X = sp.csr_matrix(A.X)
    A.layers["counts"] = A.X.copy()

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        sio.mmwrite(str(td / "c.mtx"), A.layers["counts"].T.tocoo())
        (td / "bc.txt").write_text("\n".join(A.obs_names) + "\n")
        (td / "g.txt").write_text("\n".join(A.var_names) + "\n")
        run_r("cp6_scdblfinder.R", td / "c.mtx", td / "bc.txt", td / "dbl.tsv")
        dcls, dscore = {}, {}
        with open(td / "dbl.tsv") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                dcls[row["barcode"]] = row["scDblFinder_class"]; dscore[row["barcode"]] = float(row["scDblFinder_score"])
        A.obs["doublet_class"] = [dcls.get(b, "NA") for b in A.obs_names]
        A.obs["doublet_score"] = [dscore.get(b, np.nan) for b in A.obs_names]
        run_r("cp6_deviance_hvg.R", td / "c.mtx", td / "g.txt", N_HVG, td / "hvg.txt")
        hvg = set((td / "hvg.txt").read_text().split())
    A.var["highly_variable"] = A.var_names.isin(hvg)

    sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    sc.pp.pca(A, n_comps=N_PCS, use_highly_variable=True, random_state=0)
    sc.pp.neighbors(A, n_neighbors=K, random_state=0)
    for name, res in LEIDEN.items():
        sc.tl.leiden(A, resolution=res, key_added=f"leiden_{name}", random_state=0,
                     flavor="igraph", n_iterations=2, directed=False)
    ncl = {k: int(A.obs[f"leiden_{k}"].nunique()) for k in LEIDEN}

    sc.tl.rank_genes_groups(A, "leiden_r1.0", method="wilcoxon", n_genes=50)
    mk = A.uns["rank_genes_groups"]; mrows = []
    for cl in mk["names"].dtype.names:
        for rank, g in enumerate(mk["names"][cl]):
            mrows.append((cl, rank, g, float(mk["scores"][cl][rank])))

    import celltypist
    pred = celltypist.annotate(A, model=model, majority_voting=True)
    A.obs["celltypist_label"] = pred.predicted_labels["majority_voting"].values

    A.write_h5ad(CP1 / "per_dataset_pipeline_outputs" / f"{short}_{method}.h5ad")
    obs_cols = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "doublet_class",
                "doublet_score", "celltypist_label"] + [f"leiden_{k}" for k in LEIDEN]
    A.obs[obs_cols].to_csv(outdir / f"{method}.obs.tsv", sep="\t")
    pd.DataFrame(mrows, columns=["cluster", "rank", "gene", "score"]).to_csv(outdir / f"{method}.markers.tsv", sep="\t", index=False)
    summary = dict(dataset=short, method=method, n_survivors=n0, n_postdoublet_genes=int(A.n_vars),
                   n_doublets=int((A.obs.doublet_class == "doublet").sum()), n_clusters=ncl,
                   n_celltypist_labels=int(A.obs["celltypist_label"].nunique()), celltypist_model=model,
                   runtime_s=round(time.time() - t0, 1))
    (outdir / f"{method}.summary.json").write_text(json.dumps(summary, indent=2))
    # NOTE: repro.lock registration is done SERIALLY by the orchestrator after the
    # parallel pipeline batch (concurrent register_output would race on the lock).
    print(f"[{short}/{method}] DONE {summary['runtime_s']}s clusters={ncl} labels={summary['n_celltypist_labels']}", flush=True)


if __name__ == "__main__":
    main()
