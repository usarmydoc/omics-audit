#!/usr/bin/env python
"""CP3 Deliverable C — downstream pipeline on a corrected-matrix condition.

Mirrors Audit 3 CP6 EXACTLY (same QC/doublet/HVG/norm/cluster/markers/annotation
and parameters), but loads input from a CP3 condition's corrected matrix instead
of a counting tool's native matrix. Reuses cp6_scdblfinder.R + cp6_deviance_hvg.R.

Usage: cp3_pipeline.py <dataset_id> <condition> <matrix_dir>
  matrix_dir contains: matrix.mtx (cells x genes), barcodes.txt, ensembl.txt
Outputs to cp3/<intestine|pbmc>/<condition>.{h5ad,obs.tsv,markers.tsv,summary.json}
"""
import sys, time, json, tempfile, subprocess, csv
from pathlib import Path
import numpy as np, scipy.io as sio, scipy.sparse as sp
import anndata as ad, scanpy as sc, pandas as pd

AUDIT = Path("/mnt/nvme1/omics-audit/phase2/audit3_counting")
CP6_SCRIPTS = AUDIT / "scripts"
CP3 = Path("/mnt/nvme1/omics-audit/phase2/audit_ambient_correction/cp3")
sys.path.insert(0, str(Path("/mnt/nvme1/omics-audit/phase2/scripts")))
from dge_native import register_output
sc.settings.verbosity = 1
QC_MIN_GENES, QC_MIN_CELLS, QC_MAX_PCT_MITO = 200, 3, 20.0
N_HVG, N_PCS, K_NEIGHBORS = 2000, 50, 15
LEIDEN_RES = {"r1.0": 1.0, "r0.5": 0.5, "r1.5": 1.5}
MODEL = {"gse288156_mouse_intestine_scrna": "Adult_Mouse_Gut.pkl", "10x_pbmc_5k_v3.1": "Immune_All_Low.pkl"}


def ensembl_to_symbol(ds):
    f = AUDIT / "processed" / ds / "star_default" / "Solo.out" / "Gene" / "raw" / "features.tsv"
    m = {}
    for line in open(f):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2: m[p[0].split(".")[0]] = p[1]
    return m


def run_r(script, *a):
    r = subprocess.run(["Rscript", str(CP6_SCRIPTS / script), *map(str, a)], capture_output=True, text=True)
    if r.returncode != 0: raise RuntimeError(f"{script} failed:\n{r.stdout}\n{r.stderr}")
    return r.stdout


def main():
    ds, cond, mdir = sys.argv[1], sys.argv[2], Path(sys.argv[3])
    t0 = time.time()
    sub = "intestine" if "intestine" in ds else "pbmc"
    outdir = CP3 / sub; outdir.mkdir(parents=True, exist_ok=True)
    model = MODEL[ds]
    print(f"[{ds}/{cond}] start model={model}", flush=True)

    # load condition corrected matrix
    X = sio.mmread(str(mdir / "matrix.mtx")).tocsr()                 # cells x genes
    cells = [l.strip() for l in open(mdir / "barcodes.txt")]
    ensembl = [e.split(".")[0] for e in (l.strip() for l in open(mdir / "ensembl.txt"))]
    sym = ensembl_to_symbol(ds)
    A = ad.AnnData(X=X.astype(np.float32))
    A.obs_names = cells; A.var_names = [sym.get(e, e) for e in ensembl]
    A.var["ensembl"] = ensembl; A.var_names_make_unique()
    n0 = A.n_obs
    print(f"  loaded {A.n_obs} cells x {A.n_vars} genes", flush=True)

    # QC (identical to cp6)
    A.var["mt"] = A.var_names.str.upper().str.startswith(("MT-", "MT."))
    sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], inplace=True, percent_top=None)
    sc.pp.filter_cells(A, min_genes=QC_MIN_GENES)
    sc.pp.filter_genes(A, min_cells=QC_MIN_CELLS)
    A = A[A.obs["pct_counts_mt"] < QC_MAX_PCT_MITO].copy()
    n_qc = A.n_obs
    print(f"  post-QC {n0}->{n_qc}", flush=True)
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
    sc.pp.neighbors(A, n_neighbors=K_NEIGHBORS, random_state=0)
    for name, res in LEIDEN_RES.items():
        sc.tl.leiden(A, resolution=res, key_added=f"leiden_{name}", random_state=0,
                     flavor="igraph", n_iterations=2, directed=False)
    ncl = {k: int(A.obs[f"leiden_{k}"].nunique()) for k in LEIDEN_RES}
    print(f"  clusters {ncl}", flush=True)

    sc.tl.rank_genes_groups(A, "leiden_r1.0", method="wilcoxon", n_genes=50)
    mk = A.uns["rank_genes_groups"]; mrows = []
    for cl in mk["names"].dtype.names:
        for rank, g in enumerate(mk["names"][cl]):
            mrows.append((cl, rank, g, float(mk["scores"][cl][rank])))

    import celltypist
    pred = celltypist.annotate(A, model=model, majority_voting=True)
    A.obs["celltypist_label"] = pred.predicted_labels["majority_voting"].values

    A.write_h5ad(outdir / f"{cond}.h5ad")
    obs_cols = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "doublet_class",
                "doublet_score", "celltypist_label"] + [f"leiden_{k}" for k in LEIDEN_RES]
    A.obs[obs_cols].to_csv(outdir / f"{cond}.obs.tsv", sep="\t")
    pd.DataFrame(mrows, columns=["cluster", "rank", "gene", "score"]).to_csv(outdir / f"{cond}.markers.tsv", sep="\t", index=False)
    summary = dict(dataset=ds, condition=cond, n_cells_input=n0, n_cells_postqc=n_qc,
                   n_genes=int(A.n_vars), n_doublets=int((A.obs.doublet_class == "doublet").sum()),
                   n_clusters=ncl, n_celltypist_labels=int(A.obs["celltypist_label"].nunique()),
                   celltypist_model=model, runtime_s=round(time.time() - t0, 1))
    (outdir / f"{cond}.summary.json").write_text(json.dumps(summary, indent=2))
    for p, k in [(outdir / f"{cond}.obs.tsv", f"cp3_{cond}_obs"), (outdir / f"{cond}.markers.tsv", f"cp3_{cond}_markers"),
                 (outdir / f"{cond}.summary.json", f"cp3_{cond}_summary")]:
        register_output(p, kind=k, dataset=ds)
    print(f"[{ds}/{cond}] DONE {summary['runtime_s']}s postQC={n_qc} clusters={ncl} labels={summary['n_celltypist_labels']}", flush=True)


if __name__ == "__main__":
    main()
