#!/usr/bin/env python3
"""Phase 2a A6 — NEBULA + muscat orchestrator for 21 P4 datasets.

For each P4-eligible dataset: stream from Census (matching Phase 1's P4
selection criteria), write count matrix + meta + genes as 10x-like files,
call the A6 R script. Aggregates per-tool per-dataset DEG TSVs.

This is the heavy overnight item. ~12-24hr unattended compute estimate.
"""
from __future__ import annotations

import gc
import gzip
import logging
import subprocess
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import io as sio
from scipy.sparse import csr_matrix

warnings.filterwarnings("ignore")

PHASE2_ROOT = Path("/mnt/nvme1/omics-audit/phase2")
P4_ROOT = Path("/mnt/nvme1/omics-audit/scrnaseq_audit")
sys.path.insert(0, str(P4_ROOT))
from projects.p4_pseudobulk.run_p4 import (  # noqa: E402
    CENSUS_VERSION, MAX_CELLS_PER_DATASET, MIN_CELLS_PER_GROUP,
    MIN_DONORS_PER_GROUP, MIN_GENES, TOP_N_DEGS, RNG_SEED,
    discover_eligible_datasets, fetch_expression_for_de,
)

sys.path.insert(0, str(PHASE2_ROOT / "scripts"))
from dge_native import register_output  # noqa: E402

OUT_BASE = Path("/mnt/nvme1/omics-audit/phase2a/a6_mixed_model")
LOG_FILE = Path("/mnt/nvme1/omics-audit/phase2a/logs/a6_nebula_muscat.log")
R_SCRIPT = PHASE2_ROOT / "scripts" / "a6_nebula_muscat.R"


def setup_logging():
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(LOG_FILE, mode="a")],
    )


def write_inputs_for_R(X, obs_df, gene_names, scratch: Path):
    """Write counts.mtx.gz + meta.tsv + genes.tsv for R consumption."""
    scratch.mkdir(parents=True, exist_ok=True)
    counts_path = scratch / "counts.mtx.gz"
    meta_path = scratch / "meta.tsv"
    genes_path = scratch / "genes.tsv"
    # X is cells × genes (scanpy convention); R wants genes × cells
    if hasattr(X, "tocsc"):
        X_T = X.T.tocoo()
    else:
        X_T = csr_matrix(np.asarray(X).T).tocoo()
    with gzip.open(counts_path, "wb") as gz:
        sio.mmwrite(gz, X_T, field="integer")
    meta = pd.DataFrame({
        "cell_id": [f"cell_{i}" for i in range(len(obs_df))],
        "donor_id": obs_df["donor_id"].astype(str).values,
        "group": obs_df["group"].astype(str).values,
    })
    meta.to_csv(meta_path, sep="\t", index=False)
    pd.Series(gene_names).to_csv(genes_path, sep="\t", index=False, header=False)
    return counts_path, meta_path, genes_path


def main():
    setup_logging()
    log = logging.getLogger("a6")
    OUT_BASE.mkdir(parents=True, exist_ok=True)
    import cellxgene_census

    log.info("Opening Census %s", CENSUS_VERSION)
    census = cellxgene_census.open_soma(census_version=CENSUS_VERSION)
    summary = []
    t_start = time.time()
    for organism in ["Homo sapiens", "Mus musculus"]:
        log.info("=== %s ===", organism)
        eligible = discover_eligible_datasets(census, organism)
        # Sort smallest-first so outputs are checkable early in the run.
        eligible = eligible.assign(
            _cost_proxy=eligible["n_cells_disease"].fillna(0) + eligible["n_cells_normal"].fillna(0)
        ).sort_values("_cost_proxy")
        log.info("eligible: %d (smallest-first by total cells)", len(eligible))
        for _, row in eligible.iterrows():
            ds_id = row.dataset_id
            if pd.isna(ds_id):
                continue
            short = ds_id[:8]
            ds_out = OUT_BASE / short
            # Skip if any output TSV exists (resume support)
            if list(ds_out.glob(f"{short}__*.tsv")):
                log.info("[%s] outputs exist, skipping", short)
                continue
            log.info("[%s] %s/%s, donors=%d+%d, cell_type=%s",
                     short, row.tissue, row.disease, row.n_donors_disease,
                     row.n_donors_normal, row.best_cell_type)
            try:
                fetched = fetch_expression_for_de(census, ds_id, organism, row.best_cell_type)
            except Exception as e:
                log.warning("[%s] fetch failed: %s", short, e)
                continue
            if fetched is None:
                continue
            X, obs_df, gene_names = fetched
            log.info("[%s]   %d cells, %d genes", short, len(obs_df), len(gene_names))

            ds_out.mkdir(parents=True, exist_ok=True)
            counts_p, meta_p, genes_p = write_inputs_for_R(X, obs_df, gene_names, ds_out / "_scratch")

            t0 = time.time()
            proc = subprocess.run(
                ["Rscript", "--vanilla", str(R_SCRIPT),
                 str(counts_p), str(meta_p), str(genes_p), str(ds_out), short],
                capture_output=True, text=True, timeout=14400,  # 4hr per dataset cap
            )
            elapsed = time.time() - t0
            # Always persist R stdout+stderr so silent per-tool failures are recoverable.
            (ds_out / f"{short}_R.log").write_text(
                f"=== STDOUT ===\n{proc.stdout}\n\n=== STDERR ===\n{proc.stderr}\n"
            )
            if proc.returncode != 0:
                log.error("[%s] R failed (rc=%d): %s", short, proc.returncode, proc.stderr[-1500:])
                summary.append({"dataset_short": short, "status": f"r_fail (rc={proc.returncode})",
                                "seconds": round(elapsed, 1)})
            else:
                # Hash-register written TSVs
                written = []
                for t in sorted(ds_out.glob(f"{short}__*.tsv")):
                    tool = t.stem.split("__", 1)[1]
                    df = pd.read_csv(t, sep="\t")
                    n_sig = int((df["padj"].fillna(1) < 0.05).sum())
                    register_output(t, kind="phase2a_a6_mixed_model",
                                    dataset_short=short, organism=organism,
                                    tool=tool, n_genes=int(len(df)), n_sig=n_sig)
                    written.append({"tool": tool, "n_genes": len(df), "n_sig": n_sig})
                log.info("[%s] done in %.1fs, %d tools written", short, elapsed, len(written))
                for w in written:
                    log.info("[%s]   %s: %d genes, %d sig", short, w["tool"],
                             w["n_genes"], w["n_sig"])
                summary.append({"dataset_short": short, "organism": organism,
                                "status": "ok", "seconds": round(elapsed, 1),
                                "n_tools_written": len(written)})

            # Clean scratch (mtx file is large)
            for f in (counts_p, meta_p, genes_p):
                f.unlink(missing_ok=True)
            (ds_out / "_scratch").rmdir()
            del X, obs_df, gene_names
            gc.collect()

    census.close()
    pd.DataFrame(summary).to_csv(OUT_BASE / "a6_summary.tsv", sep="\t", index=False)
    log.info("\nDONE in %.1f min", (time.time() - t_start) / 60)


if __name__ == "__main__":
    main()
