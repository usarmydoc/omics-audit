#!/usr/bin/env python3
"""Audit 1 main CP3 — E1: tool agreement on identical DEG inputs.

Runs 3 pathway tools (fgsea via Rscript, gseapy.enrich locally,
clusterProfiler ORA via Rscript) against MSigDB Hallmark on all 125
DEG TSVs from CP1. Database is held constant; tool varies.

Per-tool unified TSV output schema:
    pathway_id, pathway_name, database, tool, input_id, comparison_id,
    p_value, padj, NES, log2_odds_ratio, gene_set_size, hit_genes

Outputs:
    audit1_main/e1/runs/<input_id>__<tool>.tsv  (375 files)
    audit1_main/e1/e1_run_summary.tsv
    audit1_main/e1/logs/

Subsequent script (e1_metrics.py) computes pairwise agreement metrics
+ bootstrap CIs across this output.
"""
from __future__ import annotations
import gzip
import json
import logging
import subprocess
import sys
import time
from pathlib import Path
import pandas as pd

ROOT = Path("/mnt/nvme1/omics-audit/phase2/audit1_main")
INVENTORY = ROOT / "datasets/audit1_main_inputs.tsv"
DB_DIR = ROOT / "databases"
OUT_DIR = ROOT / "e1/runs"
LOG_DIR = ROOT / "e1/logs"
SUMMARY = ROOT / "e1/e1_run_summary.tsv"
FGSEA_SCRIPT = ROOT / "scripts" / "_e1_fgsea.R"
CP_ORA_SCRIPT = ROOT / "scripts" / "_e1_cp_ora.R"
RSCRATCH = ROOT / "e1/_rscratch"


def setup_logging():
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(LOG_DIR / "e1_run.log", mode="a"),
        ],
    )


def load_hallmark(organism: str) -> pd.DataFrame:
    org = "human" if organism in ("Hs", "human", "Homo sapiens") else "mouse"
    p = DB_DIR / f"msigdb_hallmark_{org}.tsv"
    return pd.read_csv(p, sep="\t")


def detect_organism(input_row: pd.Series) -> str:
    if input_row["gene_id_format"] in ("symbol_mouse", "ensembl_mouse"):
        return "mouse"
    return "human"  # default — TCGA, GTEx, most Census are human


def detect_id_system(input_row: pd.Series, deg: pd.DataFrame) -> str:
    """Return 'symbol', 'ensembl', or 'entrez' — what the DEG's gene_id column uses.

    Trust majority of actual file contents rather than CP1 inventory hint —
    CP1's `any()`-based heuristic mis-classified files with a few ensembl
    contaminants in an otherwise-symbol corpus.
    """
    sample = deg["gene_id"].dropna().astype(str).head(200).tolist()
    n_ens = sum(1 for g in sample if g.startswith("ENS"))
    n_total = len(sample)
    if n_total == 0:
        return "symbol"  # default fallback
    if n_ens > n_total * 0.5:
        return "ensembl"
    return "symbol"


def run_gseapy(deg: pd.DataFrame, hallmark: pd.DataFrame, id_system: str,
               input_id: str, out_path: Path) -> int:
    """Run gseapy ORA (enrichr-style) against local Hallmark."""
    import gseapy as gp

    # Significant gene list (padj < 0.05)
    import re
    # Strip Ensembl version suffix (.N) so GTEx-style IDs match bare Hallmark Ensembl IDs
    if id_system == "ensembl":
        deg = deg.copy()
        deg["gene_id"] = deg["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    sig = deg[deg["padj"].fillna(1) < 0.05].copy()
    if len(sig) < 10:
        out_path.write_text("pathway_id\tpathway_name\tdatabase\ttool\tinput_id\tcomparison_id\tp_value\tpadj\tNES\tlog2_odds_ratio\tgene_set_size\thit_genes\n")
        return 0
    gene_col_map = {"symbol": "gene_symbol", "ensembl": "gene_ensembl", "entrez": "gene_entrez"}
    gene_col = gene_col_map[id_system]
    gene_list = sig["gene_id"].dropna().astype(str).unique().tolist()
    background = deg["gene_id"].dropna().astype(str).unique().tolist()

    # Build gene_sets dict from cached Hallmark
    gene_sets = {}
    for pwy, sub in hallmark.groupby("pathway_name"):
        genes = sub[gene_col].dropna().astype(str).unique().tolist()
        if len(genes) >= 5:
            gene_sets[pwy] = genes

    try:
        enr = gp.enrich(
            gene_list=gene_list,
            gene_sets=gene_sets,
            background=background,
            outdir=None,
            verbose=False,
        )
        raw = enr.results
        # gseapy may return list (multi-library) or single DataFrame
        if isinstance(raw, list):
            if not raw:
                df = pd.DataFrame()
            else:
                df = pd.concat(raw, ignore_index=True)
        else:
            df = raw.copy() if raw is not None else pd.DataFrame()
    except Exception as e:
        out_path.write_text("pathway_id\tpathway_name\tdatabase\ttool\tinput_id\tcomparison_id\tp_value\tpadj\tNES\tlog2_odds_ratio\tgene_set_size\thit_genes\n")
        return 0

    if df.empty:
        out_path.write_text("pathway_id\tpathway_name\tdatabase\ttool\tinput_id\tcomparison_id\tp_value\tpadj\tNES\tlog2_odds_ratio\tgene_set_size\thit_genes\n")
        return 0

    # Map to unified schema
    out_rows = []
    for _, r in df.iterrows():
        out_rows.append({
            "pathway_id": r.get("Term", ""),
            "pathway_name": r.get("Term", ""),
            "database": "MSigDB_Hallmark",
            "tool": "gseapy_enrichr",
            "input_id": input_id,
            "comparison_id": input_id,
            "p_value": r.get("P-value", float("nan")),
            "padj": r.get("Adjusted P-value", float("nan")),
            "NES": float("nan"),
            "log2_odds_ratio": (
                __import__("math").log2(float(r["Odds Ratio"]))
                if r.get("Odds Ratio") and float(r["Odds Ratio"]) > 0
                else float("nan")
            ),
            "gene_set_size": len(gene_sets.get(r.get("Term", ""), [])),
            "hit_genes": r.get("Genes", ""),
        })
    pd.DataFrame(out_rows).to_csv(out_path, sep="\t", index=False)
    return len(out_rows)


def write_fgsea_helper():
    """One-time write of the fgsea Rscript helper."""
    script = '''#!/usr/bin/env Rscript
# Args: deg_tsv hallmark_tsv id_system out_tsv input_id
suppressPackageStartupMessages({library(fgsea); library(dplyr)})
args <- commandArgs(trailingOnly = TRUE)
deg <- read.table(args[1], sep="\\t", header=TRUE, quote="", stringsAsFactors=FALSE, comment.char="")
hm  <- read.table(args[2], sep="\\t", header=TRUE, quote="", stringsAsFactors=FALSE, comment.char="")
id_system <- args[3]; out_path <- args[4]; input_id <- args[5]
gene_col <- switch(id_system, "symbol"="gene_symbol", "ensembl"="gene_ensembl", "entrez"="gene_entrez")
if (id_system == "ensembl") deg$gene_id <- sub("\\\\.\\\\d+$", "", as.character(deg$gene_id))
deg$rank_stat <- sign(deg$log2FC) * -log10(pmax(deg$pvalue, 1e-300))
deg <- deg[!is.na(deg$rank_stat) & !is.na(deg$gene_id), ]
deg <- deg[order(deg$rank_stat, decreasing = TRUE), ]
ranks <- deg$rank_stat
names(ranks) <- deg$gene_id
ranks <- ranks[!duplicated(names(ranks))]
pathways <- split(hm[[gene_col]], hm$pathway_name)
pathways <- lapply(pathways, function(g) unique(g[!is.na(g) & g != ""]))
pathways <- pathways[sapply(pathways, length) >= 5]
res <- tryCatch(
  fgsea::fgsea(pathways = pathways, stats = ranks, minSize = 5, maxSize = 500),
  error = function(e) NULL
)
if (is.null(res) || nrow(res) == 0) {
  out <- data.frame(pathway_id=character(), pathway_name=character(), database=character(),
                    tool=character(), input_id=character(), comparison_id=character(),
                    p_value=numeric(), padj=numeric(), NES=numeric(),
                    log2_odds_ratio=numeric(), gene_set_size=integer(), hit_genes=character())
} else {
  out <- data.frame(
    pathway_id     = res$pathway,
    pathway_name   = res$pathway,
    database       = "MSigDB_Hallmark",
    tool           = "fgsea",
    input_id       = input_id,
    comparison_id  = input_id,
    p_value        = res$pval,
    padj           = res$padj,
    NES            = res$NES,
    log2_odds_ratio = NA_real_,
    gene_set_size  = res$size,
    hit_genes      = sapply(res$leadingEdge, function(x) paste(x, collapse=";"))
  )
}
write.table(out, out_path, sep="\\t", quote=FALSE, row.names=FALSE)
'''
    FGSEA_SCRIPT.write_text(script)
    FGSEA_SCRIPT.chmod(0o755)


def write_cp_ora_helper():
    """One-time write of the clusterProfiler ORA Rscript helper."""
    script = '''#!/usr/bin/env Rscript
# Args: deg_tsv hallmark_tsv id_system out_tsv input_id
suppressPackageStartupMessages({library(clusterProfiler); library(dplyr)})
args <- commandArgs(trailingOnly = TRUE)
deg <- read.table(args[1], sep="\\t", header=TRUE, quote="", stringsAsFactors=FALSE, comment.char="")
hm  <- read.table(args[2], sep="\\t", header=TRUE, quote="", stringsAsFactors=FALSE, comment.char="")
id_system <- args[3]; out_path <- args[4]; input_id <- args[5]
gene_col <- switch(id_system, "symbol"="gene_symbol", "ensembl"="gene_ensembl", "entrez"="gene_entrez")
if (id_system == "ensembl") deg$gene_id <- sub("\\\\.\\\\d+$", "", as.character(deg$gene_id))
sig <- deg$gene_id[!is.na(deg$padj) & deg$padj < 0.05]
bg  <- unique(deg$gene_id[!is.na(deg$gene_id)])
if (length(sig) < 10) {
  out <- data.frame(pathway_id=character(), pathway_name=character(), database=character(),
                    tool=character(), input_id=character(), comparison_id=character(),
                    p_value=numeric(), padj=numeric(), NES=numeric(),
                    log2_odds_ratio=numeric(), gene_set_size=integer(), hit_genes=character())
} else {
  t2g <- data.frame(term = hm$pathway_name, gene = hm[[gene_col]], stringsAsFactors = FALSE)
  t2g <- t2g[!is.na(t2g$gene) & t2g$gene != "", ]
  res <- tryCatch(
    clusterProfiler::enricher(gene = unique(sig), TERM2GENE = t2g, universe = bg,
                              minGSSize = 5, maxGSSize = 500,
                              pvalueCutoff = 1.0, qvalueCutoff = 1.0),
    error = function(e) NULL
  )
  if (is.null(res) || nrow(res@result) == 0) {
    out <- data.frame(pathway_id=character(), pathway_name=character(), database=character(),
                      tool=character(), input_id=character(), comparison_id=character(),
                      p_value=numeric(), padj=numeric(), NES=numeric(),
                      log2_odds_ratio=numeric(), gene_set_size=integer(), hit_genes=character())
  } else {
    r <- res@result
    out <- data.frame(
      pathway_id     = r$ID,
      pathway_name   = r$Description,
      database       = "MSigDB_Hallmark",
      tool           = "clusterProfiler_ORA",
      input_id       = input_id,
      comparison_id  = input_id,
      p_value        = r$pvalue,
      padj           = r$p.adjust,
      NES            = NA_real_,
      log2_odds_ratio = log2(sapply(strsplit(r$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])) /
                              sapply(strsplit(r$BgRatio,   "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))),
      gene_set_size  = sapply(strsplit(r$BgRatio, "/"), function(x) as.numeric(x[1])),
      hit_genes      = r$geneID
    )
  }
}
write.table(out, out_path, sep="\\t", quote=FALSE, row.names=FALSE)
'''
    CP_ORA_SCRIPT.write_text(script)
    CP_ORA_SCRIPT.chmod(0o755)


def run_r_tool(tool_script: Path, deg_path: str, hm_path: Path,
               id_system: str, out_path: Path, input_id: str) -> int:
    proc = subprocess.run(
        ["Rscript", "--vanilla", str(tool_script), deg_path, str(hm_path),
         id_system, str(out_path), input_id],
        capture_output=True, text=True, timeout=600,
    )
    if proc.returncode != 0:
        logging.warning("R tool failed for %s: rc=%d, stderr_tail=%s",
                        input_id, proc.returncode, proc.stderr[-400:])
        return 0
    if not out_path.exists():
        return 0
    return sum(1 for _ in out_path.read_text().splitlines()) - 1


def main():
    setup_logging()
    log = logging.getLogger("e1")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RSCRATCH.mkdir(parents=True, exist_ok=True)
    write_fgsea_helper()
    write_cp_ora_helper()

    inv = pd.read_csv(INVENTORY, sep="\t")
    log.info("Loaded %d inputs from inventory", len(inv))

    # Cache Hallmark tables by organism (avoid re-reading per input)
    hallmark_cache = {
        "human": (load_hallmark("human"), DB_DIR / "msigdb_hallmark_human.tsv"),
        "mouse": (load_hallmark("mouse"), DB_DIR / "msigdb_hallmark_mouse.tsv"),
    }

    summary_rows = []
    t_start = time.time()

    for idx, row in inv.iterrows():
        input_id = row["input_id"]
        deg_path = row["file_path"]
        organism = detect_organism(row)
        hm_df, hm_path = hallmark_cache[organism]
        deg = pd.read_csv(deg_path, sep="\t")
        id_system = detect_id_system(row, deg)

        for tool in ["fgsea", "gseapy_enrichr", "clusterProfiler_ORA"]:
            out_path = OUT_DIR / f"{input_id}__{tool}.tsv"
            if out_path.exists() and out_path.stat().st_size > 100:
                continue  # resume
            t0 = time.time()
            try:
                if tool == "fgsea":
                    n = run_r_tool(FGSEA_SCRIPT, deg_path, hm_path, id_system, out_path, input_id)
                elif tool == "clusterProfiler_ORA":
                    n = run_r_tool(CP_ORA_SCRIPT, deg_path, hm_path, id_system, out_path, input_id)
                elif tool == "gseapy_enrichr":
                    n = run_gseapy(deg, hm_df, id_system, input_id, out_path)
                status = "ok"
            except Exception as e:
                log.exception("tool %s failed on %s: %s", tool, input_id, e)
                n = 0
                status = "fail"
            elapsed = time.time() - t0
            summary_rows.append({
                "input_id": input_id,
                "tool": tool,
                "organism": organism,
                "id_system": id_system,
                "n_pathways_returned": n,
                "seconds": round(elapsed, 2),
                "status": status,
            })

        if (idx + 1) % 10 == 0:
            log.info("Progress: %d / %d inputs done, elapsed %.1f min",
                     idx + 1, len(inv), (time.time() - t_start) / 60)

    pd.DataFrame(summary_rows).to_csv(SUMMARY, sep="\t", index=False)
    log.info("DONE — %d (input × tool) runs in %.1f min",
             len(summary_rows), (time.time() - t_start) / 60)

    # Register summary in lock
    sys.path.insert(0, "/mnt/nvme1/omics-audit/phase2/scripts")
    from dge_native import register_output
    register_output(SUMMARY, kind="audit1_main_e1_run_summary",
                    n_runs=len(summary_rows))


if __name__ == "__main__":
    main()
