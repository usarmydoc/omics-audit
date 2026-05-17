#!/usr/bin/env python3
"""Audit 1 main CP4 — E2: database choice effect.

Tool held constant = fgsea. Database varies across 5 MSigDB collections.
All on the 125 DEG TSVs from CP1.

Caveat (per CP4 user guidance, 2026-05-16):
  E2 holds tool=fgsea. Findings are GSEA-paradigm-specific. ORA tools
  may show different database sensitivity. If findings are strong,
  consider an ORA-companion run before drafting rules (decide at CP4
  close, not now).

Output per run: e2/runs/<input_id>__<database_tag>.tsv with same
unified schema as E1.

Databases (all from cached MSigDB v25.1 via msigdbr 26.1.0):
  hallmark         50 pathways
  c2_kegg         186 pathways
  c2_reactome   1839 / 1817 pathways (Hs / Mm)
  c2_wikipathways 925 / 922
  c5_go_bp       7538 / 7402
"""
from __future__ import annotations
import logging
import subprocess
import sys
import time
from pathlib import Path
import pandas as pd

ROOT = Path("/mnt/nvme1/omics-audit/phase2/audit1_main")
INVENTORY = ROOT / "datasets/audit1_main_inputs.tsv"
DB_DIR = ROOT / "databases"
OUT_DIR = ROOT / "e2/runs"
LOG_DIR = ROOT / "e2/logs"
SUMMARY = ROOT / "e2/e2_run_summary.tsv"
FGSEA_HELPER = ROOT / "scripts" / "_e2_fgsea_dbvariant.R"

DATABASES = [
    ("hallmark",        "msigdb_hallmark_{org}.tsv"),
    ("c2_kegg",         "msigdb_c2_kegg_{org}.tsv"),
    ("c2_reactome",     "msigdb_c2_reactome_{org}.tsv"),
    ("c2_wikipathways", "msigdb_c2_wikipathways_{org}.tsv"),
    ("c5_go_bp",        "msigdb_c5_go_bp_{org}.tsv"),
]


def setup_logging():
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(LOG_DIR / "e2_run.log", mode="a"),
        ],
    )


def detect_organism(input_row: pd.Series) -> str:
    if input_row["gene_id_format"] in ("symbol_mouse", "ensembl_mouse"):
        return "mouse"
    return "human"


def detect_id_system(input_row: pd.Series, deg: pd.DataFrame) -> str:
    """Majority-vote re-detection from file (CP1 inventory hint is unreliable)."""
    sample = deg["gene_id"].dropna().astype(str).head(200).tolist()
    if not sample:
        return "symbol"
    n_ens = sum(1 for g in sample if g.startswith("ENS"))
    return "ensembl" if n_ens > len(sample) * 0.5 else "symbol"


def write_fgsea_helper():
    """Same fgsea logic as E1 but database_tag included in output for E2."""
    script = r'''#!/usr/bin/env Rscript
# Args: deg_tsv db_tsv id_system out_tsv input_id database_tag
suppressPackageStartupMessages({library(fgsea); library(dplyr)})
args <- commandArgs(trailingOnly = TRUE)
deg <- read.table(args[1], sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, comment.char="")
db  <- read.table(args[2], sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, comment.char="")
id_system <- args[3]; out_path <- args[4]; input_id <- args[5]; database_tag <- args[6]
gene_col <- switch(id_system, "symbol"="gene_symbol", "ensembl"="gene_ensembl", "entrez"="gene_entrez")
if (id_system == "ensembl") deg$gene_id <- sub("\\.\\d+$", "", as.character(deg$gene_id))
deg$rank_stat <- sign(deg$log2FC) * -log10(pmax(deg$pvalue, 1e-300))
deg <- deg[!is.na(deg$rank_stat) & !is.na(deg$gene_id), ]
deg <- deg[order(deg$rank_stat, decreasing = TRUE), ]
ranks <- deg$rank_stat
names(ranks) <- deg$gene_id
ranks <- ranks[!duplicated(names(ranks))]
pathways <- split(db[[gene_col]], db$pathway_name)
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
    database       = database_tag,
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
write.table(out, out_path, sep="\t", quote=FALSE, row.names=FALSE)
'''
    FGSEA_HELPER.write_text(script)
    FGSEA_HELPER.chmod(0o755)


def run_one(deg_path: str, db_path: Path, id_system: str,
            out_path: Path, input_id: str, database_tag: str) -> int:
    proc = subprocess.run(
        ["Rscript", "--vanilla", str(FGSEA_HELPER),
         deg_path, str(db_path), id_system, str(out_path), input_id, database_tag],
        capture_output=True, text=True, timeout=900,
    )
    if proc.returncode != 0:
        logging.warning("fgsea failed for %s × %s: rc=%d, stderr_tail=%s",
                        input_id, database_tag, proc.returncode, proc.stderr[-400:])
        return 0
    if not out_path.exists():
        return 0
    return sum(1 for _ in out_path.read_text().splitlines()) - 1


def main():
    setup_logging()
    log = logging.getLogger("e2")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    write_fgsea_helper()

    inv = pd.read_csv(INVENTORY, sep="\t")
    log.info("Loaded %d inputs from inventory", len(inv))

    summary_rows = []
    t_start = time.time()
    for idx, row in inv.iterrows():
        input_id = row["input_id"]
        deg_path = row["file_path"]
        organism = detect_organism(row)
        org_tag = "human" if organism == "human" else "mouse"
        deg = pd.read_csv(deg_path, sep="\t")
        id_system = detect_id_system(row, deg)

        for db_tag, db_name_pat in DATABASES:
            db_path = DB_DIR / db_name_pat.format(org=org_tag)
            out_path = OUT_DIR / f"{input_id}__{db_tag}.tsv"
            if out_path.exists() and out_path.stat().st_size > 100:
                continue
            t0 = time.time()
            try:
                n = run_one(deg_path, db_path, id_system, out_path, input_id, db_tag)
                status = "ok"
            except Exception as e:
                log.exception("run failed: %s × %s: %s", input_id, db_tag, e)
                n = 0
                status = "fail"
            elapsed = time.time() - t0
            summary_rows.append({
                "input_id": input_id,
                "database": db_tag,
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
    log.info("DONE — %d (input × database) runs in %.1f min",
             len(summary_rows), (time.time() - t_start) / 60)

    sys.path.insert(0, "/mnt/nvme1/omics-audit/phase2/scripts")
    from dge_native import register_output
    register_output(SUMMARY, kind="audit1_main_e2_run_summary",
                    n_runs=len(summary_rows))


if __name__ == "__main__":
    main()
