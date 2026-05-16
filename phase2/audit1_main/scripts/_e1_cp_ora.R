#!/usr/bin/env Rscript
# Args: deg_tsv hallmark_tsv id_system out_tsv input_id
suppressPackageStartupMessages({library(clusterProfiler); library(dplyr)})
args <- commandArgs(trailingOnly = TRUE)
deg <- read.table(args[1], sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, comment.char="")
hm  <- read.table(args[2], sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, comment.char="")
id_system <- args[3]; out_path <- args[4]; input_id <- args[5]
gene_col <- switch(id_system, "symbol"="gene_symbol", "ensembl"="gene_ensembl", "entrez"="gene_entrez")
if (id_system == "ensembl") deg$gene_id <- sub("\\.\\d+$", "", as.character(deg$gene_id))
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
write.table(out, out_path, sep="\t", quote=FALSE, row.names=FALSE)
