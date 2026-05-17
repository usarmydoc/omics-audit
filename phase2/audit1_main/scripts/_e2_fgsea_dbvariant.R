#!/usr/bin/env Rscript
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
