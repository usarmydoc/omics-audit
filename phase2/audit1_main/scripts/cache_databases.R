#!/usr/bin/env Rscript
# Audit 1 main CP2 — cache pathway databases to local TSVs.
#
# Databases (Hs + Mm where applicable):
#   - MSigDB Hallmark (H)
#   - MSigDB C2 curated (KEGG, Reactome, WikiPathways subsets)
#   - MSigDB C5 GO BP
#   - Reactome (via msigdbr C2:CP:REACTOME or full via reactome.db)
#   - KEGG (via msigdbr C2:CP:KEGG_LEGACY — note KEGG licensing)
#   - WikiPathways (via msigdbr C2:CP:WIKIPATHWAYS or rWikiPathways for fresh GMT)
#
# Output: one TSV per (database × organism) with columns:
#   pathway_id, pathway_name, organism, source, gene_symbol, gene_ensembl
#
# All TSVs hash-registered in phase2/repro.lock.

suppressPackageStartupMessages({
  library(msigdbr)
  library(dplyr)
})

OUT_DIR <- "/mnt/nvme1/omics-audit/phase2/audit1_main/databases"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("msigdbr version:", as.character(packageVersion("msigdbr")), "\n")
cat("Available species:", paste(head(msigdbr_species()$species_name, 5), collapse = ", "), "...\n\n")

# msigdbr v26 uses gs_collection / gs_subcollection
collections <- list(
  hallmark = list(collection = "H",  subcollection = NULL),
  c2_kegg  = list(collection = "C2", subcollection = "CP:KEGG_LEGACY"),
  c2_reactome = list(collection = "C2", subcollection = "CP:REACTOME"),
  c2_wikipathways = list(collection = "C2", subcollection = "CP:WIKIPATHWAYS"),
  c5_go_bp = list(collection = "C5", subcollection = "GO:BP")
)
organisms <- list(human = "Homo sapiens", mouse = "Mus musculus")

cache_msigdbr <- function() {
  for (db_key in names(collections)) {
    cfg <- collections[[db_key]]
    for (org_key in names(organisms)) {
      org <- organisms[[org_key]]
      out_file <- file.path(OUT_DIR, sprintf("msigdb_%s_%s.tsv", db_key, org_key))
      cat(sprintf("  fetching %s × %s ...\n", db_key, org_key))
      args <- list(species = org, collection = cfg$collection)
      if (!is.null(cfg$subcollection)) args$subcollection <- cfg$subcollection
      df <- do.call(msigdbr, args)
      out <- df %>% transmute(
        pathway_id     = gs_id,
        pathway_name   = gs_name,
        collection     = gs_collection,
        subcollection  = gs_subcollection,
        organism       = !!org,
        gene_symbol    = gene_symbol,
        gene_ensembl   = ensembl_gene,
        gene_entrez    = ncbi_gene      # msigdbr v26 renamed entrez_gene -> ncbi_gene
      )
      write.table(out, file = out_file, sep = "\t",
                  quote = FALSE, row.names = FALSE)
      cat(sprintf("    -> %s (%d gene-pathway rows, %d unique pathways)\n",
                  basename(out_file), nrow(out), length(unique(out$pathway_id))))
    }
  }
}

cache_wikipathways_gmt <- function() {
  # rWikiPathways downloads the latest GMT directly
  suppressPackageStartupMessages(library(rWikiPathways))
  for (org_key in names(organisms)) {
    org <- organisms[[org_key]]
    out_gmt <- file.path(OUT_DIR, sprintf("wikipathways_%s_latest.gmt", org_key))
    cat(sprintf("  fetching WikiPathways GMT for %s ...\n", org))
    tryCatch({
      gmt <- downloadPathwayArchive(organism = org, format = "gmt", destpath = OUT_DIR)
      cat(sprintf("    -> %s\n", gmt))
    }, error = function(e) {
      cat(sprintf("    [warn] WikiPathways GMT fetch failed for %s: %s\n", org, conditionMessage(e)))
      cat("    (Falls back on msigdbr C2:CP:WIKIPATHWAYS subset for the audit.)\n")
    })
  }
}

cat("=== Caching MSigDB collections (Hallmark, C2 subsets, C5 GO BP) ===\n")
cache_msigdbr()

cat("\n=== Caching fresh WikiPathways GMTs ===\n")
cache_wikipathways_gmt()

cat("\n=== Caching reactome.db full table ===\n")
tryCatch({
  suppressPackageStartupMessages(library(reactome.db))
  cat("  reactome.db version:", as.character(packageVersion("reactome.db")), "\n")
  # Extract pathway names + their gene memberships
  pathways <- as.list(reactomePATHID2EXTID)
  pathway_names <- as.list(reactomePATHID2NAME)
  ent <- do.call(rbind, lapply(names(pathways), function(pid) {
    if (length(pathways[[pid]]) == 0) return(NULL)
    data.frame(
      pathway_id = pid,
      pathway_name = if (!is.null(pathway_names[[pid]])) pathway_names[[pid]][1] else NA,
      gene_entrez = pathways[[pid]],
      stringsAsFactors = FALSE
    )
  }))
  # Split by organism via name prefix ("Homo sapiens: ..." / "Mus musculus: ...")
  for (org_key in names(organisms)) {
    org <- organisms[[org_key]]
    sub <- ent[grepl(paste0("^", org, ":"), ent$pathway_name), ]
    if (nrow(sub) == 0) {
      cat(sprintf("  [info] no Reactome rows for %s after filter — skipping\n", org))
      next
    }
    out_file <- file.path(OUT_DIR, sprintf("reactome_db_%s.tsv", org_key))
    write.table(sub, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("    -> %s (%d rows, %d pathways)\n",
                basename(out_file), nrow(sub), length(unique(sub$pathway_id))))
  }
}, error = function(e) {
  cat(sprintf("  [skip] reactome.db not installed or fetch failed: %s\n", conditionMessage(e)))
  cat("  (Reactome coverage via msigdbr C2:CP:REACTOME suffices for E2.)\n")
})

cat("\n=== Capturing R environment for env_id stamp ===\n")
env_file <- file.path(OUT_DIR, "R_session_info.txt")
sink(env_file); print(sessionInfo()); sink()
cat(sprintf("R sessionInfo written to %s\n", env_file))

cat("\n=== KEGG licensing note ===\n")
cat("KEGG REST API is free for academic/non-commercial use.\n")
cat("KEGG_LEGACY collection in MSigDB C2 is licensed for research use under\n")
cat("Bioconductor's data redistribution terms. For this audit:\n")
cat("  - msigdbr C2:CP:KEGG_LEGACY: research-use OK, cached as TSV\n")
cat("  - Direct KEGGREST API calls: rate-limited (3 req/sec), used only\n")
cat("    if msigdbr coverage is insufficient. Document in INTEGRATION_NOTES.md.\n")

cat("\n=== Cache complete ===\n")
