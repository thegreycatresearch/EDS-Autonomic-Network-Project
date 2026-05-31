# ============================================================
# 04_network_analysis.R
# Multi-cohort WGCNA + consensus network biology
# EDS + dysautonomia systems biology framework
# ============================================================

# ============================
# 0. PACKAGES
# ============================

suppressPackageStartupMessages({

  library(WGCNA)
  library(dplyr)
  library(igraph)

})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# ============================
# 1. PATHS
# ============================

PROCESSED_DIR <- "data/processed"
RESULTS_DIR <- "results"
NET_DIR <- file.path(RESULTS_DIR, "networks")

dir.create(NET_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================
# 2. LOAD DATA
# ============================

files <- list.files(PROCESSED_DIR,
                    pattern = "_processed.rds",
                    full.names = TRUE)

datasets <- lapply(files, readRDS)
names(datasets) <- gsub("_processed.rds", "", basename(files))

# ============================
# 3. LOAD DEG CORE (FROM STEP 3)
# ============================

eds_core <- read.csv(file.path(RESULTS_DIR, "DEG/EDS_core_genes.csv"))[,1]
pots_core <- read.csv(file.path(RESULTS_DIR, "DEG/POTS_core_genes.csv"))[,1]
shared_core <- read.csv(file.path(RESULTS_DIR, "DEG/shared_core_genes.csv"))[,1]

all_core <- unique(c(eds_core, pots_core, shared_core))

# ============================
# 4. FILTER EXPRESSION (IMPORTANT)
# ============================

filter_to_core <- function(expr, genes) {

  common <- intersect(rownames(expr), genes)

  expr[common, , drop = FALSE]

}

# ============================
# 5. WGCNA PER DATASET
# ============================

run_wgcna <- function(expr, name) {

  datExpr <- t(expr)

  # QC samples
  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

  # soft threshold
  powers <- c(1:10, seq(12, 30, 2))

  sft <- pickSoftThreshold(datExpr,
                           powerVector = powers,
                           verbose = 5)

  softPower <- ifelse(is.na(sft$powerEstimate), 6, sft$powerEstimate)

  # adjacency
  adj <- adjacency(datExpr, power = softPower)

  TOM <- TOMsimilarity(adj)
  dissTOM <- 1 - TOM

  geneTree <- hclust(as.dist(dissTOM), method = "average")

  dynamicMods <- cutreeDynamic(
    dendro = geneTree,
    distM = dissTOM,
    deepSplit = 3,
    pamRespectsDendro = FALSE,
    minClusterSize = 30
  )

  colors <- labels2colors(dynamicMods)

  MEs <- moduleEigengenes(datExpr, colors)$eigengenes

  list(
    colors = colors,
    MEs = MEs,
    TOM = TOM,
    genes = colnames(datExpr),
    name = name
  )

}

# ============================
# 6. RUN PER DATASET
# ============================

wgcna_results <- list()

for (name in names(datasets)) {

  cat("WGCNA:", name, "\n")

  expr <- datasets[[name]]$expression

  expr <- filter_to_core(expr, all_core)

  if (nrow(expr) < 50) next

  res <- run_wgcna(expr, name)

  wgcna_results[[name]] <- res

}

# ============================
# 7. CONSENSUS MODULE ANALYSIS
# ============================

module_overlap <- function(res_list) {

  all_modules <- lapply(res_list, function(x) x$colors)

  gene_sets <- lapply(res_list, function(x) x$genes)

  common_genes <- Reduce(intersect, gene_sets)

  consensus <- data.frame(gene = common_genes)

  for (i in seq_along(res_list)) {

    res <- res_list[[i]]

    idx <- match(common_genes, res$genes)

    consensus[[names(res_list)[i]]] <- res$colors[idx]

  }

  consensus

}

consensus_modules <- module_overlap(wgcna_results)

# ============================
# 8. HUB GENE DETECTION (STRICT VERSION)
# ============================

detect_hubs <- function(expr, MEs, genes) {

  datExpr <- t(expr)

  kME <- cor(datExpr, MEs, use = "p")

  hub_scores <- apply(kME, 2, function(x) abs(x))

  hubs <- apply(hub_scores, 2, function(x) {

    names(sort(x, decreasing = TRUE))[1:20]

  })

  hubs

}

hub_results <- list()

for (name in names(datasets)) {

  expr <- datasets[[name]]$expression

  expr <- filter_to_core(expr, all_core)

  if (nrow(expr) < 50) next

  hub_results[[name]] <- detect_hubs(
    expr,
    wgcna_results[[name]]$MEs,
    rownames(expr)
  )

}

# ============================
# 9. MODULE–DEG INTEGRATION
# ============================

module_deg_overlap <- data.frame()

for (name in names(wgcna_results)) {

  colors <- wgcna_results[[name]]$colors
  genes <- wgcna_results[[name]]$genes

  eds_hits <- intersect(eds_core, genes)
  pots_hits <- intersect(pots_core, genes)

  module_deg_overlap <- rbind(module_deg_overlap,
                              data.frame(
                                dataset = name,
                                eds_overlap = length(eds_hits),
                                pots_overlap = length(pots_hits)
                              ))

}

# ============================
# 10. EXPORT NETWORK
# ============================

for (name in names(wgcna_results)) {

  res <- wgcna_results[[name]]

  net <- data.frame(
    gene = res$genes,
    module = res$colors
  )

  write.csv(net,
            file.path(NET_DIR, paste0(name, "_modules.csv")),
            row.names = FALSE)

}

# ============================
# 11. SAVE CONSENSUS
# ============================

write.csv(consensus_modules,
          file.path(NET_DIR, "consensus_modules.csv"),
          row.names = FALSE)

write.csv(module_deg_overlap,
          file.path(NET_DIR, "module_deg_overlap.csv"),
          row.names = FALSE)

message("04_network_analysis completed")
