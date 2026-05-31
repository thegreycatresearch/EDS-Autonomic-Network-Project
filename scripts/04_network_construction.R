# ============================================================
# 04_network_analysis.R
# Advanced WGCNA + network biology integration
# EDS + POTS systems biology framework
# ============================================================

# ============================
# 0. PACKAGES
# ============================

suppressPackageStartupMessages({

  library(WGCNA)
  library(igraph)
  library(dplyr)
  library(Matrix)

})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# ============================
# 1. PATHS
# ============================

RESULTS_DIR <- "results"
NETWORK_DIR <- file.path(RESULTS_DIR, "networks")
WGCNA_DIR <- file.path(NETWORK_DIR, "wgcna")

dir.create(NETWORK_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(WGCNA_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================
# 2. LOAD DEG DATA
# ============================

eds_degs <- read.csv(file.path(RESULTS_DIR, "EDS_DEGs.csv"), row.names = 1)
pots_degs <- read.csv(file.path(RESULTS_DIR, "POTS_DEGs.csv"), row.names = 1)

eds_sig <- rownames(subset(eds_degs, adj.P.Val < 0.05))
pots_sig <- rownames(subset(pots_degs, adj.P.Val < 0.05))

all_genes <- unique(c(eds_sig, pots_sig))

# ============================
# 3. LOAD EXPRESSION MATRIX
# ============================

load_expr <- function() {

  files <- list.files("data/processed", full.names = TRUE)

  expr_list <- list()

  for (f in files) {

    if (grepl("_processed.rds", f)) {

      data <- readRDS(f)
      expr_list[[f]] <- data$expression
    }
  }

  common_genes <- Reduce(intersect, lapply(expr_list, rownames))
  expr_list <- lapply(expr_list, function(x) x[common_genes, ])

  merged <- do.call(cbind, expr_list)

  return(merged)
}

expr <- load_expr()

# keep only biologically relevant genes
expr <- expr[intersect(rownames(expr), all_genes), ]

# ============================
# 4. VARIANCE FILTERING (IMPORTANT UPGRADE)
# ============================

gene_var <- apply(expr, 1, var)

top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(5000, length(gene_var))]

expr <- expr[top_genes, ]

# ============================
# 5. TRANSPOSE FOR WGCNA
# ============================

datExpr <- t(expr)

# ============================
# 6. SAMPLE QC
# ============================

gsg <- goodSamplesGenes(datExpr, verbose = 3)

if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# ============================
# 7. SAMPLE CLUSTERING (OUTLIER DETECTION)
# ============================

sampleTree <- hclust(dist(datExpr), method = "average")

# optional: visualize later
saveRDS(sampleTree, file = file.path(WGCNA_DIR, "sample_tree.rds"))

# ============================
# 8. SOFT THRESHOLD SELECTION (IMPROVED)
# ============================

powers <- c(seq(1, 10, 1), seq(12, 30, 2))

sft <- pickSoftThreshold(datExpr,
                         powerVector = powers,
                         verbose = 5)

softPower <- ifelse(is.na(sft$powerEstimate), 6, sft$powerEstimate)

message(paste("Selected soft power:", softPower))

# ============================
# 9. ADJACENCY + TOM
# ============================

adjacency_matrix <- adjacency(datExpr, power = softPower)

TOM <- TOMsimilarity(adjacency_matrix)

dissTOM <- 1 - TOM

geneTree <- hclust(as.dist(dissTOM), method = "average")

# ============================
# 10. DYNAMIC MODULE DETECTION
# ============================

dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = 3,
  pamRespectsDendro = FALSE,
  minClusterSize = 30
)

moduleColors <- labels2colors(dynamicMods)

# ============================
# 11. MODULE MERGING (IMPORTANT UPGRADE)
# ============================

MEs <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes

MEDiss <- 1 - cor(MEs)

METree <- hclust(as.dist(MEDiss), method = "average")

merge <- mergeCloseModules(datExpr,
                           moduleColors,
                           cutHeight = 0.25,
                           verbose = 3)

moduleColors <- merge$colors
MEs <- merge$newMEs

# ============================
# 12. TRAIT DESIGN (EDS vs POTS)
# ============================

nSamples <- nrow(datExpr)

traitData <- data.frame(

  EDS  = rep(c(1, 0), length.out = nSamples),
  POTS = rep(c(0, 1), length.out = nSamples)

)

rownames(traitData) <- rownames(datExpr)

# ============================
# 13. MODULE–TRAIT CORRELATION
# ============================

moduleTraitCor <- cor(MEs, traitData, use = "p")

moduleTraitP <- corPvalueStudent(moduleTraitCor, nSamples)

write.csv(moduleTraitCor,
          file = file.path(WGCNA_DIR, "module_trait_correlation.csv"))

write.csv(moduleTraitP,
          file = file.path(WGCNA_DIR, "module_trait_pvalues.csv"))

# ============================
# 14. HUB GENE DETECTION (UPGRADED)
# ============================

kME <- signedKME(datExpr, MEs)

hub_scores <- apply(kME, 2, function(x) sort(x, decreasing = TRUE))

top_hubs <- lapply(hub_scores, function(x) names(x)[1:15])

write.csv(as.data.frame(top_hubs),
          file = file.path(WGCNA_DIR, "hub_genes.csv"))

# ============================
# 15. MODULE MEMBERSHIP EXPORT (Cytoscape-ready)
# ============================

module_assignment <- data.frame(
  gene = colnames(datExpr),
  module = moduleColors
)

write.csv(module_assignment,
          file = file.path(WGCNA_DIR, "module_assignment.csv"),
          row.names = FALSE)

# ============================
# 16. NETWORK EXPORT (FOR VISUALIZATION)
# ============================

# export TOM subset (top connections only)
TOM_export <- TOM

saveRDS(TOM_export,
        file = file.path(WGCNA_DIR, "TOM_network.rds"))

# ============================
# 17. BIOLOGICAL CORE GENES
# ============================

shared_core <- intersect(eds_sig, pots_sig)

write.csv(data.frame(shared_core),
          file = file.path(WGCNA_DIR, "shared_core_genes.csv"))

# ============================
# 18. FINAL SUMMARY
# ============================

summary <- data.frame(
  modules = length(unique(moduleColors)),
  genes_used = ncol(datExpr),
  samples = nrow(datExpr)
)

write.csv(summary,
          file = file.path(WGCNA_DIR, "network_summary.csv"))

message("Advanced network analysis completed")
