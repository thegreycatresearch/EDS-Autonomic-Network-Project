# ============================================================
# 05_enrichment_integration.R
# Multi-layer functional interpretation + systems biology model
# EDS + dysautonomia network framework
# ============================================================

# ============================
# 0. PACKAGES
# ============================

suppressPackageStartupMessages({

  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(tidyr)
  library(ggplot2)

})

# ============================
# 1. PATHS
# ============================

RESULTS_DIR <- "results"
DE_DIR <- file.path(RESULTS_DIR, "DEG")
NET_DIR <- file.path(RESULTS_DIR, "networks")
ENR_DIR <- file.path(RESULTS_DIR, "enrichment")

dir.create(ENR_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================
# 2. LOAD CORE DATA
# ============================

eds_core <- read.csv(file.path(DE_DIR, "EDS_core_genes.csv"))[,1]
pots_core <- read.csv(file.path(DE_DIR, "POTS_core_genes.csv"))[,1]
shared_core <- read.csv(file.path(DE_DIR, "shared_core_genes.csv"))[,1]

modules <- read.csv(file.path(NET_DIR, "consensus_modules.csv"))

all_genes <- unique(c(eds_core, pots_core, shared_core))

# ============================
# 3. GENE CONVERSION
# ============================

gene_map <- bitr(all_genes,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)

entrez <- gene_map$ENTREZID

# ============================
# 4. GLOBAL ENRICHMENT (BASE)
# ============================

ego <- enrichGO(
  gene = entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  readable = TRUE,
  pAdjustMethod = "BH"
)

ekegg <- enrichKEGG(
  gene = entrez,
  organism = "hsa"
)

# ============================
# 5. REDUNDANCY REDUCTION (IMPORTANT UPGRADE)
# ============================

simplified_go <- simplify(ego, cutoff = 0.7, by = "p.adjust")

# ============================
# 6. MODULE-LEVEL ENRICHMENT (KEY STEP)
# ============================

module_enrichment <- list()

for (mod in unique(modules$EDS)) {

  mod_genes <- modules$gene[modules$EDS == mod]

  mod_map <- bitr(mod_genes,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)

  if (nrow(mod_map) < 10) next

  ego_mod <- enrichGO(
    gene = mod_map$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    readable = TRUE
  )

  module_enrichment[[mod]] <- ego_mod

}

# ============================
# 7. NEUROVASCULAR AXIS EXTRACTION
# ============================

keywords <- c(
  "vascular",
  "autonomic",
  "adrenergic",
  "extracellular matrix",
  "collagen",
  "smooth muscle",
  "calcium",
  "ion transport",
  "endothelial"
)

go_df <- as.data.frame(ego)

neuro_axis <- go_df %>%
  filter(grepl(paste(keywords, collapse="|"),
               Description,
               ignore.case = TRUE))

# ============================
# 8. PATHWAY SCORING (NEW IMPORTANT STEP)
# ============================

score_pathways <- function(df) {

  df$score <- -log10(df$p.adjust) * df$Count

  df <- df %>%
    arrange(desc(score))

  df

}

scored_go <- score_pathways(go_df)

scored_kegg <- score_pathways(as.data.frame(ekegg))

# ============================
# 9. BIOLOGICAL AXES CONSTRUCTION
# ============================

biological_axes <- list(

  ECM_axis = c("collagen", "extracellular matrix", "fibrosis"),

  autonomic_axis = c("adrenergic", "nervous system", "synapse"),

  vascular_axis = c("endothelial", "vascular", "smooth muscle"),

  ion_axis = c("ion transport", "calcium", "channel")

)

# ============================
# 10. AXIS MAPPING FUNCTION
# ============================

map_to_axis <- function(df, axis_terms) {

  df$axis <- "other"

  for (axis in names(axis_terms)) {

    terms <- axis_terms[[axis]]

    idx <- grepl(paste(terms, collapse="|"),
                 df$Description,
                 ignore.case = TRUE)

    df$axis[idx] <- axis
  }

  df

}

go_axis <- map_to_axis(go_df, biological_axes)

# ============================
# 11. FINAL MODEL CONSTRUCTION
# ============================

final_model <- data.frame(

  gene_count = c(
    length(eds_core),
    length(pots_core),
    length(shared_core)
  ),

  category = c("EDS", "POTS", "Shared")

)

# ============================
# 12. EXPORT RESULTS
# ============================

write.csv(go_df,
          file.path(ENR_DIR, "GO_full.csv"))

write.csv(neuro_axis,
          file.path(ENR_DIR, "neurovascular_axis.csv"))

write.csv(scored_go,
          file.path(ENR_DIR, "GO_scored.csv"))

write.csv(scored_kegg,
          file.path(ENR_DIR, "KEGG_scored.csv"))

write.csv(go_axis,
          file.path(ENR_DIR, "GO_axis_mapped.csv"))

write.csv(final_model,
          file.path(ENR_DIR, "biological_summary.csv"))

message("05_enrichment_integration completed")
