# ============================================================
# 05_enrichment_integration.R
# Final integrative biological interpretation
# EDS + POTS systems biology model
# ============================================================

# ============================
# 0. PACKAGES
# ============================

suppressPackageStartupMessages({

  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(ggplot2)

})

# ============================
# 1. PATHS
# ============================

RESULTS_DIR <- "results"
NETWORK_DIR <- file.path(RESULTS_DIR, "networks")
ENRICH_DIR <- file.path(RESULTS_DIR, "enrichment")

dir.create(ENRICH_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================
# 2. LOAD CORE DATA
# ============================

eds_degs <- read.csv(file.path(RESULTS_DIR, "EDS_DEGs.csv"), row.names = 1)
pots_degs <- read.csv(file.path(RESULTS_DIR, "POTS_DEGs.csv"), row.names = 1)

shared_genes <- read.csv(file.path(NETWORK_DIR, "shared_core_genes.csv"))[,1]

module_assign <- read.csv(file.path(NETWORK_DIR, "module_assignment.csv"))

# ============================
# 3. GENE FILTERING
# ============================

eds_sig <- rownames(subset(eds_degs, adj.P.Val < 0.05))
pots_sig <- rownames(subset(pots_degs, adj.P.Val < 0.05))

all_significant <- unique(c(eds_sig, pots_sig, shared_genes))

# ============================
# 4. GENE CONVERSION (ENTREZ ID)
# ============================

gene_conversion <- bitr(all_significant,
                        fromType = "SYMBOL",
                        toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db)

entrez_genes <- gene_conversion$ENTREZID

# ============================
# 5. GO ENRICHMENT
# ============================

go_enrich <- enrichGO(
  gene = entrez_genes,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

write.csv(as.data.frame(go_enrich),
          file = file.path(ENRICH_DIR, "GO_enrichment.csv"))

# ============================
# 6. KEGG PATHWAY ANALYSIS
# ============================

kegg_enrich <- enrichKEGG(
  gene = entrez_genes,
  organism = "hsa",
  pAdjustMethod = "BH"
)

write.csv(as.data.frame(kegg_enrich),
          file = file.path(ENRICH_DIR, "KEGG_enrichment.csv"))

# ============================
# 7. NEUROVASCULAR FILTERING (KEY STEP)
# ============================

neurovascular_terms <- c(
  "vascular",
  "nervous system",
  "autonomic",
  "adrenergic",
  "smooth muscle",
  "calcium signaling",
  "extracellular matrix"
)

go_df <- as.data.frame(go_enrich)

neuro_go <- go_df %>%
  filter(grepl(paste(neurovascular_terms, collapse="|"),
                Description,
                ignore.case = TRUE))

write.csv(neuro_go,
          file = file.path(ENRICH_DIR, "neurovascular_GO.csv"))

# ============================
# 8. MODULE BIOLOGICAL ANNOTATION
# ============================

module_summary <- module_assign %>%
  group_by(module) %>%
  summarise(
    genes = n(),
    genes_list = paste(gene, collapse = ";")
  )

write.csv(module_summary,
          file = file.path(ENRICH_DIR, "module_summary.csv"))

# ============================
# 9. CORE MODEL CONSTRUCTION
# ============================

model <- list(

  ECM = c("COL1A1", "COL3A1", "COL5A1", "ELN", "TNXB"),

  autonomic = c("SLC6A2", "DBH", "TH", "ADRB1", "ADRB2"),

  ion_channels = c("SCN9A", "CACNA1C", "TRPV1")

)

write.csv(data.frame(model),
          file = file.path(ENRICH_DIR, "biological_model_core.csv"))

# ============================
# 10. FINAL INTEGRATED SUMMARY
# ============================

summary <- data.frame(
  total_genes = length(all_significant),
  go_terms = nrow(go_df),
  neurovascular_terms = nrow(neuro_go),
  kegg_pathways = nrow(as.data.frame(kegg_enrich))
)

write.csv(summary,
          file = file.path(ENRICH_DIR, "final_summary.csv"))

message("Integration + enrichment completed")
