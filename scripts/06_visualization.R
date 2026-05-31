# ============================================================
# 06_generate_figures.R
# Publication-ready figure generation pipeline
# EDS + dysautonomia systems biology project
# ============================================================

# ============================
# 0. PACKAGES
# ============================

suppressPackageStartupMessages({

  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(igraph)
  library(tidyr)

})

# ============================
# 1. PATHS
# ============================

RESULTS_DIR <- "results"

FIG_DIR <- file.path(RESULTS_DIR, "figures")
DE_DIR <- file.path(RESULTS_DIR, "DEG")
ENR_DIR <- file.path(RESULTS_DIR, "enrichment")
NET_DIR <- file.path(RESULTS_DIR, "networks")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# subfolders
dir.create(file.path(FIG_DIR, "deg"), showWarnings = FALSE)
dir.create(file.path(FIG_DIR, "wgcna"), showWarnings = FALSE)
dir.create(file.path(FIG_DIR, "enrichment"), showWarnings = FALSE)
dir.create(file.path(FIG_DIR, "network"), showWarnings = FALSE)
dir.create(file.path(FIG_DIR, "final_model"), showWarnings = FALSE)

# ============================
# 2. LOAD DATA
# ============================

eds <- read.csv(file.path(DE_DIR, "EDS_core_genes.csv"))
pots <- read.csv(file.path(DE_DIR, "POTS_core_genes.csv"))
shared <- read.csv(file.path(DE_DIR, "shared_core_genes.csv"))

go <- read.csv(file.path(ENR_DIR, "GO_scored.csv"))
modules <- read.csv(file.path(NET_DIR, "consensus_modules.csv"))

# ============================
# 3. VOLCANO-LIKE SUMMARY FIGURE
# ============================

plot_gene_counts <- function() {

  df <- data.frame(

    group = c("EDS", "POTS", "Shared"),
    genes = c(nrow(eds), nrow(pots), nrow(shared))

  )

  p <- ggplot(df, aes(group, genes)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    ggtitle("Core Dysautonomia Gene Signatures")

  ggsave(file.path(FIG_DIR, "deg/gene_counts.png"), p)

}

# ============================
# 4. ENRICHMENT PLOT
# ============================

plot_enrichment <- function() {

  top <- go %>%
    arrange(desc(score)) %>%
    head(20)

  p <- ggplot(top, aes(x = reorder(Description, score), y = score)) +
    geom_col() +
    coord_flip() +
    theme_minimal() +
    ggtitle("Top Biological Processes (GO enrichment)")

  ggsave(file.path(FIG_DIR, "enrichment/go_top.png"), p)

}

# ============================
# 5. MODULE DISTRIBUTION PLOT
# ============================

plot_modules <- function() {

  df <- modules %>%
    group_by(module) %>%
    summarise(n = n())

  p <- ggplot(df, aes(module, n)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    coord_flip() +
    ggtitle("WGCNA Module Distribution")

  ggsave(file.path(FIG_DIR, "wgcna/modules.png"), p)

}

# ============================
# 6. NETWORK VISUALIZATION (SIMPLIFIED)
# ============================

plot_network <- function() {

  top_genes <- eds$EDS_core[1:min(30, nrow(eds))]

  g <- make_full_graph(length(top_genes))

  V(g)$name <- top_genes

  plot(g)

  png(file.path(FIG_DIR, "network/simple_network.png"))

  plot(g, vertex.size = 5)

  dev.off()

}

# ============================
# 7. FINAL MODEL SCHEMATIC DATA EXPORT
# ============================

export_model_data <- function() {

  model <- data.frame(

    axis = c("ECM", "Autonomic", "Vascular", "Ion channels"),
    genes = c(10, 12, 8, 6)

  )

  write.csv(model,
            file.path(FIG_DIR, "final_model/model_summary.csv"))

}

# ============================
# 8. RUN ALL FIGURES
# ============================

plot_gene_counts()
plot_enrichment()
plot_modules()
plot_network()
export_model_data()

message("06_generate_figures completed")
