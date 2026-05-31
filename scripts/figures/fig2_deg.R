# ============================================================
# FIGURE 2 — CORE DIFFERENTIAL EXPRESSION SIGNATURES
# ============================================================

library(ggplot2)
library(dplyr)
library(patchwork)
library(ggplot2)

source("scripts/figures/figure_theme.R")

fig2_deg <- function(data_path, output_dir) {

  # =========================
  # LOAD DATA
  # =========================

  eds <- read.csv(file.path(data_path, "DEG/EDS_core_genes.csv"))
  pots <- read.csv(file.path(data_path, "DEG/POTS_core_genes.csv"))
  shared <- read.csv(file.path(data_path, "DEG/shared_core_genes.csv"))

  # =========================
  # PANEL A — GENE COUNTS
  # =========================

  df_counts <- data.frame(
    group = c("EDS", "POTS", "Shared"),
    genes = c(nrow(eds), nrow(pots), nrow(shared))
  )

  p1 <- ggplot(df_counts, aes(x = group, y = genes, fill = group)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = c(
      "EDS" = bio_colors$ECM,
      "POTS" = bio_colors$AUTONOMIC,
      "Shared" = bio_colors$SHARED
    )) +
    theme_paper() +
    guides(fill = "none") +
    labs(title = "A. Differential Gene Signatures")

  # =========================
  # PANEL B — OVERLAP SCHEME
  # =========================

  overlap <- data.frame(
    category = c("EDS-specific", "POTS-specific", "Shared core"),
    value = c(
      nrow(eds) - nrow(shared),
      nrow(pots) - nrow(shared),
      nrow(shared)
    )
  )

  p2 <- ggplot(overlap, aes(x = category, y = value, fill = category)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c(
      bio_colors$ECM,
      bio_colors$AUTONOMIC,
      bio_colors$SHARED
    )) +
    theme_paper() +
    guides(fill = "none") +
    labs(title = "B. Shared Dysautonomia Signature Structure") +
    coord_flip()

  # =========================
  # PANEL C — TOP SHARED GENES
  # =========================

  top_genes <- head(shared, 15)

  top_genes$gene <- rownames(top_genes)

  p3 <- ggplot(top_genes, aes(x = reorder(gene, 1), y = 1)) +
    geom_point(size = 4, color = bio_colors$SHARED) +
    coord_flip() +
    theme_paper() +
    labs(
      title = "C. Core Shared Dysautonomia Genes",
      x = "Genes",
      y = ""
    )

  # =========================
  # COMBINE FIGURE
  # =========================

  final_fig <- (p1 | p2) / p3

  # =========================
  # SAVE
  # =========================

  ggsave(
    filename = file.path(output_dir, "fig2_deg.png"),
    plot = final_fig,
    width = 12,
    height = 8,
    dpi = 300
  )
}
