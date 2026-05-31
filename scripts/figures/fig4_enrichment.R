# ============================================================
# FIGURE 4 — FUNCTIONAL ENRICHMENT ANALYSIS
# ============================================================

library(ggplot2)
library(dplyr)
library(patchwork)

source("scripts/figures/figure_theme.R")

fig4_enrichment <- function(data_path, output_dir) {

  # =========================
  # LOAD DATA
  # =========================

  go <- read.csv(file.path(data_path, "enrichment/GO_scored.csv"))

  # =========================
  # PANEL A — TOP ENRICHMENT TERMS
  # =========================

  top <- go %>%
    arrange(desc(score)) %>%
    head(15)

  p1 <- ggplot(top, aes(
    x = reorder(Description, score),
    y = score
  )) +
    geom_col(fill = bio_colors$AUTONOMIC) +
    coord_flip() +
    theme_paper() +
    labs(title = "A. Functional Enrichment (GO Terms)")

  # =========================
  # PANEL B — PROCESS CLUSTERING (SIMPLIFIED)
  # =========================

  go$category <- ifelse(
    grepl("matrix|collagen|extracellular", tolower(go$Description)),
    "ECM",
    ifelse(
      grepl("vascular|endothelial|blood", tolower(go$Description)),
      "Vascular",
      ifelse(
        grepl("neuron|autonomic|synapse", tolower(go$Description)),
        "Autonomic",
        "Other"
      )
    )
  )

  df_cat <- go %>%
    group_by(category) %>%
    summarise(score = sum(score))

  p2 <- ggplot(df_cat, aes(category, score, fill = category)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = c(
      "ECM" = bio_colors$ECM,
      "Vascular" = bio_colors$VASCULAR,
      "Autonomic" = bio_colors$AUTONOMIC,
      "Other" = "gray70"
    )) +
    theme_paper() +
    guides(fill = "none") +
    labs(title = "B. Biological Process Clusters")

  # =========================
  # PANEL C — SYSTEM AXIS MODEL
  # =========================

  axis_df <- data.frame(
    axis = c("ECM Remodeling", "Vascular Regulation", "Autonomic Signaling"),
    strength = c(
      sum(grepl("matrix|collagen", tolower(go$Description))),
      sum(grepl("vascular|endothelial", tolower(go$Description))),
      sum(grepl("neuron|autonomic", tolower(go$Description)))
    )
  )

  p3 <- ggplot(axis_df, aes(axis, strength, fill = axis)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = c(
      "ECM Remodeling" = bio_colors$ECM,
      "Vascular Regulation" = bio_colors$VASCULAR,
      "Autonomic Signaling" = bio_colors$AUTONOMIC
    )) +
    theme_paper() +
    guides(fill = "none") +
    labs(title = "C. Emergent Neurovascular–ECM Axis")

  # =========================
  # COMBINE FIGURE
  # =========================

  final_fig <- (p1 | p2) / p3

  # =========================
  # SAVE
  # =========================

  ggsave(
    filename = file.path(output_dir, "fig4_enrichment.png"),
    plot = final_fig,
    width = 12,
    height = 8,
    dpi = 300
  )
}
