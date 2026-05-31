# ============================================================
# FIGURE 5 — INTEGRATED MECHANISTIC MODEL
# ============================================================

library(ggplot2)
library(patchwork)

source("scripts/figures/figure_theme.R")

fig5_model <- function(output_dir) {

  # =========================
  # PANEL A — ECM DYSREGULATION
  # =========================

  p1 <- ggplot() +
    annotate("text", x = 1, y = 3,
             label = "ECM Dysregulation\n(collagen, matrix instability)",
             size = 5, fontface = "bold") +
    geom_point(aes(x = 1, y = 2), color = bio_colors$ECM, size = 8) +
    theme_void() +
    ggtitle("A. Extracellular Matrix")

  # =========================
  # PANEL B — VASCULAR SYSTEM
  # =========================

  p2 <- ggplot() +
    annotate("text", x = 1, y = 3,
             label = "Vascular instability\n(endothelial dysfunction)",
             size = 5, fontface = "bold") +
    geom_point(aes(x = 1, y = 2), color = bio_colors$VASCULAR, size = 8) +
    theme_void() +
    ggtitle("B. Vascular Regulation")

  # =========================
  # PANEL C — AUTONOMIC SYSTEM
  # =========================

  p3 <- ggplot() +
    annotate("text", x = 1, y = 3,
             label = "Autonomic dysfunction\n(POTS-like phenotype)",
             size = 5, fontface = "bold") +
    geom_point(aes(x = 1, y = 2), color = bio_colors$AUTONOMIC, size = 8) +
    theme_void() +
    ggtitle("C. Autonomic System")

  # =========================
  # PANEL D — INTEGRATED MODEL
  # =========================

  p4 <- ggplot() +
    annotate("text", x = 1, y = 3,
             label = "ECM → Vascular → Autonomic Axis\nIntegrated dysautonomia model",
             size = 5, fontface = "bold") +
    geom_segment(aes(x = 0.8, xend = 1.2, y = 2, yend = 2),
                 arrow = arrow(length = unit(0.3, "cm")),
                 color = "black") +
    theme_void() +
    ggtitle("D. Integrated Disease Model")

  # =========================
  # COMBINE FIGURE
  # =========================

  final_fig <- (p1 | p2) / (p3 | p4)

  # =========================
  # SAVE
  # =========================

  ggsave(
    filename = file.path(output_dir, "fig5_model.png"),
    plot = final_fig,
    width = 12,
    height = 8,
    dpi = 300
  )
}
