# ============================================================
# FIGURE 1 — PIPELINE (NATURE-LEVEL SCHEMATIC)
# ============================================================

library(ggplot2)

source("scripts/figures/figure_theme.R")

fig1_pipeline <- function(output_dir) {

  # =========================
  # STRUCTURE DATA
  # =========================

  df <- data.frame(
    step = c(
      "GEO datasets\n(EDS + POTS)",
      "Quality control\n(normalization + filtering)",
      "Differential expression\n(limma / statistical modeling)",
      "Co-expression networks\n(WGCNA)",
      "Functional enrichment\n(GO / KEGG)",
      "Integrated model\n(neurovascular–ECM axis)"
    ),
    x = 1:6,
    y = c(6, 5, 4, 3, 2, 1)
  )

  # =========================
  # PLOT
  # =========================

  p <- ggplot(df, aes(x = x, y = y)) +

    # nodes
    geom_point(size = 6, color = bio_colors$SHARED) +

    # labels
    geom_text(aes(label = step),
              vjust = -1,
              size = 3.8,
              fontface = "bold") +

    # flow arrows (simple approximation)
    geom_segment(
      aes(x = x, xend = x,
          y = y - 0.3, yend = y + 0.3),
      arrow = arrow(length = unit(0.15, "cm")),
      color = "gray40"
    ) +

    theme_void() +

    ggtitle("Study Design and Analytical Workflow") +

    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )

  # =========================
  # SAVE
  # =========================

  ggsave(
    filename = file.path(output_dir, "fig1_pipeline.png"),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
  )
}
