# ============================================================
# GRAPHICAL ABSTRACT — EDS + DYSAUTONOMIA MODEL
# ============================================================

library(ggplot2)

source("scripts/figures/figure_theme.R")

graphical_abstract <- function(output_dir) {

  # =========================
  # FLOW DATA
  # =========================

  df <- data.frame(
    x = c(1, 2, 3, 4, 5),
    y = c(3, 3, 3, 3, 3),
    label = c(
      "Genetic predisposition\n(EDS / connective tissue)",
      "ECM dysfunction\n(collagen instability)",
      "Vascular instability\n(endothelial fragility)",
      "Autonomic dysregulation\n(POTS phenotype)",
      "Integrated dysautonomia model"
    ),
    color = c(
      "#999999",
      bio_colors$ECM,
      bio_colors$VASCULAR,
      bio_colors$AUTONOMIC,
      bio_colors$SHARED
    )
  )

  # =========================
  # BASE PLOT
  # =========================

  p <- ggplot(df, aes(x = x, y = y)) +

    # nodes
    geom_point(aes(color = color), size = 10) +
    scale_color_identity() +

    # labels
    geom_text(aes(label = label),
              vjust = -2,
              size = 3.8,
              fontface = "bold") +

    # flow arrows
    geom_segment(
      aes(x = x, xend = x + 1, y = y, yend = y),
      arrow = arrow(length = unit(0.25, "cm")),
      color = "gray40",
      linewidth = 0.8
    ) +

    theme_void() +

    ggtitle("Integrated Model of EDS-associated Dysautonomia") +

    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )

  # =========================
  # SAVE
  # =========================

  ggsave(
    filename = file.path(output_dir, "graphical_abstract.png"),
    plot = p,
    width = 14,
    height = 4,
    dpi = 300
  )
}
