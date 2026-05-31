bio_colors <- list(
  ECM = "#2C7BB6",
  AUTONOMIC = "#D7191C",
  VASCULAR = "#FDAE61",
  ION = "#ABD9E9",
  SHARED = "#2E2E2E"
)

theme_paper <-function() {
  theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.line = element_line(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    )
}
