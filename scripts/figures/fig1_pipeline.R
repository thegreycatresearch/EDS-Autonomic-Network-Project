library(ggplot2)

fig1_pipeline <- function(output_dir) {

  p <- ggplot() +
    annotate("text", x = 1, y = 5, label = "GEO datasets") +
    annotate("text", x = 1, y = 4, label = "Quality Control") +
    annotate("text", x = 1, y = 3, label = "DEGs Analysis") +
    annotate("text", x = 1, y = 2, label = "WGCNA Networks") +
    annotate("text", x = 1, y = 1, label = "Neurovascular-ECM Model") +
    theme_void() +
    ggtitle("Study Workflow")

  ggsave(file.path(output_dir, "fig1_pipeline.png"), p, width = 6, height = 8)

}
