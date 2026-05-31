library(ggplot2)

fig5_model <- function(output_dir) {

  p <- ggplot() +
    annotate("text", x = 1, y = 3, label = "ECM dysregulation") +
    annotate("text", x = 1, y = 2, label = "Vascular instability") +
    annotate("text", x = 1, y = 1, label = "Autonomic dysfunction") +
    theme_void() +
    ggtitle("Integrated Disease Model")

  ggsave(file.path(output_dir, "fig5_model.png"), p)

}
