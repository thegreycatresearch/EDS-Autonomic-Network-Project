library(ggplot2)

fig2_deg <- function(data_path, output_dir) {

  eds <- read.csv(file.path(data_path, "DEG/EDS_core_genes.csv"))
  pots <- read.csv(file.path(data_path, "DEG/POTS_core_genes.csv"))
  shared <- read.csv(file.path(data_path, "DEG/shared_core_genes.csv"))

  df <- data.frame(
    group = c("EDS", "POTS", "Shared"),
    genes = c(nrow(eds), nrow(pots), nrow(shared))
  )

  p <- ggplot(df, aes(group, genes)) +
    geom_col() +
    theme_minimal() +
    labs(title = "Core Dysregulated Genes")

  ggsave(file.path(output_dir, "fig2_deg.png"), p)

}
