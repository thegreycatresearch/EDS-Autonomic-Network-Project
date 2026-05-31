library(ggplot2)
library(dplyr)

fig3_wgcna <- function(data_path, output_dir) {

  modules <- read.csv(file.path(data_path, "networks/consensus_modules.csv"))

  df <- modules %>%
    group_by(module) %>%
    summarise(n = n())

  p <- ggplot(df, aes(module, n)) +
    geom_col() +
    coord_flip() +
    theme_minimal() +
    labs(title = "WGCNA Module Distribution")

  ggsave(file.path(output_dir, "fig3_wgcna.png"), p)

}
