library(ggplot2)
library(dplyr)

fig4_enrichment <- function(data_path, output_dir) {

  go <- read.csv(file.path(data_path, "enrichment/GO_scored.csv"))

  top <- go %>%
    arrange(desc(score)) %>%
    head(15)

  p <- ggplot(top, aes(x = reorder(Description, score), y = score)) +
    geom_col() +
    coord_flip() +
    theme_minimal() +
    labs(title = "Biological Enrichment (GO)")

  ggsave(file.path(output_dir, "fig4_enrichment.png"), p)

}
