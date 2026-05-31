# ============================================================
# 03_differential_expression.R
# EDS vs POTS vs Controls
# Integrative transcriptomic analysis
# ============================================================

# ============================
# 0. PACKAGES
# ============================

suppressPackageStartupMessages({

  library(Biobase)
  library(limma)
  library(dplyr)

})

# ============================
# 1. PATHS
# ============================

PROCESSED_DIR <- "data/processed"
RESULTS_DIR <- "results"
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================
# 2. LOAD PROCESSED DATASETS
# ============================

load_processed <- function(gse_id) {

  file_path <- file.path(PROCESSED_DIR, paste0(gse_id, "_processed.rds"))

  if (!file.exists(file_path)) {
    stop(paste("Missing processed dataset:", gse_id))
  }

  readRDS(file_path)
}

# ============================
# 3. DEFINE GROUPS (IMPORTANT)
# ============================

# EDIT THESE BASED ON YOUR META-DATA
EDS_DATASETS <- c("GSE58072")
POTS_DATASETS <- c("GSE190123", "GSE145410")

# ============================
# 4. BUILD DESIGN MATRIX
# ============================

build_limma_analysis <- function(expr_matrix, group_labels) {

  group_labels <- factor(group_labels)

  design <- model.matrix(~0 + group_labels)
  colnames(design) <- levels(group_labels)

  fit <- lmFit(expr_matrix, design)

  return(list(fit = fit, design = design))
}

# ============================
# 5. DIFFERENTIAL EXPRESSION FUNCTION
# ============================

run_de_analysis <- function(expr_matrix, group_labels, contrast_name) {

  model <- build_limma_analysis(expr_matrix, group_labels)

  fit <- model$fit

  design <- model$design

  # create contrast
  contrast_matrix <- makeContrasts(
    contrasts = contrast_name,
    levels = design
  )

  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)

  results <- topTable(fit2, number = Inf, adjust.method = "BH")

  return(results)
}

# ============================
# 6. LOAD ALL DATASETS
# ============================

load_all_expr <- function(dataset_ids) {

  expr_list <- list()

  for (id in dataset_ids) {

    data <- load_processed(id)

    expr_list[[id]] <- data$expression
  }

  return(expr_list)
}

# ============================
# 7. MERGE DATA (simplified approach)
# ============================

merge_datasets <- function(expr_list) {

  # intersect genes across datasets
  common_genes <- Reduce(intersect, lapply(expr_list, rownames))

  expr_list <- lapply(expr_list, function(x) x[common_genes, ])

  merged <- do.call(cbind, expr_list)

  return(list(
    expr = merged,
    genes = common_genes
  ))
}

# ============================
# 8. BUILD GROUP LABELS
# ============================

build_labels <- function(n_eds, n_control) {

  c(
    rep("EDS", n_eds),
    rep("Control", n_control)
  )
}

build_labels_pots <- function(n_pots, n_control) {

  c(
    rep("POTS", n_pots),
    rep("Control", n_control)
  )
}

# ============================
# 9. RUN EDS ANALYSIS
# ============================

eds_list <- load_all_expr(EDS_DATASETS)
eds_merged <- merge_datasets(eds_list)

eds_labels <- build_labels(
  ncol(eds_merged$expr) / 2,
  ncol(eds_merged$expr) / 2
)

eds_degs <- run_de_analysis(
  eds_merged$expr,
  eds_labels,
  "EDS - Control"
)

write.csv(eds_degs,
          file = file.path(RESULTS_DIR, "EDS_DEGs.csv"))

# ============================
# 10. RUN POTS ANALYSIS
# ============================

pots_list <- load_all_expr(POTS_DATASETS)
pots_merged <- merge_datasets(pots_list)

pots_labels <- build_labels_pots(
  ncol(pots_merged$expr) / 2,
  ncol(pots_merged$expr) / 2
)

pots_degs <- run_de_analysis(
  pots_merged$expr,
  pots_labels,
  "POTS - Control"
)

write.csv(pots_degs,
          file = file.path(RESULTS_DIR, "POTS_DEGs.csv"))

# ============================
# 11. INTERSECTION ANALYSIS
# ============================

eds_sig <- eds_degs %>%
  filter(adj.P.Val < 0.05) %>%
  rownames()

pots_sig <- pots_degs %>%
  filter(adj.P.Val < 0.05) %>%
  rownames()

shared_genes <- intersect(eds_sig, pots_sig)

write.csv(data.frame(shared_genes),
          file = file.path(RESULTS_DIR, "shared_DEGs.csv"))

# ============================
# 12. BASIC SUMMARY
# ============================

summary <- data.frame(
  dataset = c("EDS", "POTS", "Shared"),
  genes = c(
    length(eds_sig),
    length(pots_sig),
    length(shared_genes)
  )
)

write.csv(summary,
          file = file.path(RESULTS_DIR, "DEG_summary.csv"))

message("Differential expression analysis completed")
