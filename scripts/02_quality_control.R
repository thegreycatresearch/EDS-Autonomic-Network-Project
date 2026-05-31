# ============================================================
# 02_quality_control.R
# Multi-cohort QC + normalization pipeline
# EDS + dysautonomia transcriptomics project
# ============================================================

# ============================
# 0. PACKAGES
# ============================

suppressPackageStartupMessages({

  library(limma)
  library(matrixStats)
  library(ggplot2)

})

# ============================
# 1. PATHS
# ============================

RAW_DIR <- "data/raw"
PROCESSED_DIR <- "data/processed"
QC_DIR <- file.path(PROCESSED_DIR, "QC")

dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================
# 2. LOAD DATASETS
# ============================

load_dataset <- function(file) {
  readRDS(file)
}

files <- list.files(RAW_DIR, pattern = ".rds", full.names = TRUE)

datasets <- lapply(files, load_dataset)
names(datasets) <- gsub(".rds", "", basename(files))

# ============================
# 3. GENE FILTERING (LOW SIGNAL REMOVAL)
# ============================

filter_low_expression <- function(expr, min_fraction = 0.2) {

  keep <- rowMeans(expr > 1) > min_fraction

  expr[keep, , drop = FALSE]

}

# ============================
# 4. SAMPLE OUTLIER DETECTION
# ============================

detect_sample_outliers <- function(expr) {

  sample_means <- colMeans(expr, na.rm = TRUE)
  z <- scale(sample_means)

  outliers <- abs(z) > 3

  list(
    filtered = expr[, !outliers, drop = FALSE],
    removed = names(sample_means)[outliers]
  )
}

# ============================
# 5. LOG TRANSFORMATION CHECK
# ============================

log_transform <- function(expr) {

  if (max(expr, na.rm = TRUE) > 100) {
    expr <- log2(expr + 1)
  }

  return(expr)
}

# ============================
# 6. NORMALIZATION (CRITICAL STEP)
# ============================

normalize_dataset <- function(expr) {

  normalizeBetweenArrays(expr, method = "quantile")

}

# ============================
# 7. BATCH EFFECT CHECK (PCA)
# ============================

run_pca <- function(expr, group_name, qc_dir) {

  pca <- prcomp(t(expr), scale. = TRUE)

  df <- data.frame(
    PC1 = pca$x[,1],
    PC2 = pca$x[,2]
  )

  p <- ggplot(df, aes(PC1, PC2)) +
    geom_point() +
    ggtitle(paste("PCA:", group_name))

  ggsave(file.path(qc_dir, paste0(group_name, "_PCA.png")), p)

}

# ============================
# 8. GENE VARIANCE FILTERING
# ============================

filter_low_variance <- function(expr, threshold = 0.25) {

  v <- apply(expr, 1, var, na.rm = TRUE)

  keep <- v > quantile(v, threshold)

  expr[keep, ]

}

# ============================
# 9. FULL QC PIPELINE PER DATASET
# ============================

process_qc <- function(data, name) {

  expr <- data$expression

  # 1. log transform
  expr <- log_transform(expr)

  # 2. remove low expression genes
  expr <- filter_low_expression(expr)

  # 3. remove sample outliers
  out <- detect_sample_outliers(expr)
  expr <- out$filtered

  # 4. normalize
  expr <- normalize_dataset(expr)

  # 5. filter low variance genes
  expr <- filter_low_variance(expr)

  # 6. PCA QC plot
  run_pca(expr, name, QC_DIR)

  data$expression <- expr
  data$qc_removed_samples <- out$removed

  return(data)
}

# ============================
# 10. APPLY TO ALL DATASETS
# ============================

processed <- list()

qc_log <- list()

for (name in names(datasets)) {

  cat("Processing:", name, "\n")

  data <- datasets[[name]]

  data <- process_qc(data, name)

  processed[[name]] <- data

  qc_log[[name]] <- list(
    samples = ncol(data$expression),
    genes = nrow(data$expression)
  )

}

# ============================
# 11. SAVE PROCESSED DATA
# ============================

for (name in names(processed)) {

  saveRDS(
    processed[[name]],
    file = file.path(PROCESSED_DIR, paste0(name, "_processed.rds"))
  )

}

# ============================
# 12. GLOBAL QC SUMMARY
# ============================

qc_df <- data.frame(
  dataset = names(qc_log),
  samples = sapply(qc_log, function(x) x$samples),
  genes = sapply(qc_log, function(x) x$genes)
)

write.csv(qc_df,
          file.path(QC_DIR, "qc_summary.csv"),
          row.names = FALSE)

message("02_quality_control.R completed successfully")
