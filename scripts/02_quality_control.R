# ============================================================
# 02_preprocessing.R
# Integrative transcriptomic preprocessing pipeline
# EDS + Dysautonomia / POTS systems biology project
# ============================================================

# ============================
# 0. PACKAGES
# ============================

suppressPackageStartupMessages({

  library(Biobase)
  library(GEOquery)
  library(limma)
  library(matrixStats)
  library(ggplot2)

})

# ============================
# 1. PATH CONFIGURATION
# ============================

BASE_DIR <- getwd()

RAW_DIR <- file.path(BASE_DIR, "data/raw")
PROCESSED_DIR <- file.path(BASE_DIR, "data/processed")
QC_DIR <- file.path(PROCESSED_DIR, "QC_reports")

dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================
# 2. DATA LOADING FUNCTION
# ============================

load_dataset <- function(gse_id) {

  file_path <- file.path(RAW_DIR, paste0(gse_id, ".rds"))

  if (!file.exists(file_path)) {
    stop(paste("Dataset not found:", gse_id))
  }

  readRDS(file_path)
}

# ============================
# 3. INITIAL QUALITY CONTROL
# ============================

qc_report <- list()

generate_qc_metrics <- function(expr_matrix, dataset_name) {

  sample_means <- colMeans(expr_matrix, na.rm = TRUE)
  sample_sds <- apply(expr_matrix, 2, sd, na.rm = TRUE)

  qc <- data.frame(
    sample = colnames(expr_matrix),
    mean_expression = sample_means,
    sd_expression = sample_sds,
    dataset = dataset_name
  )

  return(qc)
}

# ============================
# 4. OUTLIER DETECTION (ROBUST)
# ============================

detect_outliers <- function(expr_matrix) {

  sample_means <- colMeans(expr_matrix, na.rm = TRUE)

  z_scores <- scale(sample_means)

  outliers <- abs(z_scores) > 3

  keep_samples <- !outliers

  list(
    filtered_matrix = expr_matrix[, keep_samples, drop = FALSE],
    removed_samples = colnames(expr_matrix)[outliers]
  )
}

# ============================
# 5. LOG TRANSFORMATION HANDLING
# ============================

log_transform_if_needed <- function(expr_matrix) {

  max_val <- max(expr_matrix, na.rm = TRUE)

  if (max_val > 100) {
    expr_matrix <- log2(expr_matrix + 1)
  }

  return(expr_matrix)
}

# ============================
# 6. NORMALIZATION (CROSS-SAMPLE)
# ============================

normalize_dataset <- function(expr_matrix) {

  # quantile normalization (critical for cross-platform comparability)
  expr_matrix <- normalizeBetweenArrays(expr_matrix, method = "quantile")

  return(expr_matrix)
}

# ============================
# 7. GENE FILTERING (LOW SIGNAL REMOVAL)
# ============================

filter_low_expression <- function(expr_matrix, threshold = 0.2) {

  gene_variance <- rowVars(expr_matrix)

  keep_genes <- gene_variance > quantile(gene_variance, threshold)

  expr_matrix <- expr_matrix[keep_genes, ]

  return(expr_matrix)
}

# ============================
# 8. FULL PIPELINE PER DATASET
# ============================

process_dataset <- function(gse_id) {

  message(paste("Processing dataset:", gse_id))

  data <- load_dataset(gse_id)

  expr <- data$expression

  # QC metrics BEFORE filtering
  qc_before <- generate_qc_metrics(expr, gse_id)

  # Outlier detection
  outlier_result <- detect_outliers(expr)
  expr <- outlier_result$filtered_matrix

  # log transform if needed
  expr <- log_transform_if_needed(expr)

  # normalization
  expr <- normalize_dataset(expr)

  # gene filtering
  expr <- filter_low_expression(expr)

  # QC metrics AFTER filtering
  qc_after <- generate_qc_metrics(expr, gse_id)

  # save QC report
  qc_report[[gse_id]] <<- list(
    before = qc_before,
    after = qc_after,
    removed_samples = outlier_result$removed_samples
  )

  # update dataset object
  data$expression <- expr
  data$qc <- qc_report[[gse_id]]

  # save processed dataset
  saveRDS(
    data,
    file = file.path(PROCESSED_DIR, paste0(gse_id, "_processed.rds"))
  )

  return(data)
}

# ============================
# 9. LOAD ALL RAW DATASETS
# ============================

files <- list.files(RAW_DIR, pattern = "*.rds", full.names = TRUE)

gse_ids <- gsub("_.*|\\.rds", "", basename(files))

results <- list()

# ============================
# 10. MAIN LOOP
# ============================

for (gse in gse_ids) {

  try({

    results[[gse]] <- process_dataset(gse)

  }, silent = FALSE)

}

# ============================
# 11. GLOBAL QC SUMMARY
# ============================

all_qc <- do.call(rbind, lapply(names(qc_report), function(id) {

  qc_report[[id]]$after

}))

write.csv(all_qc,
          file = file.path(QC_DIR, "qc_summary.csv"),
          row.names = FALSE)

# ============================
# 12. DATASET SUMMARY
# ============================

summary_df <- data.frame(
  dataset = names(results),
  samples = sapply(results, function(x) ncol(x$expression)),
  genes = sapply(results, function(x) nrow(x$expression))
)

write.csv(summary_df,
          file = file.path(PROCESSED_DIR, "dataset_summary.csv"),
          row.names = FALSE)

# ============================
# 13. OPTIONAL VISUAL QC PLOT
# ============================

pdf(file.path(QC_DIR, "qc_expression_distribution.pdf"))

for (id in names(qc_report)) {

  hist(qc_report[[id]]$after$mean_expression,
       main = paste("Expression distribution:", id),
       xlab = "Mean expression")

}

dev.off()

# ============================
# 14. SESSION INFO (REPRODUCIBILITY)
# ============================

writeLines(capture.output(sessionInfo()),
           file.path(PROCESSED_DIR, "session_info.txt"))

message("Preprocessing pipeline completed successfully")
