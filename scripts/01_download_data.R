# ============================================================
# 01_download_data.R
# Automated GEO acquisition + QC + standardization pipeline
# EDS + dysautonomia systems biology project
# ============================================================

# ============================
# 0. PACKAGES
# ============================

suppressPackageStartupMessages({

  library(GEOquery)
  library(Biobase)
  library(dplyr)
  library(stringr)

})

# ============================
# 1. PATHS
# ============================

BASE_DIR <- getwd()

DATA_DIR <- file.path(BASE_DIR, "data")
RAW_DIR <- file.path(DATA_DIR, "raw")
LOG_DIR <- file.path(DATA_DIR, "logs")

dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================
# 2. LOGGING
# ============================

log_file <- file.path(LOG_DIR, "download_log.txt")

write_log <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"),
      file = log_file,
      append = TRUE)
}

# ============================
# 3. DATASET STRUCTURE (EXPANDED)
# ============================

dataset_map <- list(

  EDS = c(
    "GSE58072",
    "GSE23186",
    "GSE112260"
  ),

  POTS = c(
    "GSE190123",
    "GSE145410",
    "GSE52038"
  )

)

# ============================
# 4. PLATFORM CHECK (NEW)
# ============================

check_platform <- function(gse) {

  platform <- annotation(gse[[1]])

  write_log(paste("Platform detected:", platform))

  return(platform)
}

# ============================
# 5. SAFE GEO DOWNLOAD (IMPROVED)
# ============================

download_geo_safe <- function(gse_id) {

  tryCatch({

    write_log(paste("Downloading", gse_id))

    gse <- getGEO(gse_id, GSEMatrix = TRUE)

    # handle multiple platforms (IMPORTANT FIX)
    if (length(gse) > 1) {
      write_log(paste("Multiple platforms detected in", gse_id, "- using first"))
      gse <- gse[[1]]
    } else {
      gse <- gse[[1]]
    }

    expr <- exprs(gse)
    pheno <- pData(gse)
    feature <- fData(gse)

    platform <- annotation(gse)

    list(
      expression = expr,
      phenotype = pheno,
      features = feature,
      gse_id = gse_id,
      platform = platform
    )

  }, error = function(e) {

    write_log(paste("FAILED:", gse_id, e$message))

    return(NULL)

  })

}

# ============================
# 6. PROBE → GENE CLEANING (IMPORTANT UPGRADE)
# ============================

clean_expression_matrix <- function(expr, feature_data) {

  # try to map probes to gene symbols
  if ("Gene Symbol" %in% colnames(feature_data)) {

    genes <- feature_data$`Gene Symbol`

    # remove empty mappings
    keep <- genes != "" & !is.na(genes)

    expr <- expr[keep, ]
    rownames(expr) <- genes[keep]

    # collapse duplicates (mean expression)
    expr <- aggregate(expr,
                      by = list(rownames(expr)),
                      FUN = mean)

    rownames(expr) <- expr$Group.1
    expr <- expr[, -1]

  }

  return(expr)
}

# ============================
# 7. BASIC QC (NEW ADDITION)
# ============================

qc_metrics <- function(expr) {

  list(
    NA_fraction = mean(is.na(expr)),
    gene_variance = mean(apply(expr, 1, var, na.rm = TRUE)),
    sample_variance = mean(apply(expr, 2, var, na.rm = TRUE))
  )

}

# ============================
# 8. SAMPLE CLASSIFICATION
# ============================

classify_samples <- function(pheno_data) {

  labels <- rep("Unknown", nrow(pheno_data))

  if ("characteristics_ch1" %in% colnames(pheno_data)) {

    text <- tolower(paste(pheno_data$characteristics_ch1, collapse = " "))

    if (grepl("eds|ehlers|connective", text)) {
      labels[] <- "EDS"
    }

    if (grepl("pots|orthostatic|dysautonomia", text)) {
      labels[] <- "POTS"
    }

    if (grepl("control|healthy", text)) {
      labels[] <- "Control"
    }

  }

  return(labels)
}

# ============================
# 9. SAVE DATASET
# ============================

save_dataset <- function(data, group_name) {

  file_path <- file.path(RAW_DIR, paste0(data$gse_id, "_", group_name, ".rds"))

  saveRDS(data, file_path)

  write_log(paste("Saved", data$gse_id, "as", group_name))

}

# ============================
# 10. MAIN PIPELINE (ENHANCED)
# ============================

all_datasets <- list()

qc_summary <- list()

for (group in names(dataset_map)) {

  for (gse_id in dataset_map[[group]]) {

    data <- download_geo_safe(gse_id)

    if (is.null(data)) next

    # QC before cleaning
    qc_summary[[gse_id]] <- qc_metrics(data$expression)

    # CLEAN EXPRESSION MATRIX (NEW STEP)
    data$expression <- clean_expression_matrix(
      data$expression,
      data$features
    )

    # sample classification
    data$sample_labels <- classify_samples(data$phenotype)

    data$group <- group

    # store
    all_datasets[[gse_id]] <- data

    # save
    save_dataset(data, group)

  }
}

# ============================
# 11. GLOBAL SUMMARY (ENHANCED)
# ============================

summary_df <- data.frame(

  dataset = names(all_datasets),
  samples = sapply(all_datasets, function(x) ncol(x$expression)),
  genes = sapply(all_datasets, function(x) nrow(x$expression)),
  platform = sapply(all_datasets, function(x) x$platform)

)

write.csv(summary_df,
          file.path(RAW_DIR, "dataset_summary.csv"),
          row.names = FALSE)

# ============================
# 12. QC EXPORT (NEW)
# ============================

qc_df <- do.call(rbind, lapply(names(qc_summary), function(id) {

  data.frame(
    dataset = id,
    NA_fraction = qc_summary[[id]]$NA_fraction,
    gene_variance = qc_summary[[id]]$gene_variance,
    sample_variance = qc_summary[[id]]$sample_variance
  )

}))

write.csv(qc_df,
          file.path(RAW_DIR, "qc_summary.csv"),
          row.names = FALSE)

# ============================
# 13. SESSION INFO
# ============================

writeLines(capture.output(sessionInfo()),
           file.path(LOG_DIR, "session_info.txt"))

write_log("Enhanced download + QC pipeline completed")

message("01_download_data.R finished (PRO VERSION)")
