# ============================================================
# 01_download_data.R
# Integrative transcriptomic pipeline:
# EDS + Dysautonomia / POTS systems biology project
# ============================================================

# ============================
# 0. PACKAGES
# ============================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

# ============================
# 1. PROJECT CONFIG
# ============================

BASE_DIR <- getwd()

DATA_DIR <- file.path(BASE_DIR, "data")
RAW_DIR <- file.path(DATA_DIR, "raw")
PROCESSED_DIR <- file.path(DATA_DIR, "processed")
LOG_DIR <- file.path(BASE_DIR, "logs")

dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOG_DIR, "download_log.txt")

write_log <- function(message) {
  cat(paste0(Sys.time(), " - ", message, "\n"),
      file = log_file, append = TRUE)
  message(message)
}

write_log("=== Starting GEO download pipeline ===")

# ============================
# 2. BIOLOGICAL KEYWORDS
# ============================

keywords <- list(
  EDS = c("Ehlers Danlos", "Ehlers-Danlos", "connective tissue"),
  POTS = c("Postural Orthostatic Tachycardia", "POTS"),
  Dysautonomia = c("dysautonomia", "orthostatic intolerance", "autonomic")
)

write_log("Keywords defined for search strategy")

# ============================
# 3. MANUAL DATASET SEED LIST
# (to be expanded dynamically later)
# ============================

datasets <- list(
  EDS = c(
    # placeholder for real GSE IDs
  ),
  POTS = c(
  ),
  Control = c(
  )
)

write_log("Dataset structure initialized (empty seed lists)")

# ============================
# 4. SAFE GEO DOWNLOAD FUNCTION
# ============================

download_geo <- function(gse_id) {

  write_log(paste0("Downloading dataset: ", gse_id))

  result <- tryCatch({

    gse <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = TRUE)

    # If multiple platforms, select first
    if (class(gse) == "list") {
      expr_set <- gse[[1]]
    } else {
      expr_set <- gse
    }

    # Extract expression matrix
    expr_matrix <- exprs(expr_set)

    # Extract phenotype data
    pheno <- pData(expr_set)

    # Extract feature data (gene annotation)
    feature <- fData(expr_set)

    # Build output structure
    output <- list(
      expression = expr_matrix,
      phenotype = pheno,
      features = feature,
      platform = annotation(expr_set),
      gse = gse_id
    )

    # Save raw RDS
    saveRDS(output,
            file = file.path(RAW_DIR, paste0(gse_id, ".rds")))

    write_log(paste0("SUCCESS: ", gse_id, " downloaded"))

    return(output)

  }, error = function(e) {

    write_log(paste0("ERROR downloading ", gse_id, ": ", e$message))
    return(NULL)

  })

  return(result)
}

# ============================
# 5. DATASET CLASSIFICATION FUNCTION
# ============================

classify_dataset <- function(gse_id) {

  # simple heuristic classification (can be improved later)

  gse_info <- getGEO(gse_id, GSEMatrix = FALSE)

  title <- Meta(gse_info)$title
  summary <- Meta(gse_info)$summary

  text <- paste(title, summary)

  category <- "Unknown"

  if (grepl("Ehlers|Danlos|connective", text, ignore.case = TRUE)) {
    category <- "EDS"
  }

  if (grepl("POTS|orthostatic|dysautonomia|autonomic", text, ignore.case = TRUE)) {
    category <- "POTS"
  }

  return(category)
}

# ============================
# 6. MAIN PIPELINE
# ============================

write_log("Starting dataset processing pipeline")

all_gse <- unique(unlist(datasets))

results <- list()

for (gse_id in all_gse) {

  write_log(paste0("Processing: ", gse_id))

  # classify dataset
  class <- tryCatch({
    classify_dataset(gse_id)
  }, error = function(e) {
    "Unknown"
  })

  write_log(paste0("Classified as: ", class))

  # download dataset
  dataset <- download_geo(gse_id)

  if (!is.null(dataset)) {

    results[[gse_id]] <- list(
      data = dataset,
      class = class
    )
  }
}

# ============================
# 7. SUMMARY REPORT
# ============================

summary_table <- data.frame(
  GSE = names(results),
  Class = sapply(results, function(x) x$class)
)

write.csv(summary_table,
          file = file.path(PROCESSED_DIR, "dataset_summary.csv"),
          row.names = FALSE)

write_log("Summary table saved")

# ============================
# 8. SESSION INFO (REPRODUCIBILITY)
# ============================

write_log("Session information:")
write_log(capture.output(sessionInfo()))

write_log("=== Pipeline finished ===")
