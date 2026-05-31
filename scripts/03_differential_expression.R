# ============================================================
# 03_differential_expression.R
# Multi-cohort differential expression + meta-signature
# EDS + POTS systems biology framework
# ============================================================

# ============================
# 0. PACKAGES
# ============================

suppressPackageStartupMessages({

  library(limma)
  library(dplyr)
  library(matrixStats)

})

# ============================
# 1. PATHS
# ============================

PROCESSED_DIR <- "data/processed"
RESULTS_DIR <- "results"
DE_DIR <- file.path(RESULTS_DIR, "DEG")

dir.create(DE_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================
# 2. LOAD PROCESSED DATASETS
# ============================

files <- list.files(PROCESSED_DIR,
                    pattern = "_processed.rds",
                    full.names = TRUE)

datasets <- lapply(files, readRDS)
names(datasets) <- gsub("_processed.rds", "", basename(files))

# ============================
# 3. GROUP DEFINITIONS
# ============================

EDS_IDS <- names(datasets)[grepl("EDS|GSE58072|GSE23186|GSE112260", names(datasets))]
POTS_IDS <- names(datasets)[grepl("POTS|GSE190123|GSE145410|GSE52038", names(datasets))]

# ============================
# 4. LIMMA CORE FUNCTION
# ============================

run_limma <- function(expr, group_labels) {

  group_labels <- factor(group_labels)

  design <- model.matrix(~0 + group_labels)
  colnames(design) <- levels(group_labels)

  fit <- lmFit(expr, design)

  if (length(levels(group_labels)) == 2) {

    contrast <- makeContrasts(
      contrasts = paste0(levels(group_labels)[1], "-", levels(group_labels)[2]),
      levels = design
    )

    fit <- contrasts.fit(fit, contrast)
    fit <- eBayes(fit)

  }

  topTable(fit, number = Inf, adjust.method = "BH")

}

# ============================
# 5. PROCESS SINGLE DATASET (IMPORTANT UPGRADE)
# ============================

process_single_dataset <- function(data, dataset_name) {

  expr <- data$expression

  # extract labels (IMPORTANT FIX: no fake balancing)
  labels <- data$sample_labels

  # ensure valid comparison exists
  if (length(unique(labels)) < 2) {
    warning(paste("Skipping", dataset_name, "- insufficient groups"))
    return(NULL)
  }

  # run limma
  res <- run_limma(expr, labels)

  res$dataset <- dataset_name

  return(res)

}

# ============================
# 6. RUN PER DATASET (NO WRONG MERGING)
# ============================

results <- list()

for (name in names(datasets)) {

  cat("Processing:", name, "\n")

  res <- process_single_dataset(datasets[[name]], name)

  if (!is.null(res)) {
    results[[name]] <- res
  }

}

# ============================
# 7. SPLIT RESULTS BY CONDITION
# ============================

eds_results <- results[EDS_IDS]
pots_results <- results[POTS_IDS]

# ============================
# 8. META-SIGNATURE EXTRACTION (KEY STEP)
# ============================

get_significant <- function(res_list) {

  sig_list <- lapply(res_list, function(x) {

    if (is.null(x)) return(NULL)

    rownames(x[x$adj.P.Val < 0.05, ])

  })

  Reduce(intersect, sig_list)

}

eds_core <- get_significant(eds_results)
pots_core <- get_significant(pots_results)

shared_core <- intersect(eds_core, pots_core)

# ============================
# 9. EFFECT SIZE CONSISTENCY FILTER (IMPORTANT UPGRADE)
# ============================

filter_consistent_direction <- function(res_list, genes) {

  keep <- c()

  for (g in genes) {

    effects <- c()

    for (res in res_list) {

      if (!g %in% rownames(res)) next

      effects <- c(effects, res[g, "logFC"])

    }

    if (length(effects) > 1) {

      if (all(sign(effects) == sign(effects[1]))) {
        keep <- c(keep, g)
      }

    }

  }

  keep

}

eds_consistent <- filter_consistent_direction(eds_results, eds_core)
pots_consistent <- filter_consistent_direction(pots_results, pots_core)
shared_consistent <- intersect(eds_consistent, pots_consistent)

# ============================
# 10. FINAL TABLES
# ============================

write.csv(data.frame(EDS_core = eds_consistent),
          file.path(DE_DIR, "EDS_core_genes.csv"))

write.csv(data.frame(POTS_core = pots_consistent),
          file.path(DE_DIR, "POTS_core_genes.csv"))

write.csv(data.frame(shared_core = shared_consistent),
          file.path(DE_DIR, "shared_core_genes.csv"))

# ============================
# 11. SUMMARY STATISTICS
# ============================

summary <- data.frame(

  category = c("EDS", "POTS", "Shared"),
  genes = c(
    length(eds_consistent),
    length(pots_consistent),
    length(shared_consistent)
  )

)

write.csv(summary,
          file.path(DE_DIR, "DEG_summary.csv"),
          row.names = FALSE)

message("03_differential_expression completed")
