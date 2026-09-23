# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("DropletUtils")

library(DropletUtils)
packageVersion("DropletUtils")



# V02A — GSE155882
# CELL / BARCODE RECOVERY AUDIT
#
# PURPOSE:
#   1. Inspect the raw droplet matrices from V01
#   2. Determine mitochondrial gene naming
#   3. Examine barcode UMI / feature distributions
#   4. Generate barcode-rank diagnostics
#   5. Run EmptyDrops as an INDEPENDENT diagnostic
#   6. Define a conservative candidate-cell set
#   7. DO NOT perform biological QC yet
#
# IMPORTANT:
#   - No CYP/AA genes are used for cell selection
#   - No condition-specific biological optimization
#   - No normalization
#   - No integration
#   - No clustering
#   - No differential expression
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(Seurat)
library(Matrix)

if (!requireNamespace("DropletUtils", quietly = TRUE)) {
  stop(
    paste0(
      "Package 'DropletUtils' is required for V02A.\n",
      "Install it with:\n",
      "if (!requireNamespace('BiocManager', quietly = TRUE)) ",
      "install.packages('BiocManager')\n",
      "BiocManager::install('DropletUtils')\n",
      "Then rerun V02A."
    )
  )
}


# ============================================================
# 2. PATHS
# ============================================================

project_dir <- "D:/Master/ScRNA seq/TAC/VALIDATION/GSE155882"

clean_dir <- file.path(
  project_dir,
  "CLEAN DATA"
)

results_dir <- file.path(
  project_dir,
  "RESULTS",
  "V02A_CELL_BARCODE_RECOVERY"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "V02A_CELL_BARCODE_RECOVERY"
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


input_file <- file.path(
  clean_dir,
  "V01_GSE155882_RAW_SEURAT_OBJECTS.rds"
)


# ============================================================
# 3. LOAD V01 OBJECTS
# ============================================================

if (!file.exists(input_file)) {
  stop(
    paste0(
      "V01 object not found:\n",
      input_file
    )
  )
}

raw_objects <- readRDS(
  input_file
)

cat("\n============================================\n")
cat("V02A — CELL / BARCODE RECOVERY AUDIT\n")
cat("============================================\n")

cat(
  "Objects loaded:",
  length(raw_objects),
  "\n"
)

print(
  names(raw_objects)
)


# ============================================================
# 4. VERIFY EXPECTED SAMPLE STRUCTURE
# ============================================================

expected_samples <- c(
  "Sham_Rep1",
  "Sham_Rep2",
  "TAC_Rep1",
  "TAC_Rep2"
)

if (!identical(
  sort(names(raw_objects)),
  sort(expected_samples)
)) {
  stop(
    "Unexpected sample structure in V01 object."
  )
}


# ============================================================
# 5. GENE-NAMING AUDIT
# ============================================================

genes <- rownames(
  raw_objects[[1]]
)

mito_patterns <- c(
  "^mt-",
  "^MT-",
  "^Mt-",
  "^mt\\.",
  "^MT\\.",
  "^Mt\\."
)

mito_pattern_audit <- data.frame(
  Pattern = mito_patterns,
  Number_of_Genes = sapply(
    mito_patterns,
    function(pattern_now) {
      length(
        grep(
          pattern_now,
          genes,
          value = TRUE
        )
      )
    }
  ),
  stringsAsFactors = FALSE
)

cat("\n============================================\n")
cat("MITOCHONDRIAL GENE-NAMING AUDIT\n")
cat("============================================\n")

print(
  mito_pattern_audit,
  row.names = FALSE
)

write.csv(
  mito_pattern_audit,
  file.path(
    results_dir,
    "V02A_Mitochondrial_Gene_Naming_Audit.csv"
  ),
  row.names = FALSE
)


best_mito_index <- which.max(
  mito_pattern_audit$Number_of_Genes
)

best_mito_pattern <- mito_pattern_audit$Pattern[
  best_mito_index
]

best_mito_n <- mito_pattern_audit$Number_of_Genes[
  best_mito_index
]


if (best_mito_n == 0) {
  best_mito_pattern <- NA_character_
}


cat(
  "\nSelected mitochondrial pattern:",
  ifelse(
    is.na(best_mito_pattern),
    "NONE FOUND",
    best_mito_pattern
  ),
  "\n"
)


if (!is.na(best_mito_pattern)) {
  
  mito_genes <- grep(
    best_mito_pattern,
    genes,
    value = TRUE
  )
  
  cat(
    "\nFirst mitochondrial genes:\n"
  )
  
  print(
    head(
      mito_genes,
      30
    )
  )
  
} else {
  
  mito_genes <- character(0)
  
  cat(
    "\nWARNING: No mitochondrial prefix detected.\n"
  )
}


# ============================================================
# 6. HELPER — GET RAW COUNTS
# ============================================================

get_counts <- function(obj) {
  
  GetAssayData(
    obj,
    assay = "RNA",
    layer = "counts"
  )
}


# ============================================================
# 7. BARCODE-LEVEL SUMMARY
#
# We calculate this directly from the sparse count matrix.
# ============================================================

barcode_summaries <- list()
rank_tables <- list()


for (sample_name in names(raw_objects)) {
  
  cat(
    "\nProcessing barcode summary:",
    sample_name,
    "\n"
  )
  
  counts <- get_counts(
    raw_objects[[sample_name]]
  )
  
  total_umi <- Matrix::colSums(
    counts
  )
  
  detected_genes <- Matrix::colSums(
    counts > 0
  )
  
  
  if (length(mito_genes) > 0) {
    
    mt_umi <- Matrix::colSums(
      counts[
        mito_genes,
        ,
        drop = FALSE
      ]
    )
    
    percent_mt <- ifelse(
      total_umi > 0,
      100 * mt_umi / total_umi,
      NA_real_
    )
    
  } else {
    
    percent_mt <- rep(
      NA_real_,
      length(total_umi)
    )
  }
  
  
  barcode_summary <- data.frame(
    Barcode = colnames(counts),
    Total_UMI = as.numeric(total_umi),
    Detected_Genes = as.numeric(detected_genes),
    Percent_MT = as.numeric(percent_mt),
    stringsAsFactors = FALSE
  )
  
  
  barcode_summary <- barcode_summary[
    order(
      barcode_summary$Total_UMI,
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  
  barcode_summary$UMI_Rank <- seq_len(
    nrow(barcode_summary)
  )
  
  
  barcode_summary$Log10_Rank <- log10(
    barcode_summary$UMI_Rank
  )
  
  barcode_summary$Log10_UMI <- log10(
    barcode_summary$Total_UMI + 1
  )
  
  
  barcode_summaries[[sample_name]] <- barcode_summary
  
  
  rank_tables[[sample_name]] <- barcode_summary[
    ,
    c(
      "Barcode",
      "UMI_Rank",
      "Total_UMI",
      "Detected_Genes",
      "Percent_MT",
      "Log10_Rank",
      "Log10_UMI"
    )
  ]
  
  
  write.csv(
    barcode_summary,
    file.path(
      results_dir,
      paste0(
        "V02A_",
        sample_name,
        "_Barcode_Summary.csv"
      )
    ),
    row.names = FALSE
  )
}


# ============================================================
# 8. GLOBAL BARCODE DISTRIBUTION SUMMARY
# ============================================================

distribution_summary_list <- list()


for (sample_name in names(barcode_summaries)) {
  
  x <- barcode_summaries[[sample_name]]
  
  distribution_summary_list[[sample_name]] <- data.frame(
    
    Sample_ID = sample_name,
    
    Condition = unique(
      raw_objects[[sample_name]]$condition
    ),
    
    Total_Barcodes = nrow(x),
    
    UMI_gt_0 = sum(
      x$Total_UMI > 0
    ),
    
    UMI_ge_10 = sum(
      x$Total_UMI >= 10
    ),
    
    UMI_ge_50 = sum(
      x$Total_UMI >= 50
    ),
    
    UMI_ge_100 = sum(
      x$Total_UMI >= 100
    ),
    
    UMI_ge_200 = sum(
      x$Total_UMI >= 200
    ),
    
    UMI_ge_500 = sum(
      x$Total_UMI >= 500
    ),
    
    UMI_ge_1000 = sum(
      x$Total_UMI >= 1000
    ),
    
    Genes_ge_50 = sum(
      x$Detected_Genes >= 50
    ),
    
    Genes_ge_100 = sum(
      x$Detected_Genes >= 100
    ),
    
    Genes_ge_200 = sum(
      x$Detected_Genes >= 200
    ),
    
    Genes_ge_300 = sum(
      x$Detected_Genes >= 300
    ),
    
    Genes_ge_500 = sum(
      x$Detected_Genes >= 500
    ),
    
    stringsAsFactors = FALSE
  )
}


distribution_summary <- do.call(
  rbind,
  distribution_summary_list
)

rownames(
  distribution_summary
) <- NULL


cat("\n============================================\n")
cat("BARCODE DISTRIBUTION SUMMARY\n")
cat("============================================\n")

print(
  distribution_summary,
  row.names = FALSE
)


write.csv(
  distribution_summary,
  file.path(
    results_dir,
    "V02A_Barcode_Distribution_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9. BARCODE-RANK PLOTS
# ============================================================

pdf(
  file.path(
    figures_dir,
    "V02A_Barcode_Rank_Plots.pdf"
  ),
  width = 7,
  height = 6
)


for (sample_name in names(rank_tables)) {
  
  x <- rank_tables[[sample_name]]
  
  plot(
    x$UMI_Rank,
    x$Total_UMI + 1,
    log = "xy",
    pch = 16,
    cex = 0.25,
    xlab = "Barcode rank",
    ylab = "Total UMI + 1",
    main = paste(
      sample_name,
      "Barcode Rank"
    )
  )
  
  abline(
    h = c(
      10,
      50,
      100,
      200,
      500,
      1000
    ),
    lty = 3
  )
}


dev.off()


# ============================================================
# 10. FEATURE-vs-UMI DIAGNOSTIC PLOTS
# ============================================================

pdf(
  file.path(
    figures_dir,
    "V02A_Features_vs_UMI.pdf"
  ),
  width = 7,
  height = 6
)


for (sample_name in names(barcode_summaries)) {
  
  x <- barcode_summaries[[sample_name]]
  
  positive <- x$Total_UMI > 0
  
  plot(
    x$Total_UMI[positive],
    x$Detected_Genes[positive],
    log = "xy",
    pch = 16,
    cex = 0.25,
    xlab = "Total UMI",
    ylab = "Detected genes",
    main = paste(
      sample_name,
      "Features vs UMI"
    )
  )
}


dev.off()


# ============================================================
# 11. EMPTYDROPS
#
# EmptyDrops is used here as a diagnostic/cell-calling tool.
#
# We use lower = 100 to avoid spending computation on the
# enormous ultra-low-count barcode background.
#
# Retain = 1000 automatically keeps high-count barcodes.
#
# FDR <= 0.01 is used for the statistical cell call.
# ============================================================

emptydrops_results <- list()
emptydrops_summary_list <- list()


for (sample_name in names(raw_objects)) {
  
  cat(
    "\n============================================\n"
  )
  
  cat(
    "Running EmptyDrops:",
    sample_name,
    "\n"
  )
  
  cat(
    "============================================\n"
  )
  
  
  counts <- get_counts(
    raw_objects[[sample_name]]
  )
  
  
  ed <- DropletUtils::emptyDrops(
    counts,
    lower = 100,
    retain = 1000,
    niters = 10000
  )
  
  
  ed_df <- as.data.frame(
    ed
  )
  
  
  ed_df$Barcode <- rownames(
    ed_df
  )
  
  
  barcode_info <- barcode_summaries[[sample_name]]
  
  
  match_index <- match(
    ed_df$Barcode,
    barcode_info$Barcode
  )
  
  
  ed_df$Total_UMI <- barcode_info$Total_UMI[
    match_index
  ]
  
  ed_df$Detected_Genes <- barcode_info$Detected_Genes[
    match_index
  ]
  
  ed_df$Percent_MT <- barcode_info$Percent_MT[
    match_index
  ]
  
  
  ed_df$EmptyDrops_FDR_0.01 <- (
    !is.na(ed_df$FDR) &
      ed_df$FDR <= 0.01
  )
  
  
  emptydrops_results[[sample_name]] <- ed_df
  
  
  called_n <- sum(
    ed_df$EmptyDrops_FDR_0.01,
    na.rm = TRUE
  )
  
  
  high_count_n <- sum(
    ed_df$Total_UMI >= 1000,
    na.rm = TRUE
  )
  
  
  emptydrops_summary_list[[sample_name]] <- data.frame(
    
    Sample_ID = sample_name,
    
    Condition = unique(
      raw_objects[[sample_name]]$condition
    ),
    
    Total_Barcodes = ncol(counts),
    
    Tested_or_Retained_by_EmptyDrops = sum(
      !is.na(ed_df$FDR)
    ),
    
    EmptyDrops_FDR_0.01 = called_n,
    
    UMI_ge_1000 = high_count_n,
    
    stringsAsFactors = FALSE
  )
  
  
  write.csv(
    ed_df,
    file.path(
      results_dir,
      paste0(
        "V02A_",
        sample_name,
        "_EmptyDrops.csv"
      )
    ),
    row.names = FALSE
  )
}


emptydrops_summary <- do.call(
  rbind,
  emptydrops_summary_list
)

rownames(
  emptydrops_summary
) <- NULL


cat("\n============================================\n")
cat("EMPTYDROPS SUMMARY\n")
cat("============================================\n")

print(
  emptydrops_summary,
  row.names = FALSE
)


write.csv(
  emptydrops_summary,
  file.path(
    results_dir,
    "V02A_EmptyDrops_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. DEFINE CANDIDATE-CELL BARCODE SET
#
# This is NOT final QC.
#
# Candidate cells =
#   EmptyDrops FDR <= 0.01
#   OR
#   total UMI >= 1000
#
# The UMI criterion protects high-count barcodes that are
# automatically retained / may have NA statistical values.
# ============================================================

candidate_objects <- list()
candidate_summary_list <- list()


for (sample_name in names(raw_objects)) {
  
  obj <- raw_objects[[sample_name]]
  
  counts <- get_counts(
    obj
  )
  
  ed_df <- emptydrops_results[[sample_name]]
  
  
  ed_called <- ed_df$Barcode[
    !is.na(ed_df$FDR) &
      ed_df$FDR <= 0.01
  ]
  
  
  total_umi <- Matrix::colSums(
    counts
  )
  
  
  high_count <- names(total_umi)[
    total_umi >= 1000
  ]
  
  
  candidate_barcodes <- union(
    ed_called,
    high_count
  )
  
  
  candidate_barcodes <- intersect(
    candidate_barcodes,
    colnames(obj)
  )
  
  
  if (length(candidate_barcodes) == 0) {
    stop(
      paste(
        "No candidate cells recovered for",
        sample_name
      )
    )
  }
  
  
  candidate_obj <- subset(
    obj,
    cells = candidate_barcodes
  )
  
  
  if (!is.na(best_mito_pattern)) {
    
    candidate_obj[["percent.mt"]] <- PercentageFeatureSet(
      candidate_obj,
      pattern = best_mito_pattern
    )
  }
  
  
  candidate_obj$barcode_recovery_method <-
    "EmptyDrops_FDR0.01_or_UMIge1000"
  
  
  candidate_objects[[sample_name]] <- candidate_obj
  
  
  candidate_summary_list[[sample_name]] <- data.frame(
    
    Sample_ID = sample_name,
    
    Condition = unique(
      candidate_obj$condition
    ),
    
    Replicate = unique(
      candidate_obj$replicate
    ),
    
    Raw_Barcodes = ncol(obj),
    
    Candidate_Cells = ncol(
      candidate_obj
    ),
    
    Median_nFeature_RNA = median(
      candidate_obj$nFeature_RNA
    ),
    
    Median_nCount_RNA = median(
      candidate_obj$nCount_RNA
    ),
    
    Median_percent_mt = median(
      candidate_obj$percent.mt,
      na.rm = TRUE
    ),
    
    stringsAsFactors = FALSE
  )
}


candidate_summary <- do.call(
  rbind,
  candidate_summary_list
)

rownames(
  candidate_summary
) <- NULL


cat("\n============================================\n")
cat("CANDIDATE CELL SUMMARY\n")
cat("NOT FINAL QC\n")
cat("============================================\n")

print(
  candidate_summary,
  row.names = FALSE
)


write.csv(
  candidate_summary,
  file.path(
    results_dir,
    "V02A_Candidate_Cell_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. CANDIDATE-CELL QC DISTRIBUTIONS
# ============================================================

candidate_qc_quantiles_list <- list()


for (sample_name in names(candidate_objects)) {
  
  obj <- candidate_objects[[sample_name]]
  
  
  q_features <- quantile(
    obj$nFeature_RNA,
    probs = c(
      0,
      0.01,
      0.05,
      0.10,
      0.25,
      0.50,
      0.75,
      0.90,
      0.95,
      0.99,
      1
    ),
    na.rm = TRUE
  )
  
  
  q_counts <- quantile(
    obj$nCount_RNA,
    probs = c(
      0,
      0.01,
      0.05,
      0.10,
      0.25,
      0.50,
      0.75,
      0.90,
      0.95,
      0.99,
      1
    ),
    na.rm = TRUE
  )
  
  
  q_mt <- quantile(
    obj$percent.mt,
    probs = c(
      0,
      0.01,
      0.05,
      0.10,
      0.25,
      0.50,
      0.75,
      0.90,
      0.95,
      0.99,
      1
    ),
    na.rm = TRUE
  )
  
  
  candidate_qc_quantiles_list[[sample_name]] <- data.frame(
    
    Sample_ID = sample_name,
    
    Quantile = names(
      q_features
    ),
    
    nFeature_RNA = as.numeric(
      q_features
    ),
    
    nCount_RNA = as.numeric(
      q_counts
    ),
    
    percent_mt = as.numeric(
      q_mt
    ),