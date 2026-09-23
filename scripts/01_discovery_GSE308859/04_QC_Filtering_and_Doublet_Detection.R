# ============================================================
# GSE308859 TAC scRNA-seq PROJECT
# STEP 04: QC FILTERING + DOUBLET DETECTION
#
# Publication-quality workflow
#
# Samples:
#   Sham
#   TAC_2W
#   TAC_4W
#   TAC_6W
#
# INPUT:
#   03_All_Samples_PreQC_Unfiltered_List.rds
#
# OUTPUT:
#   QC thresholds
#   QC exclusion audit
#   scDblFinder doublet calls
#   Post-QC singlet Seurat objects
#   Post-QC diagnostic figures
#   CYP QC sensitivity analysis
#
# IMPORTANT:
#   - QC is performed separately for each sample.
#   - Doublet detection is performed separately for each sample.
#   - Ribosomal and hemoglobin percentages are NOT used
#     as exclusion criteria.
#   - Raw/pre-QC objects are preserved.
# ============================================================


# ============================================================
# 0. CLEAN ENVIRONMENT
# ============================================================

rm(list = ls())
gc()

set.seed(20260917)

options(stringsAsFactors = FALSE)


# ============================================================
# 1. PROJECT DIRECTORIES
# ============================================================

project_dir <- "D:/Master/ScRNA seq/TAC"

clean_dir <- file.path(
  project_dir,
  "CLEAN DATA"
)

results_dir <- file.path(
  project_dir,
  "RESULTS",
  "04_QC"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "04_QC"
)

notes_dir <- file.path(
  project_dir,
  "NOTES"
)


dirs_to_create <- c(
  clean_dir,
  results_dir,
  figures_dir,
  notes_dir
)


for (d in dirs_to_create) {
  
  if (!dir.exists(d)) {
    
    dir.create(
      d,
      recursive = TRUE
    )
  }
}


cat("\nProject directory:\n")
cat(project_dir, "\n")


# ============================================================
# 2. INSTALL / LOAD REQUIRED PACKAGES
# ============================================================

cran_packages <- c(
  "Seurat",
  "ggplot2",
  "patchwork"
)


for (pkg in cran_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    install.packages(
      pkg,
      dependencies = TRUE
    )
  }
}


if (!requireNamespace(
  "BiocManager",
  quietly = TRUE
)) {
  
  install.packages("BiocManager")
}


bioc_packages <- c(
  "SingleCellExperiment",
  "scDblFinder"
)


for (pkg in bioc_packages) {
  
  if (!requireNamespace(
    pkg,
    quietly = TRUE
  )) {
    
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
  }
}


library(Seurat)
library(ggplot2)
library(patchwork)
library(SingleCellExperiment)
library(scDblFinder)


# ============================================================
# 3. PACKAGE VERSIONS
# ============================================================

cat("\n========================================\n")
cat("PACKAGE VERSIONS\n")
cat("========================================\n")

cat(
  "R:",
  R.version.string,
  "\n"
)

cat(
  "Seurat:",
  as.character(
    packageVersion("Seurat")
  ),
  "\n"
)

cat(
  "scDblFinder:",
  as.character(
    packageVersion("scDblFinder")
  ),
  "\n"
)

cat(
  "SingleCellExperiment:",
  as.character(
    packageVersion("SingleCellExperiment")
  ),
  "\n"
)

cat(
  "ggplot2:",
  as.character(
    packageVersion("ggplot2")
  ),
  "\n"
)


# ============================================================
# 4. LOAD STEP-03 UNFILTERED OBJECTS
# ============================================================

input_file <- file.path(
  clean_dir,
  "03_All_Samples_PreQC_Unfiltered_List.rds"
)


if (!file.exists(input_file)) {
  
  stop(
    paste0(
      "Cannot find Step-03 input file:\n",
      input_file
    )
  )
}


seurat_list <- readRDS(
  input_file
)


# ============================================================
# 5. VERIFY SAMPLES
# ============================================================

sample_order <- c(
  "Sham",
  "TAC_2W",
  "TAC_4W",
  "TAC_6W"
)


missing_samples <- setdiff(
  sample_order,
  names(seurat_list)
)


if (length(missing_samples) > 0) {
  
  stop(
    paste(
      "Missing samples:",
      paste(
        missing_samples,
        collapse = ", "
      )
    )
  )
}


seurat_list <- seurat_list[
  sample_order
]


cat("\n========================================\n")
cat("SAMPLES LOADED\n")
cat("========================================\n")

print(
  names(seurat_list)
)


for (sample_name in sample_order) {
  
  cat(
    sample_name,
    ":",
    ncol(
      seurat_list[[sample_name]]
    ),
    "cells\n"
  )
}


# ============================================================
# 6. ROBUST MAD FUNCTION
# ============================================================

get_mad_thresholds <- function(x) {
  
  x <- x[
    is.finite(x)
  ]
  
  
  med <- median(
    x,
    na.rm = TRUE
  )
  
  
  mad_value <- mad(
    x,
    center = med,
    constant = 1.4826,
    na.rm = TRUE
  )
  
  
  lower <- med - (3 * mad_value)
  
  upper <- med + (3 * mad_value)
  
  
  return(
    c(
      median = med,
      MAD = mad_value,
      lower = lower,
      upper = upper
    )
  )
}


# ============================================================
# 7. CALCULATE RAW SAMPLE-SPECIFIC MAD THRESHOLDS
# ============================================================

threshold_list <- list()


for (sample_name in sample_order) {
  
  md <- seurat_list[[sample_name]]@meta.data
  
  
  feature_thr <- get_mad_thresholds(
    md$nFeature_RNA
  )
  
  
  count_thr <- get_mad_thresholds(
    md$nCount_RNA
  )
  
  
  mt_thr <- get_mad_thresholds(
    md$percent.mt
  )
  
  
  threshold_list[[sample_name]] <- data.frame(
    
    Sample =
      sample_name,
    
    Median_Features =
      as.numeric(
        feature_thr["median"]
      ),
    
    MAD_Features =
      as.numeric(
        feature_thr["MAD"]
      ),
    
    MAD_Lower_Features =
      as.numeric(
        feature_thr["lower"]
      ),
    
    MAD_Upper_Features =
      as.numeric(
        feature_thr["upper"]
      ),
    
    Median_Counts =
      as.numeric(
        count_thr["median"]
      ),
    
    MAD_Counts =
      as.numeric(
        count_thr["MAD"]
      ),
    
    MAD_Lower_Counts =
      as.numeric(
        count_thr["lower"]
      ),
    
    MAD_Upper_Counts =
      as.numeric(
        count_thr["upper"]
      ),
    
    Median_MT =
      as.numeric(
        mt_thr["median"]
      ),
    
    MAD_MT =
      as.numeric(
        mt_thr["MAD"]
      ),
    
    MAD_Upper_MT =
      as.numeric(
        mt_thr["upper"]
      ),
    
    stringsAsFactors = FALSE
  )
}


threshold_table <- do.call(
  rbind,
  threshold_list
)


rownames(
  threshold_table
) <- NULL


cat("\n========================================\n")
cat("RAW MAD THRESHOLDS\n")
cat("========================================\n")

print(
  threshold_table
)


write.csv(
  threshold_table,
  file.path(
    results_dir,
    "04_Raw_MAD_Thresholds.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 8. DEFINE HYBRID QC PRINCIPLES
# ============================================================

# We use a hybrid approach:
#
# LOW COMPLEXITY:
#   Sample-specific MAD lower boundary,
#   but never below 200 genes
#   and never above 500 genes.
#
# HIGH COMPLEXITY:
#   Candidate extreme cells are identified using the
#   smaller of:
#     1) median + 3 MAD
#     2) sample-specific 99.5th percentile
#
# HIGH UMI:
#   Same strategy:
#     smaller of median + 3 MAD and 99.5th percentile.
#
# MITOCHONDRIAL:
#   sample-specific median + 3 MAD,
#   but minimum permitted threshold = 15%
#   and maximum permitted threshold = 30%.
#
# RIBOSOMAL:
#   diagnostic only.
#
# HEMOGLOBIN:
#   diagnostic only.
#
# NOTE:
#   Doublet detection is performed AFTER the initial
#   permissive quality screen.


absolute_min_features <- 200

maximum_allowed_min_features <- 500

absolute_mt_floor <- 15

absolute_mt_ceiling <- 30


# ============================================================
# 9. CALCULATE FINAL QC THRESHOLDS
# ============================================================

final_threshold_list <- list()


for (sample_name in sample_order) {
  
  obj <- seurat_list[[sample_name]]
  
  md <- obj@meta.data
  
  
  feature_mad <- get_mad_thresholds(
    md$nFeature_RNA
  )
  
  
  count_mad <- get_mad_thresholds(
    md$nCount_RNA
  )
  
  
  mt_mad <- get_mad_thresholds(
    md$percent.mt
  )
  
  
  # ----------------------------------------------------------
  # Lower feature threshold
  # ----------------------------------------------------------
  
  min_features <- max(
    absolute_min_features,
    as.numeric(
      feature_mad["lower"]
    )
  )
  
  
  min_features <- min(
    min_features,
    maximum_allowed_min_features
  )
  
  
  min_features <- round(
    min_features
  )
  
  
  # ----------------------------------------------------------
  # Upper feature threshold
  # ----------------------------------------------------------
  
  feature_q995 <- as.numeric(
    quantile(
      md$nFeature_RNA,
      probs = 0.995,
      na.rm = TRUE,
      names = FALSE
    )
  )
  
  
  max_features <- min(
    as.numeric(
      feature_mad["upper"]
    ),
    feature_q995
  )
  
  
  max_features <- round(
    max_features
  )
  
  
  # ----------------------------------------------------------
  # Upper count threshold
  # ----------------------------------------------------------
  
  count_q995 <- as.numeric(
    quantile(
      md$nCount_RNA,
      probs = 0.995,
      na.rm = TRUE,
      names = FALSE
    )
  )
  
  
  max_counts <- min(
    as.numeric(
      count_mad["upper"]
    ),
    count_q995
  )
  
  
  max_counts <- round(
    max_counts
  )
  
  
  # ----------------------------------------------------------
  # Mitochondrial threshold
  # ----------------------------------------------------------
  
  max_mt <- max(
    absolute_mt_floor,
    as.numeric(
      mt_mad["upper"]
    )
  )
  
  
  max_mt <- min(
    max_mt,
    absolute_mt_ceiling
  )
  
  
  max_mt <- round(
    max_mt,
    digits = 2
  )
  
  
  # ----------------------------------------------------------
  # Store
  # ----------------------------------------------------------
  
  final_threshold_list[[sample_name]] <- data.frame(
    
    Sample =
      sample_name,
    
    Min_Features =
      min_features,
    
    Max_Features =
      max_features,
    
    Max_Counts =
      max_counts,
    
    Max_Percent_MT =
      max_mt,
    
    Feature_Q99_5 =
      feature_q995,
    
    Count_Q99_5 =
      count_q995,
    
    stringsAsFactors = FALSE
  )
}


final_thresholds <- do.call(
  rbind,
  final_threshold_list
)


rownames(
  final_thresholds
) <- NULL


cat("\n========================================\n")
cat("FINAL PROPOSED QC THRESHOLDS\n")
cat("========================================\n")

print(
  final_thresholds
)


write.csv(
  final_thresholds,
  file.path(
    results_dir,
    "04_Final_QC_Thresholds.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. ADD INITIAL QC FLAGS
# ============================================================

for (sample_name in sample_order) {
  
  obj <- seurat_list[[sample_name]]
  
  
  thr <- final_thresholds[
    final_thresholds$Sample == sample_name,
    ,
    drop = FALSE
  ]
  
  
  md <- obj@meta.data
  
  
  md$QC_low_features <-
    md$nFeature_RNA <
    thr$Min_Features
  
  
  md$QC_high_features <-
    md$nFeature_RNA >
    thr$Max_Features
  
  
  md$QC_high_counts <-
    md$nCount_RNA >
    thr$Max_Counts
  
  
  md$QC_high_mt <-
    md$percent.mt >
    thr$Max_Percent_MT
  
  
  md$QC_number_of_failures <-
    rowSums(
      cbind(
        md$QC_low_features,
        md$QC_high_features,
        md$QC_high_counts,
        md$QC_high_mt
      )
    )
  
  
  md$QC_initial_fail <-
    md$QC_number_of_failures > 0
  
  
  md$QC_reason <- "Pass_initial_QC"
  
  
  only_low_features <-
    md$QC_low_features &
    md$QC_number_of_failures == 1
  
  
  only_high_features <-
    md$QC_high_features &
    md$QC_number_of_failures == 1
  
  
  only_high_counts <-
    md$QC_high_counts &
    md$QC_number_of_failures == 1
  
  
  only_high_mt <-
    md$QC_high_mt &
    md$QC_number_of_failures == 1
  
  
  md$QC_reason[
    only_low_features
  ] <- "Low_features"
  
  
  md$QC_reason[
    only_high_features
  ] <- "High_features"
  
  
  md$QC_reason[
    only_high_counts
  ] <- "High_counts"
  
  
  md$QC_reason[
    only_high_mt
  ] <- "High_mitochondrial"
  
  
  md$QC_reason[
    md$QC_number_of_failures > 1
  ] <- "Multiple_QC_failures"
  
  
  obj@meta.data <- md
  
  
  seurat_list[[sample_name]] <- obj
}


# ============================================================
# 11. INITIAL QC SUMMARY
# ============================================================

initial_qc_summary_list <- list()


for (sample_name in sample_order) {
  
  md <- seurat_list[[sample_name]]@meta.data
  
  
  starting_cells <- nrow(
    md
  )
  
  
  failed_cells <- sum(
    md$QC_initial_fail
  )
  
  
  retained_cells <- sum(
    !md$QC_initial_fail
  )
  
  
  initial_qc_summary_list[[sample_name]] <- data.frame(
    
    Sample =
      sample_name,
    
    Starting_Cells =
      starting_cells,
    
    Low_Features =
      sum(
        md$QC_low_features
      ),
    
    High_Features =
      sum(
        md$QC_high_features
      ),
    
    High_Counts =
      sum(
        md$QC_high_counts
      ),
    
    High_MT =
      sum(
        md$QC_high_mt
      ),
    
    Any_Initial_QC_Failure =
      failed_cells,
    
    Initial_QC_Pass =
      retained_cells,
    
    Percent_Retained =
      round(
        100 *
          retained_cells /
          starting_cells,
        2
      ),
    
    stringsAsFactors = FALSE
  )
}


initial_qc_summary <- do.call(
  rbind,
  initial_qc_summary_list
)


rownames(
  initial_qc_summary
) <- NULL


cat("\n========================================\n")
cat("INITIAL QC SUMMARY\n")
cat("========================================\n")

print(
  initial_qc_summary
)


write.csv(
  initial_qc_summary,
  file.path(
    results_dir,
    "04_Initial_QC_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. SAVE COMPLETE INITIAL QC AUDIT
# ============================================================

initial_qc_audit_list <- list()


for (sample_name in sample_order) {
  
  md <- seurat_list[[sample_name]]@meta.data
  
  
  md$Cell_Barcode <- rownames(
    md
  )
  
  
  initial_qc_audit_list[[sample_name]] <- md
}


initial_qc_audit <- do.call(
  rbind,
  initial_qc_audit_list
)


rownames(
  initial_qc_audit
) <- NULL


write.csv(
  initial_qc_audit,
  file.path(
    results_dir,
    "04_All_Cells_Initial_QC_Audit.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. CREATE INITIAL-QC-PASS OBJECTS
# ============================================================

qc_pass_list <- list()


for (sample_name in sample_order) {
  
  obj <- seurat_list[[sample_name]]
  
  
  md <- obj@meta.data
  
  
  cells_keep <- rownames(
    md
  )[
    !md$QC_initial_fail
  ]
  
  
  qc_pass_list[[sample_name]] <- subset(
    obj,
    cells = cells_keep
  )
  
  
  cat(
    "\n",
    sample_name,
    ": ",
    ncol(obj),
    " -> ",
    ncol(
      qc_pass_list[[sample_name]]
    ),
    " cells after initial QC",
    "\n",
    sep = ""
  )
}


# ============================================================
# 14. RUN scDblFinder SEPARATELY FOR EACH SAMPLE
# ============================================================

doublet_summary_list <- list()


for (sample_name in sample_order) {
  
  cat("\n========================================\n")