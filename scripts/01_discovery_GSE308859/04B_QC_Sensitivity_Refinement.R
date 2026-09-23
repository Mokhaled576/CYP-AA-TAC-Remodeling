# ============================================================
# GSE308859 TAC scRNA-seq PROJECT
# STEP 04B: QC SENSITIVITY / REFINEMENT
#
# PURPOSE
# Compare:
#
# STRATEGY A = original Step-04 QC
#   Low features + high features + high counts + high MT
#   followed by scDblFinder
#
# STRATEGY B = refined QC
#   Low features + high MT are exclusion criteria
#   High features and high counts are FLAGS ONLY
#   followed by scDblFinder
#
# IMPORTANT
# - Starts again from untouched Step-03 objects
# - Does NOT overwrite Step-04 outputs
# - Runs doublet detection independently per library
# - Tracks CYP-positive cells through every stage
# ============================================================


# ============================================================
# 0. CLEAN ENVIRONMENT
# ============================================================

rm(list = ls())
gc()

set.seed(20260917)

options(stringsAsFactors = FALSE)


# ============================================================
# 1. DIRECTORIES
# ============================================================

project_dir <- "D:/Master/ScRNA seq/TAC"

clean_dir <- file.path(
  project_dir,
  "CLEAN DATA"
)

results_dir <- file.path(
  project_dir,
  "RESULTS",
  "04B_QC_REFINEMENT"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "04B_QC_REFINEMENT"
)

notes_dir <- file.path(
  project_dir,
  "NOTES"
)


for (d in c(
  clean_dir,
  results_dir,
  figures_dir,
  notes_dir
)) {
  
  if (!dir.exists(d)) {
    
    dir.create(
      d,
      recursive = TRUE
    )
  }
}


# ============================================================
# 2. PACKAGES
# ============================================================

cran_packages <- c(
  "Seurat",
  "ggplot2",
  "patchwork"
)


for (pkg in cran_packages) {
  
  if (!requireNamespace(
    pkg,
    quietly = TRUE
  )) {
    
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


# ============================================================
# 4. INPUT FILES
# ============================================================

step03_file <- file.path(
  clean_dir,
  "03_All_Samples_PreQC_Unfiltered_List.rds"
)

step04_file <- file.path(
  clean_dir,
  "04_All_Samples_PostQC_Singlets_List.rds"
)


if (!file.exists(step03_file)) {
  
  stop(
    paste0(
      "Cannot find Step-03 file:\n",
      step03_file
    )
  )
}


if (!file.exists(step04_file)) {
  
  stop(
    paste0(
      "Cannot find Step-04 final singlet file:\n",
      step04_file
    )
  )
}


# ============================================================
# 5. LOAD ORIGINAL AND STRATEGY-A OBJECTS
# ============================================================

original_list <- readRDS(
  step03_file
)

strategy_A_list <- readRDS(
  step04_file
)


sample_order <- c(
  "Sham",
  "TAC_2W",
  "TAC_4W",
  "TAC_6W"
)


missing_original <- setdiff(
  sample_order,
  names(original_list)
)


missing_A <- setdiff(
  sample_order,
  names(strategy_A_list)
)


if (length(missing_original) > 0) {
  
  stop(
    paste(
      "Missing original samples:",
      paste(
        missing_original,
        collapse = ", "
      )
    )
  )
}


if (length(missing_A) > 0) {
  
  stop(
    paste(
      "Missing Strategy-A samples:",
      paste(
        missing_A,
        collapse = ", "
      )
    )
  )
}


original_list <- original_list[
  sample_order
]

strategy_A_list <- strategy_A_list[
  sample_order
]


cat("\n========================================\n")
cat("INPUT CELL NUMBERS\n")
cat("========================================\n")


for (sample_name in sample_order) {
  
  cat(
    sample_name,
    " original = ",
    ncol(
      original_list[[sample_name]]
    ),
    " | Strategy A final = ",
    ncol(
      strategy_A_list[[sample_name]]
    ),
    "\n",
    sep = ""
  )
}


# ============================================================
# 6. LOAD EXACT STEP-04 THRESHOLDS
# ============================================================

step04_threshold_file <- file.path(
  project_dir,
  "RESULTS",
  "04_QC",
  "04_Final_QC_Thresholds.csv"
)


if (!file.exists(step04_threshold_file)) {
  
  stop(
    paste0(
      "Cannot find Step-04 threshold file:\n",
      step04_threshold_file
    )
  )
}


step04_thresholds <- read.csv(
  step04_threshold_file,
  stringsAsFactors = FALSE
)


required_threshold_columns <- c(
  "Sample",
  "Min_Features",
  "Max_Features",
  "Max_Counts",
  "Max_Percent_MT"
)


missing_threshold_columns <- setdiff(
  required_threshold_columns,
  colnames(step04_thresholds)
)


if (length(missing_threshold_columns) > 0) {
  
  stop(
    paste(
      "Missing threshold columns:",
      paste(
        missing_threshold_columns,
        collapse = ", "
      )
    )
  )
}


cat("\n========================================\n")
cat("STEP-04 THRESHOLDS IMPORTED\n")
cat("========================================\n")

print(
  step04_thresholds
)


# ============================================================
# 7. DEFINE STRATEGY B
# ============================================================

# Strategy B retains the SAME:
#
#   minimum nFeature threshold
#   mitochondrial threshold
#
# as Step 04.
#
# However:
#
#   high nFeature = FLAG ONLY
#   high nCount   = FLAG ONLY
#
# They are NOT exclusion criteria.
#
# This isolates the effect of the controversial
# high-complexity filtering.


strategy_B_flagged_list <- list()


for (sample_name in sample_order) {
  
  obj <- original_list[[sample_name]]
  
  
  thr <- step04_thresholds[
    step04_thresholds$Sample == sample_name,
    ,
    drop = FALSE
  ]
  
  
  if (nrow(thr) != 1) {
    
    stop(
      paste(
        "Could not uniquely identify threshold row for",
        sample_name
      )
    )
  }
  
  
  md <- obj@meta.data
  
  
  # ----------------------------------------------------------
  # True exclusion criteria
  # ----------------------------------------------------------
  
  md$B_fail_low_features <-
    md$nFeature_RNA <
    thr$Min_Features
  
  
  md$B_fail_high_mt <-
    md$percent.mt >
    thr$Max_Percent_MT
  
  
  # ----------------------------------------------------------
  # Diagnostic flags only
  # ----------------------------------------------------------
  
  md$B_flag_high_features <-
    md$nFeature_RNA >
    thr$Max_Features
  
  
  md$B_flag_high_counts <-
    md$nCount_RNA >
    thr$Max_Counts
  
  
  # ----------------------------------------------------------
  # Initial Strategy-B QC
  # ----------------------------------------------------------
  
  md$B_initial_fail <-
    md$B_fail_low_features |
    md$B_fail_high_mt
  
  
  md$B_initial_pass <-
    !md$B_initial_fail
  
  
  # ----------------------------------------------------------
  # Exclusion reason
  # ----------------------------------------------------------
  
  md$B_QC_reason <- "Pass_initial_QC"
  
  
  low_only <-
    md$B_fail_low_features &
    !md$B_fail_high_mt
  
  
  mt_only <-
    !md$B_fail_low_features &
    md$B_fail_high_mt
  
  
  both_low_mt <-
    md$B_fail_low_features &
    md$B_fail_high_mt
  
  
  md$B_QC_reason[
    low_only
  ] <- "Low_features"
  
  
  md$B_QC_reason[
    mt_only
  ] <- "High_mitochondrial"
  
  
  md$B_QC_reason[
    both_low_mt
  ] <- "Low_features_and_high_mitochondrial"
  
  
  obj@meta.data <- md
  
  
  strategy_B_flagged_list[[sample_name]] <- obj
}


# ============================================================
# 8. STRATEGY-B INITIAL QC SUMMARY
# ============================================================

B_initial_summary_list <- list()


for (sample_name in sample_order) {
  
  md <- strategy_B_flagged_list[[sample_name]]@meta.data
  
  
  total_cells <- nrow(
    md
  )
  
  
  retained_cells <- sum(
    md$B_initial_pass
  )
  
  
  B_initial_summary_list[[sample_name]] <- data.frame(
    
    Sample =
      sample_name,
    
    Starting_Cells =
      total_cells,
    
    Low_Features_Removed =
      sum(
        md$B_fail_low_features
      ),
    
    High_MT_Removed =
      sum(
        md$B_fail_high_mt
      ),
    
    High_Feature_Flag_Only =
      sum(
        md$B_flag_high_features
      ),
    
    High_Count_Flag_Only =
      sum(
        md$B_flag_high_counts
      ),
    
    Initial_QC_Removed =
      sum(
        md$B_initial_fail
      ),
    
    Initial_QC_Pass =
      retained_cells,
    
    Percent_Initial_Retained =
      round(
        100 *
          retained_cells /
          total_cells,
        2
      ),
    
    stringsAsFactors = FALSE
  )
}


B_initial_summary <- do.call(
  rbind,
  B_initial_summary_list
)


rownames(
  B_initial_summary
) <- NULL


cat("\n========================================\n")
cat("STRATEGY B - INITIAL QC SUMMARY\n")
cat("========================================\n")

print(
  B_initial_summary
)


write.csv(
  B_initial_summary,
  file.path(
    results_dir,
    "04B_Strategy_B_Initial_QC_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9. CREATE STRATEGY-B INITIAL-QC OBJECTS
# ============================================================

strategy_B_initial_list <- list()


for (sample_name in sample_order) {
  
  obj <- strategy_B_flagged_list[[sample_name]]
  
  
  cells_keep <- rownames(
    obj@meta.data
  )[
    obj$B_initial_pass
  ]
  
  
  strategy_B_initial_list[[sample_name]] <- subset(
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
      strategy_B_initial_list[[sample_name]]
    ),
    " after Strategy-B initial QC",
    "\n",
    sep = ""
  )
}


# ============================================================
# 10. RUN scDblFinder ON STRATEGY B
# ============================================================

B_doublet_summary_list <- list()


for (sample_name in sample_order) {
  
  cat("\n========================================\n")
  cat(
    "STRATEGY B scDblFinder: ",
    sample_name,
    "\n",
    sep = ""
  )
  cat("========================================\n")
  
  
  obj <- strategy_B_initial_list[[sample_name]]
  
  
  sce <- as.SingleCellExperiment(
    obj,
    assay = "RNA"
  )
  
  
  set.seed(20260917)
  
  
  sce <- scDblFinder(
    sce
  )
  
  
  required_columns <- c(
    "scDblFinder.score",
    "scDblFinder.class"
  )
  
  
  missing_columns <- setdiff(
    required_columns,
    colnames(
      colData(sce)
    )
  )
  
  
  if (length(missing_columns) > 0) {
    
    stop(
      paste(
        "Missing scDblFinder output for",
        sample_name
      )
    )
  }
  
  
  if (!identical(
    colnames(obj),
    colnames(sce)
  )) {
    
    stop(
      paste(
        "Cell order mismatch for",
        sample_name
      )
    )
  }
  
  
  obj$B_scDblFinder_score <-
    as.numeric(
      colData(sce)$scDblFinder.score
    )
  
  
  obj$B_scDblFinder_class <-
    as.character(
      colData(sce)$scDblFinder.class
    )
  
  
  strategy_B_initial_list[[sample_name]] <- obj
  
  
  n_singlet <- sum(
    obj$B_scDblFinder_class == "singlet"
  )
  
  
  n_doublet <- sum(
    obj$B_scDblFinder_class == "doublet"
  )
  
  
  B_doublet_summary_list[[sample_name]] <- data.frame(
    
    Sample =
      sample_name,
    
    Cells_After_Initial_QC =
      ncol(obj),
    
    Singlets =
      n_singlet,
    
    Doublets =
      n_doublet,
    
    Percent_Doublets =
      round(
        100 *
          n_doublet /
          ncol(obj),
        2
      ),
    
    stringsAsFactors = FALSE
  )
  
  
  print(
    table(
      obj$B_scDblFinder_class
    )
  )
  
  
  rm(sce)
  
  gc()
}


B_doublet_summary <- do.call(
  rbind,
  B_doublet_summary_list
)


rownames(
  B_doublet_summary
) <- NULL


cat("\n========================================\n")
cat("STRATEGY B - DOUBLET SUMMARY\n")
cat("========================================\n")

print(
  B_doublet_summary
)


write.csv(
  B_doublet_summary,
  file.path(
    results_dir,
    "04B_Strategy_B_Doublet_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. CREATE FINAL STRATEGY-B SINGLET OBJECTS
# ============================================================

strategy_B_final_list <- list()


for (sample_name in sample_order) {
  
  obj <- strategy_B_initial_list[[sample_name]]
  
  
  singlet_cells <- colnames(
    obj
  )[
    obj$B_scDblFinder_class == "singlet"
  ]
  
  
  strategy_B_final_list[[sample_name]] <- subset(
    obj,
    cells = singlet_cells
  )
  
  
  cat(
    "\n",
    sample_name,
    ": ",
    ncol(obj),
    " -> ",
    ncol(
      strategy_B_final_list[[sample_name]]
    ),
    " final Strategy-B singlets",
    "\n",
    sep = ""
  )
}


# ============================================================
# 12. FINAL STRATEGY-B SUMMARY
# ============================================================

B_final_summary_list <- list()


for (sample_name in sample_order) {
  
  starting_cells <- ncol(
    original_list[[sample_name]]
  )
  
  
  initial_pass <- ncol(
    strategy_B_initial_list[[sample_name]]
  )
  
  
  final_cells <- ncol(
    strategy_B_final_list[[sample_name]]
  )
  
  
  B_final_summary_list[[sample_name]] <- data.frame(
    
    Sample =
      sample_name,
    
    Starting_Cells =
      starting_cells,
    
    Removed_LowFeature_or_HighMT =
      starting_cells -
      initial_pass,
    
    Initial_QC_Pass =
      initial_pass,
    
    Doublets_Removed =
      initial_pass -
      final_cells,
    
    Final_Singlets =
      final_cells,
    
    Total_Removed =
      starting_cells -
      final_cells,
    
    Percent_Final_Retained =
      round(
        100 *
          final_cells /
          starting_cells,
        2
      ),
    
    stringsAsFactors = FALSE
  )
}


B_final_summary <- do.call(
  rbind,
  B_final_summary_list
)


rownames(
  B_final_summary
) <- NULL


cat("\n========================================\n")
cat("STRATEGY B - FINAL SUMMARY\n")
cat("========================================\n")

print(
  B_final_summary
)


write.csv(
  B_final_summary,
  file.path(
    results_dir,
    "04B_Strategy_B_Final_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. STRATEGY A vs STRATEGY B
# ============================================================

strategy_comparison_list <- list()


for (sample_name in sample_order) {
  
  starting_cells <- ncol(
    original_list[[sample_name]]
  )
  
  
  strategy_A_cells <- ncol(
    strategy_A_list[[sample_name]]
  )
  
  
  strategy_B_cells <- ncol(
    strategy_B_final_list[[sample_name]]
  )
  
  
  A_cell_names <- colnames(
    strategy_A_list[[sample_name]]
  )
  
  
  B_cell_names <- colnames(
    strategy_B_final_list[[sample_name]]
  )
  
  
  common_cells <- length(
    intersect(
      A_cell_names,