# ============================================================
# V03 — GSE155882
# REPLICATE-SPECIFIC DOUBLET DETECTION
#
# INPUT:
#   V02B_GSE155882_FINAL_QC_OBJECTS.rds
#
# OUTPUTS:
#   1. V03_GSE155882_DOUBLET_ANNOTATED_OBJECTS.rds
#   2. V03_GSE155882_SINGLET_OBJECTS.rds
#   3. Replicate-level doublet summaries
#   4. Cell-level scDblFinder calls
#   5. Post-doublet CYP/AA feasibility audit
#   6. Individual 300-dpi PNG figures
#
# DESIGN:
#   - Each 10x library processed independently
#   - scDblFinder default doublet-rate estimation
#   - No manual doublet-rate tuning
#   - No CYP/AA-driven filtering
#   - Replicates remain separate
#   - No cross-library integration
#   - No biological clustering
#   - No differential expression
#
# NOTE:
#   scDblFinder may internally perform dimensional reduction
#   and clustering as part of doublet detection.
#   These are algorithmic operations only and are NOT used
#   for biological interpretation.
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(Seurat)
library(Matrix)


if (!requireNamespace(
  "SingleCellExperiment",
  quietly = TRUE
)) {
  
  stop(
    paste0(
      "Package 'SingleCellExperiment' is required.\n",
      "Install with:\n\n",
      "if (!requireNamespace('BiocManager', quietly = TRUE)) ",
      "install.packages('BiocManager')\n",
      "BiocManager::install('SingleCellExperiment')\n\n",
      "Then rerun V03."
    )
  )
}


if (!requireNamespace(
  "scDblFinder",
  quietly = TRUE
)) {
  
  stop(
    paste0(
      "Package 'scDblFinder' is required.\n",
      "Install with:\n\n",
      "if (!requireNamespace('BiocManager', quietly = TRUE)) ",
      "install.packages('BiocManager')\n",
      "BiocManager::install('scDblFinder')\n\n",
      "Then rerun V03."
    )
  )
}


if (!requireNamespace(
  "SummarizedExperiment",
  quietly = TRUE
)) {
  
  stop(
    "Package 'SummarizedExperiment' is required."
  )
}


# ============================================================
# 2. PROJECT PATHS
# ============================================================

project_dir <- "D:/Master/ScRNA seq/TAC/VALIDATION/GSE155882"


clean_dir <- file.path(
  project_dir,
  "CLEAN DATA"
)


results_dir <- file.path(
  project_dir,
  "RESULTS",
  "V03_DOUBLET_DETECTION"
)


figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "V03_DOUBLET_DETECTION"
)


png_dir <- file.path(
  figures_dir,
  "PNG"
)


dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  png_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


input_file <- file.path(
  clean_dir,
  "V02B_GSE155882_FINAL_QC_OBJECTS.rds"
)


# ============================================================
# 3. LOAD FROZEN V02B OBJECTS
# ============================================================

if (!file.exists(input_file)) {
  
  stop(
    paste0(
      "V02B object not found:\n",
      input_file
    )
  )
}


qc_objects <- readRDS(
  input_file
)


expected_samples <- c(
  "Sham_Rep1",
  "Sham_Rep2",
  "TAC_Rep1",
  "TAC_Rep2"
)


if (!identical(
  sort(names(qc_objects)),
  sort(expected_samples)
)) {
  
  stop(
    "Unexpected sample structure in V02B object."
  )
}


cat("\n============================================\n")
cat("V03 — REPLICATE-SPECIFIC DOUBLET DETECTION\n")
cat("============================================\n")


cat(
  "Cells entering V03:",
  sum(
    sapply(
      qc_objects,
      ncol
    )
  ),
  "\n"
)


cat("\nCells entering V03 by library:\n")


print(
  data.frame(
    
    Sample_ID = names(qc_objects),
    
    Cells = as.integer(
      sapply(
        qc_objects,
        ncol
      )
    ),
    
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)


# ============================================================
# 4. HELPER — RAW COUNTS
# ============================================================

get_counts <- function(obj) {
  
  GetAssayData(
    obj,
    assay = "RNA",
    layer = "counts"
  )
}


# ============================================================
# 5. PRE-DOUBLET AUDIT
# ============================================================

pre_doublet_summary_list <- list()


for (sample_name in names(qc_objects)) {
  
  obj <- qc_objects[[sample_name]]
  
  
  pre_doublet_summary_list[[sample_name]] <- data.frame(
    
    Sample_ID = sample_name,
    
    Condition = unique(
      obj$condition
    ),
    
    Replicate = unique(
      obj$replicate
    ),
    
    Cells = ncol(obj),
    
    Median_nFeature_RNA = median(
      obj$nFeature_RNA,
      na.rm = TRUE
    ),
    
    Median_nCount_RNA = median(
      obj$nCount_RNA,
      na.rm = TRUE
    ),
    
    Median_percent_mt = median(
      obj$percent.mt,
      na.rm = TRUE
    ),
    
    stringsAsFactors = FALSE
  )
}


pre_doublet_summary <- do.call(
  rbind,
  pre_doublet_summary_list
)


rownames(
  pre_doublet_summary
) <- NULL


cat("\n============================================\n")
cat("PRE-DOUBLET SUMMARY\n")
cat("============================================\n")


print(
  pre_doublet_summary,
  row.names = FALSE
)


write.csv(
  pre_doublet_summary,
  file.path(
    results_dir,
    "V03_Pre_Doublet_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 6. RUN scDblFinder
#
# IMPORTANT:
#   Each biological library is processed independently.
#
#   We deliberately do NOT specify a manual dbr.
#   scDblFinder therefore uses its standard expected-rate
#   estimation.
#
#   CYP/AA genes do not influence doublet calling.
# ============================================================

doublet_annotated_objects <- list()

doublet_call_tables <- list()

doublet_summary_list <- list()


for (sample_name in names(qc_objects)) {
  
  cat("\n============================================\n")
  
  cat(
    "Running scDblFinder:",
    sample_name,
    "\n"
  )
  
  cat("============================================\n")
  
  
  obj <- qc_objects[[sample_name]]
  
  
  counts <- get_counts(
    obj
  )
  
  
  # ----------------------------------------------------------
  # Construct SingleCellExperiment from raw counts
  # ----------------------------------------------------------
  
  sce <- SingleCellExperiment::SingleCellExperiment(
    
    assays = list(
      counts = counts
    )
  )
  
  
  # ----------------------------------------------------------
  # Explicitly preserve cell names
  # ----------------------------------------------------------
  
  colnames(sce) <- colnames(obj)
  
  
  # ----------------------------------------------------------
  # Reproducible scDblFinder run
  # ----------------------------------------------------------
  
  set.seed(20260920)
  
  
  sce <- scDblFinder::scDblFinder(
    sce
  )
  
  
  # ----------------------------------------------------------
  # Extract scDblFinder metadata
  # ----------------------------------------------------------
  
  sce_metadata <- as.data.frame(
    SummarizedExperiment::colData(
      sce
    )
  )
  
  
  required_columns <- c(
    "scDblFinder.score",
    "scDblFinder.class"
  )
  
  
  if (!all(
    required_columns %in%
    colnames(sce_metadata)
  )) {
    
    stop(
      paste(
        "Required scDblFinder output missing for",
        sample_name
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Match results back to Seurat object
  # ----------------------------------------------------------
  
  match_index <- match(
    colnames(obj),
    rownames(sce_metadata)
  )
  
  
  if (anyNA(match_index)) {
    
    stop(
      paste(
        "Could not match all scDblFinder cells for",
        sample_name
      )
    )
  }
  
  
  obj$scDblFinder_score <-
    sce_metadata$scDblFinder.score[
      match_index
    ]
  
  
  obj$scDblFinder_class <-
    as.character(
      sce_metadata$scDblFinder.class[
        match_index
      ]
    )
  
  
  # ----------------------------------------------------------
  # Verify classes
  # ----------------------------------------------------------
  
  allowed_classes <- c(
    "singlet",
    "doublet"
  )
  
  
  observed_classes <- unique(
    obj$scDblFinder_class
  )
  
  
  if (!all(
    observed_classes %in%
    allowed_classes
  )) {
    
    stop(
      paste(
        "Unexpected scDblFinder classification for",
        sample_name,
        ":",
        paste(
          observed_classes,
          collapse = ", "
        )
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Audit classifications
  # ----------------------------------------------------------
  
  class_table <- table(
    obj$scDblFinder_class
  )
  
  
  singlet_n <- sum(
    obj$scDblFinder_class ==
      "singlet"
  )
  
  
  doublet_n <- sum(
    obj$scDblFinder_class ==
      "doublet"
  )
  
  
  doublet_percent <- (
    100 *
      doublet_n /
      ncol(obj)
  )
  
  
  cat(
    "Total cells:",
    ncol(obj),
    "\n"
  )
  
  
  cat(
    "Singlets:",
    singlet_n,
    "\n"
  )
  
  
  cat(
    "Doublets:",
    doublet_n,
    "\n"
  )
  
  
  cat(
    "Doublet percentage:",
    round(
      doublet_percent,
      2
    ),
    "%\n"
  )
  
  
  print(
    class_table
  )
  
  
  # ----------------------------------------------------------
  # Cell-level call table
  # ----------------------------------------------------------
  
  call_table <- data.frame(
    
    Barcode = colnames(obj),
    
    Sample_ID = sample_name,
    
    Condition = unique(
      obj$condition
    ),
    
    Replicate = unique(
      obj$replicate
    ),
    
    nFeature_RNA = obj$nFeature_RNA,
    
    nCount_RNA = obj$nCount_RNA,
    
    percent_mt = obj$percent.mt,
    
    scDblFinder_score =
      obj$scDblFinder_score,
    
    scDblFinder_class =
      obj$scDblFinder_class,
    
    stringsAsFactors = FALSE
  )
  
  
  rownames(
    call_table
  ) <- NULL
  
  
  doublet_call_tables[[sample_name]] <-
    call_table
  
  
  write.csv(
    call_table,
    file.path(
      results_dir,
      paste0(
        "V03_",
        sample_name,
        "_scDblFinder_Cell_Calls.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Replicate summary
  # ----------------------------------------------------------
  
  doublet_summary_list[[sample_name]] <- data.frame(
    
    Sample_ID = sample_name,
    
    Condition = unique(
      obj$condition
    ),
    
    Replicate = unique(
      obj$replicate
    ),
    
    Total_Cells = ncol(obj),
    
    Singlets = singlet_n,
    
    Doublets = doublet_n,
    
    Doublet_Percent =
      doublet_percent,
    
    Median_Singlet_Score =
      median(
        obj$scDblFinder_score[
          obj$scDblFinder_class ==
            "singlet"
        ],
        na.rm = TRUE
      ),
    
    Median_Doublet_Score =
      median(
        obj$scDblFinder_score[
          obj$scDblFinder_class ==
            "doublet"
        ],
        na.rm = TRUE
      ),
    
    stringsAsFactors = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Store annotated Seurat object
  # ----------------------------------------------------------
  
  doublet_annotated_objects[[sample_name]] <-
    obj
  
  
  rm(
    sce,
    sce_metadata,
    counts
  )
  
  
  gc()
}


# ============================================================
# 7. COMBINED DOUBLET SUMMARY
# ============================================================

doublet_summary <- do.call(
  rbind,
  doublet_summary_list
)


rownames(
  doublet_summary
) <- NULL


cat("\n============================================\n")
cat("scDblFinder SUMMARY\n")
cat("============================================\n")


print(
  doublet_summary,
  row.names = FALSE
)


write.csv(
  doublet_summary,
  file.path(
    results_dir,
    "V03_scDblFinder_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 8. DOUBLET SCORE DISTRIBUTION PNGs
# ============================================================

for (sample_name in names(doublet_annotated_objects)) {
  
  obj <- doublet_annotated_objects[[sample_name]]
  
  
  output_png <- file.path(
    png_dir,
    paste0(
      "V03_",
      sample_name,
      "_Doublet_Score_Distribution.png"
    )
  )
  
  
  png(
    filename = output_png,
    width = 7,
    height = 6,
    units = "in",
    res = 300,
    bg = "white"
  )
  
  
  hist(
    obj$scDblFinder_score,
    breaks = 80,
    main = paste(
      sample_name,
      "scDblFinder Score"
    ),
    xlab = "scDblFinder score"
  )
  
  
  dev.off()
}


# ============================================================
# 9. nFEATURE vs DOUBLET SCORE PNGs
# ============================================================

for (sample_name in names(doublet_annotated_objects)) {
  
  obj <- doublet_annotated_objects[[sample_name]]
  
  
  output_png <- file.path(
    png_dir,
    paste0(
      "V03_",
      sample_name,
      "_nFeature_vs_DoubletScore.png"
    )
  )
  
  
  png(
    filename = output_png,
    width = 7,
    height = 6,
    units = "in",
    res = 300,
    bg = "white"
  )
  
  
  plot(
    obj$nFeature_RNA,
    obj$scDblFinder_score,
    pch = 16,
    cex = 0.35,
    xlab = "nFeature_RNA",
    ylab = "scDblFinder score",
    main = paste(
      sample_name,
      "Features vs Doublet Score"
    )
  )
  
  
  dev.off()
}


# ============================================================
# 10. nCOUNT vs DOUBLET SCORE PNGs
# ============================================================

for (sample_name in names(doublet_annotated_objects)) {
  
  obj <- doublet_annotated_objects[[sample_name]]
  
  
  output_png <- file.path(
    png_dir,
    paste0(
      "V03_",
      sample_name,
      "_nCount_vs_DoubletScore.png"
    )
  )
  
  
  png(
    filename = output_png,
    width = 7,
    height = 6,
    units = "in",
    res = 300,
    bg = "white"
  )
  
  
  plot(
    obj$nCount_RNA,
    obj$scDblFinder_score,
    pch = 16,
    cex = 0.35,
    xlab = "nCount_RNA",
    ylab = "scDblFinder score",
    main = paste(
      sample_name,
      "UMIs vs Doublet Score"
    )
  )
  
  
  dev.off()
}


# ============================================================
# 11. SINGLET vs DOUBLET QC COMPARISON
# ============================================================

qc_class_summary_list <- list()


for (sample_name in names(doublet_annotated_objects)) {
  
  obj <- doublet_annotated_objects[[sample_name]]
  
  
  classes <- sort(
    unique(
      obj$scDblFinder_class
    )
  )
  
  
  for (class_now in classes) {
    
    idx <- (
      obj$scDblFinder_class ==
        class_now
    )
    
    
    qc_class_summary_list[[
      paste(
        sample_name,
        class_now,
        sep = "_"
      )
    ]] <- data.frame(
      
      Sample_ID = sample_name,
      
      Condition = unique(
        obj$condition
      ),
      
      Classification = class_now,
      
      Cells = sum(idx),
      
      Median_nFeature_RNA =
        median(
          obj$nFeature_RNA[idx],
          na.rm = TRUE
        ),
      
      Median_nCount_RNA =
        median(
          obj$nCount_RNA[idx],
          na.rm = TRUE
        ),
      
      Median_percent_mt =
        median(
          obj$percent.mt[idx],
          na.rm = TRUE
        ),
      
      Median_scDblFinder_score =
        median(
          obj$scDblFinder_score[idx],
          na.rm = TRUE
        ),
      
      stringsAsFactors = FALSE
    )
  }
}


qc_class_summary <- do.call(
  rbind,
  qc_class_summary_list
)


rownames(
  qc_class_summary
) <- NULL


cat("\n============================================\n")
cat("SINGLET vs DOUBLET QC SUMMARY\n")
cat("============================================\n")


print(
  qc_class_summary,
  row.names = FALSE
)


write.csv(
  qc_class_summary,
  file.path(
    results_dir,
    "V03_Singlet_vs_Doublet_QC_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. DOUBLET RATE SUMMARY PNG
# ============================================================

output_png <- file.path(
  png_dir,
  "V03_Doublet_Rate_by_Library.png"
)


png(
  filename = output_png,
  width = 8,
  height = 6,
  units = "in",
  res = 300,
  bg = "white"
)


barplot(
  doublet_summary$Doublet_Percent,
  names.arg =
    doublet_summary$Sample_ID,
  ylim = c(
    0,
    max(
      doublet_summary$Doublet_Percent,
      na.rm = TRUE
    ) * 1.20
  ),
  ylab = "Predicted doublets (%)",
  xlab = "Library",
  main = "scDblFinder Doublet Rate",
  las = 2
)


dev.off()


# ============================================================
# 13. REMOVE PREDICTED DOUBLETS
# ============================================================

singlet_objects <- list()


for (sample_name in names(doublet_annotated_objects)) {
  