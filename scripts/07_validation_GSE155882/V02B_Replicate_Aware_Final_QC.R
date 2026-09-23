# ============================================================
# V02B — GSE155882
# REPLICATE-AWARE FINAL QC
#
# INPUT:
#   V02A candidate-cell objects
#
# PURPOSE:
#   1. Recalculate QC metrics
#   2. Inspect each biological library independently
#   3. Apply conservative, predefined QC rules
#   4. Remove low-complexity cells
#   5. Remove extreme mitochondrial cells
#   6. Flag/remove extreme upper-complexity cells
#   7. Preserve biological replicates separately
#   8. Freeze final QC objects
#
# IMPORTANT:
#   - NO CYP/AA genes used for QC
#   - NO normalization
#   - NO integration
#   - NO clustering
#   - NO differential expression
#   - NO biological interpretation
#
# FINAL QC RULES:
#
# Lower complexity:
#   nFeature_RNA >= 200
#   nCount_RNA   >= 300
#
# Mitochondrial:
#   percent.mt <= sample-specific
#   median + 3*MAD
#
#   BUT:
#     minimum permitted threshold = 15%
#     maximum permitted threshold = 30%
#
# Upper complexity:
#   nFeature_RNA <= median + 5*MAD
#   nCount_RNA   <= median + 5*MAD
#
#   These are conservative extreme-tail rules.
#
# NOTE:
#   Doublet detection is NOT performed here.
#   It will be handled separately after QC.
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(Seurat)
library(Matrix)


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
  "V02B_FINAL_QC"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "V02B_FINAL_QC"
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
  "V02A_GSE155882_CANDIDATE_CELL_OBJECTS.rds"
)


# ============================================================
# 3. LOAD V02A
# ============================================================

if (!file.exists(input_file)) {
  
  stop(
    paste0(
      "V02A object not found:\n",
      input_file
    )
  )
}


candidate_objects <- readRDS(
  input_file
)


expected_samples <- c(
  "Sham_Rep1",
  "Sham_Rep2",
  "TAC_Rep1",
  "TAC_Rep2"
)


if (!identical(
  sort(names(candidate_objects)),
  sort(expected_samples)
)) {
  
  stop(
    "Unexpected sample structure in V02A object."
  )
}


cat("\n============================================\n")
cat("V02B — REPLICATE-AWARE FINAL QC\n")
cat("============================================\n")

cat(
  "Candidate cells entering V02B:",
  sum(
    sapply(
      candidate_objects,
      ncol
    )
  ),
  "\n"
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
# 5. RECALCULATE QC METRICS
# ============================================================

for (sample_name in names(candidate_objects)) {
  
  obj <- candidate_objects[[sample_name]]
  
  
  obj[["percent.mt"]] <- PercentageFeatureSet(
    obj,
    pattern = "^mt-"
  )
  
  
  obj[["percent.ribo"]] <- PercentageFeatureSet(
    obj,
    pattern = "^Rp[sl]"
  )
  
  
  obj$log10GenesPerUMI <- log10(
    obj$nFeature_RNA + 1
  ) / log10(
    obj$nCount_RNA + 1
  )
  
  
  candidate_objects[[sample_name]] <- obj
}


# ============================================================
# 6. PRE-QC SUMMARY
# ============================================================

pre_qc_summary_list <- list()


for (sample_name in names(candidate_objects)) {
  
  obj <- candidate_objects[[sample_name]]
  
  
  pre_qc_summary_list[[sample_name]] <- data.frame(
    
    Sample_ID = sample_name,
    
    Condition = unique(
      obj$condition
    ),
    
    Replicate = unique(
      obj$replicate
    ),
    
    Cells = ncol(obj),
    
    Median_nFeature_RNA = median(
      obj$nFeature_RNA
    ),
    
    Median_nCount_RNA = median(
      obj$nCount_RNA
    ),
    
    Median_percent_mt = median(
      obj$percent.mt
    ),
    
    Median_percent_ribo = median(
      obj$percent.ribo
    ),
    
    Median_log10GenesPerUMI = median(
      obj$log10GenesPerUMI
    ),
    
    stringsAsFactors = FALSE
  )
}


pre_qc_summary <- do.call(
  rbind,
  pre_qc_summary_list
)

rownames(pre_qc_summary) <- NULL


cat("\n============================================\n")
cat("PRE-QC SUMMARY\n")
cat("============================================\n")

print(
  pre_qc_summary,
  row.names = FALSE
)


write.csv(
  pre_qc_summary,
  file.path(
    results_dir,
    "V02B_Pre_QC_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 7. ROBUST THRESHOLD HELPER
# ============================================================

robust_upper <- function(
    x,
    nmads = 3
) {
  
  median(
    x,
    na.rm = TRUE
  ) +
    nmads *
    mad(
      x,
      center = median(
        x,
        na.rm = TRUE
      ),
      constant = 1.4826,
      na.rm = TRUE
    )
}


# ============================================================
# 8. DEFINE QC THRESHOLDS
#
# Lower floors are identical across samples.
#
# Upper thresholds are calculated independently for
# each biological library.
# ============================================================

threshold_list <- list()


for (sample_name in names(candidate_objects)) {
  
  obj <- candidate_objects[[sample_name]]
  
  
  mt_raw_threshold <- robust_upper(
    obj$percent.mt,
    nmads = 3
  )
  
  
  mt_final_threshold <- min(
    max(
      mt_raw_threshold,
      15
    ),
    30
  )
  
  
  feature_upper <- robust_upper(
    obj$nFeature_RNA,
    nmads = 5
  )
  
  
  count_upper <- robust_upper(
    obj$nCount_RNA,
    nmads = 5
  )
  
  
  threshold_list[[sample_name]] <- data.frame(
    
    Sample_ID = sample_name,
    
    Condition = unique(
      obj$condition
    ),
    
    Replicate = unique(
      obj$replicate
    ),
    
    nFeature_Lower = 200,
    
    nCount_Lower = 300,
    
    nFeature_Upper = feature_upper,
    
    nCount_Upper = count_upper,
    
    MT_Raw_MedianPlus3MAD = mt_raw_threshold,
    
    MT_Final_Threshold = mt_final_threshold,
    
    stringsAsFactors = FALSE
  )
}


qc_thresholds <- do.call(
  rbind,
  threshold_list
)

rownames(qc_thresholds) <- NULL


cat("\n============================================\n")
cat("REPLICATE-SPECIFIC QC THRESHOLDS\n")
cat("============================================\n")

print(
  qc_thresholds,
  row.names = FALSE
)


write.csv(
  qc_thresholds,
  file.path(
    results_dir,
    "V02B_QC_Thresholds.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9. AUDIT EACH FILTER SEPARATELY
# ============================================================

filter_audit_list <- list()


for (sample_name in names(candidate_objects)) {
  
  obj <- candidate_objects[[sample_name]]
  
  
  threshold_row <- qc_thresholds[
    qc_thresholds$Sample_ID == sample_name,
    ,
    drop = FALSE
  ]
  
  
  low_feature <- (
    obj$nFeature_RNA <
      threshold_row$nFeature_Lower
  )
  
  
  low_count <- (
    obj$nCount_RNA <
      threshold_row$nCount_Lower
  )
  
  
  high_feature <- (
    obj$nFeature_RNA >
      threshold_row$nFeature_Upper
  )
  
  
  high_count <- (
    obj$nCount_RNA >
      threshold_row$nCount_Upper
  )
  
  
  high_mt <- (
    obj$percent.mt >
      threshold_row$MT_Final_Threshold
  )
  
  
  keep <- (
    !low_feature &
      !low_count &
      !high_feature &
      !high_count &
      !high_mt
  )
  
  
  filter_audit_list[[sample_name]] <- data.frame(
    
    Sample_ID = sample_name,
    
    Condition = unique(
      obj$condition
    ),
    
    Starting_Cells = ncol(obj),
    
    Fail_Low_Features = sum(
      low_feature
    ),
    
    Fail_Low_Counts = sum(
      low_count
    ),
    
    Fail_High_Features = sum(
      high_feature
    ),
    
    Fail_High_Counts = sum(
      high_count
    ),
    
    Fail_High_MT = sum(
      high_mt
    ),
    
    Fail_Any_QC = sum(
      !keep
    ),
    
    Pass_All_QC = sum(
      keep
    ),
    
    Percent_Retained = 100 * mean(
      keep
    ),
    
    stringsAsFactors = FALSE
  )
}


filter_audit <- do.call(
  rbind,
  filter_audit_list
)

rownames(filter_audit) <- NULL


cat("\n============================================\n")
cat("QC FILTER AUDIT\n")
cat("============================================\n")

print(
  filter_audit,
  row.names = FALSE
)


write.csv(
  filter_audit,
  file.path(
    results_dir,
    "V02B_Filter_Audit.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. CREATE FINAL QC OBJECTS
# ============================================================

qc_objects <- list()


for (sample_name in names(candidate_objects)) {
  
  obj <- candidate_objects[[sample_name]]
  
  
  threshold_row <- qc_thresholds[
    qc_thresholds$Sample_ID == sample_name,
    ,
    drop = FALSE
  ]
  
  
  keep <- (
    obj$nFeature_RNA >=
      threshold_row$nFeature_Lower &
      obj$nCount_RNA >=
      threshold_row$nCount_Lower &
      obj$nFeature_RNA <=
      threshold_row$nFeature_Upper &
      obj$nCount_RNA <=
      threshold_row$nCount_Upper &
      obj$percent.mt <=
      threshold_row$MT_Final_Threshold
  )
  
  
  cells_keep <- colnames(obj)[
    keep
  ]
  
  
  qc_obj <- subset(
    obj,
    cells = cells_keep
  )
  
  
  qc_obj$V02B_QC_status <- "PASS"
  
  qc_obj$V02B_QC_version <-
    "ReplicateAware_QC_v1"
  
  
  qc_objects[[sample_name]] <- qc_obj
}


# ============================================================
# 11. POST-QC SUMMARY
# ============================================================

post_qc_summary_list <- list()


for (sample_name in names(qc_objects)) {
  
  obj <- qc_objects[[sample_name]]
  
  
  post_qc_summary_list[[sample_name]] <- data.frame(
    
    Sample_ID = sample_name,
    
    Condition = unique(
      obj$condition
    ),
    
    Replicate = unique(
      obj$replicate
    ),
    
    Cells = ncol(obj),
    
    Median_nFeature_RNA = median(
      obj$nFeature_RNA
    ),
    
    Median_nCount_RNA = median(
      obj$nCount_RNA
    ),
    
    Median_percent_mt = median(
      obj$percent.mt
    ),
    
    Median_percent_ribo = median(
      obj$percent.ribo
    ),
    
    Median_log10GenesPerUMI = median(
      obj$log10GenesPerUMI
    ),
    
    stringsAsFactors = FALSE
  )
}


post_qc_summary <- do.call(
  rbind,
  post_qc_summary_list
)

rownames(post_qc_summary) <- NULL


cat("\n============================================\n")
cat("POST-QC SUMMARY\n")
cat("============================================\n")

print(
  post_qc_summary,
  row.names = FALSE
)


write.csv(
  post_qc_summary,
  file.path(
    results_dir,
    "V02B_Post_QC_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. PRE-vs-POST CELL COUNTS
# ============================================================

cell_retention <- merge(
  pre_qc_summary[
    ,
    c(
      "Sample_ID",
      "Condition",
      "Cells"
    )
  ],
  post_qc_summary[
    ,
    c(
      "Sample_ID",
      "Cells"
    )
  ],
  by = "Sample_ID",
  suffixes = c(
    "_Before",
    "_After"
  )
)


cell_retention$Removed <- (
  cell_retention$Cells_Before -
    cell_retention$Cells_After
)


cell_retention$Percent_Retained <- (
  100 *
    cell_retention$Cells_After /
    cell_retention$Cells_Before
)


cat("\n============================================\n")
cat("CELL RETENTION\n")
cat("============================================\n")

print(
  cell_retention,
  row.names = FALSE
)


write.csv(
  cell_retention,
  file.path(
    results_dir,
    "V02B_Cell_Retention.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. POST-QC QUANTILES
# ============================================================

post_qc_quantiles_list <- list()


for (sample_name in names(qc_objects)) {
  
  obj <- qc_objects[[sample_name]]
  
  
  probs_now <- c(
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
  )
  
  
  q_feature <- quantile(
    obj$nFeature_RNA,
    probs = probs_now,
    na.rm = TRUE
  )
  
  
  q_count <- quantile(
    obj$nCount_RNA,
    probs = probs_now,
    na.rm = TRUE
  )
  
  
  q_mt <- quantile(
    obj$percent.mt,
    probs = probs_now,
    na.rm = TRUE
  )
  
  
  post_qc_quantiles_list[[sample_name]] <- data.frame(
    
    Sample_ID = sample_name,
    
    Quantile = names(
      q_feature
    ),
    
    nFeature_RNA = as.numeric(
      q_feature
    ),
    
    nCount_RNA = as.numeric(
      q_count
    ),
    
    percent_mt = as.numeric(
      q_mt
    ),
    
    stringsAsFactors = FALSE
  )
}


post_qc_quantiles <- do.call(
  rbind,
  post_qc_quantiles_list
)

rownames(post_qc_quantiles) <- NULL


write.csv(
  post_qc_quantiles,
  file.path(
    results_dir,
    "V02B_Post_QC_Quantiles.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. QC PLOTS — PRE vs POST
# ============================================================

pdf(
  file.path(
    figures_dir,
    "V02B_Pre_Post_QC_Distributions.pdf"
  ),
  width = 12,
  height = 8
)


for (sample_name in names(candidate_objects)) {
  
  pre_obj <- candidate_objects[[sample_name]]
  
  post_obj <- qc_objects[[sample_name]]
  
  
  old_par <- par(
    mfrow = c(
      2,
      3
    )
  )
  
  
  hist(
    pre_obj$nFeature_RNA,
    breaks = 80,
    main = paste(
      sample_name,
      "PRE — Features"
    ),
    xlab = "nFeature_RNA"
  )
  
  
  hist(
    pre_obj$nCount_RNA,
    breaks = 80,
    main = paste(
      sample_name,
      "PRE — UMIs"
    ),
    xlab = "nCount_RNA"
  )
  
  
  hist(
    pre_obj$percent.mt,
    breaks = 80,
    main = paste(
      sample_name,
      "PRE — mt%"
    ),
    xlab = "percent.mt"
  )
  
  
  hist(
    post_obj$nFeature_RNA,
    breaks = 80,
    main = paste(
      sample_name,
      "POST — Features"
    ),
    xlab = "nFeature_RNA"
  )
  
  
  hist(
    post_obj$nCount_RNA,
    breaks = 80,
    main = paste(
      sample_name,
      "POST — UMIs"
    ),
    xlab = "nCount_RNA"
  )
  
  
  hist(
    post_obj$percent.mt,
    breaks = 80,
    main = paste(
      sample_name,
      "POST — mt%"
    ),
    xlab = "percent.mt"
  )
  
  
  par(
    old_par
  )
}


dev.off()


# ============================================================
# 15. FEATURE-vs-COUNT QC PLOTS
# ============================================================

pdf(
  file.path(
    figures_dir,
    "V02B_Post_QC_Feature_vs_Count.pdf"
  ),
  width = 7,
  height = 6
)


for (sample_name in names(qc_objects)) {
  
  obj <- qc_objects[[sample_name]]
  
  
  plot(
    obj$nCount_RNA,
    obj$nFeature_RNA,
    log = "xy",
    pch = 16,
    cex = 0.35,
    xlab = "nCount_RNA",
    ylab = "nFeature_RNA",
    main = paste(
      sample_name,
      "Post-QC"
    )
  )
}


dev.off()


# ============================================================
# 16. POST-QC CYP/AA FEASIBILITY AUDIT
#
# IMPORTANT:
# Descriptive only.
#
# These genes did NOT influence QC thresholds.
# ============================================================

priority_genes <- c(
  "Cyp1b1",
  "Cyp2j6",
  "Cyp2j9",
  "Cyp4f13",
  "Cyp4f16",
  "Cyp4f17",
  "Cyp4f18",
  "Ephx2"
)


priority_detection_list <- list()


for (sample_name in names(qc_objects)) {
  
  obj <- qc_objects[[sample_name]]
  
  counts <- get_counts(
    obj
  )
  