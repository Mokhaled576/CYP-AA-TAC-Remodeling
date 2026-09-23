# ============================================================
# STEP 11B
# CONSERVATIVE SPATIAL QC FILTERING + NORMALIZATION
# GSE308859
#
# PURPOSE
# ------------------------------------------------------------
# 1. Load the raw spatial objects created in Step 11A.
# 2. Recalculate/verify spatial QC metrics.
# 3. Evaluate several lower-tail QC strategies.
# 4. Use the 1st-percentile rule as the primary QC rule.
# 5. Require BOTH:
#       nFeature_Spatial >= sample-specific 1st percentile
#       nCount_Spatial   >= sample-specific 1st percentile
# 6. Do NOT filter using mitochondrial percentage.
# 7. Do NOT impose upper nFeature/nCount cutoffs.
# 8. Map retained/removed spots spatially.
# 9. Create filtered spatial objects.
# 10. LogNormalize each filtered spatial object.
# 11. Save normalized objects and QC audit tables.
#
# IMPORTANT
# ------------------------------------------------------------
# - No CYP/eicosanoid biology is used to determine QC.
# - No cell-type transfer/deconvolution is performed here.
# - No condition-level statistical testing is performed.
# - High mitochondrial percentage is retained as QC metadata.
# - High-complexity spots are NOT treated as doublets.
# ============================================================


# ============================================================
# 1. CLEAN WORKSPACE
# ============================================================

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
set.seed(12345)


# ============================================================
# 2. LOAD PACKAGES
# ============================================================

required_packages <- c(
  "Seurat",
  "Matrix",
  "dplyr",
  "tidyr",
  "ggplot2",
  "patchwork"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})


# ============================================================
# 3. PROJECT PATHS
# ============================================================

project_dir <- "D:/Master/ScRNA seq/TAC"

input_file <- file.path(
  project_dir,
  "CLEAN DATA",
  "11A_ALL_RAW_SPATIAL_OBJECTS.rds"
)

results_dir <- file.path(
  project_dir,
  "RESULTS",
  "11_SPATIAL",
  "11B_QC_NORMALIZATION"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "11_SPATIAL",
  "11B_QC_NORMALIZATION"
)

clean_dir <- file.path(
  project_dir,
  "CLEAN DATA"
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

dir.create(
  clean_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 4. VERIFY INPUT
# ============================================================

if (!file.exists(input_file)) {
  stop(
    paste(
      "Input file not found:",
      input_file
    )
  )
}

spatial_objects <- readRDS(input_file)

if (!is.list(spatial_objects)) {
  stop(
    "11A_ALL_RAW_SPATIAL_OBJECTS.rds is not a list."
  )
}

cat("\nLoaded spatial object list.\n")
cat("Number of objects:", length(spatial_objects), "\n")
cat("Object names:\n")
print(names(spatial_objects))


# ============================================================
# 5. STANDARDIZE SAMPLE ORDER
# ============================================================

sample_order <- c(
  "Sham",
  "TAC_2W",
  "TAC_4W",
  "TAC_6W"
)

missing_samples <- setdiff(
  sample_order,
  names(spatial_objects)
)

if (length(missing_samples) > 0) {
  stop(
    paste(
      "Missing expected spatial samples:",
      paste(missing_samples, collapse = ", ")
    )
  )
}

spatial_objects <- spatial_objects[sample_order]


# ============================================================
# 6. VERIFY SPATIAL ASSAY
# ============================================================

for (sid in sample_order) {
  
  sobj <- spatial_objects[[sid]]
  
  assay_names <- names(sobj@assays)
  
  if (!"Spatial" %in% assay_names) {
    stop(
      paste(
        "Spatial assay missing from sample:",
        sid
      )
    )
  }
  
  DefaultAssay(sobj) <- "Spatial"
  
  spatial_objects[[sid]] <- sobj
}

cat("\nSpatial assay verified for all samples.\n")


# ============================================================
# 7. VERIFY / RECALCULATE QC METRICS
# ============================================================

for (sid in sample_order) {
  
  cat(
    "\nCalculating QC metrics:",
    sid,
    "\n"
  )
  
  sobj <- spatial_objects[[sid]]
  
  DefaultAssay(sobj) <- "Spatial"
  
  sobj[["percent.mt"]] <- PercentageFeatureSet(
    sobj,
    pattern = "^mt-",
    assay = "Spatial"
  )
  
  sobj[["percent.ribo"]] <- PercentageFeatureSet(
    sobj,
    pattern = "^Rp[sl]",
    assay = "Spatial"
  )
  
  sobj[["percent.hb"]] <- PercentageFeatureSet(
    sobj,
    pattern = "^Hb[ab]",
    assay = "Spatial"
  )
  
  sobj$sample_id_11B <- sid
  
  spatial_objects[[sid]] <- sobj
}


# ============================================================
# 8. EXTRACT COMPLETE QC METADATA
# ============================================================

qc_list <- vector(
  mode = "list",
  length = length(sample_order)
)

names(qc_list) <- sample_order

for (sid in sample_order) {
  
  sobj <- spatial_objects[[sid]]
  
  md <- sobj@meta.data
  
  required_columns <- c(
    "nFeature_Spatial",
    "nCount_Spatial",
    "percent.mt",
    "percent.ribo",
    "percent.hb"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(md)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      paste(
        "Missing QC columns in",
        sid,
        ":",
        paste(missing_columns, collapse = ", ")
      )
    )
  }
  
  temp <- data.frame(
    Barcode = rownames(md),
    Sample = sid,
    nFeature_Spatial = md$nFeature_Spatial,
    nCount_Spatial = md$nCount_Spatial,
    percent.mt = md$percent.mt,
    percent.ribo = md$percent.ribo,
    percent.hb = md$percent.hb,
    stringsAsFactors = FALSE
  )
  
  if ("in_tissue" %in% colnames(md)) {
    temp$in_tissue <- md$in_tissue
  } else {
    temp$in_tissue <- 1
  }
  
  coordinate_columns <- c(
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres"
  )
  
  for (coord_col in coordinate_columns) {
    
    if (coord_col %in% colnames(md)) {
      temp[[coord_col]] <- md[[coord_col]]
    } else {
      temp[[coord_col]] <- NA_real_
    }
  }
  
  qc_list[[sid]] <- temp
}

all_qc <- bind_rows(qc_list)

all_qc$Sample <- factor(
  all_qc$Sample,
  levels = sample_order
)

cat(
  "\nTotal raw spatial spots:",
  nrow(all_qc),
  "\n"
)


# ============================================================
# 9. DEFINE SENSITIVITY THRESHOLDS
# ============================================================
#
# We compare:
#   0.5%
#   1%
#   2%
#   5%
#
# Primary analysis:
#   1%
#
# A spot passes a rule only if it satisfies BOTH:
#   nFeature >= threshold
#   nCount   >= threshold
#
# No percent.mt threshold.
# No upper threshold.
# ============================================================

tail_probabilities <- c(
  0.005,
  0.01,
  0.02,
  0.05
)

tail_names <- c(
  "P0.5",
  "P1",
  "P2",
  "P5"
)


# ============================================================
# 10. CALCULATE SAMPLE-SPECIFIC THRESHOLDS
# ============================================================

threshold_list <- list()

counter <- 1

for (sid in sample_order) {
  
  sample_qc <- all_qc[
    all_qc$Sample == sid,
    ,
    drop = FALSE
  ]
  
  for (i in seq_along(tail_probabilities)) {
    
    prob_i <- tail_probabilities[i]
    rule_i <- tail_names[i]
    
    feature_threshold <- as.numeric(
      quantile(
        sample_qc$nFeature_Spatial,
        probs = prob_i,
        na.rm = TRUE,
        names = FALSE,
        type = 7
      )
    )
    
    count_threshold <- as.numeric(
      quantile(
        sample_qc$nCount_Spatial,
        probs = prob_i,
        na.rm = TRUE,
        names = FALSE,
        type = 7
      )
    )
    
    threshold_list[[counter]] <- data.frame(
      Sample = sid,
      Rule = rule_i,
      Tail_Probability = prob_i,
      Min_nFeature = feature_threshold,
      Min_nCount = count_threshold,
      stringsAsFactors = FALSE
    )
    
    counter <- counter + 1
  }
}

threshold_table <- bind_rows(
  threshold_list
)

threshold_table$Sample <- factor(
  threshold_table$Sample,
  levels = sample_order
)

threshold_table$Rule <- factor(
  threshold_table$Rule,
  levels = tail_names
)

write.csv(
  threshold_table,
  file.path(
    results_dir,
    "11B_Sample_Specific_QC_Thresholds.csv"
  ),
  row.names = FALSE
)

cat("\nQC threshold table:\n")
print(
  as.data.frame(threshold_table),
  row.names = FALSE
)


# ============================================================
# 11. SENSITIVITY ANALYSIS
# ============================================================

sensitivity_list <- list()

counter <- 1

for (sid in sample_order) {
  
  sample_qc <- all_qc[
    all_qc$Sample == sid,
    ,
    drop = FALSE
  ]
  
  for (rule_i in tail_names) {
    
    threshold_row <- threshold_table[
      threshold_table$Sample == sid &
        threshold_table$Rule == rule_i,
      ,
      drop = FALSE
    ]
    
    min_feature <- threshold_row$Min_nFeature[1]
    min_count <- threshold_row$Min_nCount[1]
    
    pass_feature <- (
      sample_qc$nFeature_Spatial >= min_feature
    )
    
    pass_count <- (
      sample_qc$nCount_Spatial >= min_count
    )
    
    pass_both <- (
      pass_feature &
        pass_count
    )
    
    total_spots <- nrow(sample_qc)
    retained_spots <- sum(pass_both, na.rm = TRUE)
    removed_spots <- total_spots - retained_spots
    
    sensitivity_list[[counter]] <- data.frame(
      Sample = sid,
      Rule = rule_i,
      Total_Spots = total_spots,
      Retained_Spots = retained_spots,
      Removed_Spots = removed_spots,
      Percent_Retained = 100 * retained_spots / total_spots,
      Percent_Removed = 100 * removed_spots / total_spots,
      Min_nFeature = min_feature,
      Min_nCount = min_count,
      stringsAsFactors = FALSE
    )
    
    counter <- counter + 1
  }
}

sensitivity_table <- bind_rows(
  sensitivity_list
)

sensitivity_table$Sample <- factor(
  sensitivity_table$Sample,
  levels = sample_order
)

sensitivity_table$Rule <- factor(
  sensitivity_table$Rule,
  levels = tail_names
)

write.csv(
  sensitivity_table,
  file.path(
    results_dir,
    "11B_QC_Sensitivity_Analysis.csv"
  ),
  row.names = FALSE
)

cat("\nQC sensitivity analysis:\n")
print(
  as.data.frame(sensitivity_table),
  row.names = FALSE
)


# ============================================================
# 12. SELECT PRIMARY QC RULE
# ============================================================

primary_rule <- "P1"

primary_thresholds <- threshold_table[
  threshold_table$Rule == primary_rule,
  ,
  drop = FALSE
]

cat("\nPrimary QC rule:", primary_rule, "\n")
cat("1st-percentile lower-tail filtering.\n")
cat("Both nFeature and nCount must pass.\n")
cat("No mitochondrial cutoff.\n")
cat("No upper complexity cutoff.\n")


# ============================================================
# 13. APPLY PRIMARY QC FLAGS
# ============================================================

primary_qc_list <- list()

for (sid in sample_order) {
  
  sample_qc <- all_qc[
    all_qc$Sample == sid,
    ,
    drop = FALSE
  ]
  
  threshold_row <- primary_thresholds[
    primary_thresholds$Sample == sid,
    ,
    drop = FALSE
  ]
  
  min_feature <- threshold_row$Min_nFeature[1]
  min_count <- threshold_row$Min_nCount[1]
  
  sample_qc$Pass_nFeature <- (
    sample_qc$nFeature_Spatial >= min_feature
  )
  
  sample_qc$Pass_nCount <- (
    sample_qc$nCount_Spatial >= min_count
  )
  
  sample_qc$Pass_In_Tissue <- (
    is.na(sample_qc$in_tissue) |
      sample_qc$in_tissue == 1
  )
  
  sample_qc$QC_Pass <- (
    sample_qc$Pass_nFeature &
      sample_qc$Pass_nCount &
      sample_qc$Pass_In_Tissue
  )
  
  sample_qc$QC_Status <- ifelse(
    sample_qc$QC_Pass,
    "Retained",
    "Removed"
  )
  
  sample_qc$QC_Failure_Reason <- "Pass"
  
  sample_qc$QC_Failure_Reason[
    !sample_qc$Pass_nFeature &
      sample_qc$Pass_nCount
  ] <- "Low_nFeature"
  
  sample_qc$QC_Failure_Reason[
    sample_qc$Pass_nFeature &
      !sample_qc$Pass_nCount
  ] <- "Low_nCount"
  
  sample_qc$QC_Failure_Reason[
    !sample_qc$Pass_nFeature &
      !sample_qc$Pass_nCount
  ] <- "Low_nFeature_and_nCount"
  
  sample_qc$QC_Failure_Reason[
    !sample_qc$Pass_In_Tissue
  ] <- "Not_in_tissue"
  
  primary_qc_list[[sid]] <- sample_qc
}

primary_qc <- bind_rows(
  primary_qc_list
)

primary_qc$Sample <- factor(
  primary_qc$Sample,
  levels = sample_order
)

write.csv(
  primary_qc,
  file.path(
    results_dir,
    "11B_All_Spots_Primary_QC_Status.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. PRIMARY FILTERING SUMMARY
# ============================================================

primary_summary <- primary_qc %>%
  group_by(Sample) %>%
  summarise(
    Raw_Spots = n(),
    Retained_Spots = sum(QC_Pass),
    Removed_Spots = sum(!QC_Pass),
    Percent_Retained = 100 * mean(QC_Pass),
    Percent_Removed = 100 * mean(!QC_Pass),
    Median_nFeature_Before = median(
      nFeature_Spatial,
      na.rm = TRUE
    ),
    Median_nFeature_After = median(
      nFeature_Spatial[QC_Pass],
      na.rm = TRUE
    ),
    Median_nCount_Before = median(
      nCount_Spatial,
      na.rm = TRUE
    ),
    Median_nCount_After = median(
      nCount_Spatial[QC_Pass],
      na.rm = TRUE
    ),
    Median_percent_mt_Before = median(
      percent.mt,
      na.rm = TRUE
    ),
    Median_percent_mt_After = median(
      percent.mt[QC_Pass],
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write.csv(
  primary_summary,
  file.path(
    results_dir,
    "11B_Primary_QC_Filtering_Summary.csv"
  ),
  row.names = FALSE
)

cat("\nPrimary filtering summary:\n")
print(
  as.data.frame(primary_summary),
  row.names = FALSE
)


# ============================================================
# 15. FAILURE REASON SUMMARY
# ============================================================

failure_summary <- primary_qc %>%
  count(
    Sample,
    QC_Failure_Reason,
    name = "Spots"
  ) %>%
  group_by(Sample) %>%
  mutate(
    Percent = 100 * Spots / sum(Spots)
  ) %>%
  ungroup()

write.csv(
  failure_summary,
  file.path(
    results_dir,
    "11B_QC_Failure_Reasons.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. FIGURE - SENSITIVITY ANALYSIS
# ============================================================

p_sensitivity <- ggplot(
  sensitivity_table,
  aes(
    x = Rule,
    y = Percent_Removed,
    group = Sample
  )
) +
  geom_line(
    linewidth = 0.8
  ) +
  geom_point(
    size = 2.5
  ) +
  facet_wrap(
    ~ Sample,
    nrow = 1
  ) +
  labs(
    title = "Sensitivity of spatial QC to lower-tail threshold",
    subtitle = paste0(
      "Spot must pass both nFeature and nCount; ",
      "mitochondrial percentage is not used for filtering"
    ),
    x = "Sample-specific lower-tail rule",
    y = "Spots removed (%)"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    strip.text = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

ggsave(
  filename = file.path(
    figures_dir,
    "11B_QC_Sensitivity_Analysis.png"
  ),
  plot = p_sensitivity,
  width = 12,
  height = 4.5,
  dpi = 400
)


# ============================================================
# 17. FIGURE - QC STATUS IN COMPLEXITY SPACE
# ============================================================

p_complexity_status <- ggplot(
  primary_qc,
  aes(
    x = nCount_Spatial,
    y = nFeature_Spatial,
    shape = QC_Status
  )
) +
  geom_point(
    alpha = 0.55,
    size = 1
  ) +
  facet_wrap(
    ~ Sample,
    scales = "free",
    ncol = 2
  ) +
  labs(
    title = "Spatial QC classification in library-complexity space",
    subtitle = paste0(
      "Primary rule: sample-specific 1st percentile ",
      "for both UMI and detected genes"
    ),
    x = "UMI count",
    y = "Detected genes",
    shape = "QC status"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    strip.text = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

ggsave(
  filename = file.path(
    figures_dir,
    "11B_QC_Status_Complexity_Space.png"
  ),
  plot = p_complexity_status,
  width = 10,
  height = 8,
  dpi = 400
)


# ============================================================
# 18. FIGURE - SPATIAL LOCATION OF REMOVED SPOTS
# ============================================================
#
# IMPORTANT:
# Do NOT use coord_fixed() together with free facet scales.
# This avoids the error encountered in Step 11A.
# ============================================================

coordinate_qc <- primary_qc %>%
  filter(
    !is.na(pxl_col_in_fullres),
    !is.na(pxl_row_in_fullres)
  )

if (nrow(coordinate_qc) > 0) {
  
  p_spatial_qc <- ggplot(
    coordinate_qc,
    aes(
      x = pxl_col_in_fullres,
      y = -pxl_row_in_fullres
    )
  ) +
    geom_point(
      aes(
        shape = QC_Status
      ),
      size = 1.25,
      alpha = 0.8
    ) +
    facet_wrap(
      ~ Sample,
      scales = "free",
      ncol = 2
    ) +
    labs(
      title = "Spatial distribution of QC-retained and removed spots",
      subtitle = paste0(
        "Primary 1st-percentile lower-complexity rule; ",
        "coordinates shown without histology image"
      ),
      x = "Full-resolution pixel X",
      y = "Full-resolution pixel Y",
      shape = "QC status"
    ) +
    theme_void(
      base_size = 11
    ) +
    theme(
      strip.text = element_text(
        face = "bold"
      ),
      plot.title = element_text(
        face = "bold"
      ),
      legend.position = "right"
    )
  
  ggsave(
    filename = file.path(
      figures_dir,
      "11B_Spatial_Distribution_QC_Status.png"
    ),
    plot = p_spatial_qc,
    width = 10,
    height = 10,
    dpi = 400
  )
}


# ============================================================
# 19. FIGURE - REMOVED SPOTS ONLY
# ============================================================

removed_coordinate_qc <- coordinate_qc %>%
  filter(
    QC_Status == "Removed"
  )

if (nrow(removed_coordinate_qc) > 0) {
  
  p_removed_only <- ggplot(
    coordinate_qc,
    aes(
      x = pxl_col_in_fullres,
      y = -pxl_row_in_fullres
    )
  ) +
    geom_point(
      size = 0.75,
      alpha = 0.15
    ) +
    geom_point(
      data = removed_coordinate_qc,
      size = 2,
      alpha = 0.95
    ) +
    facet_wrap(
      ~ Sample,
      scales = "free",
      ncol = 2
    ) +
    labs(
      title = "Spatial location of low-complexity spots removed by QC",
      subtitle = paste0(
        "Background shows all spots; emphasized points are ",
        "spots failing the primary QC rule"
      ),
      x = "Full-resolution pixel X",
      y = "Full-resolution pixel Y"
    ) +
    theme_void(
      base_size = 11
    ) +
    theme(
      strip.text = element_text(
        face = "bold"
      ),
      plot.title = element_text(
        face = "bold"
      )
    )
  
  ggsave(
    filename = file.path(
      figures_dir,
      "11B_Removed_Spots_Spatial_Map.png"
    ),
    plot = p_removed_only,
    width = 10,
    height = 10,
    dpi = 400
  )
}


# ============================================================
# 20. FIGURE - BEFORE / AFTER nFeature
# ============================================================

feature_plot_data <- bind_rows(
  primary_qc %>%
    transmute(
      Sample = Sample,
      Stage = "Before QC",
      Value = nFeature_Spatial
    ),
  primary_qc %>%
    filter(QC_Pass) %>%
    transmute(
      Sample = Sample,