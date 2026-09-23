# ============================================================
# V05 — GSE155882
# UNBIASED GRAPH CLUSTERING + RESOLUTION AUDIT +
# CLUSTER MARKERS + ANNOTATION EVIDENCE
#
# INPUT:
#   V04_GSE155882_NORMALIZED_PCA_UNINTEGRATED_ATLAS.rds
#
# PURPOSE:
#   1. Preserve the frozen V04 unintegrated representation
#   2. Construct SNN graph using V04 PCs 1:32
#   3. Perform predefined resolution sweep
#   4. Audit cluster number and minimum cluster size
#   5. Audit sample/condition representation in each cluster
#   6. Select a primary resolution using structural criteria
#   7. Calculate unbiased cluster markers
#   8. Generate canonical-marker evidence for annotation
#   9. Save 300-dpi PNG figures
#
# IMPORTANT:
#   - NO cells removed
#   - NO new normalization
#   - NO new HVG selection
#   - NO batch integration
#   - NO CYP/AA-driven clustering
#   - NO CYP/AA-driven annotation
#   - NO final cell-type labels assigned
#   - NO Sham-vs-TAC differential expression
#   - NO condition-level inference
#
# V05 clusters are exploratory transcriptional groups.
# Final biological annotation occurs only after review.
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(Seurat)
library(Matrix)
library(ggplot2)


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
  "V05_UNBIASED_CLUSTERING_MARKER_AUDIT"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "V05_UNBIASED_CLUSTERING_MARKER_AUDIT"
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
  "V04_GSE155882_NORMALIZED_PCA_UNINTEGRATED_ATLAS.rds"
)

output_file <- file.path(
  clean_dir,
  "V05_GSE155882_UNBIASED_CLUSTERED_ATLAS.rds"
)


# ============================================================
# 3. LOAD FROZEN V04 ATLAS
# ============================================================

if (!file.exists(input_file)) {
  stop(
    paste0(
      "V04 object not found:\n",
      input_file
    )
  )
}

obj <- readRDS(
  input_file
)


cat("\n============================================\n")
cat("V05 — UNBIASED CLUSTERING + MARKER AUDIT\n")
cat("============================================\n")

cat(
  "Cells entering V05:",
  ncol(obj),
  "\n"
)

cat(
  "Genes:",
  nrow(obj),
  "\n"
)


# ============================================================
# 4. PRE-V05 INTEGRITY CHECK
# ============================================================

required_metadata <- c(
  "sample_id",
  "condition",
  "replicate",
  "nFeature_RNA",
  "nCount_RNA",
  "percent.mt",
  "scDblFinder_score",
  "scDblFinder_class",
  "V04_integration_status",
  "V04_selected_npcs"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(obj@meta.data)
)

if (length(missing_metadata) > 0) {
  
  stop(
    paste(
      "Missing required metadata:",
      paste(
        missing_metadata,
        collapse = ", "
      )
    )
  )
}


if (!"pca" %in% names(obj@reductions)) {
  stop(
    "PCA reduction is missing."
  )
}


if (!"umap.unintegrated" %in% names(obj@reductions)) {
  stop(
    "V04 unintegrated UMAP is missing."
  )
}


if (!all(
  as.character(
    obj$V04_integration_status
  ) == "UNINTEGRATED"
)) {
  
  stop(
    "V04 object is not marked UNINTEGRATED."
  )
}


if (!all(
  obj$scDblFinder_class ==
  "singlet"
)) {
  
  stop(
    "Non-singlet cells detected in V05 input."
  )
}


input_cells <- colnames(obj)

input_n <- ncol(obj)

input_sample_counts <- table(
  obj$sample_id
)


# ============================================================
# 5. LOCK V04 PC COUNT
# ============================================================

selected_npcs_values <- unique(
  obj$V04_selected_npcs
)

selected_npcs_values <- as.numeric(
  selected_npcs_values
)


if (length(selected_npcs_values) != 1) {
  stop(
    "V04 selected PC count is not unique."
  )
}


selected_npcs <- selected_npcs_values[1]


cat(
  "Using frozen V04 PC count:",
  selected_npcs,
  "\n"
)


if (selected_npcs != 32) {
  
  warning(
    paste(
      "Expected 32 PCs from reviewed V04, but object reports",
      selected_npcs
    )
  )
}


# ============================================================
# 6. GRAPH CONSTRUCTION
#
# Uses the frozen UNINTEGRATED PCA.
#
# No sample/condition/CYP information enters graph
# construction.
# ============================================================

set.seed(20260920)


obj <- FindNeighbors(
  object = obj,
  reduction = "pca",
  dims = 1:selected_npcs,
  k.param = 20,
  graph.name = c(
    "V05_nn",
    "V05_snn"
  ),
  verbose = FALSE
)


# ============================================================
# 7. PREDEFINED RESOLUTION SWEEP
#
# We deliberately examine a broad but reasonable range.
#
# No resolution is selected from CYP/AA behavior.
# ============================================================

resolution_values <- c(
  0.1,
  0.2,
  0.3,
  0.4,
  0.5,
  0.6,
  0.8,
  1.0,
  1.2
)


for (resolution_now in resolution_values) {
  
  cat(
    "\nClustering resolution:",
    resolution_now,
    "\n"
  )
  
  set.seed(20260920)
  
  obj <- FindClusters(
    object = obj,
    graph.name = "V05_snn",
    resolution = resolution_now,
    algorithm = 1,
    random.seed = 20260920,
    verbose = FALSE
  )
  
  source_column <- "seurat_clusters"
  
  target_column <- paste0(
    "V05_res_",
    gsub(
      "\\.",
      "_",
      as.character(
        resolution_now
      )
    )
  )
  
  obj[[target_column]] <- as.character(
    obj[[source_column]][, 1]
  )
}


# ============================================================
# 8. RESOLUTION STRUCTURE AUDIT
#
# For each resolution:
#   - number of clusters
#   - smallest cluster
#   - median cluster size
#   - largest cluster
#   - number of tiny clusters (<50 cells)
#   - number of small clusters (<100 cells)
# ============================================================

resolution_audit_list <- list()


for (resolution_now in resolution_values) {
  
  cluster_column <- paste0(
    "V05_res_",
    gsub(
      "\\.",
      "_",
      as.character(
        resolution_now
      )
    )
  )
  
  cluster_vector <- as.character(
    obj[[cluster_column]][, 1]
  )
  
  cluster_sizes <- table(
    cluster_vector
  )
  
  resolution_audit_list[[
    as.character(
      resolution_now
    )
  ]] <- data.frame(
    
    Resolution = resolution_now,
    
    Number_of_Clusters =
      length(
        cluster_sizes
      ),
    
    Minimum_Cluster_Size =
      min(
        cluster_sizes
      ),
    
    Median_Cluster_Size =
      median(
        as.numeric(
          cluster_sizes
        )
      ),
    
    Maximum_Cluster_Size =
      max(
        cluster_sizes
      ),
    
    Clusters_Under_50_Cells =
      sum(
        cluster_sizes < 50
      ),
    
    Clusters_Under_100_Cells =
      sum(
        cluster_sizes < 100
      ),
    
    stringsAsFactors = FALSE
  )
}


resolution_audit <- do.call(
  rbind,
  resolution_audit_list
)

rownames(
  resolution_audit
) <- NULL


cat("\n============================================\n")
cat("RESOLUTION STRUCTURE AUDIT\n")
cat("============================================\n")

print(
  resolution_audit,
  row.names = FALSE
)


write.csv(
  resolution_audit,
  file.path(
    results_dir,
    "V05_Resolution_Structure_Audit.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9. SAMPLE REPRESENTATION AUDIT AT EACH RESOLUTION
#
# A cluster is considered "sample-restricted" for this
# structural audit if >90% of its cells come from ONE library.
#
# This is NOT automatically a bad cluster.
# Rare biology may genuinely be sample enriched.
#
# We only use this metric as one safeguard against choosing
# a heavily sample-driven global resolution.
# ============================================================

sample_representation_list <- list()


for (resolution_now in resolution_values) {
  
  cluster_column <- paste0(
    "V05_res_",
    gsub(
      "\\.",
      "_",
      as.character(
        resolution_now
      )
    )
  )
  
  cluster_vector <- as.character(
    obj[[cluster_column]][, 1]
  )
  
  sample_vector <- as.character(
    obj$sample_id
  )
  
  tab <- table(
    cluster_vector,
    sample_vector
  )
  
  proportions <- prop.table(
    tab,
    margin = 1
  )
  
  for (cluster_now in rownames(tab)) {
    
    cluster_counts <- tab[
      cluster_now,
      ,
      drop = TRUE
    ]
    
    cluster_props <- proportions[
      cluster_now,
      ,
      drop = TRUE
    ]
    
    dominant_sample <- names(
      which.max(
        cluster_props
      )
    )
    
    dominant_fraction <- max(
      cluster_props
    )
    
    represented_samples <- sum(
      cluster_counts > 0
    )
    
    represented_samples_1pct <- sum(
      cluster_props >= 0.01
    )
    
    sample_representation_list[[
      paste(
        resolution_now,
        cluster_now,
        sep = "_"
      )
    ]] <- data.frame(
      
      Resolution =
        resolution_now,
      
      Cluster =
        cluster_now,
      
      Cluster_Size =
        sum(
          cluster_counts
        ),
      
      Dominant_Sample =
        dominant_sample,
      
      Dominant_Sample_Fraction =
        dominant_fraction,
      
      Samples_With_Any_Cells =
        represented_samples,
      
      Samples_With_At_Least_1pct =
        represented_samples_1pct,
      
      Sample_Restricted_90pct =
        dominant_fraction > 0.90,
      
      stringsAsFactors = FALSE
    )
  }
}


sample_representation <- do.call(
  rbind,
  sample_representation_list
)

rownames(
  sample_representation
) <- NULL


write.csv(
  sample_representation,
  file.path(
    results_dir,
    "V05_Resolution_Sample_Representation.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. SUMMARIZE SAMPLE RESTRICTION BY RESOLUTION
# ============================================================

sample_restriction_summary_list <- list()


for (resolution_now in resolution_values) {
  
  x <- sample_representation[
    sample_representation$Resolution ==
      resolution_now,
    ,
    drop = FALSE
  ]
  
  sample_restriction_summary_list[[
    as.character(
      resolution_now
    )
  ]] <- data.frame(
    
    Resolution =
      resolution_now,
    
    Number_of_Clusters =
      nrow(x),
    
    Sample_Restricted_Clusters =
      sum(
        x$Sample_Restricted_90pct
      ),
    
    Fraction_Sample_Restricted =
      mean(
        x$Sample_Restricted_90pct
      ),
    
    Median_Dominant_Sample_Fraction =
      median(
        x$Dominant_Sample_Fraction
      ),
    
    stringsAsFactors = FALSE
  )
}


sample_restriction_summary <- do.call(
  rbind,
  sample_restriction_summary_list
)

rownames(
  sample_restriction_summary
) <- NULL


cat("\n============================================\n")
cat("SAMPLE-RESTRICTION SUMMARY\n")
cat("============================================\n")

print(
  sample_restriction_summary,
  row.names = FALSE
)


write.csv(
  sample_restriction_summary,
  file.path(
    results_dir,
    "V05_Resolution_Sample_Restriction_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. PREDEFINED PRIMARY-RESOLUTION RULE
#
# Structural criteria:
#
# Candidate resolution must:
#   1. Have no cluster <50 cells
#   2. Have <=10% sample-restricted clusters
#   3. Have at least 8 clusters
#
# Among passing resolutions:
#   choose the HIGHEST resolution.
#
# Why highest?
#   We want sufficient granularity for major stromal,
#   vascular and immune compartments without choosing a
#   resolution that generates tiny/sample-specific groups.
#
# IMPORTANT:
#   CYP/AA expression does NOT enter this decision.
#
# If no resolution passes:
#   fall back to 0.4 and flag for manual review.
# ============================================================

resolution_decision <- merge(
  resolution_audit,
  sample_restriction_summary[
    ,
    c(
      "Resolution",
      "Sample_Restricted_Clusters",
      "Fraction_Sample_Restricted",
      "Median_Dominant_Sample_Fraction"
    )
  ],
  by = "Resolution",
  all.x = TRUE
)


resolution_decision$Pass_Min_Size <- (
  resolution_decision$Minimum_Cluster_Size >=
    50
)


resolution_decision$Pass_Sample_Restriction <- (
  resolution_decision$Fraction_Sample_Restricted <=
    0.10
)


resolution_decision$Pass_Min_Cluster_Number <- (
  resolution_decision$Number_of_Clusters >=
    8
)


resolution_decision$Pass_All <- (
  resolution_decision$Pass_Min_Size &
    resolution_decision$Pass_Sample_Restriction &
    resolution_decision$Pass_Min_Cluster_Number
)


passing_resolutions <- resolution_decision$Resolution[
  resolution_decision$Pass_All
]


if (length(
  passing_resolutions
) > 0) {
  
  primary_resolution <- max(
    passing_resolutions
  )
  
  primary_resolution_rule <-
    "highest_resolution_passing_predefined_structural_criteria"
  
} else {
  
  primary_resolution <- 0.4
  
  primary_resolution_rule <-
    "fallback_0.4_no_resolution_passed_all_predefined_criteria"
}


cat("\n============================================\n")
cat("PRIMARY RESOLUTION DECISION\n")
cat("============================================\n")

print(
  resolution_decision,
  row.names = FALSE
)

cat(
  "\nSelected primary resolution:",
  primary_resolution,
  "\n"
)

cat(
  "Selection rule:",
  primary_resolution_rule,
  "\n"
)


write.csv(
  resolution_decision,
  file.path(
    results_dir,
    "V05_Primary_Resolution_Decision.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. LOCK PRIMARY CLUSTERS
# ============================================================

primary_cluster_column <- paste0(
  "V05_res_",
  gsub(
    "\\.",
    "_",
    as.character(
      primary_resolution
    )
  )
)


obj$V05_primary_cluster <- as.character(
  obj[[primary_cluster_column]][, 1]
)


cluster_levels_numeric <- sort(
  unique(
    as.numeric(
      obj$V05_primary_cluster
    )
  )
)


obj$V05_primary_cluster <- factor(
  obj$V05_primary_cluster,
  levels = as.character(
    cluster_levels_numeric
  )
)


Idents(
  obj
) <- "V05_primary_cluster"


# ============================================================
# 13. PRIMARY CLUSTER SIZE TABLE
# ============================================================

primary_cluster_sizes <- as.data.frame(
  table(
    Cluster = obj$V05_primary_cluster
  )
)


colnames(
  primary_cluster_sizes
) <- c(
  "Cluster",
  "Cells"
)


primary_cluster_sizes$Percent_of_Atlas <- (
  100 *
    primary_cluster_sizes$Cells /
    ncol(obj)
)


cat("\n============================================\n")
cat("PRIMARY CLUSTER SIZES\n")
cat("============================================\n")

print(
  primary_cluster_sizes,
  row.names = FALSE
)


write.csv(
  primary_cluster_sizes,
  file.path(
    results_dir,
    "V05_Primary_Cluster_Sizes.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. PRIMARY CLUSTER × SAMPLE COUNTS
# ============================================================

cluster_sample_counts <- as.data.frame.matrix(
  table(
    obj$V05_primary_cluster,
    obj$sample_id
  )
)


cluster_sample_counts$Cluster <- rownames(
  cluster_sample_counts
)


cluster_sample_counts <- cluster_sample_counts[
  ,
  c(
    "Cluster",
    setdiff(
      colnames(
        cluster_sample_counts
      ),
      "Cluster"
    )
  ),
  drop = FALSE
]


rownames(
  cluster_sample_counts
) <- NULL


write.csv(
  cluster_sample_counts,
  file.path(
    results_dir,
    "V05_Primary_Cluster_by_Sample_Counts.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. PRIMARY CLUSTER × SAMPLE PROPORTIONS
#
# ROW proportions:
#   composition of each cluster across libraries.
# ============================================================

cluster_sample_table <- table(
  obj$V05_primary_cluster,
  obj$sample_id
)


cluster_sample_rowprop <- prop.table(
  cluster_sample_table,
  margin = 1
)


write.csv(
  as.data.frame.matrix(
    cluster_sample_rowprop
  ),
  file.path(
    results_dir,
    "V05_Primary_Cluster_by_Sample_RowProportions.csv"
  ),
  row.names = TRUE
)


# ============================================================
# 16. PRIMARY CLUSTER × CONDITION COUNTS
#
# DESCRIPTIVE ONLY.
# No statistical testing.
# ============================================================

cluster_condition_table <- table(
  obj$V05_primary_cluster,
  obj$condition
)


write.csv(
  as.data.frame.matrix(
    cluster_condition_table
  ),
  file.path(
    results_dir,
    "V05_Primary_Cluster_by_Condition_Counts.csv"
  ),
  row.names = TRUE
)


# ============================================================
# 17. WITHIN-LIBRARY CLUSTER PROPORTIONS
#
# This is the useful descriptive composition view:
# percentage of each library belonging to each cluster.
#
# We keep individual replicates separate.
# ============================================================

cluster_by_sample_prop <- prop.table(
  cluster_sample_table,
  margin = 2
)


write.csv(
  as.data.frame.matrix(
    cluster_by_sample_prop
  ),