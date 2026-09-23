# ============================================================
# V07 — FINAL ANNOTATION FREEZE
# GSE155882 INDEPENDENT VALIDATION ATLAS
# ============================================================
#
# PURPOSE
# -------
# Freeze final biological annotations for the 22 primary
# resolution-0.6 clusters after V06 + V06B evidence review.
#
# PRIMARY ANNOTATION:
#   resolution 0.6
#
# FINE-STATE METADATA:
#   resolution 1.2
#
# IMPORTANT:
#   NO cell filtering
#   NO reclustering
#   NO integration
#   NO graph rebuilding
#   NO PCA recalculation
#   NO UMAP recalculation
#   NO Sham-vs-TAC differential expression
#   NO CYP/AA validation
#
# Cluster 10 is intentionally retained as an ambiguous
# stromal-myeloid mixed state and is NOT forced into either
# fibroblast or myeloid lineage.
#
# ============================================================


# ============================================================
# 1. CLEAN WORKSPACE
# ============================================================

rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  max.print = 100000
)

set.seed(20260921)


# ============================================================
# 2. PACKAGES
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(Matrix)
  library(patchwork)
})


# ============================================================
# 3. PATHS
# ============================================================

project_root <- "D:/Master/ScRNA seq/TAC/VALIDATION/GSE155882"

input_file <- file.path(
  project_root,
  "CLEAN DATA",
  "V06_GSE155882_ANNOTATION_EVIDENCE_ATLAS.rds"
)

output_file <- file.path(
  project_root,
  "CLEAN DATA",
  "V07_GSE155882_FINAL_ANNOTATED_ATLAS.rds"
)

results_dir <- file.path(
  project_root,
  "RESULTS",
  "V07_FINAL_ANNOTATION_FREEZE"
)

figures_dir <- file.path(
  project_root,
  "FIGURES",
  "V07_FINAL_ANNOTATION_FREEZE"
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


# ============================================================
# 4. LOAD V06 ATLAS
# ============================================================

if (!file.exists(input_file)) {
  stop("V06 annotation-evidence atlas not found.")
}

obj <- readRDS(input_file)

input_n <- ncol(obj)
input_cells <- colnames(obj)

cat("\n============================================\n")
cat("V07 — FINAL ANNOTATION FREEZE\n")
cat("============================================\n")

cat("\nInput cells:", input_n, "\n")
cat("Input genes:", nrow(obj), "\n")


# ============================================================
# 5. INPUT INTEGRITY
# ============================================================

required_metadata <- c(
  "sample_id",
  "condition",
  "replicate",
  "scDblFinder_class",
  "V06_cluster_0.6",
  "V06_fine_cluster_1.2"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(obj@meta.data)
)

if (length(missing_metadata) > 0) {
  stop(
    paste(
      "Missing required metadata:",
      paste(missing_metadata, collapse = ", ")
    )
  )
}

if (!"RNA" %in% Assays(obj)) {
  stop("RNA assay missing.")
}

if (!"pca" %in% Reductions(obj)) {
  stop("Frozen PCA missing.")
}

if (!"umap.unintegrated" %in% Reductions(obj)) {
  stop("Frozen unintegrated UMAP missing.")
}

if (input_n != 28388) {
  stop(
    paste(
      "Unexpected cell number:",
      input_n,
      "Expected: 28388"
    )
  )
}

if (!all(as.character(obj$scDblFinder_class) == "singlet")) {
  stop("Non-singlet cells detected in V06 input.")
}


# ============================================================
# 6. VERIFY PRIMARY / FINE CLUSTERS
# ============================================================

primary_clusters <- sort(
  unique(
    as.character(obj$V06_cluster_0.6)
  )
)

fine_clusters <- sort(
  unique(
    as.character(obj$V06_fine_cluster_1.2)
  )
)

cat("\nPrimary resolution-0.6 clusters:\n")
print(primary_clusters)

cat("\nNumber of primary clusters:",
    length(primary_clusters), "\n")

cat("\nNumber of fine resolution-1.2 clusters:",
    length(fine_clusters), "\n")

if (length(primary_clusters) != 22) {
  stop("Expected exactly 22 primary clusters.")
}

if (length(fine_clusters) != 33) {
  stop("Expected exactly 33 fine clusters.")
}


# ============================================================
# 7. FINAL ANNOTATION MAP
# ============================================================
#
# Based on:
#   V05 unbiased marker discovery
#   V06 canonical annotation audit
#   V06 resolution 0.6 / 1.2 crosswalk
#   V06B targeted review of clusters 8, 10, 16
#
# CYP/AA genes were NOT used to determine these labels.
#
# ============================================================

annotation_map <- data.frame(
  
  Cluster = as.character(0:21),
  
  Broad_Cell_Type = c(
    
    # 0
    "Fibroblasts",
    
    # 1
    "Macrophages",
    
    # 2
    "Endothelial",
    
    # 3
    "Myeloid",
    
    # 4
    "Fibroblasts",
    
    # 5
    "Endothelial",
    
    # 6
    "Fibroblasts",
    
    # 7
    "Fibroblasts",
    
    # 8
    "Endothelial",
    
    # 9
    "Endothelial",
    
    # 10
    "Ambiguous/Mixed",
    
    # 11
    "Vascular mural",
    
    # 12
    "Endothelial",
    
    # 13
    "Pericytes",
    
    # 14
    "Fibroblasts",
    
    # 15
    "Fibroblasts",
    
    # 16
    "Mesothelial/Epicardial",
    
    # 17
    "Fibroblasts",
    
    # 18
    "Interferon-response",
    
    # 19
    "Endothelial",
    
    # 20
    "Myeloid",
    
    # 21
    "Schwann/Neural"
  ),
  
  Refined_Cell_Type = c(
    
    # 0
    "Fibroblasts",
    
    # 1
    "Resident macrophages",
    
    # 2
    "Capillary endothelial",
    
    # 3
    "Antigen-presenting myeloid",
    
    # 4
    "Activated fibroblasts",
    
    # 5
    "Capillary endothelial",
    
    # 6
    "Activated fibroblasts",
    
    # 7
    "Fibroblasts",
    
    # 8
    "Capillary endothelial",
    
    # 9
    "Endothelial",
    
    # 10
    "Stromal-myeloid mixed state",
    
    # 11
    "Vascular mural cells",
    
    # 12
    "Arterial endothelial",
    
    # 13
    "Pericytes",
    
    # 14
    "Stromal/mural fibroblasts",
    
    # 15
    "Activated fibroblasts",
    
    # 16
    "Mesothelial cells",
    
    # 17
    "Activated fibroblasts",
    
    # 18
    "Interferon-response cells",
    
    # 19
    "Lymphatic endothelial",
    
    # 20
    "Myeloid cells",
    
    # 21
    "Schwann/neural cells"
  ),
  
  Confidence = c(
    "High",       # 0
    "High",       # 1
    "High",       # 2
    "High",       # 3
    "High",       # 4
    "High",       # 5
    "High",       # 6
    "High",       # 7
    "High",       # 8
    "Moderate",   # 9
    "Low",        # 10
    "Moderate",   # 11
    "High",       # 12
    "High",       # 13
    "Moderate",   # 14
    "High",       # 15
    "Very high",  # 16
    "Moderate",   # 17
    "High",       # 18
    "High",       # 19
    "Moderate",   # 20
    "High"        # 21
  ),
  
  Annotation_Status = c(
    rep("Resolved", 10),
    "Ambiguous - retained separately",
    rep("Resolved", 11)
  ),
  
  stringsAsFactors = FALSE
)


# ============================================================
# 8. VERIFY ANNOTATION MAP
# ============================================================

if (nrow(annotation_map) != 22) {
  stop("Annotation map does not contain 22 clusters.")
}

if (anyDuplicated(annotation_map$Cluster) > 0) {
  stop("Duplicated cluster IDs in annotation map.")
}

missing_annotation_clusters <- setdiff(
  primary_clusters,
  annotation_map$Cluster
)

extra_annotation_clusters <- setdiff(
  annotation_map$Cluster,
  primary_clusters
)

if (length(missing_annotation_clusters) > 0) {
  stop(
    paste(
      "Missing annotation clusters:",
      paste(
        missing_annotation_clusters,
        collapse = ", "
      )
    )
  )
}

if (length(extra_annotation_clusters) > 0) {
  stop(
    paste(
      "Annotation map contains unexpected clusters:",
      paste(
        extra_annotation_clusters,
        collapse = ", "
      )
    )
  )
}

write.csv(
  annotation_map,
  file.path(
    results_dir,
    "V07_Final_Cluster_Annotation_Map.csv"
  ),
  row.names = FALSE
)

cat("\nFinal annotation map:\n")
print(annotation_map, row.names = FALSE)


# ============================================================
# 9. MAP ANNOTATIONS TO CELLS
# ============================================================

cluster_lookup_broad <- setNames(
  annotation_map$Broad_Cell_Type,
  annotation_map$Cluster
)

cluster_lookup_refined <- setNames(
  annotation_map$Refined_Cell_Type,
  annotation_map$Cluster
)

cluster_lookup_confidence <- setNames(
  annotation_map$Confidence,
  annotation_map$Cluster
)

cluster_lookup_status <- setNames(
  annotation_map$Annotation_Status,
  annotation_map$Cluster
)

cell_cluster <- as.character(
  obj$V06_cluster_0.6
)

obj$V07_primary_cluster <- cell_cluster

obj$V07_broad_cell_type <- unname(
  cluster_lookup_broad[cell_cluster]
)

obj$V07_refined_cell_type <- unname(
  cluster_lookup_refined[cell_cluster]
)

obj$V07_annotation_confidence <- unname(
  cluster_lookup_confidence[cell_cluster]
)

obj$V07_annotation_status <- unname(
  cluster_lookup_status[cell_cluster]
)

obj$V07_fine_cluster <- as.character(
  obj$V06_fine_cluster_1.2
)


# ============================================================
# 10. FREEZE VERSION METADATA
# ============================================================

obj$V07_annotation_version <-
  "GSE155882_Validation_Atlas_v1_Final"

obj$V07_annotation_step <-
  "V07"

obj$V07_annotation_date <-
  "2026-09-21"

obj$V07_primary_resolution <-
  "0.6"

obj$V07_fine_state_resolution <-
  "1.2"

obj$V07_annotation_basis <-
  paste0(
    "V05 unbiased clustering + ",
    "V06 annotation evidence audit + ",
    "V06B targeted identity adjudication"
  )


# ============================================================
# 11. CHECK FOR UNMAPPED CELLS
# ============================================================

annotation_fields <- c(
  "V07_broad_cell_type",
  "V07_refined_cell_type",
  "V07_annotation_confidence",
  "V07_annotation_status"
)

for (field in annotation_fields) {
  
  if (any(is.na(obj@meta.data[[field]]))) {
    
    stop(
      paste(
        "NA annotations detected in:",
        field
      )
    )
  }
}


# ============================================================
# 12. BROAD CELL-TYPE COUNTS
# ============================================================

broad_counts <- as.data.frame(
  table(
    Broad_Cell_Type =
      obj$V07_broad_cell_type
  )
)

colnames(broad_counts)[2] <- "Cells"

broad_counts$Percent_of_Atlas <-
  100 * broad_counts$Cells / ncol(obj)

broad_counts <- broad_counts %>%
  arrange(
    desc(Cells)
  )

write.csv(
  broad_counts,
  file.path(
    results_dir,
    "V07_Broad_Cell_Type_Counts.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. REFINED CELL-TYPE COUNTS
# ============================================================

refined_counts <- as.data.frame(
  table(
    Refined_Cell_Type =
      obj$V07_refined_cell_type
  )
)

colnames(refined_counts)[2] <- "Cells"

refined_counts$Percent_of_Atlas <-
  100 * refined_counts$Cells / ncol(obj)

refined_counts <- refined_counts %>%
  arrange(
    desc(Cells)
  )

write.csv(
  refined_counts,
  file.path(
    results_dir,
    "V07_Refined_Cell_Type_Counts.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. CONFIDENCE COUNTS
# ============================================================

confidence_counts <- as.data.frame(
  table(
    Confidence =
      obj$V07_annotation_confidence
  )
)

colnames(confidence_counts)[2] <- "Cells"

confidence_counts$Percent_of_Atlas <-
  100 * confidence_counts$Cells / ncol(obj)

write.csv(
  confidence_counts,
  file.path(
    results_dir,
    "V07_Annotation_Confidence_Counts.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. CLUSTER × FINAL ANNOTATION
# ============================================================

cluster_annotation_summary <-
  obj@meta.data %>%
  mutate(
    Cluster = as.character(
      V07_primary_cluster
    )
  ) %>%
  count(
    Cluster,
    V07_broad_cell_type,
    V07_refined_cell_type,
    V07_annotation_confidence,
    V07_annotation_status,
    name = "Cells"
  ) %>%
  arrange(
    as.numeric(Cluster)
  )

write.csv(
  cluster_annotation_summary,
  file.path(
    results_dir,
    "V07_Cluster_Annotation_Cell_Counts.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. SAMPLE × BROAD CELL TYPE
# ============================================================

sample_broad_counts <- as.data.frame(
  table(
    Sample = obj$sample_id,
    Broad_Cell_Type =
      obj$V07_broad_cell_type
  )
)

sample_broad_counts <- sample_broad_counts[
  sample_broad_counts$Freq > 0,
]

sample_broad_counts <- sample_broad_counts %>%
  group_by(
    Sample
  ) %>%
  mutate(
    Fraction_of_Sample =
      Freq / sum(Freq),
    Percent_of_Sample =
      100 * Fraction_of_Sample
  ) %>%
  ungroup()

write.csv(
  sample_broad_counts,
  file.path(
    results_dir,
    "V07_Sample_by_Broad_Cell_Type.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. SAMPLE × REFINED CELL TYPE
# ============================================================

sample_refined_counts <- as.data.frame(
  table(
    Sample = obj$sample_id,
    Refined_Cell_Type =
      obj$V07_refined_cell_type
  )
)

sample_refined_counts <- sample_refined_counts[
  sample_refined_counts$Freq > 0,
]

sample_refined_counts <- sample_refined_counts %>%
  group_by(
    Sample
  ) %>%
  mutate(
    Fraction_of_Sample =
      Freq / sum(Freq),
    Percent_of_Sample =
      100 * Fraction_of_Sample
  ) %>%
  ungroup()

write.csv(
  sample_refined_counts,
  file.path(
    results_dir,
    "V07_Sample_by_Refined_Cell_Type.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. CONDITION × BROAD CELL TYPE
# ============================================================
#
# DESCRIPTIVE ONLY.
# These are not biological replicate-level statistical tests.
#
# ============================================================

condition_broad_counts <- as.data.frame(
  table(
    Condition = obj$condition,
    Broad_Cell_Type =
      obj$V07_broad_cell_type
  )
)

condition_broad_counts <- condition_broad_counts[
  condition_broad_counts$Freq > 0,
]

condition_broad_counts <- condition_broad_counts %>%
  group_by(
    Condition
  ) %>%
  mutate(
    Fraction_of_Condition =
      Freq / sum(Freq),
    Percent_of_Condition =
      100 * Fraction_of_Condition
  ) %>%
  ungroup()

write.csv(
  condition_broad_counts,
  file.path(
    results_dir,
    "V07_Condition_by_Broad_Cell_Type_DESCRIPTIVE.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. CONDITION × REFINED CELL TYPE
# ============================================================

condition_refined_counts <- as.data.frame(
  table(
    Condition = obj$condition,
    Refined_Cell_Type =
      obj$V07_refined_cell_type
  )
)

condition_refined_counts <- condition_refined_counts[
  condition_refined_counts$Freq > 0,
]

condition_refined_counts <- condition_refined_counts %>%
  group_by(
    Condition
  ) %>%
  mutate(
    Fraction_of_Condition =
      Freq / sum(Freq),
    Percent_of_Condition =
      100 * Fraction_of_Condition
  ) %>%
  ungroup()

write.csv(
  condition_refined_counts,
  file.path(
    results_dir,
    "V07_Condition_by_Refined_Cell_Type_DESCRIPTIVE.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 20. FINE-STATE CROSSWALK
# ============================================================

fine_crosswalk <- as.data.frame(
  table(
    Primary_Cluster =
      obj$V07_primary_cluster,
    Fine_Cluster =
      obj$V07_fine_cluster
  )
)

fine_crosswalk <- fine_crosswalk[
  fine_crosswalk$Freq > 0,
]

fine_crosswalk <- fine_crosswalk %>%
  group_by(
    Primary_Cluster
  ) %>%
  mutate(
    Fraction_of_Primary_Cluster =
      Freq / sum(Freq)
  ) %>%
  ungroup()

write.csv(
  fine_crosswalk,
  file.path(
    results_dir,
    "V07_Primary_to_Fine_Cluster_Crosswalk.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 21. FINE CLUSTER → FINAL ANNOTATION CROSSWALK
# ============================================================

fine_annotation_crosswalk <-
  obj@meta.data %>%
  transmute(
    Fine_Cluster =
      as.character(V07_fine_cluster),
    Primary_Cluster =
      as.character(V07_primary_cluster),
    Broad_Cell_Type =
      V07_broad_cell_type,
    Refined_Cell_Type =
      V07_refined_cell_type
  ) %>%
  count(
    Fine_Cluster,
    Primary_Cluster,
    Broad_Cell_Type,
    Refined_Cell_Type,
    name = "Cells"
  ) %>%
  group_by(
    Fine_Cluster
  ) %>%
  mutate(
    Fraction_of_Fine_Cluster =
      Cells / sum(Cells)
  ) %>%
  ungroup()

write.csv(
  fine_annotation_crosswalk,
  file.path(
    results_dir,
    "V07_Fine_Cluster_Annotation_Crosswalk.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 22. CLUSTER-10 SAFETY AUDIT
# ============================================================
#
# Cluster 10 must remain separate and must not accidentally
# enter fibroblast- or myeloid-specific analyses later.
#
# ============================================================

cluster10_cells <- colnames(obj)[
  obj$V07_primary_cluster == "10"
]

cluster10_audit <- data.frame(
  Metric = c(
    "Cluster_10_cells",
    "Broad_annotation_unique_count",
    "Refined_annotation_unique_count",
    "Broad_annotation_is_Ambiguous_Mixed",
    "Refined_annotation_is_Stromal_myeloid_mixed_state"
  ),
  Value = c(
    length(cluster10_cells),
    length(
      unique(
        obj$V07_broad_cell_type[
          obj$V07_primary_cluster == "10"
        ]
      )
    ),
    length(
      unique(
        obj$V07_refined_cell_type[
          obj$V07_primary_cluster == "10"
        ]
      )
    ),
    all(
      obj$V07_broad_cell_type[
        obj$V07_primary_cluster == "10"
      ] == "Ambiguous/Mixed"
    ),
    all(
      obj$V07_refined_cell_type[
        obj$V07_primary_cluster == "10"
      ] == "Stromal-myeloid mixed state"
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  cluster10_audit,
  file.path(
    results_dir,
    "V07_Cluster10_Ambiguous_State_Safety_Audit.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 23. CREATE FACTORS FOR PLOTTING
# ============================================================

obj$V07_primary_cluster <- factor(
  obj$V07_primary_cluster,
  levels = as.character(0:21)
)
