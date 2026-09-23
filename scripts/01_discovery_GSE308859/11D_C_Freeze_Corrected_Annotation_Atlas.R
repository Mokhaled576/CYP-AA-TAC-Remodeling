# ============================================================
# STEP 11D-C
# FREEZE CORRECTED ANNOTATION ATLAS
#
# PURPOSE:
# Write the manually validated 29-cluster annotation into NEW
# metadata columns and save a new corrected frozen atlas.
#
# IMPORTANT:
# - NO reclustering
# - NO cell filtering
# - NO expression data modification
# - NO original annotations overwritten
# - Original broad_cell_type/refined_cell_type retained
# - Corrected labels stored in NEW metadata columns
#
# INPUT:
# 07C_FINAL_FROZEN_ANNOTATED_ATLAS.rds
#
# OUTPUT:
# 11D_C_FINAL_CORRECTED_ANNOTATED_ATLAS.rds
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)


# ============================================================
# 2. PATHS
# ============================================================

project_dir <- "D:/Master/ScRNA seq/TAC"

input_file <- file.path(
  project_dir,
  "CLEAN DATA",
  "07C_FINAL_FROZEN_ANNOTATED_ATLAS.rds"
)

output_file <- file.path(
  project_dir,
  "CLEAN DATA",
  "11D_C_FINAL_CORRECTED_ANNOTATED_ATLAS.rds"
)

results_dir <- file.path(
  project_dir,
  "RESULTS",
  "11_SPATIAL",
  "11D_C_FINAL_CORRECTED_ANNOTATION"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "11_SPATIAL",
  "11D_C_FINAL_CORRECTED_ANNOTATION"
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


# ============================================================
# 3. LOAD ORIGINAL FROZEN ATLAS
# ============================================================

cat("\n============================================\n")
cat("STEP 11D-C\n")
cat("FREEZE CORRECTED ANNOTATION ATLAS\n")
cat("============================================\n")

scrna <- readRDS(input_file)

DefaultAssay(scrna) <- "RNA"

cat("Input atlas:", input_file, "\n")
cat("Cells:", ncol(scrna), "\n")
cat("Genes:", nrow(scrna), "\n")


# ============================================================
# 4. VERIFY REQUIRED METADATA
# ============================================================

required_metadata <- c(
  "seurat_clusters",
  "broad_cell_type",
  "refined_cell_type",
  "sample_id"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(scrna@meta.data)
)

if (length(missing_metadata) > 0) {
  
  stop(
    paste(
      "Missing required metadata:",
      paste(missing_metadata, collapse = ", ")
    )
  )
}

cat("Required metadata verified.\n")


# ============================================================
# 5. VERIFY EXACT CLUSTER STRUCTURE
# ============================================================

cluster_ids <- sort(
  unique(
    as.numeric(
      as.character(
        scrna$seurat_clusters
      )
    )
  )
)

expected_clusters <- as.numeric(0:28)

if (!setequal(cluster_ids, expected_clusters)) {
  
  stop(
    paste(
      "Unexpected cluster IDs:",
      paste(cluster_ids, collapse = ", ")
    )
  )
}

if (length(cluster_ids) != 29) {
  
  stop(
    paste(
      "Expected 29 clusters but found",
      length(cluster_ids)
    )
  )
}

if (!all(cluster_ids == expected_clusters)) {
  
  stop(
    "Cluster IDs are not exactly 0-28."
  )
}

cat("Verified 29 frozen clusters: 0-28.\n")


# ============================================================
# 6. RECORD ORIGINAL ATLAS STATE
#
# These will be checked after annotation.
# ============================================================

original_cells <- colnames(scrna)

original_genes <- rownames(scrna)

original_clusters <- as.character(
  scrna$seurat_clusters
)

original_broad <- as.character(
  scrna$broad_cell_type
)

original_refined <- as.character(
  scrna$refined_cell_type
)

original_ncells <- ncol(scrna)

original_ngenes <- nrow(scrna)


# ============================================================
# 7. FINAL VALIDATED CLUSTER ANNOTATION MAP
#
# Based on:
# - 11D-A annotation audit
# - 11D-B1 reconstruction evidence
# - 11D-B2 final annotation validation
# - 11D-B3 targeted resolution
#
# Cluster 10 corrected to Fibroblasts after B3.
# Cluster 16 retained conservatively as Cycling myeloid cells.
# Cluster 19 retained as Dendritic cells.
# Cluster 24 confirmed Vascular smooth muscle cells.
# ============================================================

final_annotation <- data.frame(
  
  Cluster = as.character(0:28),
  
  Corrected_Broad = c(
    
    "Fibroblasts",             # 0
    "Macrophages",             # 1
    "B cells",                 # 2
    "Fibroblasts",             # 3
    "Fibroblasts",             # 4
    "Macrophages",             # 5
    "Fibroblasts",             # 6
    "Neutrophils",             # 7
    "Endothelial cells",       # 8
    "B cells",                 # 9
    
    "Fibroblasts",             # 10  <-- B3 corrected
    
    "Fibroblasts",             # 11
    "Platelets",               # 12
    "T cells",                 # 13
    "NK cells",                # 14
    "Fibroblasts",             # 15
    "Myeloid cells",           # 16
    "Monocytes",               # 17
    "Endothelial cells",       # 18
    "Dendritic cells",         # 19
    "Neutrophils",             # 20
    "Endothelial cells",       # 21
    "Cardiomyocytes",          # 22
    "Erythroid cells",         # 23
    "Smooth muscle cells",     # 24
    "T cells",                 # 25
    "T cells",                 # 26
    "Pericytes",               # 27
    "B cells"                  # 28
  ),
  
  Corrected_Refined = c(
    
    "Fibroblasts",                     # 0
    "Macrophages",                     # 1
    "B cells",                         # 2
    "Fibroblasts",                     # 3
    "Fibroblasts",                     # 4
    "Resident macrophages",            # 5
    "Activated fibroblasts",           # 6
    "Inflammatory neutrophils",        # 7
    "Capillary endothelial cells",     # 8
    "B cells",                         # 9
    
    "Fibroblasts",                     # 10 <-- B3 corrected
    
    "Fibroblasts",                     # 11
    "Platelets/megakaryocytes",        # 12
    "T cells",                         # 13
    "NK/cytotoxic lymphocytes",        # 14
    "Fibroblasts",                     # 15
    "Cycling myeloid cells",           # 16
    "Monocytes",                       # 17
    "Endothelial cells",               # 18
    "Dendritic cells",                 # 19
    "Neutrophils",                     # 20
    "Arterial endothelial cells",      # 21
    "Cardiomyocytes",                  # 22
    "Erythroid cells",                 # 23
    "Vascular smooth muscle cells",    # 24
    "Gamma-delta T cells",             # 25
    "Immature T cells",                # 26
    "Pericytes",                       # 27
    "Immature B cells"                 # 28
  ),
  
  Annotation_Confidence = c(
    
    "High",       # 0
    "High",       # 1
    "Very high",  # 2
    "High",       # 3
    "High",       # 4
    "Very high",  # 5
    "Very high",  # 6
    "Very high",  # 7
    "Very high",  # 8
    "High",       # 9
    
    "Very high",  # 10 after B3
    
    "Very high",  # 11
    "Very high",  # 12
    "Very high",  # 13
    "Very high",  # 14
    "High",       # 15
    
    "High",       # 16 cycling lineage
    
    "Very high",  # 17
    "High",       # 18
    "High",       # 19
    "Very high",  # 20
    "Very high",  # 21
    "Very high",  # 22
    "Very high",  # 23
    "Very high",  # 24
    "Very high",  # 25
    "Very high",  # 26
    "Very high",  # 27
    "High"        # 28
  ),
  
  Validation_Source = c(
    
    rep(
      "11D-A + 11D-B1 + 11D-B2",
      10
    ),
    
    "11D-A + 11D-B1 + 11D-B2 + 11D-B3",
    
    rep(
      "11D-A + 11D-B1 + 11D-B2",
      5
    ),
    
    "11D-A + 11D-B1 + 11D-B2 + 11D-B3", # 16
    
    rep(
      "11D-A + 11D-B1 + 11D-B2",
      2
    ),
    
    "11D-A + 11D-B1 + 11D-B2 + 11D-B3", # 19
    
    rep(
      "11D-A + 11D-B1 + 11D-B2",
      4
    ),
    
    "11D-A + 11D-B1 + 11D-B2 + 11D-B3", # 24
    
    rep(
      "11D-A + 11D-B1 + 11D-B2",
      4
    )
  ),
  
  stringsAsFactors = FALSE
)


# ============================================================
# 8. VALIDATE ANNOTATION TABLE
# ============================================================

if (nrow(final_annotation) != 29) {
  
  stop(
    "Final annotation table does not contain exactly 29 rows."
  )
}

if (anyDuplicated(final_annotation$Cluster) > 0) {
  
  stop(
    "Duplicate cluster IDs detected in final annotation table."
  )
}

if (!setequal(
  final_annotation$Cluster,
  as.character(0:28)
)) {
  
  stop(
    "Final annotation table does not contain exactly clusters 0-28."
  )
}

if (any(is.na(final_annotation$Corrected_Broad))) {
  
  stop(
    "NA broad annotation detected."
  )
}

if (any(is.na(final_annotation$Corrected_Refined))) {
  
  stop(
    "NA refined annotation detected."
  )
}

cat("Final annotation map passed validation.\n")


# ============================================================
# 9. SAVE FINAL CLUSTER-LEVEL ANNOTATION TABLE
# ============================================================

write.csv(
  final_annotation,
  file.path(
    results_dir,
    "11D_C_FINAL_CLUSTER_ANNOTATION_MAP.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. CREATE LOOKUP VECTORS
# ============================================================

broad_lookup <- setNames(
  final_annotation$Corrected_Broad,
  final_annotation$Cluster
)

refined_lookup <- setNames(
  final_annotation$Corrected_Refined,
  final_annotation$Cluster
)

confidence_lookup <- setNames(
  final_annotation$Annotation_Confidence,
  final_annotation$Cluster
)

validation_lookup <- setNames(
  final_annotation$Validation_Source,
  final_annotation$Cluster
)


# ============================================================
# 11. MAP CORRECTED ANNOTATIONS TO EVERY CELL
# ============================================================

cell_clusters <- as.character(
  scrna$seurat_clusters
)

corrected_broad <- unname(
  broad_lookup[
    cell_clusters
  ]
)

corrected_refined <- unname(
  refined_lookup[
    cell_clusters
  ]
)

corrected_confidence <- unname(
  confidence_lookup[
    cell_clusters
  ]
)

corrected_validation <- unname(
  validation_lookup[
    cell_clusters
  ]
)


# ============================================================
# 12. STRICT MAPPING CHECK
# ============================================================

if (any(is.na(corrected_broad))) {
  
  stop(
    "Some cells failed broad annotation mapping."
  )
}

if (any(is.na(corrected_refined))) {
  
  stop(
    "Some cells failed refined annotation mapping."
  )
}

if (any(is.na(corrected_confidence))) {
  
  stop(
    "Some cells failed confidence mapping."
  )
}

if (length(corrected_broad) != ncol(scrna)) {
  
  stop(
    "Corrected annotation vector length does not equal cell count."
  )
}

cat(
  "All",
  ncol(scrna),
  "cells mapped successfully.\n"
)


# ============================================================
# 13. WRITE NEW METADATA COLUMNS
#
# ORIGINAL COLUMNS ARE NOT TOUCHED.
# ============================================================

scrna$corrected_broad_cell_type <-
  corrected_broad

scrna$corrected_refined_cell_type <-
  corrected_refined

scrna$corrected_annotation_confidence <-
  corrected_confidence

scrna$corrected_annotation_validation <-
  corrected_validation

scrna$corrected_annotation_version <-
  "GSE308859_TAC_Atlas_v2_Corrected"

scrna$corrected_annotation_step <-
  "11D-C"

scrna$corrected_annotation_date <-
  as.character(Sys.Date())


# ============================================================
# 14. VERIFY ORIGINAL ANNOTATIONS WERE NOT ALTERED
# ============================================================

if (!identical(
  original_broad,
  as.character(
    scrna$broad_cell_type
  )
)) {
  
  stop(
    "ERROR: Original broad_cell_type was altered."
  )
}

if (!identical(
  original_refined,
  as.character(
    scrna$refined_cell_type
  )
)) {
  
  stop(
    "ERROR: Original refined_cell_type was altered."
  )
}

cat(
  "Original annotation columns verified unchanged.\n"
)


# ============================================================
# 15. VERIFY CELLS, GENES AND CLUSTERS
# ============================================================

if (!identical(
  original_cells,
  colnames(scrna)
)) {
  
  stop(
    "Cell identities/order changed."
  )
}

if (!identical(
  original_genes,
  rownames(scrna)
)) {
  
  stop(
    "Gene identities/order changed."
  )
}

if (!identical(
  original_clusters,
  as.character(
    scrna$seurat_clusters
  )
)) {
  
  stop(
    "Cluster assignments changed."
  )
}

if (ncol(scrna) != original_ncells) {
  
  stop(
    "Cell count changed."
  )
}

if (nrow(scrna) != original_ngenes) {
  
  stop(
    "Gene count changed."
  )
}

cat(
  "Cells, genes and cluster assignments verified unchanged.\n"
)


# ============================================================
# 16. CLUSTER-LEVEL VERIFICATION
# ============================================================

cluster_verification <- scrna@meta.data %>%
  
  mutate(
    Cluster =
      as.character(
        seurat_clusters
      )
  ) %>%
  
  group_by(
    Cluster
  ) %>%
  
  summarise(
    
    Cells = n(),
    
    Corrected_Broad =
      paste(
        unique(
          corrected_broad_cell_type
        ),
        collapse = "; "
      ),
    
    Corrected_Refined =
      paste(
        unique(
          corrected_refined_cell_type
        ),
        collapse = "; "
      ),
    
    Confidence =
      paste(
        unique(
          corrected_annotation_confidence
        ),
        collapse = "; "
      ),
    
    .groups = "drop"
  ) %>%
  
  arrange(
    as.numeric(Cluster)
  )

write.csv(
  cluster_verification,
  file.path(
    results_dir,
    "11D_C_CLUSTER_VERIFICATION.csv"
  ),
  row.names = FALSE
)

cat("\n============================================\n")
cat("FINAL CLUSTER ANNOTATION\n")
cat("============================================\n")

print(
  as.data.frame(
    cluster_verification
  ),
  row.names = FALSE
)


# ============================================================
# 17. ENSURE EXACTLY ONE ANNOTATION PER CLUSTER
# ============================================================

annotation_uniqueness <- scrna@meta.data %>%
  
  mutate(
    Cluster =
      as.character(
        seurat_clusters
      )
  ) %>%
  
  group_by(
    Cluster
  ) %>%
  
  summarise(
    
    N_Broad =
      n_distinct(
        corrected_broad_cell_type
      ),
    
    N_Refined =
      n_distinct(
        corrected_refined_cell_type
      ),
    
    N_Confidence =
      n_distinct(
        corrected_annotation_confidence
      ),
    
    .groups = "drop"
  )

if (any(annotation_uniqueness$N_Broad != 1)) {
  
  stop(
    "A cluster contains multiple corrected broad labels."
  )
}

if (any(annotation_uniqueness$N_Refined != 1)) {
  
  stop(
    "A cluster contains multiple corrected refined labels."
  )
}

if (any(annotation_uniqueness$N_Confidence != 1)) {
  
  stop(
    "A cluster contains multiple confidence labels."
  )
}

cat(
  "\nEach cluster has exactly one corrected annotation.\n"
)


# ============================================================
# 18. BROAD CELL-TYPE COMPOSITION
# ============================================================

broad_composition <- scrna@meta.data %>%
  
  count(
    corrected_broad_cell_type,
    name = "Cells"
  ) %>%
  
  mutate(
    Percent_of_Atlas =
      100 * Cells / sum(Cells)
  ) %>%
  
  arrange(
    desc(Cells)
  )

write.csv(
  broad_composition,
  file.path(
    results_dir,
    "11D_C_CORRECTED_BROAD_COMPOSITION.csv"
  ),
  row.names = FALSE
)

cat("\n============================================\n")
cat("CORRECTED BROAD CELL-TYPE COMPOSITION\n")
cat("============================================\n")

print(
  as.data.frame(
    broad_composition
  ),
  row.names = FALSE
)


# ============================================================
# 19. REFINED CELL-TYPE COMPOSITION
# ============================================================

refined_composition <- scrna@meta.data %>%
  
  count(
    corrected_refined_cell_type,
    name = "Cells"
  ) %>%
  
  mutate(
    Percent_of_Atlas =
      100 * Cells / sum(Cells)
  ) %>%
  
  arrange(
    desc(Cells)
  )

write.csv(
  refined_composition,
  file.path(
    results_dir,
    "11D_C_CORRECTED_REFINED_COMPOSITION.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 20. SAMPLE × BROAD CELL-TYPE TABLE
#
# DESCRIPTIVE ONLY.
# No replicate-level inference.
# ============================================================

sample_broad_counts <- scrna@meta.data %>%
  
  count(
    sample_id,
    corrected_broad_cell_type,
    name = "Cells"
  ) %>%
  
  group_by(
    sample_id
  ) %>%
  
  mutate(
    Percent_Within_Sample =
      100 * Cells / sum(Cells)
  ) %>%
  
  ungroup()

write.csv(
  sample_broad_counts,
  file.path(
    results_dir,
    "11D_C_SAMPLE_BY_BROAD_CELL_TYPE.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 21. ORIGINAL VS CORRECTED BROAD ANNOTATION
#
# This provides a permanent audit trail of what changed.
# ============================================================

old_vs_new_broad <- scrna@meta.data %>%
  
  count(
    broad_cell_type,
    corrected_broad_cell_type,
    name = "Cells"
  ) %>%
  
  arrange(
    desc(Cells)
  )

write.csv(
  old_vs_new_broad,
  file.path(
    results_dir,
    "11D_C_ORIGINAL_VS_CORRECTED_BROAD.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 22. ORIGINAL VS CORRECTED REFINED ANNOTATION
# ============================================================

old_vs_new_refined <- scrna@meta.data %>%
  
  count(
    refined_cell_type,
    corrected_refined_cell_type,
    name = "Cells"
  ) %>%
  
  arrange(
    desc(Cells)
  )

write.csv(
  old_vs_new_refined,
  file.path(
    results_dir,
    "11D_C_ORIGINAL_VS_CORRECTED_REFINED.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 23. UMAP — CORRECTED BROAD ANNOTATION
# ============================================================

if ("umap" %in% Reductions(scrna)) {
  
  p_broad <- DimPlot(
    scrna,
    reduction = "umap",
    group.by = "corrected_broad_cell_type",
    label = TRUE,
    repel = TRUE,
    raster = FALSE
  ) +
    labs(
      title = "Corrected broad cell-type annotation",
      subtitle = "GSE308859 TAC scRNA-seq atlas"
    )
  
  ggsave(
    file.path(
      figures_dir,
      "11D_C_Corrected_Broad_Annotation_UMAP.png"
    ),
    p_broad,
    width = 13,
    height = 9,
    dpi = 400
  )
}


# ============================================================
# 24. UMAP — CORRECTED REFINED ANNOTATION
# ============================================================

if ("umap" %in% Reductions(scrna)) {
  
  p_refined <- DimPlot(
    scrna,
    reduction = "umap",
    group.by = "corrected_refined_cell_type",
    label = TRUE,
    repel = TRUE,
    raster = FALSE
  ) +
    labs(
      title = "Corrected refined cell-type annotation",
      subtitle = "Validated 29-cluster annotation"
    )
  
  ggsave(
    file.path(
      figures_dir,
      "11D_C_Corrected_Refined_Annotation_UMAP.png"
    ),
    p_refined,
    width = 15,
    height = 10,
    dpi = 400
  )
}


# ============================================================
# 25. UMAP — CLUSTER NUMBERS
#