# ============================================================
# GSE308859 TAC scRNA-seq PROJECT
# STEP 07C: FINAL CELL-TYPE ANNOTATION AND ATLAS FREEZE
#
# PURPOSE:
#   1. Freeze resolution-0.4 atlas
#   2. Preserve original atlas_cluster
#   3. Assign broad cell types
#   4. Assign refined cell types/states
#   5. Generate final annotated UMAPs
#   6. Generate annotation and composition tables
#   7. Perform final marker sanity checks
#   8. Save locked annotated Seurat object
#
# IMPORTANT:
#   - NO reclustering
#   - NO reintegration
#   - NO resolution changes
#   - NO CYP/pathway analysis
#   - NO statistical comparison of condition proportions
#
# CONDITION COMPOSITION IS DESCRIPTIVE ONLY:
#   one scRNA-seq library per condition.
# ============================================================


# ============================================================
# 0. CLEAN ENVIRONMENT
# ============================================================

rm(list = ls())
gc()

set.seed(20260918)

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
  "07C_FINAL_CELL_ATLAS"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "07C_FINAL_CELL_ATLAS"
)

notes_dir <- file.path(
  project_dir,
  "NOTES"
)


for (d in c(
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

required_packages <- c(
  "Seurat",
  "ggplot2",
  "patchwork",
  "Matrix"
)


for (pkg in required_packages) {
  
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


library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)


# ============================================================
# 3. LOAD STEP 07B OBJECT
# ============================================================

input_file <- file.path(
  clean_dir,
  "07B_Targeted_Validation_Diagnostic_Atlas.rds"
)


if (!file.exists(input_file)) {
  
  stop(
    paste0(
      "Cannot find Step 07B object:\n",
      input_file
    )
  )
}


obj <- readRDS(
  input_file
)


DefaultAssay(obj) <- "RNA"


cat("\n========================================\n")
cat("STEP 07C: FINAL CELL ATLAS\n")
cat("========================================\n")

cat(
  "Input cells:",
  ncol(obj),
  "\n"
)

cat(
  "Input genes:",
  nrow(obj),
  "\n"
)


# ============================================================
# 4. VERIFY REQUIRED METADATA
# ============================================================

required_metadata <- c(
  "atlas_cluster",
  "sample_id"
)


missing_metadata <- setdiff(
  required_metadata,
  colnames(obj@meta.data)
)


if (length(missing_metadata) > 0) {
  
  stop(
    paste(
      "Missing metadata:",
      paste(
        missing_metadata,
        collapse = ", "
      )
    )
  )
}


cluster_ids <- sort(
  unique(
    as.integer(
      as.character(
        obj$atlas_cluster
      )
    )
  )
)


expected_clusters <- 0:21


if (!identical(
  cluster_ids,
  expected_clusters
)) {
  
  stop(
    paste0(
      "Expected clusters 0:21 but found: ",
      paste(
        cluster_ids,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 5. FINAL ANNOTATION MAP
# ============================================================
#
# These labels freeze the annotation decisions from
# Steps 07 and 07B.
#
# Cluster 11:
#   Cycling myeloid cells
#
# Cluster 12:
#   Monocytes
#
# Cluster 16:
#   Activated dendritic cells
#
# Cluster 20:
#   Immature T cells
#
# Cluster 21:
#   Immature B cells
#
# ============================================================

annotation_map <- data.frame(
  
  atlas_cluster = as.character(0:21),
  
  broad_cell_type = c(
    
    # 0
    "Fibroblasts",
    
    # 1
    "Fibroblasts",
    
    # 2
    "Endothelial cells",
    
    # 3
    "Macrophages",
    
    # 4
    "Macrophages",
    
    # 5
    "T cells",
    
    # 6
    "Macrophages",
    
    # 7
    "Endothelial cells",
    
    # 8
    "Pericytes",
    
    # 9
    "Cardiomyocytes",
    
    # 10
    "Fibroblasts",
    
    # 11
    "Myeloid cells",
    
    # 12
    "Myeloid cells",
    
    # 13
    "Macrophages",
    
    # 14
    "Fibroblasts",
    
    # 15
    "T cells",
    
    # 16
    "Dendritic cells",
    
    # 17
    "Cardiomyocytes",
    
    # 18
    "Erythroid cells",
    
    # 19
    "B cells",
    
    # 20
    "T cells",
    
    # 21
    "B cells"
  ),
  
  refined_cell_type = c(
    
    # 0
    "Fibroblasts",
    
    # 1
    "Fibroblasts",
    
    # 2
    "Endothelial cells",
    
    # 3
    "Macrophages",
    
    # 4
    "Macrophages",
    
    # 5
    "T cells",
    
    # 6
    "Macrophages",
    
    # 7
    "Endothelial cells",
    
    # 8
    "Pericytes",
    
    # 9
    "Cardiomyocytes",
    
    # 10
    "Fibroblasts",
    
    # 11
    "Cycling myeloid cells",
    
    # 12
    "Monocytes",
    
    # 13
    "Macrophages",
    
    # 14
    "Fibroblasts",
    
    # 15
    "T cells",
    
    # 16
    "Activated dendritic cells",
    
    # 17
    "Cardiomyocytes",
    
    # 18
    "Erythroid cells",
    
    # 19
    "B cells",
    
    # 20
    "Immature T cells",
    
    # 21
    "Immature B cells"
  ),
  
  annotation_confidence = c(
    
    "High",   # 0
    "High",   # 1
    "High",   # 2
    "High",   # 3
    "High",   # 4
    "High",   # 5
    "High",   # 6
    "High",   # 7
    "High",   # 8
    "High",   # 9
    "High",   # 10
    "High",   # 11
    "High",   # 12
    "High",   # 13
    "High",   # 14
    "High",   # 15
    "High",   # 16
    "High",   # 17
    "High",   # 18
    "High",   # 19
    "High",   # 20
    "High"    # 21
  ),
  
  stringsAsFactors = FALSE
)


# ============================================================
# 6. VERIFY ANNOTATION MAP
# ============================================================

if (anyDuplicated(
  annotation_map$atlas_cluster
) > 0) {
  
  stop(
    "Duplicate cluster IDs detected in annotation map."
  )
}


if (!setequal(
  annotation_map$atlas_cluster,
  as.character(0:21)
)) {
  
  stop(
    "Annotation map does not contain exactly clusters 0-21."
  )
}


write.csv(
  annotation_map,
  file.path(
    results_dir,
    "07C_Final_Cluster_Annotation_Map.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 7. ADD FINAL ANNOTATIONS TO METADATA
# ============================================================

cluster_character <- as.character(
  obj$atlas_cluster
)


map_index <- match(
  cluster_character,
  annotation_map$atlas_cluster
)


if (any(is.na(map_index))) {
  
  stop(
    "At least one cell could not be matched to annotation map."
  )
}


obj$broad_cell_type <-
  annotation_map$broad_cell_type[
    map_index
  ]


obj$refined_cell_type <-
  annotation_map$refined_cell_type[
    map_index
  ]


obj$annotation_confidence <-
  annotation_map$annotation_confidence[
    map_index
  ]


# ============================================================
# 8. DEFINE FACTOR ORDER
# ============================================================

broad_order <- c(
  "Cardiomyocytes",
  "Fibroblasts",
  "Endothelial cells",
  "Pericytes",
  "Macrophages",
  "Myeloid cells",
  "Dendritic cells",
  "T cells",
  "B cells",
  "Erythroid cells"
)


refined_order <- c(
  "Cardiomyocytes",
  "Fibroblasts",
  "Endothelial cells",
  "Pericytes",
  "Macrophages",
  "Cycling myeloid cells",
  "Monocytes",
  "Activated dendritic cells",
  "T cells",
  "Immature T cells",
  "B cells",
  "Immature B cells",
  "Erythroid cells"
)


obj$broad_cell_type <- factor(
  obj$broad_cell_type,
  levels = broad_order
)


obj$refined_cell_type <- factor(
  obj$refined_cell_type,
  levels = refined_order
)


obj$atlas_cluster <- factor(
  as.character(
    obj$atlas_cluster
  ),
  levels = as.character(0:21)
)


# ============================================================
# 9. VERIFY NO MISSING FINAL ANNOTATIONS
# ============================================================

if (any(is.na(
  obj$broad_cell_type
))) {
  
  stop(
    "NA broad cell-type annotations detected."
  )
}


if (any(is.na(
  obj$refined_cell_type
))) {
  
  stop(
    "NA refined cell-type annotations detected."
  )
}


cat("\nAll cells received final annotations.\n")


# ============================================================
# 10. CELL NUMBER BY CLUSTER
# ============================================================

cluster_counts <- as.data.frame(
  table(
    obj$atlas_cluster
  ),
  stringsAsFactors = FALSE
)


colnames(
  cluster_counts
) <- c(
  "atlas_cluster",
  "Cells"
)


cluster_counts$Percent_of_Atlas <-
  100 *
  cluster_counts$Cells /
  sum(
    cluster_counts$Cells
  )


cluster_counts <- merge(
  cluster_counts,
  annotation_map,
  by = "atlas_cluster",
  all.x = TRUE,
  sort = FALSE
)


cluster_counts <- cluster_counts[
  match(
    as.character(0:21),
    cluster_counts$atlas_cluster
  ),
  ,
  drop = FALSE
]


rownames(
  cluster_counts
) <- NULL


write.csv(
  cluster_counts,
  file.path(
    results_dir,
    "07C_Final_Cell_Numbers_By_Cluster.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. BROAD CELL-TYPE COUNTS
# ============================================================

broad_counts <- as.data.frame(
  table(
    obj$broad_cell_type
  ),
  stringsAsFactors = FALSE
)


colnames(
  broad_counts
) <- c(
  "Broad_Cell_Type",
  "Cells"
)


broad_counts$Percent_of_Atlas <-
  100 *
  broad_counts$Cells /
  sum(
    broad_counts$Cells
  )


write.csv(
  broad_counts,
  file.path(
    results_dir,
    "07C_Broad_Cell_Type_Counts.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. REFINED CELL-TYPE COUNTS
# ============================================================

refined_counts <- as.data.frame(
  table(
    obj$refined_cell_type
  ),
  stringsAsFactors = FALSE
)


colnames(
  refined_counts
) <- c(
  "Refined_Cell_Type",
  "Cells"
)


refined_counts$Percent_of_Atlas <-
  100 *
  refined_counts$Cells /
  sum(
    refined_counts$Cells
  )


write.csv(
  refined_counts,
  file.path(
    results_dir,
    "07C_Refined_Cell_Type_Counts.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. SAMPLE ORDER
# ============================================================

sample_order <- c(
  "Sham",
  "TAC_2W",
  "TAC_4W",
  "TAC_6W"
)


obj$sample_id <- factor(
  as.character(
    obj$sample_id
  ),
  levels = sample_order
)


# ============================================================
# 14. CLUSTER × SAMPLE COMPOSITION
# ============================================================

cluster_sample_table <- as.data.frame(
  table(
    Cluster =
      obj$atlas_cluster,
    Sample =
      obj$sample_id
  ),
  stringsAsFactors = FALSE
)


cluster_sample_totals <- aggregate(
  Freq ~ Sample,
  data = cluster_sample_table,
  FUN = sum
)


colnames(
  cluster_sample_totals
)[2] <- "Sample_Total"


cluster_sample_table <- merge(
  cluster_sample_table,
  cluster_sample_totals,
  by = "Sample",
  all.x = TRUE,
  sort = FALSE
)


cluster_sample_table$Percent_of_Sample <-
  100 *
  cluster_sample_table$Freq /
  cluster_sample_table$Sample_Total


write.csv(
  cluster_sample_table,
  file.path(
    results_dir,
    "07C_Cluster_By_Sample_Composition.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. BROAD CELL TYPE × SAMPLE COMPOSITION
# ============================================================

broad_sample_table <- as.data.frame(
  table(
    Broad_Cell_Type =
      obj$broad_cell_type,
    Sample =
      obj$sample_id
  ),
  stringsAsFactors = FALSE
)


broad_sample_totals <- aggregate(
  Freq ~ Sample,
  data = broad_sample_table,
  FUN = sum
)


colnames(
  broad_sample_totals
)[2] <- "Sample_Total"


broad_sample_table <- merge(
  broad_sample_table,
  broad_sample_totals,
  by = "Sample",
  all.x = TRUE,
  sort = FALSE
)


broad_sample_table$Percent_of_Sample <-
  100 *
  broad_sample_table$Freq /
  broad_sample_table$Sample_Total


write.csv(
  broad_sample_table,
  file.path(
    results_dir,
    "07C_Broad_Cell_Type_By_Sample_Composition.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. REFINED CELL TYPE × SAMPLE COMPOSITION
# ============================================================

refined_sample_table <- as.data.frame(
  table(
    Refined_Cell_Type =
      obj$refined_cell_type,
    Sample =
      obj$sample_id
  ),
  stringsAsFactors = FALSE
)


refined_sample_totals <- aggregate(
  Freq ~ Sample,
  data = refined_sample_table,
  FUN = sum
)


colnames(
  refined_sample_totals
)[2] <- "Sample_Total"


refined_sample_table <- merge(
  refined_sample_table,
  refined_sample_totals,
  by = "Sample",
  all.x = TRUE,
  sort = FALSE
)


refined_sample_table$Percent_of_Sample <-
  100 *
  refined_sample_table$Freq /
  refined_sample_table$Sample_Total


write.csv(
  refined_sample_table,
  file.path(
    results_dir,
    "07C_Refined_Cell_Type_By_Sample_Composition.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. FINAL BROAD UMAP
# ============================================================

p_broad <- DimPlot(
  object = obj,
  reduction = "umap.unintegrated",
  group.by = "broad_cell_type",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.25
) +
  labs(
    title =
      "Final broad cardiac cell atlas"
  ) +
  theme_classic(
    base_size = 13
  )


ggsave(
  filename = file.path(
    figures_dir,
    "07C_Final_Broad_Cell_Type_UMAP.png"
  ),
  plot = p_broad,
  width = 12,
  height = 9,
  dpi = 600
)


# ============================================================
# 18. FINAL REFINED UMAP
# ============================================================

p_refined <- DimPlot(
  object = obj,
  reduction = "umap.unintegrated",
  group.by = "refined_cell_type",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.25
) +
  labs(
    title =
      "Final refined cardiac cell atlas"
  ) +
  theme_classic(
    base_size = 13
  )


ggsave(
  filename = file.path(
    figures_dir,
    "07C_Final_Refined_Cell_Type_UMAP.png"
  ),
  plot = p_refined,
  width = 13,
  height = 9,
  dpi = 600
)


# ============================================================
# 19. ORIGINAL CLUSTER UMAP
# ============================================================

p_cluster <- DimPlot(
  object = obj,
  reduction = "umap.unintegrated",
  group.by = "atlas_cluster",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.2
) +
  labs(
    title =
      "Frozen atlas clusters - resolution 0.4"
  ) +
  theme_classic(
    base_size = 13
  )


ggsave(
  filename = file.path(
    figures_dir,
    "07C_Frozen_Atlas_Clusters_UMAP.png"
  ),
  plot = p_cluster,
  width = 12,
  height = 9,
  dpi = 600
)


# ============================================================
# 20. REFINED ATLAS SPLIT BY SAMPLE
# ============================================================

p_refined_split <- DimPlot(
  object = obj,
  reduction = "umap.unintegrated",
  group.by = "refined_cell_type",
  split.by = "sample_id",
  pt.size = 0.15,