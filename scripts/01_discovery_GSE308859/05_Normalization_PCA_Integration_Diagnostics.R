# ============================================================
# GSE308859 TAC scRNA-seq PROJECT
# STEP 05: NORMALIZATION, PCA & INTEGRATION DIAGNOSTICS
#
# FINAL QC INPUT:
#   04B_All_Samples_RefinedQC_Singlets_List.rds
#
# PURPOSE:
#   1. Normalize each library independently
#   2. Identify variable features independently
#   3. Select shared integration features
#   4. Construct an UNINTEGRATED merged object
#   5. Perform PCA
#   6. Examine PC significance / variance
#   7. Construct a diagnostic unintegrated UMAP
#   8. Examine sample mixing and biological structure
#   9. Examine CYP genes in the unintegrated space
#
# IMPORTANT:
#   NO integration is performed in this script.
#   NO final clustering is performed in this script.
#
# We first determine whether integration is appropriate.
# ============================================================


# ============================================================
# 0. CLEAN ENVIRONMENT
# ============================================================

rm(list = ls())
gc()

set.seed(20260918)

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
  "05_NORMALIZATION_PCA"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "05_NORMALIZATION_PCA"
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

required_packages <- c(
  "Seurat",
  "ggplot2",
  "patchwork"
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
  "ggplot2:",
  as.character(
    packageVersion("ggplot2")
  ),
  "\n"
)


# ============================================================
# 4. LOAD FINAL REFINED-QC OBJECTS
# ============================================================

input_file <- file.path(
  clean_dir,
  "04B_All_Samples_RefinedQC_Singlets_List.rds"
)


if (!file.exists(input_file)) {
  
  stop(
    paste0(
      "Cannot find final Step-04B object:\n",
      input_file
    )
  )
}


seurat_list <- readRDS(
  input_file
)


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


# ============================================================
# 5. VERIFY FINAL CELL NUMBERS
# ============================================================

cat("\n========================================\n")
cat("FINAL QC INPUT CELL NUMBERS\n")
cat("========================================\n")


cell_number_list <- list()


for (sample_name in sample_order) {
  
  obj <- seurat_list[[sample_name]]
  
  
  cell_number_list[[sample_name]] <- data.frame(
    
    Sample =
      sample_name,
    
    Cells =
      ncol(obj),
    
    Genes =
      nrow(obj),
    
    stringsAsFactors = FALSE
  )
}


cell_numbers <- do.call(
  rbind,
  cell_number_list
)


rownames(
  cell_numbers
) <- NULL


print(
  cell_numbers
)


cat(
  "\nTotal cells:",
  sum(
    cell_numbers$Cells
  ),
  "\n"
)


write.csv(
  cell_numbers,
  file.path(
    results_dir,
    "05_Input_Cell_Numbers.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 6. VERIFY METADATA
# ============================================================

cat("\n========================================\n")
cat("METADATA CHECK\n")
cat("========================================\n")


for (sample_name in sample_order) {
  
  obj <- seurat_list[[sample_name]]
  
  
  cat(
    "\n--- ",
    sample_name,
    " ---\n",
    sep = ""
  )
  
  
  print(
    colnames(
      obj@meta.data
    )
  )
  
  
  if (!"sample_id" %in% colnames(obj@meta.data)) {
    
    obj$sample_id <- sample_name
  }
  
  
  if (!"condition" %in% colnames(obj@meta.data)) {
    
    obj$condition <- sample_name
  }
  
  
  seurat_list[[sample_name]] <- obj
}


# ============================================================
# 7. NORMALIZE EACH SAMPLE INDEPENDENTLY
# ============================================================

cat("\n========================================\n")
cat("NORMALIZING EACH SAMPLE\n")
cat("========================================\n")


for (sample_name in sample_order) {
  
  cat(
    "\nNormalizing ",
    sample_name,
    "...\n",
    sep = ""
  )
  
  
  obj <- seurat_list[[sample_name]]
  
  
  DefaultAssay(obj) <- "RNA"
  
  
  obj <- NormalizeData(
    object = obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  
  
  seurat_list[[sample_name]] <- obj
}


# ============================================================
# 8. IDENTIFY VARIABLE FEATURES PER SAMPLE
# ============================================================

cat("\n========================================\n")
cat("IDENTIFYING VARIABLE FEATURES\n")
cat("========================================\n")


for (sample_name in sample_order) {
  
  cat(
    "\nVariable features: ",
    sample_name,
    "\n",
    sep = ""
  )
  
  
  obj <- seurat_list[[sample_name]]
  
  
  obj <- FindVariableFeatures(
    object = obj,
    selection.method = "vst",
    nfeatures = 3000,
    verbose = FALSE
  )
  
  
  cat(
    "Variable features identified:",
    length(
      VariableFeatures(obj)
    ),
    "\n"
  )
  
  
  seurat_list[[sample_name]] <- obj
}


# ============================================================
# 9. SAVE VARIABLE FEATURES FOR EACH SAMPLE
# ============================================================

for (sample_name in sample_order) {
  
  variable_gene_table <- data.frame(
    
    Gene =
      VariableFeatures(
        seurat_list[[sample_name]]
      ),
    
    stringsAsFactors = FALSE
  )
  
  
  write.csv(
    variable_gene_table,
    file.path(
      results_dir,
      paste0(
        "05_",
        sample_name,
        "_Variable_Features.csv"
      )
    ),
    row.names = FALSE
  )
}


# ============================================================
# 10. SELECT SHARED FEATURES
# ============================================================

shared_features <- SelectIntegrationFeatures(
  object.list = seurat_list,
  nfeatures = 3000
)


cat("\n========================================\n")
cat("SHARED FEATURES\n")
cat("========================================\n")

cat(
  "Selected shared features:",
  length(shared_features),
  "\n"
)


write.csv(
  data.frame(
    Gene = shared_features
  ),
  file.path(
    results_dir,
    "05_Shared_Features_3000.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. VARIABLE-FEATURE OVERLAP
# ============================================================

variable_sets <- lapply(
  seurat_list,
  VariableFeatures
)


overlap_matrix <- matrix(
  NA_integer_,
  nrow = length(sample_order),
  ncol = length(sample_order),
  dimnames = list(
    sample_order,
    sample_order
  )
)


for (i in seq_along(sample_order)) {
  
  for (j in seq_along(sample_order)) {
    
    overlap_matrix[i, j] <- length(
      intersect(
        variable_sets[[sample_order[i]]],
        variable_sets[[sample_order[j]]]
      )
    )
  }
}


cat("\n========================================\n")
cat("VARIABLE FEATURE OVERLAP\n")
cat("========================================\n")

print(
  overlap_matrix
)


write.csv(
  overlap_matrix,
  file.path(
    results_dir,
    "05_Variable_Feature_Overlap.csv"
  ),
  row.names = TRUE
)


# ============================================================
# 12. SAVE NORMALIZED INDIVIDUAL OBJECTS
# ============================================================

saveRDS(
  seurat_list,
  file = file.path(
    clean_dir,
    "05_Normalized_Individual_Samples_List.rds"
  ),
  compress = TRUE
)


# ============================================================
# 13. MERGE NORMALIZED SAMPLES — UNINTEGRATED
# ============================================================

cat("\n========================================\n")
cat("MERGING SAMPLES - UNINTEGRATED\n")
cat("========================================\n")


merged_obj <- merge(
  x = seurat_list[["Sham"]],
  y = list(
    seurat_list[["TAC_2W"]],
    seurat_list[["TAC_4W"]],
    seurat_list[["TAC_6W"]]
  ),
  add.cell.ids = c(
    "Sham",
    "TAC2W",
    "TAC4W",
    "TAC6W"
  ),
  project = "GSE308859_TAC"
)


DefaultAssay(
  merged_obj
) <- "RNA"


cat(
  "Merged cells:",
  ncol(merged_obj),
  "\n"
)

cat(
  "Merged genes:",
  nrow(merged_obj),
  "\n"
)


# ============================================================
# 14. JOIN RNA LAYERS IF REQUIRED BY SEURAT v5
# ============================================================

# merge() in Seurat v5 may retain separate layers
# for each original sample.
#
# JoinLayers() consolidates them for downstream analysis.

if (
  "JoinLayers" %in%
  getNamespaceExports("SeuratObject")
) {
  
  merged_obj <- JoinLayers(
    merged_obj,
    assay = "RNA"
  )
}


# ============================================================
# 15. SET VARIABLE FEATURES OF MERGED OBJECT
# ============================================================

VariableFeatures(
  merged_obj
) <- shared_features


cat(
  "Merged object variable features:",
  length(
    VariableFeatures(merged_obj)
  ),
  "\n"
)


# ============================================================
# 16. SCALE SHARED FEATURES
# ============================================================

cat("\n========================================\n")
cat("SCALING SHARED FEATURES\n")
cat("========================================\n")


merged_obj <- ScaleData(
  object = merged_obj,
  features = shared_features,
  verbose = FALSE
)


# ============================================================
# 17. PCA
# ============================================================

cat("\n========================================\n")
cat("RUNNING PCA\n")
cat("========================================\n")


merged_obj <- RunPCA(
  object = merged_obj,
  features = shared_features,
  npcs = 50,
  verbose = FALSE
)


cat(
  "PCA dimensions generated:",
  ncol(
    Embeddings(
      merged_obj,
      reduction = "pca"
    )
  ),
  "\n"
)


# ============================================================
# 18. PCA STANDARD DEVIATIONS AND VARIANCE
# ============================================================

pca_stdev <- Stdev(
  merged_obj,
  reduction = "pca"
)


variance <- pca_stdev^2


percent_variance <- 100 *
  variance /
  sum(variance)


cumulative_variance <- cumsum(
  percent_variance
)


pca_variance_table <- data.frame(
  
  PC =
    seq_along(
      pca_stdev
    ),
  
  Standard_Deviation =
    pca_stdev,
  
  Percent_Variance =
    percent_variance,
  
  Cumulative_Percent_Variance =
    cumulative_variance
)


cat("\n========================================\n")
cat("PCA VARIANCE - FIRST 30 PCs\n")
cat("========================================\n")


print(
  head(
    pca_variance_table,
    30
  )
)


write.csv(
  pca_variance_table,
  file.path(
    results_dir,
    "05_PCA_Variance.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. ELBOW PLOT
# ============================================================

p_elbow <- ElbowPlot(
  merged_obj,
  ndims = 50
) +
  labs(
    title =
      "PCA elbow plot - unintegrated data"
  ) +
  theme_classic(
    base_size = 13
  )


ggsave(
  filename = file.path(
    figures_dir,
    "05_Unintegrated_PCA_Elbow.png"
  ),
  plot = p_elbow,
  width = 8,
  height = 6,
  dpi = 600
)


# ============================================================
# 20. PCA BY SAMPLE
# ============================================================

p_pca_sample <- DimPlot(
  merged_obj,
  reduction = "pca",
  group.by = "sample_id",
  pt.size = 0.15
) +
  labs(
    title =
      "Unintegrated PCA by sample"
  ) +
  theme_classic(
    base_size = 13
  )


ggsave(
  filename = file.path(
    figures_dir,
    "05_Unintegrated_PCA_By_Sample.png"
  ),
  plot = p_pca_sample,
  width = 9,
  height = 7,
  dpi = 600
)


# ============================================================
# 21. PCA SPLIT BY SAMPLE
# ============================================================

p_pca_split <- DimPlot(
  merged_obj,
  reduction = "pca",
  group.by = "sample_id",
  split.by = "sample_id",
  pt.size = 0.10,
  ncol = 2
) +
  labs(
    title =
      "Unintegrated PCA - sample-specific distributions"
  )


ggsave(
  filename = file.path(
    figures_dir,
    "05_Unintegrated_PCA_Split_By_Sample.png"
  ),
  plot = p_pca_split,
  width = 12,
  height = 10,
  dpi = 600
)


# ============================================================
# 22. PCA LOADINGS
# ============================================================

pca_loadings <- Loadings(
  merged_obj,
  reduction = "pca"
)


top_loading_list <- list()

loading_counter <- 1


for (pc_number in 1:20) {
  
  pc_name <- paste0(
    "PC_",
    pc_number
  )
  
  
  loadings_vector <- pca_loadings[
    ,
    pc_name
  ]
  
  
  top_positive <- names(
    sort(
      loadings_vector,
      decreasing = TRUE
    )
  )[1:20]
  
  
  top_negative <- names(
    sort(
      loadings_vector,
      decreasing = FALSE
    )
  )[1:20]
  
  
  top_loading_list[[loading_counter]] <- data.frame(
    
    PC =
      pc_number,
    
    Direction =
      "Positive",
    
    Gene =
      top_positive,
    
    Loading =
      loadings_vector[
        top_positive
      ],
    
    stringsAsFactors = FALSE
  )
  
  
  loading_counter <- loading_counter + 1
  
  
  top_loading_list[[loading_counter]] <- data.frame(
    
    PC =
      pc_number,
    
    Direction =
      "Negative",
    
    Gene =
      top_negative,
    
    Loading =
      loadings_vector[
        top_negative
      ],
    
    stringsAsFactors = FALSE
  )
  
  
  loading_counter <- loading_counter + 1
}


top_loadings <- do.call(
  rbind,
  top_loading_list
)


rownames(
  top_loadings
) <- NULL


write.csv(
  top_loadings,
  file.path(
    results_dir,
    "05_Top_PCA_Loadings_PC1_to_PC20.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 23. HEATMAP OF EARLY PCs
# ============================================================

# Downsample cells for visualization only.
# This does NOT alter the analytical object.

set.seed(20260918)


heatmap_cells <- sample(
  colnames(merged_obj),
  size = min(
    3000,
    ncol(merged_obj)
  )
)


p_pca_heatmap <- DimHeatmap(
  merged_obj,
  dims = 1:12,
  cells = heatmap_cells,
  balanced = TRUE,
  fast = FALSE
)


ggsave(
  filename = file.path(
    figures_dir,
    "05_PCA_Heatmap_PC1_to_PC12.png"
  ),
  plot = p_pca_heatmap,
  width = 14,
  height = 18,
  dpi = 400
)


# ============================================================
# 24. RUN DIAGNOSTIC UNINTEGRATED UMAP
# ============================================================

# This is NOT the final manuscript UMAP.
#
# We use 30 PCs provisionally so that we can examine
# whether sample/library separation dominates the data.
#
# The final PC choice will be made after reviewing
# the PCA diagnostics.

cat("\n========================================\n")
cat("RUNNING DIAGNOSTIC UNINTEGRATED UMAP\n")
cat("========================================\n")


merged_obj <- RunUMAP(
  object = merged_obj,
  reduction = "pca",
  dims = 1:30,
  n.neighbors = 30,
  min.dist = 0.3,
  metric = "cosine",
  reduction.name = "umap.unintegrated",
  reduction.key = "UMAPUNINT_",
  seed.use = 20260918,
  verbose = FALSE
)


# ============================================================
# 25. UNINTEGRATED UMAP BY SAMPLE
# ============================================================