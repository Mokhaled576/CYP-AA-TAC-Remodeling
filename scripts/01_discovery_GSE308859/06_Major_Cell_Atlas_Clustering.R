# ============================================================
# GSE308859 TAC scRNA-seq PROJECT
# STEP 06: MAJOR CELL ATLAS & CLUSTERING STABILITY
#
# INPUT:
#   05_Unintegrated_Normalized_PCA_UMAP.rds
#
# PURPOSE:
#   1. Use the accepted unintegrated PCA representation
#   2. Construct graph using PC1-PC30
#   3. Test multiple clustering resolutions
#   4. Quantify cluster stability / granularity
#   5. Select a provisional broad-atlas resolution
#   6. Generate a clean UMAP
#   7. Examine cluster composition by condition
#   8. Find broad cluster markers
#   9. Examine canonical cardiac lineage markers
#  10. Produce annotation diagnostics
#
# IMPORTANT:
#   - NO integration
#   - NO CYP differential expression
#   - NO final biological labels yet
#   - NO cells are removed
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
  "06_MAJOR_CELL_ATLAS"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "06_MAJOR_CELL_ATLAS"
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


# ============================================================
# 4. LOAD STEP-05 OBJECT
# ============================================================

input_file <- file.path(
  clean_dir,
  "05_Unintegrated_Normalized_PCA_UMAP.rds"
)


if (!file.exists(input_file)) {
  
  stop(
    paste0(
      "Cannot find Step-05 object:\n",
      input_file
    )
  )
}


obj <- readRDS(
  input_file
)


DefaultAssay(obj) <- "RNA"


cat("\n========================================\n")
cat("INPUT OBJECT\n")
cat("========================================\n")

cat(
  "Cells:",
  ncol(obj),
  "\n"
)

cat(
  "Genes:",
  nrow(obj),
  "\n"
)


if (!"pca" %in% names(obj@reductions)) {
  
  stop(
    "PCA reduction is missing from Step-05 object."
  )
}


if (!"umap.unintegrated" %in% names(obj@reductions)) {
  
  stop(
    "Unintegrated UMAP is missing from Step-05 object."
  )
}


# ============================================================
# 5. VERIFY SAMPLE METADATA
# ============================================================

sample_order <- c(
  "Sham",
  "TAC_2W",
  "TAC_4W",
  "TAC_6W"
)


if (!"sample_id" %in% colnames(obj@meta.data)) {
  
  stop(
    "sample_id metadata is missing."
  )
}


obj$sample_id <- factor(
  as.character(obj$sample_id),
  levels = sample_order
)


cat("\nSample cell numbers:\n")

print(
  table(
    obj$sample_id
  )
)


# ============================================================
# 6. REMOVE OLD DIAGNOSTIC GRAPH IF PRESENT
# ============================================================

# We construct a dedicated Step-06 graph.
# Existing diagnostic graphs from Step 05 are not used
# for final clustering.


# ============================================================
# 7. BUILD STEP-06 NEIGHBOR GRAPH
# ============================================================

cat("\n========================================\n")
cat("BUILDING STEP-06 NEIGHBOR GRAPH\n")
cat("PCs 1-30\n")
cat("========================================\n")


obj <- FindNeighbors(
  object = obj,
  reduction = "pca",
  dims = 1:30,
  k.param = 30,
  graph.name = c(
    "atlas_nn",
    "atlas_snn"
  ),
  verbose = FALSE
)


# ============================================================
# 8. TEST MULTIPLE CLUSTERING RESOLUTIONS
# ============================================================

# We do not assume that one arbitrary resolution is correct.
#
# These resolutions span broad to moderately granular
# clustering appropriate for major-cell-type annotation.


resolutions <- c(
  0.1,
  0.2,
  0.3,
  0.4,
  0.5,
  0.6,
  0.8,
  1.0
)


cat("\n========================================\n")
cat("TESTING CLUSTERING RESOLUTIONS\n")
cat("========================================\n")


for (res in resolutions) {
  
  cat(
    "Resolution:",
    res,
    "\n"
  )
  
  
  obj <- FindClusters(
    object = obj,
    graph.name = "atlas_snn",
    resolution = res,
    algorithm = 1,
    random.seed = 20260918,
    verbose = FALSE
  )
  
  
  generated_column <- paste0(
    "atlas_snn_res.",
    res
  )
  
  
  desired_column <- paste0(
    "atlas_res_",
    gsub(
      "\\.",
      "_",
      as.character(res)
    )
  )
  
  
  if (!generated_column %in% colnames(obj@meta.data)) {
    
    stop(
      paste(
        "Expected clustering column not found:",
        generated_column
      )
    )
  }
  
  
  obj[[desired_column]] <- obj[[generated_column]]
}


# ============================================================
# 9. COUNT CLUSTERS AT EACH RESOLUTION
# ============================================================

resolution_summary_list <- list()


for (res in resolutions) {
  
  column_name <- paste0(
    "atlas_res_",
    gsub(
      "\\.",
      "_",
      as.character(res)
    )
  )
  
  
  cluster_vector <- as.character(
    obj@meta.data[[column_name]]
  )
  
  
  cluster_sizes <- table(
    cluster_vector
  )
  
  
  resolution_summary_list[[as.character(res)]] <- data.frame(
    
    Resolution =
      res,
    
    Number_of_Clusters =
      length(
        cluster_sizes
      ),
    
    Smallest_Cluster =
      min(
        cluster_sizes
      ),
    
    Median_Cluster_Size =
      median(
        as.numeric(
          cluster_sizes
        )
      ),
    
    Largest_Cluster =
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


resolution_summary <- do.call(
  rbind,
  resolution_summary_list
)


rownames(
  resolution_summary
) <- NULL


cat("\n========================================\n")
cat("RESOLUTION SUMMARY\n")
cat("========================================\n")

print(
  resolution_summary
)


write.csv(
  resolution_summary,
  file.path(
    results_dir,
    "06_Clustering_Resolution_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. SAVE CLUSTER SIZE TABLES FOR EVERY RESOLUTION
# ============================================================

cluster_size_all_list <- list()

cluster_size_counter <- 1


for (res in resolutions) {
  
  column_name <- paste0(
    "atlas_res_",
    gsub(
      "\\.",
      "_",
      as.character(res)
    )
  )
  
  
  cluster_sizes <- table(
    obj@meta.data[[column_name]]
  )
  
  
  temp_df <- data.frame(
    
    Resolution =
      res,
    
    Cluster =
      names(
        cluster_sizes
      ),
    
    Cells =
      as.numeric(
        cluster_sizes
      ),
    
    Percent_of_All_Cells =
      100 *
      as.numeric(
        cluster_sizes
      ) /
      ncol(obj),
    
    stringsAsFactors = FALSE
  )
  
  
  cluster_size_all_list[[cluster_size_counter]] <- temp_df
  
  cluster_size_counter <- cluster_size_counter + 1
}


cluster_size_all <- do.call(
  rbind,
  cluster_size_all_list
)


rownames(
  cluster_size_all
) <- NULL


write.csv(
  cluster_size_all,
  file.path(
    results_dir,
    "06_All_Resolution_Cluster_Sizes.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. CREATE UMAP PANELS FOR ALL RESOLUTIONS
# ============================================================

resolution_plot_list <- list()


for (res in resolutions) {
  
  column_name <- paste0(
    "atlas_res_",
    gsub(
      "\\.",
      "_",
      as.character(res)
    )
  )
  
  
  p <- DimPlot(
    obj,
    reduction = "umap.unintegrated",
    group.by = column_name,
    label = TRUE,
    repel = TRUE,
    pt.size = 0.08
  ) +
    labs(
      title = paste0(
        "Resolution ",
        res
      )
    ) +
    NoLegend()
  
  
  resolution_plot_list[[as.character(res)]] <- p
}


resolution_panel <- wrap_plots(
  resolution_plot_list,
  ncol = 2
)


ggsave(
  filename = file.path(
    figures_dir,
    "06_Clustering_Resolution_Comparison.png"
  ),
  plot = resolution_panel,
  width = 16,
  height = 24,
  dpi = 500
)


# ============================================================
# 12. CLUSTER TRANSITION TABLES
# ============================================================

# These tables show how clusters split when resolution
# increases.


transition_pairs <- list(
  c(0.1, 0.2),
  c(0.2, 0.3),
  c(0.3, 0.4),
  c(0.4, 0.5),
  c(0.5, 0.6),
  c(0.6, 0.8),
  c(0.8, 1.0)
)


for (pair in transition_pairs) {
  
  res1 <- pair[1]
  res2 <- pair[2]
  
  
  col1 <- paste0(
    "atlas_res_",
    gsub(
      "\\.",
      "_",
      as.character(res1)
    )
  )
  
  
  col2 <- paste0(
    "atlas_res_",
    gsub(
      "\\.",
      "_",
      as.character(res2)
    )
  )
  
  
  transition_table <- as.data.frame.matrix(
    table(
      obj@meta.data[[col1]],
      obj@meta.data[[col2]]
    )
  )
  
  
  transition_table$Cluster_From <-
    rownames(
      transition_table
    )
  
  
  transition_table <- transition_table[
    ,
    c(
      "Cluster_From",
      setdiff(
        colnames(transition_table),
        "Cluster_From"
      )
    ),
    drop = FALSE
  ]
  
  
  rownames(
    transition_table
  ) <- NULL
  
  
  write.csv(
    transition_table,
    file.path(
      results_dir,
      paste0(
        "06_Cluster_Transition_",
        res1,
        "_to_",
        res2,
        ".csv"
      )
    ),
    row.names = FALSE
  )
}


# ============================================================
# 13. PROVISIONAL BROAD-ATLAS RESOLUTION
# ============================================================

# IMPORTANT:
#
# Resolution 0.4 is selected only as a PROVISIONAL
# broad-atlas working resolution.
#
# We will inspect:
#   - resolution comparison
#   - marker specificity
#   - lineage coherence
#   - sample composition
#
# before accepting it.
#
# If the diagnostics show over- or under-clustering,
# we will change it.


provisional_resolution <- 0.4

provisional_column <- "atlas_res_0_4"


if (!provisional_column %in% colnames(obj@meta.data)) {
  
  stop(
    "Provisional clustering column atlas_res_0_4 is missing."
  )
}


obj$atlas_cluster <- factor(
  as.character(
    obj@meta.data[[provisional_column]]
  )
)


Idents(obj) <- "atlas_cluster"


cat("\n========================================\n")
cat("PROVISIONAL BROAD ATLAS\n")
cat("Resolution = 0.4\n")
cat("========================================\n")


print(
  table(
    Idents(obj)
  )
)


# ============================================================
# 14. PROVISIONAL CLUSTER SIZE TABLE
# ============================================================

provisional_sizes <- as.data.frame(
  table(
    Cluster =
      obj$atlas_cluster
  )
)


colnames(
  provisional_sizes
) <- c(
  "Cluster",
  "Cells"
)


provisional_sizes$Percent_of_All_Cells <-
  100 *
  provisional_sizes$Cells /
  ncol(obj)


provisional_sizes <- provisional_sizes[
  order(
    as.numeric(
      as.character(
        provisional_sizes$Cluster
      )
    )
  ),
  ,
  drop = FALSE
]


rownames(
  provisional_sizes
) <- NULL


cat("\nProvisional cluster sizes:\n")

print(
  provisional_sizes
)


write.csv(
  provisional_sizes,
  file.path(
    results_dir,
    "06_Provisional_Cluster_Sizes_Res0.4.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. PROVISIONAL ATLAS UMAP
# ============================================================

p_atlas <- DimPlot(
  obj,
  reduction = "umap.unintegrated",
  group.by = "atlas_cluster",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.12
) +
  labs(
    title =
      "Provisional major-cell atlas - resolution 0.4"
  ) +
  theme_classic(
    base_size = 13
  )


ggsave(
  filename = file.path(
    figures_dir,
    "06_Provisional_Atlas_UMAP_Res0.4.png"
  ),
  plot = p_atlas,
  width = 10,
  height = 8,
  dpi = 600
)


# ============================================================
# 16. PROVISIONAL ATLAS SPLIT BY SAMPLE
# ============================================================

p_atlas_split <- DimPlot(
  obj,
  reduction = "umap.unintegrated",
  group.by = "atlas_cluster",
  split.by = "sample_id",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.08,
  ncol = 2
)


ggsave(
  filename = file.path(
    figures_dir,
    "06_Provisional_Atlas_Split_By_Sample.png"
  ),
  plot = p_atlas_split,
  width = 16,
  height = 12,
  dpi = 600
)


# ============================================================
# 17. CLUSTER x SAMPLE CELL COUNTS
# ============================================================

cluster_sample_counts <- as.data.frame.matrix(
  table(
    obj$atlas_cluster,
    obj$sample_id
  )
)


cluster_sample_counts$Cluster <-
  rownames(
    cluster_sample_counts
  )


cluster_sample_counts <- cluster_sample_counts[
  ,
  c(
    "Cluster",
    sample_order
  ),
  drop = FALSE
]


rownames(
  cluster_sample_counts
) <- NULL


cat("\n========================================\n")
cat("CLUSTER x SAMPLE COUNTS\n")
cat("========================================\n")

print(
  cluster_sample_counts
)


write.csv(
  cluster_sample_counts,
  file.path(
    results_dir,
    "06_Cluster_By_Sample_Counts.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. WITHIN-CLUSTER SAMPLE PROPORTIONS
# ============================================================

cluster_sample_matrix <- table(
  obj$atlas_cluster,
  obj$sample_id
)


cluster_sample_proportions <- prop.table(
  cluster_sample_matrix,
  margin = 1
)


cluster_sample_prop_df <- as.data.frame.matrix(
  cluster_sample_proportions
)


cluster_sample_prop_df$Cluster <-
  rownames(
    cluster_sample_prop_df
  )


cluster_sample_prop_df <- cluster_sample_prop_df[
  ,
  c(
    "Cluster",
    sample_order
  ),
  drop = FALSE
]


rownames(
  cluster_sample_prop_df
) <- NULL


for (sample_name in sample_order) {
  
  cluster_sample_prop_df[[sample_name]] <-
    100 *
    cluster_sample_prop_df[[sample_name]]
}


write.csv(
  cluster_sample_prop_df,
  file.path(
    results_dir,
    "06_Within_Cluster_Sample_Percentages.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. CLUSTER COMPOSITION BARPLOT
# ============================================================

composition_long <- as.data.frame(
  cluster_sample_matrix
)


colnames(
  composition_long
) <- c(
  "Cluster",
  "Sample",
  "Cells"
)


composition_long$Cluster <- factor(
  composition_long$Cluster,
  levels = levels(
    obj$atlas_cluster
  )
)


composition_long$Sample <- factor(
  composition_long$Sample,
  levels = sample_order
)

