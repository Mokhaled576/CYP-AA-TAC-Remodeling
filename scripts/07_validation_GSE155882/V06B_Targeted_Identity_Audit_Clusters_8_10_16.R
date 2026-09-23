# ============================================================
# V06B — TARGETED IDENTITY AUDIT
# CLUSTERS 8, 10, AND 16
# ============================================================
#
# PURPOSE
# -------
# Resolve the three remaining ambiguous resolution-0.6 clusters:
#
#   Cluster 8
#   Cluster 10
#   Cluster 16
#
# IMPORTANT
# ---------
# This is an annotation evidence audit ONLY.
#
# NO:
#   - cell filtering
#   - reclustering
#   - integration
#   - graph rebuilding
#   - UMAP recalculation
#   - Sham-vs-TAC DE
#   - CYP/AA validation
#   - final annotation freezing
#
# Existing V06 clusters and frozen V04/V05 reductions are used.
#
# CYP/AA pathway genes are excluded from marker evidence used
# for biological annotation.
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
# 2. LOAD PACKAGES
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

results_dir <- file.path(
  project_root,
  "RESULTS",
  "V06B_TARGETED_IDENTITY_AUDIT"
)

figures_dir <- file.path(
  project_root,
  "FIGURES",
  "V06B_TARGETED_IDENTITY_AUDIT"
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
# 4. LOAD FROZEN V06 OBJECT
# ============================================================

if (!file.exists(input_file)) {
  stop("V06 input object does not exist.")
}

obj <- readRDS(input_file)

cat("\nV06 object loaded successfully.\n")
cat("Cells:", ncol(obj), "\n")
cat("Genes:", nrow(obj), "\n")


# ============================================================
# 5. INPUT INTEGRITY AUDIT
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

input_n <- ncol(obj)
input_cells <- colnames(obj)

cat("\nAvailable reductions:\n")
print(Reductions(obj))


# ============================================================
# 6. TARGET CLUSTERS
# ============================================================

target_clusters <- c(
  "8",
  "10",
  "16"
)

available_clusters <- unique(
  as.character(obj$V06_cluster_0.6)
)

missing_targets <- setdiff(
  target_clusters,
  available_clusters
)

if (length(missing_targets) > 0) {
  stop(
    paste(
      "Target clusters missing:",
      paste(missing_targets, collapse = ", ")
    )
  )
}

cat("\nTarget clusters confirmed:\n")
print(target_clusters)


# ============================================================
# 7. TARGET CLUSTER CELL COUNTS
# ============================================================

target_counts <- as.data.frame(
  table(
    Cluster = obj$V06_cluster_0.6
  )
)

target_counts$Cluster <- as.character(
  target_counts$Cluster
)

target_counts <- target_counts[
  target_counts$Cluster %in% target_clusters,
]

colnames(target_counts)[2] <- "Cells"

write.csv(
  target_counts,
  file.path(
    results_dir,
    "V06B_Target_Cluster_Cell_Counts.csv"
  ),
  row.names = FALSE
)

cat("\nTarget cluster sizes:\n")
print(target_counts, row.names = FALSE)


# ============================================================
# 8. CYP/AA ANNOTATION EXCLUSION SET
# ============================================================

all_genes <- rownames(obj)

strict_cyp_genes <- grep(
  "^Cyp[0-9]",
  all_genes,
  value = TRUE
)

aa_exclusion_genes <- unique(
  c(
    strict_cyp_genes,
    "Ephx1",
    "Ephx2",
    "Pla2g4a",
    "Pla2g4b",
    "Pla2g4c",
    "Pla2g6",
    "Pla2g2a",
    "Pla2g5",
    "Ptgs1",
    "Ptgs2",
    "Ptges",
    "Ptges2",
    "Ptges3",
    "Hpgds",
    "Tbxas1",
    "Alox5",
    "Alox5ap",
    "Alox12",
    "Alox12b",
    "Alox15",
    "Alox15b",
    "Ptger1",
    "Ptger2",
    "Ptger3",
    "Ptger4",
    "Ptgdr",
    "Ptgdr2",
    "Tbxa2r",
    "Ltb4r1",
    "Ltb4r2",
    "Cysltr1",
    "Cysltr2"
  )
)

aa_exclusion_genes <- intersect(
  aa_exclusion_genes,
  all_genes
)

cat(
  "\nCYP/AA genes excluded from annotation evidence:",
  length(aa_exclusion_genes),
  "\n"
)


# ============================================================
# 9. RNA ASSAY
# ============================================================

DefaultAssay(obj) <- "RNA"

obj <- JoinLayers(
  obj,
  assay = "RNA"
)


# ============================================================
# 10. DEFINE COMPARISON NEIGHBORHOODS
# ============================================================
#
# We compare each ambiguous cluster against biologically
# plausible neighboring identities.
#
# This does NOT change cluster assignments.
#
# ============================================================

comparison_clusters <- list(
  
  Cluster8 = c(
    "0", "5", "6", "7", "8", "9", "11"
  ),
  
  Cluster10 = c(
    "1", "2", "3", "4", "10", "13", "14", "19", "20"
  ),
  
  Cluster16 = c(
    "1", "2", "11", "13", "16"
  )
)

writeLines(
  capture.output(
    print(comparison_clusters)
  ),
  file.path(
    results_dir,
    "V06B_Comparison_Cluster_Definitions.txt"
  )
)


# ============================================================
# 11. TARGETED CANONICAL MARKER PANEL
# ============================================================

marker_sets <- list(
  
  Fibroblast = c(
    "Col1a1",
    "Col1a2",
    "Col3a1",
    "Dcn",
    "Lum",
    "Pdgfra",
    "Pcolce",
    "Col6a1",
    "Col6a2",
    "Col6a3"
  ),
  
  Activated_Fibroblast = c(
    "Postn",
    "Cthrc1",
    "Comp",
    "Cilp",
    "Thbs2",
    "Tnc",
    "Fn1",
    "Ltbp2",
    "Ctgf",
    "Mfap4"
  ),
  
  Endothelial = c(
    "Pecam1",
    "Cdh5",
    "Kdr",
    "Emcn",
    "Esam",
    "Ramp2",
    "Eng",
    "Egfl7",
    "Klf2",
    "Klf4"
  ),
  
  Capillary_Endothelial = c(
    "Car4",
    "Gpihbp1",
    "Rgcc",
    "Aplnr",
    "Emcn",
    "Kdr",
    "Rbp7"
  ),
  
  Arterial_Endothelial = c(
    "Gja5",
    "Efnb2",
    "Sox17",
    "Bmx",
    "Hey1",
    "Nrp1"
  ),
  
  Venous_Endothelial = c(
    "Nr2f2",
    "Vcam1",
    "Ackr1"
  ),
  
  Pericyte = c(
    "Rgs5",
    "Cspg4",
    "Pdgfrb",
    "Kcnj8",
    "Abcc9",
    "Notch3",
    "Des",
    "Vtn"
  ),
  
  Smooth_Muscle = c(
    "Acta2",
    "Tagln",
    "Myh11",
    "Cnn1",
    "Actg2",
    "Lmod1",
    "Myocd",
    "Kcnmb1"
  ),
  
  Macrophage = c(
    "Adgre1",
    "C1qa",
    "C1qb",
    "C1qc",
    "Cd68",
    "Mrc1",
    "Fcgr1",
    "Ms4a7"
  ),
  
  Resident_Macrophage = c(
    "Timd4",
    "Folr2",
    "Vsig4",
    "Cd163",
    "Mrc1",
    "Lyve1"
  ),
  
  Monocyte = c(
    "Lyz2",
    "Ly6c2",
    "Ccr2",
    "Plac8",
    "Ms4a8a",
    "Ctss"
  ),
  
  Dendritic = c(
    "Flt3",
    "Xcr1",
    "Clec9a",
    "Cd209a",
    "Cd74",
    "H2-Aa",
    "H2-Ab1"
  ),
  
  Neutrophil = c(
    "S100a8",
    "S100a9",
    "Ly6g",
    "Csf3r",
    "Cxcr2",
    "Ngp",
    "Camp",
    "Ltf"
  ),
  
  Cycling = c(
    "Mki67",
    "Top2a",
    "Cenpf",
    "Ube2c",
    "Birc5",
    "Ccna2",
    "Ccnb1",
    "Plk1"
  ),
  
  Mesothelial_Epicardial = c(
    "Wt1",
    "Tbx18",
    "Upk3b",
    "Msln",
    "Krt8",
    "Krt18",
    "Krt19",
    "Aldh1a2",
    "Lrrn4",
    "Muc16"
  ),
  
  Schwann_Neural = c(
    "Sox10",
    "S100b",
    "Plp1",
    "Mpz",
    "Gfra3",
    "Gjc3"
  ),
  
  Cardiomyocyte = c(
    "Tnnt2",
    "Tnni3",
    "Myh6",
    "Myh7",
    "Actc1",
    "Myl2",
    "Myl3"
  ),
  
  Erythroid = c(
    "Hba-a1",
    "Hba-a2",
    "Hbb-bs",
    "Hbb-bt",
    "Alas2",
    "Gypa"
  )
)


# ============================================================
# 12. FILTER MARKER PANEL
# ============================================================

marker_sets <- lapply(
  marker_sets,
  function(x) {
    intersect(
      setdiff(
        x,
        aa_exclusion_genes
      ),
      rownames(obj)
    )
  }
)

marker_sets <- marker_sets[
  lengths(marker_sets) > 0
]

marker_features <- unique(
  unlist(
    marker_sets,
    use.names = FALSE
  )
)

cat(
  "\nTargeted canonical marker genes available:",
  length(marker_features),
  "\n"
)


# ============================================================
# 13. TARGET CLUSTER SAMPLE/CONDITION AUDIT
# ============================================================

target_metadata <- obj@meta.data %>%
  mutate(
    Cluster = as.character(
      V06_cluster_0.6
    )
  ) %>%
  filter(
    Cluster %in% target_clusters
  )

target_sample_counts <- as.data.frame(
  table(
    Cluster = target_metadata$Cluster,
    Sample = target_metadata$sample_id
  )
)

target_sample_counts <- target_sample_counts[
  target_sample_counts$Freq > 0,
]

target_sample_counts <- target_sample_counts %>%
  group_by(Cluster) %>%
  mutate(
    Fraction_of_Cluster =
      Freq / sum(Freq)
  ) %>%
  ungroup()

write.csv(
  target_sample_counts,
  file.path(
    results_dir,
    "V06B_Target_Cluster_Sample_Composition.csv"
  ),
  row.names = FALSE
)


target_condition_counts <- as.data.frame(
  table(
    Cluster = target_metadata$Cluster,
    Condition = target_metadata$condition
  )
)

target_condition_counts <- target_condition_counts[
  target_condition_counts$Freq > 0,
]

target_condition_counts <- target_condition_counts %>%
  group_by(Cluster) %>%
  mutate(
    Fraction_of_Cluster =
      Freq / sum(Freq)
  ) %>%
  ungroup()

write.csv(
  target_condition_counts,
  file.path(
    results_dir,
    "V06B_Target_Cluster_Condition_Composition.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. RESOLUTION 1.2 SUBSTRUCTURE
# ============================================================

fine_crosswalk <- as.data.frame(
  table(
    Cluster_0.6 = obj$V06_cluster_0.6,
    Cluster_1.2 = obj$V06_fine_cluster_1.2
  )
)

fine_crosswalk$Cluster_0.6 <- as.character(
  fine_crosswalk$Cluster_0.6
)

fine_crosswalk$Cluster_1.2 <- as.character(
  fine_crosswalk$Cluster_1.2
)

fine_crosswalk <- fine_crosswalk[
  fine_crosswalk$Cluster_0.6 %in% target_clusters &
    fine_crosswalk$Freq > 0,
]

fine_crosswalk <- fine_crosswalk %>%
  group_by(
    Cluster_0.6
  ) %>%
  mutate(
    Fraction_of_0.6_Cluster =
      Freq / sum(Freq)
  ) %>%
  ungroup()

write.csv(
  fine_crosswalk,
  file.path(
    results_dir,
    "V06B_Target_Clusters_Resolution_1.2_Substructure.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. RETRIEVE V06 UNBIASED MARKERS
# ============================================================
#
# IMPORTANT:
# We do NOT rerun FindAllMarkers().
# We use the marker evidence already generated in V06.
#
# ============================================================

v06_marker_file <- file.path(
  project_root,
  "RESULTS",
  "V06_ANNOTATION_EVIDENCE_AUDIT",
  "V06_Annotation_Markers_CYP_AA_Excluded.csv"
)

if (!file.exists(v06_marker_file)) {
  stop(
    "V06 annotation marker table not found."
  )
}

v06_markers <- read.csv(
  v06_marker_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(c("cluster", "gene") %in% colnames(v06_markers))) {
  stop("Unexpected V06 marker table format.")
}

v06_markers$cluster <- as.character(
  v06_markers$cluster
)

target_markers <- v06_markers[
  v06_markers$cluster %in% target_clusters,
]

target_markers <- target_markers[
  !target_markers$gene %in% aa_exclusion_genes,
]

target_markers <- target_markers %>%
  group_by(
    cluster
  ) %>%
  arrange(
    desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  mutate(
    Marker_Rank = row_number()
  ) %>%
  ungroup()

top50_target_markers <- target_markers %>%
  filter(
    Marker_Rank <= 50
  )

write.csv(
  top50_target_markers,
  file.path(
    results_dir,
    "V06B_Top50_Markers_Clusters_8_10_16.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. TARGET CLUSTER MARKER DETECTION / EXPRESSION
# ============================================================

rna_data <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "data"
)

cluster_vector <- as.character(
  obj$V06_cluster_0.6
)

canonical_summary_list <- list()

counter <- 1

for (cluster_id in target_clusters) {
  
  target_cells <- colnames(obj)[
    cluster_vector == cluster_id
  ]
  
  for (lineage in names(marker_sets)) {
    
    genes <- marker_sets[[lineage]]
    
    mat <- rna_data[
      genes,
      target_cells,
      drop = FALSE
    ]
    
    pct <- Matrix::rowMeans(
      mat > 0
    ) * 100
    
    avg <- Matrix::rowMeans(
      mat
    )
    
    canonical_summary_list[[counter]] <- data.frame(
      Cluster = cluster_id,
      Program = lineage,
      Number_of_Markers = length(genes),
      Mean_Percent_Detected = mean(pct),
      Median_Percent_Detected = median(pct),
      Mean_Normalized_Expression = mean(avg),
      stringsAsFactors = FALSE
    )
    
    counter <- counter + 1
  }
}

canonical_summary <- bind_rows(
  canonical_summary_list
)

canonical_summary <- canonical_summary %>%
  group_by(
    Cluster
  ) %>%
  arrange(
    desc(Mean_Percent_Detected),
    desc(Mean_Normalized_Expression),
    .by_group = TRUE
  ) %>%
  mutate(
    Rank = row_number()
  ) %>%
  ungroup()

write.csv(
  canonical_summary,
  file.path(
    results_dir,
    "V06B_Target_Cluster_Canonical_Evidence.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. INDIVIDUAL MARKER SUMMARY
# ============================================================

individual_marker_summary_list <- list()

counter <- 1

for (cluster_id in target_clusters) {
  
  cells_cluster <- colnames(obj)[
    cluster_vector == cluster_id
  ]
  
  mat <- rna_data[
    marker_features,
    cells_cluster,
    drop = FALSE
  ]
  
  pct <- Matrix::rowMeans(
    mat > 0
  ) * 100
  
  avg <- Matrix::rowMeans(
    mat
  )
  
  individual_marker_summary_list[[counter]] <- data.frame(
    Cluster = cluster_id,
    Gene = marker_features,
    Percent_Detected = as.numeric(pct),
    Mean_Normalized_Expression = as.numeric(avg),
    stringsAsFactors = FALSE
  )
  
  counter <- counter + 1
}

individual_marker_summary <- bind_rows(
  individual_marker_summary_list
)

write.csv(
  individual_marker_summary,
  file.path(
    results_dir,
    "V06B_Individual_Canonical_Marker_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. TARGET VS BIOLOGICAL NEIGHBORS
# ============================================================

neighbor_summary_list <- list()

counter <- 1

for (target_name in names(comparison_clusters)) {
  
  cluster_set <- comparison_clusters[[target_name]]
  
  for (cluster_id in cluster_set) {
    
    cells_cluster <- colnames(obj)[
      cluster_vector == cluster_id
    ]
    
    for (lineage in names(marker_sets)) {
      
      genes <- marker_sets[[lineage]]
      
      mat <- rna_data[
        genes,
        cells_cluster,
        drop = FALSE
      ]
      
      pct <- Matrix::rowMeans(
        mat > 0
      ) * 100
      
      avg <- Matrix::rowMeans(
        mat
      )
      
      neighbor_summary_list[[counter]] <- data.frame(
        Comparison_Set = target_name,
        Cluster = cluster_id,
        Program = lineage,
        Mean_Percent_Detected = mean(pct),
        Mean_Normalized_Expression = mean(avg),
        stringsAsFactors = FALSE
      )
      
      counter <- counter + 1
    }
  }
}

neighbor_summary <- bind_rows(
  neighbor_summary_list
)

write.csv(
  neighbor_summary,
  file.path(
    results_dir,
    "V06B_Target_vs_Neighbor_Lineage_Evidence.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. UMAP — TARGET CLUSTERS
# ============================================================

obj$V06B_target_status <- ifelse(
  as.character(obj$V06_cluster_0.6) %in% target_clusters,
  paste0(
    "Cluster ",
    as.character(obj$V06_cluster_0.6)
  ),
  "Other"
)

target_levels <- c(
  "Other",
  "Cluster 8",
  "Cluster 10",
  "Cluster 16"
)

obj$V06B_target_status <- factor(
  obj$V06B_target_status,
  levels = target_levels
)

p_target_umap <- DimPlot(
  obj,
  reduction = "umap.unintegrated",
  group.by = "V06B_target_status",
  raster = TRUE,
  shuffle = FALSE
) +
  ggtitle(
    "Targeted Identity Audit — Clusters 8, 10, and 16"
  ) +