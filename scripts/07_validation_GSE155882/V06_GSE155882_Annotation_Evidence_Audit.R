# ============================================================
# V06 — GSE155882 ANNOTATION EVIDENCE AUDIT
# CORRECTED VERSION
# ============================================================
#
# PURPOSE
# -------
# Build an independent, publication-grade annotation evidence
# framework for the GSE155882 validation atlas.
#
# PRIMARY ANNOTATION RESOLUTION:
#   V05 resolution 0.6
#
# SUPPORTING FINE-STATE RESOLUTION:
#   V05 resolution 1.2
#
# IMPORTANT:
#   - NO CYP/AA genes used to assign cell identity
#   - NO cells removed
#   - NO reclustering
#   - NO batch integration
#   - NO Sham-vs-TAC DE
#   - NO final biological annotation frozen here
#
# V06 produces:
#   1. 22-cluster resolution-0.6 atlas
#   2. unbiased markers for resolution-0.6 clusters
#   3. canonical lineage-marker evidence
#   4. cluster/sample/condition composition
#   5. relationship between resolution 0.6 and 1.2
#   6. annotation-support tables
#   7. high-resolution PNG figures
#
# CORRECTION:
#   Frozen V04/V05 UMAP reduction is named:
#       "umap.unintegrated"
#   NOT:
#       "umap_unintegrated"
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
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(Matrix)
})


# ============================================================
# 3. PATHS
# ============================================================

project_root <- "D:/Master/ScRNA seq/TAC/VALIDATION/GSE155882"

input_file <- file.path(
  project_root,
  "CLEAN DATA",
  "V05_GSE155882_UNBIASED_CLUSTERED_ATLAS.rds"
)

output_file <- file.path(
  project_root,
  "CLEAN DATA",
  "V06_GSE155882_ANNOTATION_EVIDENCE_ATLAS.rds"
)

results_dir <- file.path(
  project_root,
  "RESULTS",
  "V06_ANNOTATION_EVIDENCE_AUDIT"
)

figures_dir <- file.path(
  project_root,
  "FIGURES",
  "V06_ANNOTATION_EVIDENCE_AUDIT"
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
# 4. LOAD V05 OBJECT
# ============================================================

stopifnot(
  file.exists(input_file)
)

obj <- readRDS(
  input_file
)

cat("\nV05 object loaded successfully.\n")
cat("Cells:", ncol(obj), "\n")
cat("Genes:", nrow(obj), "\n")

cat("\nAvailable reductions:\n")
print(
  Reductions(obj)
)


# ============================================================
# 5. BASIC INPUT AUDIT
# ============================================================

input_n <- ncol(obj)

input_cell_names <- colnames(obj)

required_metadata <- c(
  "sample_id",
  "condition",
  "replicate",
  "scDblFinder_score",
  "scDblFinder_class"
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

if (!"RNA" %in% Assays(obj)) {
  stop("RNA assay is missing.")
}

if (!"pca" %in% Reductions(obj)) {
  stop("Frozen V04 PCA is missing.")
}

if (!"umap.unintegrated" %in% Reductions(obj)) {
  stop("Frozen unintegrated UMAP is missing.")
}

cat("\nFrozen reductions confirmed:\n")
cat(
  "PCA:",
  nrow(
    Embeddings(
      obj,
      reduction = "pca"
    )
  ),
  "cells x",
  ncol(
    Embeddings(
      obj,
      reduction = "pca"
    )
  ),
  "dimensions\n"
)

cat(
  "UMAP:",
  nrow(
    Embeddings(
      obj,
      reduction = "umap.unintegrated"
    )
  ),
  "cells x",
  ncol(
    Embeddings(
      obj,
      reduction = "umap.unintegrated"
    )
  ),
  "dimensions\n"
)


# ============================================================
# 6. IDENTIFY RESOLUTION METADATA COLUMNS
# ============================================================

resolution_columns <- grep(
  "snn_res",
  colnames(obj@meta.data),
  value = TRUE
)

cat("\nAvailable clustering columns:\n")
print(
  resolution_columns
)


find_resolution_column <- function(
    metadata_columns,
    resolution
) {
  
  candidates <- metadata_columns[
    grepl(
      paste0(
        "res\\.",
        resolution,
        "$"
      ),
      metadata_columns
    )
  ]
  
  if (length(candidates) == 0) {
    
    candidates <- metadata_columns[
      grepl(
        paste0(
          "res\\.",
          resolution
        ),
        metadata_columns
      )
    ]
    
  }
  
  if (length(candidates) != 1) {
    
    stop(
      paste0(
        "Could not uniquely identify resolution ",
        resolution,
        " column. Candidates: ",
        paste(
          candidates,
          collapse = ", "
        )
      )
    )
    
  }
  
  return(candidates)
}


resolution_06_col <- find_resolution_column(
  resolution_columns,
  "0.6"
)

resolution_12_col <- find_resolution_column(
  resolution_columns,
  "1.2"
)

cat(
  "\nResolution 0.6 column:",
  resolution_06_col,
  "\n"
)

cat(
  "Resolution 1.2 column:",
  resolution_12_col,
  "\n"
)


# ============================================================
# 7. FREEZE V06 CLUSTER IDENTITIES
# ============================================================

obj$V06_cluster_0.6 <- factor(
  obj@meta.data[[resolution_06_col]]
)

obj$V06_fine_cluster_1.2 <- factor(
  obj@meta.data[[resolution_12_col]]
)

Idents(obj) <- "V06_cluster_0.6"

number_clusters_06 <- length(
  levels(
    Idents(obj)
  )
)

number_clusters_12 <- length(
  unique(
    obj$V06_fine_cluster_1.2
  )
)

cat(
  "\nResolution 0.6 clusters:",
  number_clusters_06,
  "\n"
)

cat(
  "Resolution 1.2 clusters:",
  number_clusters_12,
  "\n"
)

if (number_clusters_06 != 22) {
  
  warning(
    paste(
      "Expected 22 clusters at resolution 0.6 but found",
      number_clusters_06
    )
  )
  
}

if (number_clusters_12 != 33) {
  
  warning(
    paste(
      "Expected 33 clusters at resolution 1.2 but found",
      number_clusters_12
    )
  )
  
}


# ============================================================
# 8. EXCLUDE CYP/AA GENES FROM ANNOTATION EVIDENCE
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
    
    # Epoxide hydrolases
    "Ephx1",
    "Ephx2",
    
    # PLA2 / AA liberation
    "Pla2g4a",
    "Pla2g4b",
    "Pla2g4c",
    "Pla2g6",
    "Pla2g2a",
    "Pla2g5",
    
    # COX arm
    "Ptgs1",
    "Ptgs2",
    "Ptges",
    "Ptges2",
    "Ptges3",
    "Hpgds",
    "Tbxas1",
    
    # LOX arm
    "Alox5",
    "Alox5ap",
    "Alox12",
    "Alox12b",
    "Alox15",
    "Alox15b",
    
    # Selected eicosanoid receptors
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
  "\nStrict CYP genes excluded from annotation evidence:",
  length(strict_cyp_genes),
  "\n"
)

cat(
  "Total CYP/AA genes excluded from annotation evidence:",
  length(aa_exclusion_genes),
  "\n"
)

write.csv(
  data.frame(
    Gene = aa_exclusion_genes,
    stringsAsFactors = FALSE
  ),
  file.path(
    results_dir,
    "V06_CYP_AA_Genes_Excluded_From_Annotation.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9. CLUSTER SIZE AUDIT
# ============================================================

cluster_sizes <- as.data.frame(
  table(
    Cluster = obj$V06_cluster_0.6
  )
)

colnames(cluster_sizes) <- c(
  "Cluster",
  "Cells"
)

cluster_sizes$Cluster <- as.character(
  cluster_sizes$Cluster
)

cluster_sizes$Percent_of_Atlas <- (
  cluster_sizes$Cells /
    ncol(obj)
) * 100

cluster_sizes <- cluster_sizes[
  order(
    as.numeric(
      cluster_sizes$Cluster
    )
  ),
]

rownames(cluster_sizes) <- NULL

write.csv(
  cluster_sizes,
  file.path(
    results_dir,
    "V06_Cluster_Sizes_Resolution_0.6.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. CLUSTER × SAMPLE COUNTS
# ============================================================

cluster_sample_counts <- as.data.frame.matrix(
  table(
    obj$V06_cluster_0.6,
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
      colnames(cluster_sample_counts),
      "Cluster"
    )
  )
]

rownames(cluster_sample_counts) <- NULL

write.csv(
  cluster_sample_counts,
  file.path(
    results_dir,
    "V06_Cluster_by_Sample_Counts.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. CLUSTER × CONDITION COUNTS
# ============================================================

cluster_condition_counts <- as.data.frame.matrix(
  table(
    obj$V06_cluster_0.6,
    obj$condition
  )
)

cluster_condition_counts$Cluster <- rownames(
  cluster_condition_counts
)

cluster_condition_counts <- cluster_condition_counts[
  ,
  c(
    "Cluster",
    setdiff(
      colnames(cluster_condition_counts),
      "Cluster"
    )
  )
]

rownames(cluster_condition_counts) <- NULL

write.csv(
  cluster_condition_counts,
  file.path(
    results_dir,
    "V06_Cluster_by_Condition_Counts.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. WITHIN-CLUSTER SAMPLE FRACTIONS
# ============================================================

sample_fraction_table <- as.data.frame(
  prop.table(
    table(
      obj$V06_cluster_0.6,
      obj$sample_id
    ),
    margin = 1
  )
)

colnames(sample_fraction_table) <- c(
  "Cluster",
  "Sample",
  "Fraction"
)

sample_fraction_table$Cluster <- as.character(
  sample_fraction_table$Cluster
)

sample_fraction_table$Sample <- as.character(
  sample_fraction_table$Sample
)

sample_fraction_table$Percent <- (
  sample_fraction_table$Fraction *
    100
)

write.csv(
  sample_fraction_table,
  file.path(
    results_dir,
    "V06_Within_Cluster_Sample_Fractions.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. RESOLUTION 0.6 ↔ 1.2 CROSSWALK
# ============================================================

cluster_crosswalk <- as.data.frame(
  table(
    Cluster_0.6 = obj$V06_cluster_0.6,
    Cluster_1.2 = obj$V06_fine_cluster_1.2
  )
)

cluster_crosswalk$Cluster_0.6 <- as.character(
  cluster_crosswalk$Cluster_0.6
)

cluster_crosswalk$Cluster_1.2 <- as.character(
  cluster_crosswalk$Cluster_1.2
)

cluster_crosswalk <- cluster_crosswalk[
  cluster_crosswalk$Freq > 0,
]

cluster_crosswalk <- cluster_crosswalk %>%
  group_by(
    Cluster_0.6
  ) %>%
  mutate(
    Fraction_of_0.6_Cluster =
      Freq /
      sum(Freq)
  ) %>%
  ungroup()

write.csv(
  cluster_crosswalk,
  file.path(
    results_dir,
    "V06_Resolution_0.6_to_1.2_Crosswalk.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. SET RNA ASSAY AND JOIN LAYERS
# ============================================================

DefaultAssay(obj) <- "RNA"

obj <- JoinLayers(
  obj,
  assay = "RNA"
)

Idents(obj) <- "V06_cluster_0.6"


# ============================================================
# 15. UNBIASED MARKERS — RESOLUTION 0.6
# ============================================================

cat(
  "\nRunning unbiased marker analysis for resolution 0.6...\n"
)

cluster_markers_all <- FindAllMarkers(
  object = obj,
  assay = "RNA",
  slot = "data",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  return.thresh = 0.05,
  verbose = TRUE
)

if (nrow(cluster_markers_all) == 0) {
  stop("FindAllMarkers returned zero marker rows.")
}

if (!"gene" %in% colnames(cluster_markers_all)) {
  
  cluster_markers_all$gene <- rownames(
    cluster_markers_all
  )
  
}

rownames(cluster_markers_all) <- NULL

write.csv(
  cluster_markers_all,
  file.path(
    results_dir,
    "V06_All_Unbiased_Markers_Resolution_0.6.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. REMOVE CYP/AA GENES FROM ANNOTATION MARKER TABLE
# ============================================================

annotation_markers <- cluster_markers_all[
  !cluster_markers_all$gene %in%
    aa_exclusion_genes,
]

if (nrow(annotation_markers) == 0) {
  stop("No annotation markers remained after CYP/AA exclusion.")
}

write.csv(
  annotation_markers,
  file.path(
    results_dir,
    "V06_Annotation_Markers_CYP_AA_Excluded.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. TOP 20 MARKERS PER CLUSTER
# ============================================================

top20_markers <- annotation_markers %>%
  group_by(
    cluster
  ) %>%
  arrange(
    desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  slice_head(
    n = 20
  ) %>%
  ungroup()

write.csv(
  top20_markers,
  file.path(
    results_dir,
    "V06_Top20_Annotation_Markers_Per_Cluster.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. CANONICAL ANNOTATION MARKER SETS
# ============================================================

canonical_marker_sets <- list(
  
  Fibroblast = c(
    "Col1a1",
    "Col1a2",
    "Dcn",
    "Lum",
    "Pdgfra",
    "Col3a1",
    "Col6a1",
    "Col6a2",
    "Col6a3",
    "Pcolce"
  ),
  
  Activated_Fibroblast = c(
    "Postn",
    "Cthrc1",
    "Comp",
    "Thbs2",
    "Tnc",
    "Fn1",
    "Ltbp2",
    "Crlf1"
  ),
  
  Endothelial = c(
    "Pecam1",
    "Cdh5",
    "Kdr",
    "Emcn",
    "Klf2",
    "Klf4",
    "Esam",
    "Ramp2",
    "Eng",
    "Egfl7"
  ),
  
  Capillary_Endothelial = c(
    "Car4",
    "Gpihbp1",
    "Rgcc",
    "Kdr",
    "Emcn"
  ),
  
  Arterial_Endothelial = c(
    "Gja5",
    "Efnb2",
    "Sox17",
    "Bmx",
    "Hey1"
  ),
  
  Venous_Endothelial = c(
    "Nr2f2",
    "Vcam1",
    "Ackr1"
  ),
  
  Lymphatic_Endothelial = c(
    "Ccl21a",
    "Lyve1",
    "Prox1",
    "Flt4",
    "Pdpn",
    "Mmrn1"
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
    "Lyve1",
    "Vsig4",
    "Cd163",
    "Mrc1"
  ),
  
  Monocyte = c(
    "Ly6c2",
    "Ccr2",
    "Plac8",
    "Ms4a8a",
    "Lyz2",
    "Ctss"
  ),
  
  Dendritic = c(
    "Flt3",
    "Xcr1",
    "Clec9a",
    "Cd209a",
    "H2-Ab1",
    "Cd74"
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
  
  T_cell = c(
    "Cd3d",
    "Cd3e",
    "Cd3g",
    "Trac",
    "Cd247",
    "Lck"
  ),
  
  NK_cell = c(
    "Ncr1",
    "Klrb1c",
    "Prf1",
    "Nkg7",
    "Gzma",
    "Eomes"
  ),
  
  B_cell = c(
    "Cd79a",
    "Cd79b",
    "Ms4a1",
    "Cd19",
    "Cd37",
    "H2-DMb1"
  ),
  
  Plasma_cell = c(
    "Jchain",
    "Mzb1",
    "Sdc1",
    "Xbp1",
    "Derl3"
  ),
  
  Pericyte = c(
    "Rgs5",
    "Cspg4",
    "Kcnj8",
    "Abcc9",
    "Pdgfrb",
    "Des"
  ),
  
  Smooth_Muscle = c(
    "Acta2",
    "Tagln",
    "Myh11",
    "Cnn1",
    "Actg2",
    "Lmod1",
    "Myocd"
  ),
  
  Schwann_Neural = c(
    "S100b",
    "Plp1",
    "Sox10",
    "Mpz",
    "Sox2"
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
  ),
  
  Platelet = c(
    "Pf4",
    "Ppbp",
    "Gp9",
    "Itga2b",
    "Tubb1",
    "Clec1b"
  ),
  
  Cycling = c(
    "Mki67",
    "Top2a",
    "Cenpf",
    "Ube2c",