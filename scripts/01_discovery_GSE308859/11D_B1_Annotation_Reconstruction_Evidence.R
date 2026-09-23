# ============================================================
# STEP 11D-B1
# ANNOTATION RECONSTRUCTION — EVIDENCE GENERATION
#
# PURPOSE:
# Independently reconstruct the biological identity of all
# frozen clusters BEFORE changing any annotation.
#
# IMPORTANT:
#   - NO reclustering
#   - NO cell filtering
#   - NO corrected labels assigned
#   - NO CYP genes used for annotation
#   - NO original metadata overwritten
#   - NO corrected atlas saved
#
# OUTPUT:
# Publication-grade evidence tables and diagnostic figures
# for manual/biological review of clusters 0–28.
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)


# ============================================================
# 2. PATHS
# ============================================================

project_dir <- "D:/Master/ScRNA seq/TAC"

input_file <- file.path(
  project_dir,
  "CLEAN DATA",
  "07C_FINAL_FROZEN_ANNOTATED_ATLAS.rds"
)

audit_marker_file <- file.path(
  project_dir,
  "RESULTS",
  "11_SPATIAL",
  "11D_ANNOTATION_AUDIT",
  "11D_A_All_Cluster_Audit_Markers.csv"
)

results_dir <- file.path(
  project_dir,
  "RESULTS",
  "11_SPATIAL",
  "11D_B1_ANNOTATION_RECONSTRUCTION"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "11_SPATIAL",
  "11D_B1_ANNOTATION_RECONSTRUCTION"
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
cat("STEP 11D-B1\n")
cat("ANNOTATION RECONSTRUCTION — EVIDENCE ONLY\n")
cat("============================================\n")

scrna <- readRDS(input_file)

DefaultAssay(scrna) <- "RNA"

cat("\nFrozen atlas loaded.\n")
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
      paste(
        missing_metadata,
        collapse = ", "
      )
    )
  )
}

cat("Required metadata present.\n")


# ============================================================
# 5. VERIFY CLUSTERS
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

cat(
  "Number of frozen clusters:",
  length(cluster_ids),
  "\n"
)

cat(
  "Cluster IDs:",
  paste(
    cluster_ids,
    collapse = ", "
  ),
  "\n"
)


# ============================================================
# 6. RECORD ORIGINAL STATE
# ============================================================

original_cells <- colnames(scrna)

original_clusters <- as.character(
  scrna$seurat_clusters
)

original_broad <- as.character(
  scrna$broad_cell_type
)

original_refined <- as.character(
  scrna$refined_cell_type
)


# ============================================================
# 7. LOAD INDEPENDENT CLUSTER MARKERS FROM 11D-A
# ============================================================

if (!file.exists(audit_marker_file)) {
  
  stop(
    paste(
      "Cannot find 11D-A marker file:",
      audit_marker_file
    )
  )
}

cluster_markers <- read.csv(
  audit_marker_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "\nLoaded",
  nrow(cluster_markers),
  "cluster-marker rows from 11D-A.\n"
)

required_marker_columns <- c(
  "gene",
  "cluster",
  "pct.1",
  "pct.2",
  "p_val_adj"
)

missing_marker_columns <- setdiff(
  required_marker_columns,
  colnames(cluster_markers)
)

if (length(missing_marker_columns) > 0) {
  
  stop(
    paste(
      "Missing marker columns:",
      paste(
        missing_marker_columns,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 8. IDENTIFY FOLD-CHANGE COLUMN
# ============================================================

possible_fc_columns <- c(
  "avg_log2FC",
  "avg_logFC"
)

fc_column <- intersect(
  possible_fc_columns,
  colnames(cluster_markers)
)

if (length(fc_column) == 0) {
  
  stop(
    paste(
      "Fold-change column not found. Available columns:",
      paste(
        colnames(cluster_markers),
        collapse = ", "
      )
    )
  )
}

fc_column <- fc_column[1]

cat(
  "Using fold-change column:",
  fc_column,
  "\n"
)


# ============================================================
# 9. DEFINE BROAD CANONICAL IDENTITY PANELS
#
# These panels intentionally avoid CYP genes.
#
# A biological identity will NOT be assigned simply because
# one panel has the highest score.
# ============================================================

marker_panels <- list(
  
  Cardiomyocyte = c(
    "Tnnt2",
    "Tnni3",
    "Myh6",
    "Myh7",
    "Actc1",
    "Myl2",
    "Myl3",
    "Pln",
    "Ryr2",
    "Atp2a2",
    "Csrp3",
    "Mybpc3"
  ),
  
  Fibroblast = c(
    "Col1a1",
    "Col1a2",
    "Col3a1",
    "Dcn",
    "Lum",
    "Pdgfra",
    "Tcf21",
    "Dpt",
    "Col5a1",
    "Pcolce2",
    "Fbln1",
    "Clec3b"
  ),
  
  Endothelial = c(
    "Pecam1",
    "Cdh5",
    "Kdr",
    "Emcn",
    "Esam",
    "Egfl7",
    "Ptprb",
    "Tek",
    "Vwf",
    "Klf2",
    "Klf4"
  ),
  
  Macrophage = c(
    "Adgre1",
    "Csf1r",
    "Cd68",
    "C1qa",
    "C1qb",
    "C1qc",
    "Fcgr1",
    "Ms4a7",
    "Mertk",
    "Aif1",
    "Tyrobp"
  ),
  
  Monocyte = c(
    "Ccr2",
    "Ly6c2",
    "Lyz2",
    "Ctss",
    "Fcgr3",
    "Plac8",
    "Tyrobp",
    "Aif1"
  ),
  
  Neutrophil = c(
    "S100a8",
    "S100a9",
    "Csf3r",
    "Cxcr2",
    "Retnlg",
    "Mmp8",
    "Mmp9",
    "Lcn2",
    "Ly6g",
    "Camp",
    "Ngp"
  ),
  
  Dendritic = c(
    "Flt3",
    "Itgax",
    "Clec9a",
    "Batf3",
    "Cd209a",
    "Ciita",
    "Wdfy4",
    "H2-DMa",
    "H2-DMb1"
  ),
  
  T_cell = c(
    "Cd3d",
    "Cd3e",
    "Cd3g",
    "Trac",
    "Lck",
    "Lat",
    "Cd247",
    "Il7r",
    "Itk"
  ),
  
  NK_cell = c(
    "Nkg7",
    "Ncr1",
    "Prf1",
    "Gzmb",
    "Ccl5",
    "Klrd1",
    "Klrk1",
    "Eomes"
  ),
  
  B_cell = c(
    "Cd79a",
    "Cd79b",
    "Ms4a1",
    "Cd19",
    "Cd22",
    "Pax5",
    "Fcmr",
    "Ighm",
    "Ighd",
    "Tnfrsf13c"
  ),
  
  Pericyte = c(
    "Rgs5",
    "Pdgfrb",
    "Cspg4",
    "Kcnj8",
    "Abcc9",
    "Rbp1",
    "Notch3"
  ),
  
  Smooth_muscle = c(
    "Acta2",
    "Tagln",
    "Myh11",
    "Cnn1",
    "Myocd",
    "Mylk",
    "Kcnmb1"
  ),
  
  Platelet_Megakaryocyte = c(
    "Pf4",
    "Ppbp",
    "Itga2b",
    "Gp9",
    "Gp5",
    "Gp6",
    "Treml1",
    "Clec1b",
    "Gp1ba"
  ),
  
  Erythroid = c(
    "Gypa",
    "Slc4a1",
    "Alas2",
    "Bpgm",
    "Fech",
    "Tspo2",
    "Hba-a1",
    "Hba-a2"
  ),
  
  Cycling = c(
    "Mki67",
    "Top2a",
    "Cdk1",
    "Ccnb1",
    "Ccnb2",
    "Birc5",
    "Pclaf",
    "Tpx2",
    "Pbk"
  )
)


# ============================================================
# 10. DEFINE SUBTYPE / STATE PANELS
#
# These are supportive evidence only.
# They are NOT automatically converted into labels.
# ============================================================

subtype_panels <- list(
  
  Resident_Macrophage = c(
    "Timd4",
    "Lyve1",
    "Folr2",
    "Cd163",
    "Mrc1",
    "Vsig4",
    "F13a1"
  ),
  
  Inflammatory_Myeloid = c(
    "Il1b",
    "S100a8",
    "S100a9",
    "Trem1",
    "Ccr2",
    "Ptgs2"
  ),
  
  Capillary_Endothelial = c(
    "Car4",
    "Aplnr",
    "Gpihbp1",
    "Aqp1",
    "Rgcc",
    "Kdr"
  ),
  
  Arterial_Endothelial = c(
    "Gja5",
    "Efnb2",
    "Sox17",
    "Bmx",
    "Fbln5"
  ),
  
  Venous_Endothelial = c(
    "Nr2f2",
    "Ackr1"
  ),
  
  Activated_Fibroblast = c(
    "Postn",
    "Cthrc1",
    "Ccn2",
    "Col8a1",
    "Lox",
    "Loxl2"
  ),
  
  ECM_Fibroblast = c(
    "Col1a1",
    "Col1a2",
    "Col3a1",
    "Dcn",
    "Lum",
    "Fmod"
  ),
  
  GammaDelta_T = c(
    "Trdc",
    "Trgc1",
    "Trgc2"
  ),
  
  Immature_Lymphoid = c(
    "Rag1",
    "Rag2",
    "Dntt"
  )
)


# ============================================================
# 11. KEEP ONLY GENES PRESENT IN DATASET
# ============================================================

filter_panel <- function(panel_list, genes_present) {
  
  output <- lapply(
    panel_list,
    function(x) {
      intersect(
        x,
        genes_present
      )
    }
  )
  
  output[
    lengths(output) > 0
  ]
}

marker_panels_present <- filter_panel(
  marker_panels,
  rownames(scrna)
)

subtype_panels_present <- filter_panel(
  subtype_panels,
  rownames(scrna)
)


# ============================================================
# 12. PANEL AVAILABILITY TABLE
# ============================================================

panel_availability <- bind_rows(
  
  lapply(
    names(marker_panels),
    function(x) {
      
      requested <- marker_panels[[x]]
      
      present <- intersect(
        requested,
        rownames(scrna)
      )
      
      data.frame(
        Panel_Type = "Broad identity",
        Panel = x,
        Requested_Genes = length(requested),
        Present_Genes = length(present),
        Genes_Present = paste(
          present,
          collapse = "; "
        ),
        stringsAsFactors = FALSE
      )
    }
  ),
  
  lapply(
    names(subtype_panels),
    function(x) {
      
      requested <- subtype_panels[[x]]
      
      present <- intersect(
        requested,
        rownames(scrna)
      )
      
      data.frame(
        Panel_Type = "Subtype/state",
        Panel = x,
        Requested_Genes = length(requested),
        Present_Genes = length(present),
        Genes_Present = paste(
          present,
          collapse = "; "
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  panel_availability,
  file.path(
    results_dir,
    "11D_B1_Marker_Panel_Availability.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. TOP 30 MARKERS PER CLUSTER
# ============================================================

cluster_markers <- cluster_markers %>%
  mutate(
    cluster = as.character(cluster),
    Detection_Difference = pct.1 - pct.2
  )

top30_markers <- cluster_markers %>%
  arrange(
    cluster,
    p_val_adj,
    desc(.data[[fc_column]]),
    desc(Detection_Difference)
  ) %>%
  group_by(cluster) %>%
  slice_head(n = 30) %>%
  ungroup()

write.csv(
  top30_markers,
  file.path(
    results_dir,
    "11D_B1_Top30_Markers_Per_Cluster.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. COMPACT TOP-MARKER TABLE
# ============================================================

top_marker_summary <- top30_markers %>%
  group_by(cluster) %>%
  summarise(
    Top30_Markers = paste(
      gene,
      collapse = ", "
    ),
    .groups = "drop"
  )

write.csv(
  top_marker_summary,
  file.path(
    results_dir,
    "11D_B1_Top30_Marker_Summary.csv"
  ),
  row.names = FALSE
)

cat("\n============================================\n")
cat("TOP 30 MARKERS PER CLUSTER\n")
cat("============================================\n")

print(
  as.data.frame(top_marker_summary),
  row.names = FALSE
)


# ============================================================
# 15. CALCULATE GENE-LEVEL EXPRESSION STATISTICS
#
# For each canonical gene and each frozen cluster:
#
#   - Mean normalized expression
#   - Percent expressing
#
# This is more transparent than relying only on module scores.
# ============================================================

all_evidence_genes <- unique(
  c(
    unlist(
      marker_panels_present,
      use.names = FALSE
    ),
    unlist(
      subtype_panels_present,
      use.names = FALSE
    )
  )
)

all_evidence_genes <- intersect(
  all_evidence_genes,
  rownames(scrna)
)

cat(
  "\nCanonical/subtype genes available:",
  length(all_evidence_genes),
  "\n"
)

Idents(scrna) <- scrna$seurat_clusters


# ============================================================
# 16. AVERAGE EXPRESSION BY CLUSTER
# ============================================================

average_expression <- AverageExpression(
  scrna,
  assays = "RNA",
  features = all_evidence_genes,
  group.by = "seurat_clusters",
  slot = "data",
  verbose = FALSE
)$RNA

average_expression_df <- as.data.frame(
  average_expression
)

average_expression_df$Gene <- rownames(
  average_expression_df
)

average_expression_long <- average_expression_df %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Cluster",
    values_to = "Mean_Normalized_Expression"
  )

average_expression_long$Cluster <- sub(
  "^g",
  "",
  average_expression_long$Cluster
)


# ============================================================
# 17. PERCENT EXPRESSION BY CLUSTER
# ============================================================

rna_data <- GetAssayData(
  scrna,
  assay = "RNA",
  layer = "data"
)

rna_data <- rna_data[
  all_evidence_genes,
  ,
  drop = FALSE
]

cluster_vector <- as.character(
  scrna$seurat_clusters
)

percent_expression_list <- vector(
  mode = "list",
  length = length(cluster_ids)
)

names(percent_expression_list) <- as.character(
  cluster_ids
)

for (cluster_now in as.character(cluster_ids)) {
  
  cells_now <- colnames(scrna)[
    cluster_vector == cluster_now
  ]
  
  pct_now <- Matrix::rowMeans(
    rna_data[
      ,
      cells_now,
      drop = FALSE
    ] > 0
  ) * 100
  
  percent_expression_list[[cluster_now]] <- data.frame(
    Gene = names(pct_now),
    Cluster = cluster_now,
    Percent_Expressing = as.numeric(pct_now),
    stringsAsFactors = FALSE
  )
}

percent_expression_long <- bind_rows(
  percent_expression_list
)


# ============================================================
# 18. COMBINE GENE-LEVEL EVIDENCE
# ============================================================

gene_evidence <- average_expression_long %>%
  left_join(
    percent_expression_long,
    by = c(
      "Gene",
      "Cluster"
    )
  )

write.csv(
  gene_evidence,
  file.path(
    results_dir,
    "11D_B1_Canonical_Gene_Evidence_By_Cluster.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. CREATE PANEL MEMBERSHIP TABLE
# ============================================================

broad_membership <- bind_rows(
  lapply(
    names(marker_panels_present),
    function(panel_name) {
      
      data.frame(
        Gene = marker_panels_present[[panel_name]],
        Panel = panel_name,
        Panel_Type = "Broad identity",
        stringsAsFactors = FALSE
      )
    }
  )
)

subtype_membership <- bind_rows(
  lapply(
    names(subtype_panels_present),
    function(panel_name) {
      
      data.frame(
        Gene = subtype_panels_present[[panel_name]],
        Panel = panel_name,
        Panel_Type = "Subtype/state",
        stringsAsFactors = FALSE
      )
    }
  )
)

panel_membership <- bind_rows(
  broad_membership,
  subtype_membership
)

write.csv(
  panel_membership,
  file.path(
    results_dir,
    "11D_B1_Gene_to_Panel_Membership.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 20. PANEL-LEVEL EVIDENCE SUMMARY
#
# IMPORTANT:
# This is descriptive evidence.
# It is NOT an automatic annotation algorithm.
# ============================================================

panel_evidence <- gene_evidence %>%
  inner_join(
    panel_membership,
    by = "Gene"
  ) %>%
  group_by(
    Cluster,
    Panel_Type,
    Panel
  ) %>%
  summarise(
    Number_of_Genes = n_distinct(Gene),
    
    Mean_Expression =
      mean(
        Mean_Normalized_Expression,
        na.rm = TRUE
      ),
    
    Mean_Percent_Expressing =
      mean(
        Percent_Expressing,
        na.rm = TRUE
      ),
    
    Median_Percent_Expressing =
      median(
        Percent_Expressing,
        na.rm = TRUE
      ),
    
    Genes_Detected_5pct =
      sum(
        Percent_Expressing >= 5,
        na.rm = TRUE
      ),
    
    Genes_Detected_10pct =
      sum(
        Percent_Expressing >= 10,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

write.csv(
  panel_evidence,
  file.path(
    results_dir,
    "11D_B1_Panel_Evidence_By_Cluster.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 21. CURRENT ANNOTATION PER CLUSTER
# ============================================================

current_annotation <- scrna@meta.data %>%
  mutate(
    Cluster = as.character(seurat_clusters),
    Broad = as.character(broad_cell_type),
    Refined = as.character(refined_cell_type)
  ) %>%
  count(
    Cluster,
    Broad,
    Refined,
    name = "Cells"
  ) %>%
  group_by(Cluster) %>%
  arrange(
    desc(Cells),
    .by_group = TRUE
  ) %>%
  slice_head(n = 1) %>%
  ungroup()