# ============================================================
# STEP 11D-B2
# FINAL ANNOTATION VALIDATION
#
# PURPOSE:
# Validate the proposed biological identity of all 29 frozen
# clusters using positive lineage markers AND contradictory
# lineage markers before writing corrected annotations.
#
# IMPORTANT:
# - NO reclustering
# - NO filtering
# - NO CYP genes used for annotation
# - NO original annotations overwritten
# - NO corrected atlas saved
# - Proposed labels are tested, not automatically accepted
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

results_dir <- file.path(
  project_dir,
  "RESULTS",
  "11_SPATIAL",
  "11D_B2_FINAL_ANNOTATION_VALIDATION"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "11_SPATIAL",
  "11D_B2_FINAL_ANNOTATION_VALIDATION"
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
cat("STEP 11D-B2\n")
cat("FINAL ANNOTATION VALIDATION\n")
cat("============================================\n")

scrna <- readRDS(input_file)

DefaultAssay(scrna) <- "RNA"

cat("Cells:", ncol(scrna), "\n")
cat("Genes:", nrow(scrna), "\n")


# ============================================================
# 4. VERIFY METADATA
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
      "Missing metadata:",
      paste(missing_metadata, collapse = ", ")
    )
  )
}


# ============================================================
# 5. VERIFY EXPECTED 29 CLUSTERS
# ============================================================
# ============================================================
# 5. VERIFY EXPECTED 29 CLUSTERS
# CORRECTED — robust to integer/numeric type differences
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

# Additional strict checks
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
    paste(
      "Cluster IDs are not in the expected 0-28 sequence:",
      paste(cluster_ids, collapse = ", ")
    )
  )
}

cat(
  "Verified 29 frozen clusters: 0-28.\n"
)

cat(
  "Cluster IDs:",
  paste(cluster_ids, collapse = ", "),
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
# 7. PROPOSED ANNOTATION MAP
#
# IMPORTANT:
# This is a hypothesis being VALIDATED here.
# It is NOT written into the Seurat metadata.
# ============================================================

proposed_annotation <- data.frame(
  
  Cluster = as.character(0:28),
  
  Proposed_Broad = c(
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
    "Endothelial cells",       # 10
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
  
  Proposed_Refined = c(
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
    "Venous endothelial cells",        # 10
    "Fibroblasts",                     # 11
    "Platelets/megakaryocytes",        # 12
    "T cells",                         # 13
    "NK/cytotoxic lymphocytes",        # 14
    "Fibroblasts",                     # 15
    "Cycling myeloid cells",           # 16
    "Monocytes",                       # 17
    "Endothelial cells",               # 18
    "cDC-like dendritic cells",        # 19
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
  
  stringsAsFactors = FALSE
)

write.csv(
  proposed_annotation,
  file.path(
    results_dir,
    "11D_B2_PROPOSED_ANNOTATION_MAP.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 8. DEFINE CORE LINEAGE MARKER PANELS
#
# Deliberately conservative.
# ============================================================

lineage_markers <- list(
  
  Fibroblast = c(
    "Col1a1",
    "Col1a2",
    "Col3a1",
    "Dcn",
    "Lum",
    "Pdgfra",
    "Tcf21",
    "Dpt"
  ),
  
  Macrophage = c(
    "Adgre1",
    "Csf1r",
    "Cd68",
    "C1qa",
    "C1qb",
    "C1qc",
    "Fcgr1",
    "Mertk"
  ),
  
  Monocyte = c(
    "Ccr2",
    "Ly6c2",
    "Lyz2",
    "Plac8",
    "Ctss",
    "Fcgr3"
  ),
  
  Neutrophil = c(
    "S100a8",
    "S100a9",
    "Csf3r",
    "Cxcr2",
    "Ly6g",
    "Camp",
    "Ngp",
    "Lcn2"
  ),
  
  Endothelial = c(
    "Pecam1",
    "Cdh5",
    "Kdr",
    "Esam",
    "Egfl7",
    "Ptprb",
    "Tek"
  ),
  
  Cardiomyocyte = c(
    "Tnnt2",
    "Tnni3",
    "Myh6",
    "Actc1",
    "Myl2",
    "Pln",
    "Ryr2",
    "Mybpc3"
  ),
  
  T_cell = c(
    "Cd3d",
    "Cd3e",
    "Cd3g",
    "Trac",
    "Lck",
    "Cd247",
    "Itk"
  ),
  
  NK_cell = c(
    "Nkg7",
    "Ncr1",
    "Prf1",
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
    "Pax5",
    "Fcmr",
    "Ighm",
    "Ighd"
  ),
  
  Dendritic = c(
    "Flt3",
    "Itgax",
    "Clec9a",
    "Batf3",
    "Cd209a",
    "Wdfy4"
  ),
  
  Pericyte = c(
    "Rgs5",
    "Pdgfrb",
    "Cspg4",
    "Kcnj8",
    "Abcc9",
    "Notch3"
  ),
  
  Smooth_muscle = c(
    "Acta2",
    "Tagln",
    "Myh11",
    "Cnn1",
    "Myocd",
    "Kcnmb1"
  ),
  
  Platelet = c(
    "Pf4",
    "Ppbp",
    "Itga2b",
    "Gp9",
    "Gp5",
    "Gp6",
    "Clec1b"
  ),
  
  Erythroid = c(
    "Gypa",
    "Slc4a1",
    "Alas2",
    "Bpgm",
    "Fech",
    "Tspo2"
  )
)


# ============================================================
# 9. DEFINE VALIDATION / SUBTYPE PANELS
# ============================================================

validation_panels <- list(
  
  Resident_Macrophage = c(
    "Timd4",
    "Lyve1",
    "Folr2",
    "Cd163",
    "Mrc1",
    "Vsig4",
    "F13a1"
  ),
  
  Activated_Fibroblast = c(
    "Postn",
    "Ccn2",
    "Col8a1",
    "Loxl2",
    "Meox1"
  ),
  
  Capillary_Endothelial = c(
    "Car4",
    "Aplnr",
    "Aqp1",
    "Gpihbp1",
    "Rgcc"
  ),
  
  Arterial_Endothelial = c(
    "Gja5",
    "Sema3g",
    "Efnb2",
    "Sox17",
    "Bmx"
  ),
  
  Venous_Endothelial = c(
    "Ackr1",
    "Ackr2",
    "Nr2f2"
  ),
  
  Cycling = c(
    "Mki67",
    "Top2a",
    "Cdk1",
    "Ccnb1",
    "Birc5",
    "Pclaf",
    "Pbk"
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
# 10. DEFINE BROAD EXCLUSION MARKERS
#
# These are used to identify obvious lineage contradictions.
#
# They are NOT used as hard automatic rejection thresholds.
# ============================================================

exclusion_panels <- list(
  
  Hematopoietic = c(
    "Ptprc"
  ),
  
  Stromal_ECM = c(
    "Col1a1",
    "Col1a2",
    "Dcn",
    "Pdgfra"
  ),
  
  Endothelial = c(
    "Pecam1",
    "Cdh5",
    "Kdr"
  ),
  
  Cardiomyocyte = c(
    "Tnnt2",
    "Tnni3",
    "Myh6"
  ),
  
  T_cell = c(
    "Cd3d",
    "Cd3e",
    "Trac"
  ),
  
  B_cell = c(
    "Cd79a",
    "Ms4a1",
    "Cd19"
  ),
  
  Myeloid = c(
    "Csf1r",
    "Adgre1",
    "Lyz2"
  ),
  
  Smooth_muscle = c(
    "Myh11",
    "Tagln",
    "Acta2"
  )
)


# ============================================================
# 11. FILTER TO GENES PRESENT
# ============================================================

filter_panels <- function(panel_list, genes_present) {
  
  result <- lapply(
    panel_list,
    function(x) {
      intersect(
        x,
        genes_present
      )
    }
  )
  
  result[
    lengths(result) > 0
  ]
}

lineage_present <- filter_panels(
  lineage_markers,
  rownames(scrna)
)

validation_present <- filter_panels(
  validation_panels,
  rownames(scrna)
)

exclusion_present <- filter_panels(
  exclusion_panels,
  rownames(scrna)
)


# ============================================================
# 12. REPORT AVAILABLE MARKERS
# ============================================================

availability_table <- bind_rows(
  
  lapply(
    names(lineage_markers),
    function(panel) {
      
      present <- intersect(
        lineage_markers[[panel]],
        rownames(scrna)
      )
      
      data.frame(
        Category = "Core lineage",
        Panel = panel,
        Requested = length(lineage_markers[[panel]]),
        Present = length(present),
        Genes = paste(present, collapse = "; "),
        stringsAsFactors = FALSE
      )
    }
  ),
  
  lapply(
    names(validation_panels),
    function(panel) {
      
      present <- intersect(
        validation_panels[[panel]],
        rownames(scrna)
      )
      
      data.frame(
        Category = "Subtype/state",
        Panel = panel,
        Requested = length(validation_panels[[panel]]),
        Present = length(present),
        Genes = paste(present, collapse = "; "),
        stringsAsFactors = FALSE
      )
    }
  ),
  
  lapply(
    names(exclusion_panels),
    function(panel) {
      
      present <- intersect(
        exclusion_panels[[panel]],
        rownames(scrna)
      )
      
      data.frame(
        Category = "Exclusion",
        Panel = panel,
        Requested = length(exclusion_panels[[panel]]),
        Present = length(present),
        Genes = paste(present, collapse = "; "),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  availability_table,
  file.path(
    results_dir,
    "11D_B2_Marker_Availability.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. GET NORMALIZED EXPRESSION MATRIX
# ============================================================

all_validation_genes <- unique(
  c(
    unlist(lineage_present, use.names = FALSE),
    unlist(validation_present, use.names = FALSE),
    unlist(exclusion_present, use.names = FALSE)
  )
)

all_validation_genes <- intersect(
  all_validation_genes,
  rownames(scrna)
)

rna_data <- GetAssayData(
  scrna,
  assay = "RNA",
  layer = "data"
)

rna_validation <- rna_data[
  all_validation_genes,
  ,
  drop = FALSE
]

cluster_vector <- as.character(
  scrna$seurat_clusters
)


# ============================================================
# 14. GENE-LEVEL EXPRESSION AND DETECTION BY CLUSTER
# ============================================================

gene_stats_list <- vector(
  "list",
  length(cluster_ids)
)

names(gene_stats_list) <- as.character(
  cluster_ids
)

for (cluster_now in as.character(cluster_ids)) {
  
  cells_now <- colnames(scrna)[
    cluster_vector == cluster_now
  ]
  
  matrix_now <- rna_validation[
    ,
    cells_now,
    drop = FALSE
  ]
  
  mean_expr <- Matrix::rowMeans(
    matrix_now
  )
  
  pct_expr <- Matrix::rowMeans(
    matrix_now > 0
  ) * 100
  
  gene_stats_list[[cluster_now]] <- data.frame(
    Cluster = cluster_now,
    Gene = rownames(matrix_now),
    Mean_Expression = as.numeric(mean_expr),
    Percent_Expressing = as.numeric(pct_expr),
    stringsAsFactors = FALSE
  )
}

gene_stats <- bind_rows(
  gene_stats_list
)

write.csv(
  gene_stats,
  file.path(
    results_dir,
    "11D_B2_Gene_Level_Validation.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. CREATE PANEL MEMBERSHIP TABLE
# ============================================================

make_membership <- function(panel_list, category) {
  
  bind_rows(
    lapply(
      names(panel_list),
      function(panel) {
        
        data.frame(
          Gene = panel_list[[panel]],
          Panel = panel,
          Category = category,
          stringsAsFactors = FALSE
        )
      }
    )
  )
}

panel_membership <- bind_rows(
  
  make_membership(
    lineage_present,
    "Core lineage"
  ),
  
  make_membership(
    validation_present,
    "Subtype/state"
  ),
  
  make_membership(
    exclusion_present,
    "Exclusion"
  )
)

write.csv(
  panel_membership,
  file.path(
    results_dir,
    "11D_B2_Panel_Membership.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. PANEL-LEVEL SUMMARY
# ============================================================

panel_summary <- gene_stats %>%
  inner_join(
    panel_membership,
    by = "Gene"
  ) %>%
  group_by(
    Cluster,
    Category,
    Panel
  ) %>%
  summarise(
    
    Number_of_Genes = n_distinct(Gene),
    
    Genes_Expressed_5pct =
      sum(
        Percent_Expressing >= 5,
        na.rm = TRUE
      ),
    
    Genes_Expressed_10pct =
      sum(
        Percent_Expressing >= 10,
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
    
    Mean_Normalized_Expression =
      mean(
        Mean_Expression,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

write.csv(
  panel_summary,
  file.path(
    results_dir,
    "11D_B2_Panel_Validation_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. CORE LINEAGE RANKING
#
# This is descriptive.
# Highest ranking does NOT automatically become annotation.
# ============================================================

lineage_ranking <- panel_summary %>%
  filter(
    Category == "Core lineage"
  ) %>%
  group_by(Cluster) %>%
  arrange(
    desc(Genes_Expressed_10pct),
    desc(Mean_Percent_Expressing),
    desc(Mean_Normalized_Expression),
    .by_group = TRUE
  ) %>%
  mutate(
    Evidence_Rank = row_number()
  ) %>%
  ungroup()

write.csv(
  lineage_ranking,
  file.path(
    results_dir,
    "11D_B2_Core_Lineage_Ranking.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. TOP THREE LINEAGE PROGRAMS
# ============================================================

top3_lineages <- lineage_ranking %>%
  filter(
    Evidence_Rank <= 3
  )

write.csv(
  top3_lineages,
  file.path(
    results_dir,
    "11D_B2_Top3_Lineages_Per_Cluster.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. PROPOSED BROAD LABEL -> EXPECTED PANEL
# ============================================================

broad_to_panel <- c(
  
  "Fibroblasts" =
    "Fibroblast",
  
  "Macrophages" =
    "Macrophage",
  
  "B cells" =
    "B_cell",
  
  "Neutrophils" =
    "Neutrophil",
  
  "Endothelial cells" =
    "Endothelial",
  
  "Platelets" =
    "Platelet",
  
  "T cells" =
    "T_cell",
  
  "NK cells" =
    "NK_cell",
  
  "Myeloid cells" =
    "Monocyte",
  
  "Monocytes" =
    "Monocyte",
  
  "Dendritic cells" =
    "Dendritic",
  
  "Cardiomyocytes" =
    "Cardiomyocyte",
  
  "Erythroid cells" =
    "Erythroid",
  
  "Smooth muscle cells" =
    "Smooth_muscle",
  
  "Pericytes" =
    "Pericyte"
)


# ============================================================
# 20. EXTRACT EVIDENCE FOR EACH PROPOSED IDENTITY
# ============================================================

expected_panel_table <- proposed_annotation %>%
  mutate(
    Expected_Core_Panel =
      unname(
        broad_to_panel[
          Proposed_Broad
        ]
      )
  )

proposed_evidence <- expected_panel_table %>%
  left_join(
    lineage_ranking %>%
      select(
        Cluster,