# ============================================================
# STEP 11D-A
# FROZEN ATLAS ANNOTATION AUDIT
#
# PURPOSE:
# Diagnose the unexpected mismatch between several broad
# cell-type labels and their marker signatures BEFORE spatial
# projection.
#
# IMPORTANT:
# - NO reclustering
# - NO filtering
# - NO annotation changes
# - NO CYP-based annotation
# - NO object overwrite
# - Diagnostic only
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
  "11D_ANNOTATION_AUDIT"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "11_SPATIAL",
  "11D_ANNOTATION_AUDIT"
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
# 3. LOAD FROZEN ATLAS
# ============================================================

cat("\n============================================\n")
cat("LOADING FROZEN ATLAS\n")
cat("============================================\n")

scrna <- readRDS(input_file)

cat(
  "Cells:",
  ncol(scrna),
  "\n"
)

cat(
  "Genes:",
  nrow(scrna),
  "\n"
)

cat(
  "Active assay:",
  DefaultAssay(scrna),
  "\n"
)


# ============================================================
# 4. VERIFY REQUIRED METADATA
# ============================================================

required_metadata <- c(
  "seurat_clusters",
  "broad_cell_type",
  "refined_cell_type"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(scrna@meta.data)
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

cat("\nRequired metadata fields found.\n")


# ============================================================
# 5. FREEZE SAFETY INFORMATION
# ============================================================

original_cells <- colnames(scrna)

original_cluster <- as.character(
  scrna$seurat_clusters
)

original_broad <- as.character(
  scrna$broad_cell_type
)

original_refined <- as.character(
  scrna$refined_cell_type
)


# ============================================================
# 6. CLUSTER × BROAD CELL-TYPE TABLE
# ============================================================

cat("\n============================================\n")
cat("CLUSTER x BROAD CELL TYPE\n")
cat("============================================\n")

cluster_broad_table <- table(
  Cluster = scrna$seurat_clusters,
  Broad_Cell_Type = scrna$broad_cell_type
)

print(
  cluster_broad_table
)

write.csv(
  as.data.frame.matrix(cluster_broad_table),
  file.path(
    results_dir,
    "11D_A_Cluster_vs_Broad_CellType.csv"
  )
)


# ============================================================
# 7. CLUSTER × REFINED CELL-TYPE TABLE
# ============================================================

cat("\n============================================\n")
cat("CLUSTER x REFINED CELL TYPE\n")
cat("============================================\n")

cluster_refined_table <- table(
  Cluster = scrna$seurat_clusters,
  Refined_Cell_Type = scrna$refined_cell_type
)

print(
  cluster_refined_table
)

write.csv(
  as.data.frame.matrix(cluster_refined_table),
  file.path(
    results_dir,
    "11D_A_Cluster_vs_Refined_CellType.csv"
  )
)


# ============================================================
# 8. BROAD × REFINED CELL-TYPE TABLE
# ============================================================

cat("\n============================================\n")
cat("BROAD x REFINED CELL TYPE\n")
cat("============================================\n")

broad_refined_table <- table(
  Broad_Cell_Type = scrna$broad_cell_type,
  Refined_Cell_Type = scrna$refined_cell_type
)

print(
  broad_refined_table
)

write.csv(
  as.data.frame.matrix(broad_refined_table),
  file.path(
    results_dir,
    "11D_A_Broad_vs_Refined_CellType.csv"
  )
)


# ============================================================
# 9. DOMINANT LABEL PER CLUSTER
# ============================================================

cluster_annotation_summary <- scrna@meta.data %>%
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
  group_by(
    Cluster
  ) %>%
  mutate(
    Cluster_Total = sum(Cells),
    Percent_of_Cluster =
      100 * Cells / Cluster_Total
  ) %>%
  ungroup() %>%
  arrange(
    as.numeric(Cluster),
    desc(Cells)
  )

write.csv(
  cluster_annotation_summary,
  file.path(
    results_dir,
    "11D_A_Cluster_Annotation_Summary.csv"
  ),
  row.names = FALSE
)

cat("\n============================================\n")
cat("DOMINANT ANNOTATION PER CLUSTER\n")
cat("============================================\n")

dominant_cluster_annotation <-
  cluster_annotation_summary %>%
  group_by(
    Cluster
  ) %>%
  slice_max(
    order_by = Cells,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

print(
  as.data.frame(dominant_cluster_annotation),
  row.names = FALSE
)


# ============================================================
# 10. DEFINE INDEPENDENT CANONICAL MARKER PANEL
#
# These genes are used ONLY to audit biological identity.
# They do not change annotations.
# ============================================================

canonical_markers <- list(
  
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
    "Atp2a2"
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
    "Pcolce2"
  ),
  
  Endothelial = c(
    "Pecam1",
    "Cdh5",
    "Kdr",
    "Emcn",
    "Esam",
    "Eng",
    "Egfl7",
    "Ptprb",
    "Tek",
    "Vwf"
  ),
  
  Pericyte_Mural = c(
    "Rgs5",
    "Pdgfrb",
    "Cspg4",
    "Des",
    "Notch3",
    "Kcnj8",
    "Abcc9",
    "Cnn1",
    "Acta2",
    "Rbp1"
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
    "Lyz2",
    "Mertk"
  ),
  
  Monocyte_Myeloid = c(
    "Lyz2",
    "Ccr2",
    "Ly6c2",
    "S100a8",
    "S100a9",
    "Ctss",
    "Fcgr3",
    "Tyrobp",
    "Aif1",
    "Csf1r"
  ),
  
  Dendritic = c(
    "Flt3",
    "Itgax",
    "Clec9a",
    "Batf3",
    "Cd209a",
    "H2-DMa",
    "H2-DMb1",
    "Ciita",
    "Wdfy4"
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
    "Ccr7"
  ),
  
  NK_Cytotoxic = c(
    "Nkg7",
    "Ncr1",
    "Prf1",
    "Gzmb",
    "Ccl5",
    "Klrc1",
    "Klrd1",
    "Klrk1"
  ),
  
  B_cell = c(
    "Cd79a",
    "Cd79b",
    "Ms4a1",
    "Cd19",
    "Cd22",
    "Cd37",
    "H2-Ob",
    "Cd74",
    "Ighm"
  ),
  
  Neutrophil_Granulocyte = c(
    "S100a8",
    "S100a9",
    "Csf3r",
    "Cxcr2",
    "Retnlg",
    "Mmp8",
    "Mmp9",
    "Lcn2",
    "Trem1"
  ),
  
  Platelet_Megakaryocyte = c(
    "Pf4",
    "Ppbp",
    "Gp9",
    "Gp5",
    "Gp6",
    "Itga2b",
    "Treml1",
    "Clec1b"
  ),
  
  Erythroid = c(
    "Gypa",
    "Slc4a1",
    "Alas2",
    "Bpgm",
    "Fech",
    "Slc25a37",
    "Tspo2"
  ),
  
  Cycling = c(
    "Mki67",
    "Top2a",
    "Cdk1",
    "Ccnb1",
    "Ccnb2",
    "Birc5",
    "Pclaf",
    "Tpx2"
  )
)


# ============================================================
# 11. CHECK MARKER AVAILABILITY
# ============================================================

all_marker_genes <- unique(
  unlist(
    canonical_markers,
    use.names = FALSE
  )
)

present_markers <- intersect(
  all_marker_genes,
  rownames(scrna)
)

missing_markers <- setdiff(
  all_marker_genes,
  rownames(scrna)
)

cat("\n============================================\n")
cat("CANONICAL MARKER AVAILABILITY\n")
cat("============================================\n")

cat(
  "Requested:",
  length(all_marker_genes),
  "\n"
)

cat(
  "Present:",
  length(present_markers),
  "\n"
)

cat(
  "Missing:",
  length(missing_markers),
  "\n"
)

if (length(missing_markers) > 0) {
  cat(
    "\nMissing genes:\n",
    paste(
      missing_markers,
      collapse = ", "
    ),
    "\n"
  )
}


# ============================================================
# 12. CANONICAL MARKER DOTPLOT BY CLUSTER
# ============================================================
# ============================================================
# 12. CANONICAL MARKER DOTPLOT BY CLUSTER
# CORRECTED: removes duplicate genes for plotting only
# ============================================================

DefaultAssay(scrna) <- "RNA"

Idents(scrna) <- scrna$seurat_clusters

canonical_markers_present <- lapply(
  canonical_markers,
  function(x) {
    intersect(
      x,
      rownames(scrna)
    )
  }
)

canonical_markers_present <-
  canonical_markers_present[
    lengths(canonical_markers_present) > 0
  ]

# ------------------------------------------------------------
# IMPORTANT:
# Some canonical genes intentionally occur in more than one
# biological program.
#
# DotPlot cannot safely plot duplicated feature names when a
# named feature list is supplied.
#
# Therefore:
# - keep the ORIGINAL overlapping lists for module scoring
# - create a UNIQUE vector only for visualization
# ------------------------------------------------------------

canonical_plot_genes <- unique(
  unlist(
    canonical_markers_present,
    use.names = FALSE
  )
)

cat(
  "\nCanonical genes available for scoring:",
  sum(lengths(canonical_markers_present)),
  "\n"
)

cat(
  "Unique canonical genes used in DotPlot:",
  length(canonical_plot_genes),
  "\n"
)

duplicated_plot_genes <- names(
  which(
    table(
      unlist(
        canonical_markers_present,
        use.names = FALSE
      )
    ) > 1
  )
)

if (length(duplicated_plot_genes) > 0) {
  
  cat(
    "\nGenes shared between canonical programs:\n"
  )
  
  cat(
    paste(
      duplicated_plot_genes,
      collapse = ", "
    ),
    "\n"
  )
}

p_cluster_dot <- DotPlot(
  object = scrna,
  features = canonical_plot_genes,
  assay = "RNA",
  dot.scale = 5
) +
  RotatedAxis() +
  labs(
    title = "Independent canonical-marker audit by frozen cluster",
    x = NULL,
    y = "Frozen Seurat cluster"
  ) +
  theme(
    axis.text.x = element_text(
      size = 7
    ),
    axis.text.y = element_text(
      size = 9
    )
  )

ggsave(
  filename = file.path(
    figures_dir,
    "11D_A_Canonical_Markers_By_Cluster.png"
  ),
  plot = p_cluster_dot,
  width = 22,
  height = 10,
  dpi = 400
)

cat("\nCluster canonical-marker DotPlot saved successfully.\n")

# ============================================================
# 13. CANONICAL MARKER DOTPLOT BY EXISTING BROAD LABEL
# ============================================================
# ============================================================
# 13. CANONICAL MARKER DOTPLOT BY EXISTING BROAD LABEL
# CORRECTED
# ============================================================

Idents(scrna) <- scrna$broad_cell_type

p_broad_dot <- DotPlot(
  object = scrna,
  features = canonical_plot_genes,
  assay = "RNA",
  dot.scale = 5
) +
  RotatedAxis() +
  labs(
    title = "Canonical-marker audit by existing broad annotation",
    x = NULL,
    y = "Existing broad cell type"
  ) +
  theme(
    axis.text.x = element_text(
      size = 7
    ),
    axis.text.y = element_text(
      size = 9
    )
  )

ggsave(
  filename = file.path(
    figures_dir,
    "11D_A_Canonical_Markers_By_Broad_Label.png"
  ),
  plot = p_broad_dot,
  width = 22,
  height = 8,
  dpi = 400
)

cat("\nBroad-label canonical-marker DotPlot saved successfully.\n")
# ============================================================
# 14. CANONICAL MARKER DOTPLOT BY REFINED LABEL
# ============================================================

# ============================================================
# 14. CANONICAL MARKER DOTPLOT BY REFINED LABEL
# CORRECTED
# ============================================================

Idents(scrna) <- scrna$refined_cell_type

p_refined_dot <- DotPlot(
  object = scrna,
  features = canonical_plot_genes,
  assay = "RNA",
  dot.scale = 5
) +
  RotatedAxis() +
  labs(
    title = "Canonical-marker audit by existing refined annotation",
    x = NULL,
    y = "Existing refined cell type"
  ) +
  theme(
    axis.text.x = element_text(
      size = 7
    ),
    axis.text.y = element_text(
      size = 8
    )
  )

ggsave(
  filename = file.path(
    figures_dir,
    "11D_A_Canonical_Markers_By_Refined_Label.png"
  ),
  plot = p_refined_dot,
  width = 22,
  height = 9,
  dpi = 400
)

cat("\nRefined-label canonical-marker DotPlot saved successfully.\n")

# ============================================================
# 15. COMPUTE CANONICAL PROGRAM SCORES
#
# These scores provide an independent numerical audit.
#
# They are diagnostic only.
# ============================================================

cat("\n============================================\n")
cat("CALCULATING CANONICAL IDENTITY SCORES\n")
cat("============================================\n")

audit_score_names <- character(0)

for (program_name in names(canonical_markers_present)) {
  
  genes_now <- canonical_markers_present[[program_name]]
  
  if (length(genes_now) < 2) {
    next
  }
  
  old_columns <- colnames(
    scrna@meta.data
  )
  
  scrna <- AddModuleScore(
    object = scrna,
    features = list(genes_now),
    name = paste0(
      "AUDIT_",
      program_name,
      "_"
    ),
    assay = "RNA",
    search = FALSE
  )
  
  new_column <- setdiff(
    colnames(scrna@meta.data),
    old_columns
  )
  
  if (length(new_column) != 1) {
    stop(
      paste(
        "Could not uniquely identify score column for",
        program_name
      )
    )
  }
  
  desired_name <- paste0(
    "AUDIT_",
    program_name
  )
  
  colnames(scrna@meta.data)[
    colnames(scrna@meta.data) == new_column
  ] <- desired_name
  
  audit_score_names <- c(
    audit_score_names,
    desired_name
  )
}


# ============================================================
# 16. MEAN CANONICAL SCORES BY CLUSTER
# ============================================================

audit_metadata <- scrna@meta.data %>%
  mutate(
    Cluster = as.character(seurat_clusters),
    Broad_Cell_Type =
      as.character(broad_cell_type),
    Refined_Cell_Type =
      as.character(refined_cell_type)
  )

cluster_score_summary <- audit_metadata %>%
  group_by(
    Cluster
  ) %>%
  summarise(
    across(
      all_of(audit_score_names),
      ~ mean(
        .x,
        na.rm = TRUE
      )
    ),
    .groups = "drop"
  )

write.csv(
  cluster_score_summary,
  file.path(
    results_dir,
    "11D_A_Canonical_Scores_By_Cluster.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. DOMINANT CANONICAL PROGRAM PER CLUSTER
# ============================================================

cluster_score_long <- cluster_score_summary %>%
  pivot_longer(
    cols = all_of(audit_score_names),
    names_to = "Canonical_Program",
    values_to = "Mean_Score"
  ) %>%
  mutate(
    Canonical_Program =
      sub(
        "^AUDIT_",
        "",
        Canonical_Program
      )
  )

dominant_canonical_program <- cluster_score_long %>%
  group_by(
    Cluster
  ) %>%
  slice_max(
    order_by = Mean_Score,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(
    as.numeric(Cluster)
  )

dominant_canonical_program <-
  dominant_canonical_program %>%
  left_join(
    dominant_cluster_annotation %>%
      select(
        Cluster,
        Broad,
        Refined,
        Cells,
        Cluster_Total,
        Percent_of_Cluster
      ),
    by = "Cluster"
  )

write.csv(
  dominant_canonical_program,
  file.path(
    results_dir,
    "11D_A_Dominant_Canonical_Program_Per_Cluster.csv"
  ),
  row.names = FALSE
)

cat("\n============================================\n")
cat("DOMINANT CANONICAL PROGRAM PER CLUSTER\n")
cat("============================================\n")

print(
  as.data.frame(dominant_canonical_program),
  row.names = FALSE
)


# ============================================================
# 18. HEATMAP OF CANONICAL SCORES
#
# ggplot heatmap avoids modifying the object or requiring
# additional heatmap packages.
# ============================================================

heatmap_df <- cluster_score_long %>%
  mutate(
    Cluster = factor(
      Cluster,
      levels = as.character(
        sort(
          unique(
            as.numeric(Cluster)
          )
        )
      )
    )
  )

p_heatmap <- ggplot(
  heatmap_df,
  aes(
    x = Canonical_Program,
    y = Cluster,
    fill = Mean_Score
  )
) +
  geom_tile() +
  labs(
    title = "Canonical cell-identity scores across frozen clusters",
    x = NULL,
    y = "Frozen Seurat cluster",
    fill = "Mean\nscore"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    panel.grid = element_blank()
  )

ggsave(
  filename = file.path(
    figures_dir,
    "11D_A_Canonical_Program_Heatmap_By_Cluster.png"
  ),
  plot = p_heatmap,
  width = 12,
  height = 9,
  dpi = 400
)


# ============================================================
# 19. TOP POSITIVE MARKERS PER CLUSTER
#
# This is independent of broad labels.
# We use frozen cluster IDs directly.
# ============================================================

cat("\n============================================\n")
cat("DERIVING CLUSTER-LEVEL AUDIT MARKERS\n")
cat("============================================\n")

DefaultAssay(scrna) <- "RNA"

Idents(scrna) <- scrna$seurat_clusters

cluster_markers <- FindAllMarkers(
  object = scrna,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.20,
  logfc.threshold = 0.50,
  test.use = "wilcox",
  return.thresh = 0.05,
  verbose = TRUE
)

if (nrow(cluster_markers) == 0) {
  stop(
    "Cluster-level marker analysis returned zero markers."
  )
}

write.csv(
  cluster_markers,
  file.path(
    results_dir,
    "11D_A_All_Cluster_Audit_Markers.csv"
  ),
  row.names = FALSE
)
