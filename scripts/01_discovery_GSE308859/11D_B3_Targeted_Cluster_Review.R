# ============================================================
# STEP 11D-B3
# TARGETED RESOLUTION OF AMBIGUOUS CLUSTERS
#
# Clusters:
# 10 = Endothelial vs Fibroblast
# 16 = Cycling macrophage/monocyte/myeloid
# 19 = Dendritic vs Macrophage
# 24 = Smooth muscle vs Fibroblast
#
# NO annotation changes are made.
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

results_dir <- file.path(
  project_dir,
  "RESULTS",
  "11_SPATIAL",
  "11D_B3_TARGETED_CLUSTER_REVIEW"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "11_SPATIAL",
  "11D_B3_TARGETED_CLUSTER_REVIEW"
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

scrna <- readRDS(input_file)

DefaultAssay(scrna) <- "RNA"

Idents(scrna) <- scrna$seurat_clusters

target_clusters <- c(
  "10",
  "16",
  "19",
  "24"
)

cat("\n============================================\n")
cat("STEP 11D-B3\n")
cat("TARGETED CLUSTER REVIEW\n")
cat("============================================\n")

cat("Cells:", ncol(scrna), "\n")
cat(
  "Target clusters:",
  paste(target_clusters, collapse = ", "),
  "\n"
)


# ============================================================
# 4. TARGETED MARKER PANELS
# ============================================================

target_panels <- list(
  
  # ----------------------------------------------------------
  # CLUSTER 10
  # ----------------------------------------------------------
  
  Endothelial_Core = c(
    "Pecam1",
    "Cdh5",
    "Kdr",
    "Ptprb",
    "Esam",
    "Tek",
    "Egfl7",
    "Emcn",
    "Klf2",
    "Klf4"
  ),
  
  Venous_Endothelial = c(
    "Ackr1",
    "Ackr2",
    "Nr2f2",
    "Vwf",
    "Selp",
    "Sele"
  ),
  
  Lymphatic_Endothelial = c(
    "Prox1",
    "Pdpn",
    "Flt4",
    "Ccl21a",
    "Lyve1"
  ),
  
  Fibroblast_Core = c(
    "Col1a1",
    "Col1a2",
    "Col3a1",
    "Dcn",
    "Lum",
    "Pdgfra",
    "Tcf21",
    "Dpt",
    "Col5a1",
    "Col5a2"
  ),
  
  # ----------------------------------------------------------
  # CLUSTER 16
  # ----------------------------------------------------------
  
  Macrophage_Core = c(
    "Adgre1",
    "Csf1r",
    "Cd68",
    "C1qa",
    "C1qb",
    "C1qc",
    "Fcgr1",
    "Mertk",
    "Ms4a7"
  ),
  
  Monocyte_Core = c(
    "Ccr2",
    "Ly6c2",
    "Plac8",
    "Ctss",
    "Fcgr3",
    "Lyz2",
    "Ms4a8a"
  ),
  
  Cycling = c(
    "Mki67",
    "Top2a",
    "Cdk1",
    "Ccnb1",
    "Ccnb2",
    "Birc5",
    "Pclaf",
    "Pbk",
    "Ncapg",
    "Tpx2"
  ),
  
  # ----------------------------------------------------------
  # CLUSTER 19
  # ----------------------------------------------------------
  
  Dendritic_Core = c(
    "Flt3",
    "Itgax",
    "Clec9a",
    "Batf3",
    "Wdfy4",
    "Cd209a",
    "Zbtb46"
  ),
  
  cDC1 = c(
    "Clec9a",
    "Xcr1",
    "Batf3",
    "Wdfy4",
    "Cadm1"
  ),
  
  cDC2 = c(
    "Cd209a",
    "Sirpa",
    "Clec10a",
    "Itgam"
  ),
  
  # ----------------------------------------------------------
  # CLUSTER 24
  # ----------------------------------------------------------
  
  Smooth_Muscle_Core = c(
    "Myh11",
    "Acta2",
    "Tagln",
    "Cnn1",
    "Myocd",
    "Lmod1",
    "Kcnmb1",
    "Mylk",
    "Actg2"
  ),
  
  Pericyte_Core = c(
    "Rgs5",
    "Pdgfrb",
    "Cspg4",
    "Kcnj8",
    "Abcc9",
    "Notch3"
  )
)


# ============================================================
# 5. KEEP ONLY GENES PRESENT
# ============================================================

target_panels_present <- lapply(
  target_panels,
  function(x) {
    intersect(
      x,
      rownames(scrna)
    )
  }
)

target_panels_present <-
  target_panels_present[
    lengths(target_panels_present) > 0
  ]


# ============================================================
# 6. SAVE MARKER AVAILABILITY
# ============================================================

marker_availability <- bind_rows(
  lapply(
    names(target_panels),
    function(panel) {
      
      present <- intersect(
        target_panels[[panel]],
        rownames(scrna)
      )
      
      data.frame(
        Panel = panel,
        Requested = length(target_panels[[panel]]),
        Present = length(present),
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
  marker_availability,
  file.path(
    results_dir,
    "11D_B3_Marker_Availability.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 7. TARGET CLUSTER CELL COUNTS
# ============================================================

target_counts <- as.data.frame(
  table(
    Cluster = as.character(
      scrna$seurat_clusters
    )
  )
) %>%
  filter(
    Cluster %in% target_clusters
  )

write.csv(
  target_counts,
  file.path(
    results_dir,
    "11D_B3_Target_Cluster_Cell_Counts.csv"
  ),
  row.names = FALSE
)

print(
  target_counts,
  row.names = FALSE
)


# ============================================================
# 8. DIFFERENTIAL MARKERS FOR EACH TARGET CLUSTER
#
# Independent cluster-vs-rest evidence.
# ============================================================

target_markers_list <- list()

for (cl in target_clusters) {
  
  cat(
    "\nFinding markers for cluster",
    cl,
    "...\n"
  )
  
  markers_now <- FindMarkers(
    object = scrna,
    ident.1 = cl,
    assay = "RNA",
    slot = "data",
    only.pos = TRUE,
    min.pct = 0.05,
    logfc.threshold = 0.10
  )
  
  markers_now$Gene <- rownames(
    markers_now
  )
  
  markers_now$Cluster <- cl
  
  markers_now <- markers_now %>%
    arrange(
      desc(avg_log2FC)
    )
  
  target_markers_list[[cl]] <-
    markers_now
  
  write.csv(
    markers_now,
    file.path(
      results_dir,
      paste0(
        "11D_B3_Cluster_",
        cl,
        "_Markers.csv"
      )
    ),
    row.names = FALSE
  )
}

target_markers <- bind_rows(
  target_markers_list
)

write.csv(
  target_markers,
  file.path(
    results_dir,
    "11D_B3_All_Target_Cluster_Markers.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9. TOP 30 MARKERS FOR EACH TARGET CLUSTER
# ============================================================

top30_markers <- target_markers %>%
  group_by(Cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 30,
    with_ties = FALSE
  ) %>%
  ungroup()

write.csv(
  top30_markers,
  file.path(
    results_dir,
    "11D_B3_Top30_Markers_Target_Clusters.csv"
  ),
  row.names = FALSE
)

cat("\n============================================\n")
cat("TOP 30 MARKERS\n")
cat("============================================\n")

print(
  as.data.frame(
    top30_markers %>%
      select(
        Cluster,
        Gene,
        avg_log2FC,
        pct.1,
        pct.2,
        p_val_adj
      )
  ),
  row.names = FALSE
)


# ============================================================
# 10. NORMALIZED EXPRESSION MATRIX
# ============================================================

all_target_genes <- unique(
  unlist(
    target_panels_present,
    use.names = FALSE
  )
)

rna_data <- GetAssayData(
  scrna,
  assay = "RNA",
  layer = "data"
)

rna_target <- rna_data[
  all_target_genes,
  ,
  drop = FALSE
]


# ============================================================
# 11. GENE-LEVEL TARGET STATISTICS
# ============================================================

target_gene_stats_list <- list()

for (cl in target_clusters) {
  
  cells_now <- colnames(scrna)[
    as.character(
      scrna$seurat_clusters
    ) == cl
  ]
  
  mat_now <- rna_target[
    ,
    cells_now,
    drop = FALSE
  ]
  
  target_gene_stats_list[[cl]] <-
    data.frame(
      
      Cluster = cl,
      
      Gene = rownames(mat_now),
      
      Mean_Expression =
        as.numeric(
          Matrix::rowMeans(
            mat_now
          )
        ),
      
      Percent_Expressing =
        as.numeric(
          Matrix::rowMeans(
            mat_now > 0
          ) * 100
        ),
      
      stringsAsFactors = FALSE
    )
}

target_gene_stats <- bind_rows(
  target_gene_stats_list
)

write.csv(
  target_gene_stats,
  file.path(
    results_dir,
    "11D_B3_Target_Gene_Statistics.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. PANEL MEMBERSHIP
# ============================================================

target_membership <- bind_rows(
  lapply(
    names(target_panels_present),
    function(panel) {
      
      data.frame(
        Gene =
          target_panels_present[[panel]],
        
        Panel =
          panel,
        
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  target_membership,
  file.path(
    results_dir,
    "11D_B3_Target_Panel_Membership.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. PANEL SUMMARY
#
# relationship = "many-to-many" is intentional because:
# - each gene occurs in multiple clusters
# - some genes belong to multiple biological panels
# ============================================================

target_panel_summary <-
  target_gene_stats %>%
  
  inner_join(
    target_membership,
    by = "Gene",
    relationship = "many-to-many"
  ) %>%
  
  group_by(
    Cluster,
    Panel
  ) %>%
  
  summarise(
    
    Number_of_Genes =
      n_distinct(Gene),
    
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
  target_panel_summary,
  file.path(
    results_dir,
    "11D_B3_Target_Panel_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. PRINT ALL TARGET PANEL EVIDENCE
# ============================================================

cat("\n============================================\n")
cat("TARGET PANEL EVIDENCE\n")
cat("============================================\n")

print(
  as.data.frame(
    target_panel_summary %>%
      arrange(
        as.numeric(Cluster),
        desc(Mean_Percent_Expressing)
      )
  ),
  row.names = FALSE
)


# ============================================================
# 15. CLUSTER 10 SPECIFIC REVIEW
# ============================================================

cluster10_review <- target_panel_summary %>%
  filter(
    Cluster == "10",
    Panel %in% c(
      "Endothelial_Core",
      "Venous_Endothelial",
      "Lymphatic_Endothelial",
      "Fibroblast_Core"
    )
  ) %>%
  arrange(
    desc(Mean_Percent_Expressing)
  )

write.csv(
  cluster10_review,
  file.path(
    results_dir,
    "11D_B3_Cluster10_Endothelial_vs_Fibroblast.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. CLUSTER 16 SPECIFIC REVIEW
# ============================================================

cluster16_review <- target_panel_summary %>%
  filter(
    Cluster == "16",
    Panel %in% c(
      "Macrophage_Core",
      "Monocyte_Core",
      "Cycling"
    )
  ) %>%
  arrange(
    desc(Mean_Percent_Expressing)
  )

write.csv(
  cluster16_review,
  file.path(
    results_dir,
    "11D_B3_Cluster16_Myeloid_Review.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. CLUSTER 19 SPECIFIC REVIEW
# ============================================================

cluster19_review <- target_panel_summary %>%
  filter(
    Cluster == "19",
    Panel %in% c(
      "Dendritic_Core",
      "cDC1",
      "cDC2",
      "Macrophage_Core",
      "Monocyte_Core"
    )
  ) %>%
  arrange(
    desc(Mean_Percent_Expressing)
  )

write.csv(
  cluster19_review,
  file.path(
    results_dir,
    "11D_B3_Cluster19_DC_vs_Myeloid.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. CLUSTER 24 SPECIFIC REVIEW
# ============================================================

cluster24_review <- target_panel_summary %>%
  filter(
    Cluster == "24",
    Panel %in% c(
      "Smooth_Muscle_Core",
      "Pericyte_Core",
      "Fibroblast_Core"
    )
  ) %>%
  arrange(
    desc(Mean_Percent_Expressing)
  )

write.csv(
  cluster24_review,
  file.path(
    results_dir,
    "11D_B3_Cluster24_Mural_vs_Fibroblast.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. PRINT FOUR FOCUSED REVIEWS
# ============================================================

cat("\n\n============================================\n")
cat("CLUSTER 10 REVIEW\n")
cat("============================================\n")

print(
  as.data.frame(cluster10_review),
  row.names = FALSE
)

cat("\n============================================\n")
cat("CLUSTER 16 REVIEW\n")
cat("============================================\n")

print(
  as.data.frame(cluster16_review),
  row.names = FALSE
)

cat("\n============================================\n")
cat("CLUSTER 19 REVIEW\n")
cat("============================================\n")

print(
  as.data.frame(cluster19_review),
  row.names = FALSE
)

cat("\n============================================\n")
cat("CLUSTER 24 REVIEW\n")
cat("============================================\n")

print(
  as.data.frame(cluster24_review),
  row.names = FALSE
)


# ============================================================
# 20. TARGETED DOTPLOT
# ============================================================

target_plot_genes <- unique(
  unlist(
    target_panels_present,
    use.names = FALSE
  )
)

target_cells <- WhichCells(
  scrna,
  idents = target_clusters
)

target_object <- subset(
  scrna,
  cells = target_cells
)

Idents(target_object) <-
  target_object$seurat_clusters

p_target_dot <- DotPlot(
  target_object,
  features = target_plot_genes,
  assay = "RNA",
  dot.scale = 6
) +
  RotatedAxis() +
  labs(
    title = "Targeted biological review of clusters 10, 16, 19 and 24",
    subtitle = "Canonical lineage and subtype markers",
    x = NULL,
    y = "Frozen cluster"
  ) +
  theme(
    axis.text.x = element_text(
      size = 7
    ),
    axis.text.y = element_text(
      size = 10
    )
  )

ggsave(
  file.path(
    figures_dir,
    "11D_B3_Targeted_Cluster_DotPlot.png"
  ),
  p_target_dot,
  width = 25,
  height = 7,
  dpi = 400
)


# ============================================================
# 21. TARGETED PANEL HEATMAP
# ============================================================

target_panel_summary$Cluster <- factor(
  target_panel_summary$Cluster,
  levels = target_clusters
)

p_target_heatmap <- ggplot(
  target_panel_summary,
  aes(
    x = Panel,
    y = Cluster,
    fill = Mean_Percent_Expressing
  )
) +
  geom_tile() +
  labs(
    title = "Targeted lineage evidence",
    subtitle = "Mean percentage of cells expressing each marker panel",
    x = NULL,
    y = "Frozen cluster",
    fill = "Mean %\nexpressing"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid = element_blank()
  )

ggsave(
  file.path(
    figures_dir,
    "11D_B3_Targeted_Panel_Heatmap.png"
  ),
  p_target_heatmap,
  width = 14,
  height = 5,
  dpi = 400
)


# ============================================================
# 22. SAMPLE DISTRIBUTION OF TARGET CLUSTERS
#
# Descriptive only.
# NOT biological replicate inference.
# ============================================================

target_sample_distribution <-
  scrna@meta.data %>%
  
  mutate(
    Cluster =
      as.character(
        seurat_clusters
      )
  ) %>%
  
  filter(
    Cluster %in% target_clusters
  ) %>%
  
  count(
    Cluster,
    sample_id,
    name = "Cells"
  ) %>%
  
  group_by(
    Cluster
  ) %>%
  
  mutate(
    Percent_Within_Cluster =
      100 * Cells / sum(Cells)
  ) %>%
  
  ungroup()

write.csv(
  target_sample_distribution,
  file.path(
    results_dir,
    "11D_B3_Target_Cluster_Sample_Distribution.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 23. INTEGRITY CHECK
# ============================================================

integrity_check <- data.frame(
  
  Check = c(
    "Cell count unchanged",
    "Cluster count remains 29",
    "Original broad annotation present",
    "Original refined annotation present"
  ),
  
  Passed = c(
    
    ncol(scrna) == 25186,
    
    length(
      unique(
        as.character(
          scrna$seurat_clusters
        )
      )
    ) == 29,
    
    "broad_cell_type" %in%
      colnames(scrna@meta.data),
    
    "refined_cell_type" %in%
      colnames(scrna@meta.data)
  )
)

write.csv(
  integrity_check,
  file.path(
    results_dir,
    "11D_B3_Integrity_Check.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 24. SESSION INFO
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(
    results_dir,
    "11D_B3_sessionInfo.txt"
  )
)


# ============================================================
# 25. FINAL STATUS
# ============================================================