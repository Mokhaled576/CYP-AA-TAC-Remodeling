# ============================================================
# GSE308859 TAC scRNA-seq PROJECT
# STEP 07B: TARGETED ANNOTATION VALIDATION
#
# TARGET CLUSTERS:
#   11, 12, 16, 20, 21
#
# PURPOSE:
#   Cluster 11:
#     Determine lineage beneath strong cycling program
#
#   Cluster 12:
#     Resolve lineage / APC-myeloid identity
#
#   Cluster 16:
#     Resolve myeloid / dendritic / APC identity
#
#   Clusters 20 and 21:
#     Validate unusual immature lymphoid signatures and
#     investigate QC / doublet-like characteristics
#
# IMPORTANT:
#   - NO reclustering
#   - NO integration
#   - NO change to resolution 0.4
#   - NO CYP analysis
#   - NO automatic final annotation
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
  "07B_TARGETED_ANNOTATION_VALIDATION"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "07B_TARGETED_ANNOTATION_VALIDATION"
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
# 3. LOAD STEP-07 DIAGNOSTIC OBJECT
# ============================================================

input_file <- file.path(
  clean_dir,
  "07_Annotation_Diagnostic_Atlas.rds"
)


if (!file.exists(input_file)) {
  
  stop(
    paste0(
      "Cannot find Step-07 object:\n",
      input_file
    )
  )
}


obj <- readRDS(
  input_file
)


DefaultAssay(obj) <- "RNA"


# ============================================================
# 4. VERIFY CLUSTERS
# ============================================================

if (!"atlas_cluster" %in% colnames(obj@meta.data)) {
  
  stop(
    "atlas_cluster metadata is missing."
  )
}


obj$atlas_cluster <- factor(
  as.character(obj$atlas_cluster),
  levels = as.character(0:21)
)


Idents(obj) <- "atlas_cluster"


target_clusters <- c(
  "11",
  "12",
  "16",
  "20",
  "21"
)


if (!all(
  target_clusters %in%
  levels(obj$atlas_cluster)
)) {
  
  stop(
    "One or more target clusters are missing."
  )
}


cat("\n========================================\n")
cat("STEP 07B INPUT\n")
cat("========================================\n")

cat(
  "Total cells:",
  ncol(obj),
  "\n"
)

cat(
  "Total genes:",
  nrow(obj),
  "\n"
)

cat(
  "Target clusters:",
  paste(
    target_clusters,
    collapse = ", "
  ),
  "\n"
)


# ============================================================
# 5. TARGET CLUSTER CELL NUMBERS
# ============================================================

target_cluster_sizes <- data.frame(
  
  Cluster =
    target_clusters,
  
  Cells =
    sapply(
      target_clusters,
      function(cl) {
        
        sum(
          as.character(
            obj$atlas_cluster
          ) == cl
        )
      }
    ),
  
  stringsAsFactors = FALSE
)


write.csv(
  target_cluster_sizes,
  file.path(
    results_dir,
    "07B_Target_Cluster_Cell_Numbers.csv"
  ),
  row.names = FALSE
)


cat("\nTARGET CLUSTER CELL NUMBERS:\n")

print(
  target_cluster_sizes
)


# ============================================================
# 6. DEFINE TARGETED LINEAGE PANELS
# ============================================================

marker_panels <- list(
  
  Cardiomyocyte = c(
    "Tnnt2",
    "Tnni3",
    "Actc1",
    "Myh6",
    "Myh7",
    "Myl2",
    "Ryr2"
  ),
  
  Fibroblast = c(
    "Col1a1",
    "Col1a2",
    "Col3a1",
    "Dcn",
    "Lum",
    "Pdgfra"
  ),
  
  Endothelial = c(
    "Pecam1",
    "Cdh5",
    "Kdr",
    "Emcn",
    "Esam",
    "Egfl7"
  ),
  
  Pericyte = c(
    "Rgs5",
    "Pdgfrb",
    "Cspg4",
    "Abcc9",
    "Kcnj8"
  ),
  
  Smooth_Muscle = c(
    "Acta2",
    "Tagln",
    "Myh11",
    "Cnn1"
  ),
  
  Pan_Leukocyte = c(
    "Ptprc",
    "Tyrobp",
    "Laptm5"
  ),
  
  Macrophage = c(
    "Lyz2",
    "Csf1r",
    "Adgre1",
    "Cd68",
    "Fcgr1",
    "C1qa",
    "C1qb",
    "C1qc"
  ),
  
  Monocyte = c(
    "Ly6c2",
    "Ccr2",
    "Plac8",
    "Ctss",
    "Lyz2"
  ),
  
  Dendritic_APC = c(
    "Cd74",
    "H2-Aa",
    "H2-Ab1",
    "H2-Eb1",
    "Itgax",
    "Flt3"
  ),
  
  cDC1 = c(
    "Xcr1",
    "Clec9a",
    "Batf3",
    "Irf8"
  ),
  
  cDC2 = c(
    "Cd209a",
    "Sirpa",
    "Clec10a",
    "Irf4"
  ),
  
  Activated_DC = c(
    "Ccl17",
    "Ccl22",
    "Cd80",
    "Cd86"
  ),
  
  Neutrophil = c(
    "S100a8",
    "S100a9",
    "Ly6g",
    "Retnlg",
    "Camp",
    "Ngp",
    "Ltf"
  ),
  
  T_Cell = c(
    "Cd3d",
    "Cd3e",
    "Cd3g",
    "Trbc1",
    "Trbc2",
    "Lck"
  ),
  
  NK_Cell = c(
    "Nkg7",
    "Klrd1",
    "Klrk1",
    "Prf1",
    "Gzmb"
  ),
  
  Mature_B = c(
    "Cd79a",
    "Cd79b",
    "Ms4a1",
    "Cd37",
    "Cd74",
    "Cd19",
    "Pax5"
  ),
  
  Plasma_Cell = c(
    "Jchain",
    "Mzb1",
    "Sdc1",
    "Xbp1",
    "Prdm1"
  ),
  
  Immature_Lymphoid = c(
    "Rag1",
    "Rag2",
    "Dntt"
  ),
  
  Immature_T = c(
    "Rag1",
    "Dntt",
    "Ccr9",
    "Notch1",
    "Tcf7"
  ),
  
  Immature_B = c(
    "Rag1",
    "Rag2",
    "Dntt",
    "Igll1",
    "Vpreb1",
    "Vpreb3",
    "Myb"
  ),
  
  Megakaryocyte_Platelet = c(
    "Ppbp",
    "Pf4",
    "Tubb1",
    "Gp5",
    "Gp9",
    "Gp1ba",
    "Gp6",
    "Clec1b",
    "Treml1"
  ),
  
  Erythroid = c(
    "Hba-a1",
    "Hba-a2",
    "Hbb-bs",
    "Hbb-bt",
    "Alas2",
    "Gypa",
    "Slc4a1"
  ),
  
  Cycling = c(
    "Mki67",
    "Top2a",
    "Pbk",
    "Ccna2",
    "Ccnb1",
    "Ccnb2",
    "Ncapg",
    "Cdca5"
  )
)


# ============================================================
# 7. FILTER PANELS TO PRESENT GENES
# ============================================================

present_panels <- list()


for (panel_name in names(marker_panels)) {
  
  genes <- marker_panels[[panel_name]]
  
  present_genes <- genes[
    genes %in% rownames(obj)
  ]
  
  
  if (length(present_genes) > 0) {
    
    present_panels[[panel_name]] <-
      present_genes
  }
}


marker_availability <- data.frame()


for (panel_name in names(marker_panels)) {
  
  genes <- marker_panels[[panel_name]]
  
  
  temp <- data.frame(
    
    Panel =
      panel_name,
    
    Gene =
      genes,
    
    Present =
      genes %in% rownames(obj),
    
    stringsAsFactors = FALSE
  )
  
  
  marker_availability <- rbind(
    marker_availability,
    temp
  )
}


write.csv(
  marker_availability,
  file.path(
    results_dir,
    "07B_Targeted_Marker_Availability.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 8. CREATE TARGET-CLUSTER OBJECT
# ============================================================

target_cells <- colnames(obj)[
  as.character(
    obj$atlas_cluster
  ) %in%
    target_clusters
]


target_obj <- subset(
  obj,
  cells = target_cells
)


target_obj$atlas_cluster <- factor(
  as.character(
    target_obj$atlas_cluster
  ),
  levels = target_clusters
)


Idents(target_obj) <- "atlas_cluster"


cat(
  "\nTarget object cells:",
  ncol(target_obj),
  "\n"
)


# ============================================================
# 9. TARGETED DOTPLOT
# ============================================================

dot_genes <- unique(
  unlist(
    present_panels,
    use.names = FALSE
  )
)


p_target_dot <- DotPlot(
  object = target_obj,
  features = dot_genes,
  group.by = "atlas_cluster",
  assay = "RNA",
  dot.scale = 7
) +
  RotatedAxis() +
  labs(
    title =
      "Targeted lineage validation of ambiguous clusters",
    x =
      "Marker",
    y =
      "Cluster"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    axis.text.x = element_text(
      size = 7,
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )


ggsave(
  filename = file.path(
    figures_dir,
    "07B_Targeted_Lineage_DotPlot.png"
  ),
  plot = p_target_dot,
  width = 26,
  height = 7,
  dpi = 600,
  limitsize = FALSE
)


# ============================================================
# 10. RAW DETECTION MATRIX
# ============================================================

rna_counts <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "counts"
)


# ============================================================
# 11. PANEL DETECTION BY TARGET CLUSTER
# ============================================================

panel_detection_list <- list()

counter <- 1


for (cluster_id in target_clusters) {
  
  cluster_cells <- colnames(obj)[
    as.character(
      obj$atlas_cluster
    ) == cluster_id
  ]
  
  
  for (panel_name in names(present_panels)) {
    
    genes <- present_panels[[panel_name]]
    
    
    sub_matrix <- rna_counts[
      genes,
      cluster_cells,
      drop = FALSE
    ]
    
    
    detected_per_cell <- Matrix::colSums(
      sub_matrix > 0
    )
    
    
    panel_detection_list[[counter]] <-
      data.frame(
        
        Cluster =
          cluster_id,
        
        Panel =
          panel_name,
        
        Cells =
          length(cluster_cells),
        
        Genes_Available =
          length(genes),
        
        Mean_Number_Markers_Detected =
          mean(
            detected_per_cell
          ),
        
        Cells_With_Any_Marker_Percent =
          100 *
          mean(
            detected_per_cell >= 1
          ),
        
        Cells_With_2plus_Markers_Percent =
          100 *
          mean(
            detected_per_cell >= 2
          ),
        
        Cells_With_3plus_Markers_Percent =
          100 *
          mean(
            detected_per_cell >= 3
          ),
        
        stringsAsFactors = FALSE
      )
    
    
    counter <- counter + 1
  }
}


panel_detection <- do.call(
  rbind,
  panel_detection_list
)


rownames(
  panel_detection
) <- NULL


write.csv(
  panel_detection,
  file.path(
    results_dir,
    "07B_Target_Cluster_Panel_Detection.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. GENE-LEVEL DETECTION BY TARGET CLUSTER
# ============================================================

target_genes <- unique(
  unlist(
    present_panels,
    use.names = FALSE
  )
)


gene_detection_list <- list()

counter <- 1


for (cluster_id in target_clusters) {
  
  cluster_cells <- colnames(obj)[
    as.character(
      obj$atlas_cluster
    ) == cluster_id
  ]
  
  
  for (gene in target_genes) {
    
    values <- rna_counts[
      gene,
      cluster_cells,
      drop = TRUE
    ]
    
    
    gene_detection_list[[counter]] <-
      data.frame(
        
        Cluster =
          cluster_id,
        
        Gene =
          gene,
        
        Cells =
          length(cluster_cells),
        
        Positive_Cells =
          sum(
            values > 0
          ),
        
        Percent_Positive =
          100 *
          mean(
            values > 0
          ),
        
        Mean_Raw_Count =
          mean(
            values
          ),
        
        Total_UMI =
          sum(
            values
          ),
        
        stringsAsFactors = FALSE
      )
    
    
    counter <- counter + 1
  }
}


gene_detection <- do.call(
  rbind,
  gene_detection_list
)


rownames(
  gene_detection
) <- NULL


write.csv(
  gene_detection,
  file.path(
    results_dir,
    "07B_Target_Cluster_Gene_Detection.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. LOAD STEP-07 TOP MARKERS
# ============================================================

marker_file <- file.path(
  project_dir,
  "RESULTS",
  "07_CELL_TYPE_ANNOTATION",
  "07_Top30_Markers_Per_Cluster.csv"
)


if (!file.exists(marker_file)) {
  
  stop(
    paste0(
      "Cannot find Step-07 top-marker file:\n",
      marker_file
    )
  )
}


top_markers <- read.csv(
  marker_file,
  check.names = FALSE
)


target_top_markers <- top_markers[
  as.character(
    top_markers$cluster
  ) %in%
    target_clusters,
  ,
  drop = FALSE
]


write.csv(
  target_top_markers,
  file.path(
    results_dir,
    "07B_Target_Cluster_Top30_Markers.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. CREATE TOP-20 MARKER SUMMARY
# ============================================================

marker_summary_list <- list()


for (cluster_id in target_clusters) {
  
  temp <- target_top_markers[
    as.character(
      target_top_markers$cluster
    ) == cluster_id,
    ,
    drop = FALSE
  ]
  
  
  genes <- head(
    temp$gene,
    20
  )
  
  
  marker_summary_list[[cluster_id]] <-
    data.frame(
      
      Cluster =
        cluster_id,
      
      Top20_Markers =
        paste(
          genes,
          collapse = "; "
        ),
      
      stringsAsFactors = FALSE
    )
}


marker_summary <- do.call(
  rbind,
  marker_summary_list
)


rownames(
  marker_summary
) <- NULL


# ============================================================
# 15. QC METADATA AUDIT
# ============================================================

qc_columns <- c(
  "nFeature_RNA",
  "nCount_RNA",
  "percent.mt",
  "percent.ribo",
  "percent.hb",
  "scDblFinder.score",
  "scDblFinder.class"
)


available_qc_columns <- qc_columns[
  qc_columns %in%
    colnames(obj@meta.data)
]


cat("\nAvailable QC columns:\n")

print(
  available_qc_columns
)


# ============================================================
# 16. TARGET CLUSTER QC SUMMARY
# ============================================================

qc_summary_list <- list()


for (cluster_id in target_clusters) {
  
  cluster_cells <- colnames(obj)[
    as.character(
      obj$atlas_cluster
    ) == cluster_id
  ]
  
  
  md <- obj@meta.data[
    cluster_cells,
    ,
    drop = FALSE
  ]
  
  
  temp <- data.frame(
    
    Cluster =
      cluster_id,
    
    Cells =
      nrow(md),
    
    Median_nFeature_RNA =
      median(
        md$nFeature_RNA,
        na.rm = TRUE
      ),
    
    Median_nCount_RNA =