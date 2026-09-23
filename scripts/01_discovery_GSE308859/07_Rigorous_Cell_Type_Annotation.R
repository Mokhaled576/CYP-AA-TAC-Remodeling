# ============================================================
# GSE308859 TAC scRNA-seq PROJECT
# STEP 07: RIGOROUS CELL-TYPE ANNOTATION AND ATLAS VALIDATION
#
# INPUT:
#   06_Provisional_Major_Cell_Atlas.rds
#
# PURPOSE:
#   1. Preserve frozen resolution-0.4 clustering
#   2. Evaluate broad cardiac lineage programs
#   3. Evaluate refined immune / vascular populations
#   4. Identify cycling and stress-associated states
#   5. Quantify positive and conflicting marker evidence
#   6. Generate cluster-level annotation evidence tables
#   7. Investigate small / ambiguous clusters
#   8. Produce diagnostic UMAPs and DotPlots
#
# IMPORTANT:
#   - NO integration
#   - NO reclustering
#   - NO CYP comparisons
#   - NO condition DE
#   - NO automatic final cell-type assignment
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
  "07_CELL_TYPE_ANNOTATION"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "07_CELL_TYPE_ANNOTATION"
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
# 3. LOAD STEP-06 OBJECT
# ============================================================

input_file <- file.path(
  clean_dir,
  "06_Provisional_Major_Cell_Atlas.rds"
)


if (!file.exists(input_file)) {
  
  stop(
    paste0(
      "Cannot find Step-06 object:\n",
      input_file
    )
  )
}


obj <- readRDS(
  input_file
)


DefaultAssay(obj) <- "RNA"


cat("\n========================================\n")
cat("STEP 07 INPUT OBJECT\n")
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


# ============================================================
# 4. VERIFY FROZEN CLUSTERING
# ============================================================

if (!"atlas_cluster" %in% colnames(obj@meta.data)) {
  
  stop(
    "atlas_cluster metadata is missing."
  )
}


obj$atlas_cluster <- factor(
  as.character(obj$atlas_cluster),
  levels = as.character(
    0:21
  )
)


Idents(obj) <- "atlas_cluster"


cluster_sizes <- table(
  obj$atlas_cluster
)


cat("\n========================================\n")
cat("FROZEN RESOLUTION-0.4 CLUSTERS\n")
cat("========================================\n")

print(
  cluster_sizes
)


if (sum(cluster_sizes) != ncol(obj)) {
  
  stop(
    "Cluster cell numbers do not equal total cell number."
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


# ============================================================
# 6. DEFINE BROAD LINEAGE MARKER PANELS
# ============================================================

# Marker panels are deliberately multi-gene.
# No cluster will be annotated from one marker alone.


broad_marker_panels <- list(
  
  Cardiomyocyte = c(
    "Tnnt2",
    "Tnni3",
    "Actc1",
    "Myh6",
    "Myh7",
    "Myl2",
    "Myl3",
    "Pln",
    "Ryr2"
  ),
  
  Fibroblast = c(
    "Col1a1",
    "Col1a2",
    "Col3a1",
    "Dcn",
    "Lum",
    "Pdgfra",
    "Col6a1",
    "Col6a2"
  ),
  
  Endothelial = c(
    "Pecam1",
    "Cdh5",
    "Kdr",
    "Emcn",
    "Esam",
    "Egfl7",
    "Klf2",
    "Tek"
  ),
  
  Pericyte = c(
    "Rgs5",
    "Pdgfrb",
    "Cspg4",
    "Des",
    "Abcc9",
    "Kcnj8"
  ),
  
  Smooth_Muscle = c(
    "Acta2",
    "Tagln",
    "Myh11",
    "Cnn1",
    "Actg2",
    "Kcnmb1"
  ),
  
  Myeloid = c(
    "Ptprc",
    "Lyz2",
    "Csf1r",
    "Adgre1",
    "Cd68",
    "Fcgr1",
    "Tyrobp",
    "Laptm5"
  ),
  
  T_Cell = c(
    "Ptprc",
    "Cd3d",
    "Cd3e",
    "Cd3g",
    "Trbc1",
    "Trbc2",
    "Lck"
  ),
  
  B_Cell = c(
    "Ptprc",
    "Cd79a",
    "Cd79b",
    "Ms4a1",
    "Cd37",
    "Cd74",
    "H2-Aa"
  ),
  
  NK_Cell = c(
    "Ptprc",
    "Nkg7",
    "Klrd1",
    "Klrk1",
    "Prf1",
    "Gzmb"
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
  
  Mast_Cell = c(
    "Kit",
    "Cpa3",
    "Mcpt4",
    "Ms4a2",
    "Fcer1a"
  ),
  
  Erythroid = c(
    "Hba-a1",
    "Hba-a2",
    "Hbb-bs",
    "Hbb-bt",
    "Alas2",
    "Gypa",
    "Slc4a1"
  )
)


# ============================================================
# 7. DEFINE REFINED MARKER PANELS
# ============================================================

refined_marker_panels <- list(
  
  Resident_Macrophage = c(
    "Timd4",
    "Folr2",
    "Lyve1",
    "Cd163",
    "Vsig4",
    "C1qa",
    "C1qb",
    "C1qc"
  ),
  
  Inflammatory_Monocyte = c(
    "Ly6c2",
    "Ccr2",
    "Plac8",
    "S100a8",
    "S100a9",
    "Lcn2"
  ),
  
  Dendritic_APC = c(
    "Cd74",
    "H2-Aa",
    "H2-Ab1",
    "H2-Eb1",
    "Flt3",
    "Itgax"
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
  
  Plasma_Cell = c(
    "Jchain",
    "Mzb1",
    "Sdc1",
    "Xbp1",
    "Prdm1"
  ),
  
  Capillary_Endothelial = c(
    "Car4",
    "Rgcc",
    "Kdr",
    "Emcn",
    "Gpihbp1",
    "Klf2"
  ),
  
  Arterial_Endothelial = c(
    "Efnb2",
    "Gja5",
    "Sox17",
    "Cxcl12",
    "Fbln5"
  ),
  
  Venous_Endothelial = c(
    "Nr2f2",
    "Ackr1",
    "Vwf",
    "Selp"
  ),
  
  Lymphatic_Endothelial = c(
    "Prox1",
    "Pdpn",
    "Flt4",
    "Lyve1",
    "Ccl21a"
  ),
  
  Activated_Fibroblast = c(
    "Postn",
    "Cthrc1",
    "Thbs4",
    "Comp",
    "Fn1"
  ),
  
  Myofibroblast = c(
    "Acta2",
    "Tagln",
    "Postn",
    "Cnn1",
    "Tpm2"
  )
)


# ============================================================
# 8. DEFINE CELL-STATE PANELS
# ============================================================

state_marker_panels <- list(
  
  Cycling = c(
    "Mki67",
    "Top2a",
    "Pbk",
    "Ccna2",
    "Ccnb1",
    "Ccnb2",
    "Ncapg",
    "Cdca5"
  ),
  
  Immediate_Early = c(
    "Fos",
    "Jun",
    "Junb",
    "Egr1",
    "Atf3"
  ),
  
  Heat_Shock = c(
    "Hspa1a",
    "Hspa1b",
    "Hsp90aa1",
    "Dnajb1"
  ),
  
  Interferon_Response = c(
    "Isg15",
    "Ifit1",
    "Ifit2",
    "Ifit3",
    "Irf7",
    "Oas1a"
  )
)


# ============================================================
# 9. FUNCTION: FILTER PANEL TO PRESENT GENES
# ============================================================

filter_marker_panel <- function(panel_list, object) {
  
  output <- list()
  
  
  for (panel_name in names(panel_list)) {
    
    genes <- panel_list[[panel_name]]
    
    present_genes <- genes[
      genes %in% rownames(object)
    ]
    
    
    if (length(present_genes) > 0) {
      
      output[[panel_name]] <- present_genes
    }
  }
  
  
  return(output)
}


broad_present <- filter_marker_panel(
  broad_marker_panels,
  obj
)


refined_present <- filter_marker_panel(
  refined_marker_panels,
  obj
)


state_present <- filter_marker_panel(
  state_marker_panels,
  obj
)


# ============================================================
# 10. MARKER AVAILABILITY AUDIT
# ============================================================

all_panel_types <- list(
  Broad = broad_marker_panels,
  Refined = refined_marker_panels,
  State = state_marker_panels
)


availability_list <- list()

availability_counter <- 1


for (panel_type in names(all_panel_types)) {
  
  panel_collection <- all_panel_types[[panel_type]]
  
  
  for (panel_name in names(panel_collection)) {
    
    genes <- panel_collection[[panel_name]]
    
    
    for (gene in genes) {
      
      availability_list[[availability_counter]] <- data.frame(
        
        Panel_Type =
          panel_type,
        
        Panel =
          panel_name,
        
        Gene =
          gene,
        
        Present =
          gene %in% rownames(obj),
        
        stringsAsFactors = FALSE
      )
      
      
      availability_counter <-
        availability_counter + 1
    }
  }
}


marker_availability <- do.call(
  rbind,
  availability_list
)


rownames(
  marker_availability
) <- NULL


write.csv(
  marker_availability,
  file.path(
    results_dir,
    "07_Marker_Panel_Availability.csv"
  ),
  row.names = FALSE
)


cat("\n========================================\n")
cat("MARKER PANEL AVAILABILITY\n")
cat("========================================\n")


availability_summary <- aggregate(
  Present ~ Panel_Type + Panel,
  data = marker_availability,
  FUN = function(x) {
    paste0(
      sum(x),
      "/",
      length(x)
    )
  }
)


print(
  availability_summary
)


# ============================================================
# 11. UNIQUE GENES FOR DOTPLOTS
# ============================================================

# DotPlot requires unique feature names.
# Shared genes such as Ptprc must appear only once.


broad_dot_genes <- unique(
  unlist(
    broad_present,
    use.names = FALSE
  )
)


refined_dot_genes <- unique(
  unlist(
    refined_present,
    use.names = FALSE
  )
)


state_dot_genes <- unique(
  unlist(
    state_present,
    use.names = FALSE
  )
)


# ============================================================
# 12. BROAD-LINEAGE DOTPLOT
# ============================================================

p_broad_dot <- DotPlot(
  object = obj,
  features = broad_dot_genes,
  group.by = "atlas_cluster",
  assay = "RNA",
  dot.scale = 6
) +
  RotatedAxis() +
  labs(
    title =
      "Broad cardiac lineage evidence",
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
    "07_Broad_Lineage_DotPlot.png"
  ),
  plot = p_broad_dot,
  width = 22,
  height = 10,
  dpi = 600
)


# ============================================================
# 13. REFINED-LINEAGE DOTPLOT
# ============================================================

p_refined_dot <- DotPlot(
  object = obj,
  features = refined_dot_genes,
  group.by = "atlas_cluster",
  assay = "RNA",
  dot.scale = 6
) +
  RotatedAxis() +
  labs(
    title =
      "Refined cell-type and subtype evidence",
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
    "07_Refined_Lineage_DotPlot.png"
  ),
  plot = p_refined_dot,
  width = 22,
  height = 10,
  dpi = 600
)


# ============================================================
# 14. CELL-STATE DOTPLOT
# ============================================================

p_state_dot <- DotPlot(
  object = obj,
  features = state_dot_genes,
  group.by = "atlas_cluster",
  assay = "RNA",
  dot.scale = 6
) +
  RotatedAxis() +
  labs(
    title =
      "Cell-state programs by cluster",
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
      size = 8,
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )


ggsave(
  filename = file.path(
    figures_dir,
    "07_Cell_State_DotPlot.png"
  ),
  plot = p_state_dot,
  width = 16,
  height = 10,
  dpi = 600
)


# ============================================================
# 15. FUNCTION: CALCULATE PANEL DETECTION BY CLUSTER
# ============================================================

rna_counts <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "counts"
)


calculate_panel_detection <- function(
    object,
    count_matrix,
    panel_list,
    panel_type
) {
  
  result_list <- list()
  
  counter <- 1
  
  
  cluster_ids <- levels(
    object$atlas_cluster
  )
  
  
  for (cluster_id in cluster_ids) {
    
    cluster_cells <- colnames(object)[
      object$atlas_cluster == cluster_id
    ]
    
    
    for (panel_name in names(panel_list)) {
      
      genes <- panel_list[[panel_name]]
      
      genes <- genes[
        genes %in% rownames(count_matrix)
      ]
      
      
      if (length(genes) == 0) {
        
        next
      }
      
      
      sub_matrix <- count_matrix[
        genes,
        cluster_cells,
        drop = FALSE
      ]
      
      
      gene_detection <- Matrix::rowMeans(
        sub_matrix > 0
      )
      
      
      cells_with_any_marker <- Matrix::colSums(
        sub_matrix > 0
      ) > 0
      
      
      cells_with_two_or_more <- Matrix::colSums(
        sub_matrix > 0
      ) >= 2
      
      
      result_list[[counter]] <- data.frame(
        
        Cluster =
          cluster_id,
        
        Panel_Type =
          panel_type,
        
        Panel =
          panel_name,
        
        Genes_Available =
          length(genes),
        
        Mean_Gene_Detection_Percent =
          100 *
          mean(
            gene_detection
          ),
        
        Cells_With_Any_Marker_Percent =
          100 *
          mean(
            cells_with_any_marker
          ),
        
        Cells_With_2plus_Markers_Percent =
          100 *
          mean(
            cells_with_two_or_more
          ),
        
        stringsAsFactors = FALSE
      )
      
      
      counter <- counter + 1
    }
  }
  
  
  output <- do.call(
    rbind,
    result_list
  )
  
  
  rownames(output) <- NULL
  
  return(output)
}


# ============================================================
# 16. CALCULATE PANEL DETECTION
# ============================================================

broad_detection <- calculate_panel_detection(
  object = obj,
  count_matrix = rna_counts,
  panel_list = broad_present,
  panel_type = "Broad"
)


refined_detection <- calculate_panel_detection(
  object = obj,
  count_matrix = rna_counts,
  panel_list = refined_present,
  panel_type = "Refined"
)


state_detection <- calculate_panel_detection(
  object = obj,
  count_matrix = rna_counts,
  panel_list = state_present,
  panel_type = "State"
)


all_panel_detection <- rbind(
  broad_detection,
  refined_detection,
  state_detection
)


write.csv(
  all_panel_detection,
  file.path(
    results_dir,
    "07_Panel_Detection_By_Cluster.csv"
  ),