# ============================================================
# GSE308859 TAC scRNA-seq PROJECT
# STEP 08: CYP–ARACHIDONIC ACID / EICOSANOID LANDSCAPE
#
# INPUT:
#   07C_FINAL_FROZEN_ANNOTATED_ATLAS.rds
#
# PURPOSE:
#   1. Identify every CYP gene detected in the frozen atlas
#   2. Quantify CYP expression by cell type and cluster
#   3. Quantify CYP expression by condition/timepoint
#   4. Map CYP-AA epoxygenase/hydroxylase machinery
#   5. Map soluble epoxide hydrolase (Ephx2)
#   6. Map phospholipase / AA-release machinery
#   7. Map COX/prostaglandin machinery
#   8. Map LOX/leukotriene machinery
#   9. Map selected eicosanoid receptors
#  10. Separate robust signals from sparse observations
#  11. Generate publication-oriented discovery figures
#
# IMPORTANT:
#   - Frozen atlas is NOT modified
#   - No reclustering
#   - No reintegration
#   - No cell-level hypothesis tests across conditions
#   - Condition comparisons are DESCRIPTIVE
#   - Cells are NOT biological replicates
# ============================================================


# ============================================================
# 0. CLEAN ENVIRONMENT
# ============================================================

rm(list = ls())
gc()

set.seed(20260919)

options(stringsAsFactors = FALSE)


# ============================================================
# 1. PROJECT DIRECTORIES
# ============================================================

project_dir <- "D:/Master/ScRNA seq/TAC"

clean_dir <- file.path(
  project_dir,
  "CLEAN DATA"
)

results_dir <- file.path(
  project_dir,
  "RESULTS",
  "08_CYP_AA_EICOSANOID"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "08_CYP_AA_EICOSANOID"
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
  "Matrix",
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
library(Matrix)
library(ggplot2)
library(patchwork)


# ============================================================
# 3. LOAD FROZEN ATLAS
# ============================================================

input_file <- file.path(
  clean_dir,
  "07C_FINAL_FROZEN_ANNOTATED_ATLAS.rds"
)


if (!file.exists(input_file)) {
  
  stop(
    paste0(
      "Cannot find frozen atlas:\n",
      input_file
    )
  )
}


obj <- readRDS(
  input_file
)


DefaultAssay(obj) <- "RNA"


cat("\n========================================\n")
cat("STEP 08\n")
cat("CYP-AA / EICOSANOID LANDSCAPE\n")
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
# 4. VERIFY FROZEN ATLAS
# ============================================================

required_metadata <- c(
  "atlas_cluster",
  "sample_id",
  "broad_cell_type",
  "refined_cell_type",
  "atlas_frozen"
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


if (!all(obj$atlas_frozen)) {
  
  stop(
    "Atlas is not marked as frozen."
  )
}


if (ncol(obj) != 25186) {
  
  stop(
    paste0(
      "Expected 25,186 frozen atlas cells, but found ",
      ncol(obj),
      "."
    )
  )
}


cat(
  "Frozen atlas verification: PASSED\n"
)


# ============================================================
# 5. DEFINE ORDERS
# ============================================================

sample_order <- c(
  "Sham",
  "TAC_2W",
  "TAC_4W",
  "TAC_6W"
)


refined_order <- c(
  "Cardiomyocytes",
  "Fibroblasts",
  "Endothelial cells",
  "Pericytes",
  "Macrophages",
  "Cycling myeloid cells",
  "Monocytes",
  "Activated dendritic cells",
  "T cells",
  "Immature T cells",
  "B cells",
  "Immature B cells",
  "Erythroid cells"
)


cluster_order <- as.character(
  0:21
)


obj$sample_id <- factor(
  as.character(obj$sample_id),
  levels = sample_order
)


obj$refined_cell_type <- factor(
  as.character(obj$refined_cell_type),
  levels = refined_order
)


obj$atlas_cluster <- factor(
  as.character(obj$atlas_cluster),
  levels = cluster_order
)


# ============================================================
# 6. GET RAW COUNTS
# ============================================================

rna_counts <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "counts"
)


cat(
  "\nRaw count matrix:",
  nrow(rna_counts),
  "genes x",
  ncol(rna_counts),
  "cells\n"
)


# ============================================================
# 7. IDENTIFY ALL CYP GENES
# ============================================================
#
# Strict CYP naming rule:
#   Cyp followed by a number.
#
# This avoids false matches such as unrelated "Cypt" genes.
# ============================================================

all_genes <- rownames(
  rna_counts
)


all_cyp_genes <- grep(
  "^Cyp[0-9]",
  all_genes,
  value = TRUE
)


all_cyp_genes <- sort(
  unique(all_cyp_genes)
)


cat(
  "\nTotal CYP genes present in matrix:",
  length(all_cyp_genes),
  "\n"
)


# ============================================================
# 8. GLOBAL RAW EXPRESSION FOR ALL CYP GENES
# ============================================================

global_cyp_list <- vector(
  "list",
  length(all_cyp_genes)
)


for (i in seq_along(all_cyp_genes)) {
  
  gene <- all_cyp_genes[i]
  
  values <- rna_counts[
    gene,
    ,
    drop = TRUE
  ]
  
  
  global_cyp_list[[i]] <- data.frame(
    
    Gene = gene,
    
    Positive_Cells =
      sum(values > 0),
    
    Percent_Positive =
      100 * mean(values > 0),
    
    Mean_Raw_Count =
      mean(values),
    
    Total_UMI =
      sum(values),
    
    Max_Raw_Count =
      max(values),
    
    stringsAsFactors = FALSE
  )
}


global_cyp <- do.call(
  rbind,
  global_cyp_list
)


global_cyp <- global_cyp[
  order(
    -global_cyp$Positive_Cells,
    -global_cyp$Total_UMI
  ),
  ,
  drop = FALSE
]


rownames(global_cyp) <- NULL


# ============================================================
# 9. CYP DETECTION CLASSIFICATION
# ============================================================
#
# These categories are descriptive signal-support categories,
# NOT statistical significance categories.
#
# Robust:
#   >= 1% of atlas OR >= 250 positive cells
#
# Moderate:
#   >= 0.1% OR >= 25 positive cells
#
# Sparse:
#   > 0 but below Moderate
#
# Undetected:
#   zero counts
#
# ============================================================

global_cyp$Detection_Class <- "Undetected"


global_cyp$Detection_Class[
  global_cyp$Positive_Cells > 0
] <- "Sparse"


global_cyp$Detection_Class[
  global_cyp$Positive_Cells >= 25 |
    global_cyp$Percent_Positive >= 0.1
] <- "Moderate"


global_cyp$Detection_Class[
  global_cyp$Positive_Cells >= 250 |
    global_cyp$Percent_Positive >= 1
] <- "Robust"


global_cyp$Detection_Class <- factor(
  global_cyp$Detection_Class,
  levels = c(
    "Robust",
    "Moderate",
    "Sparse",
    "Undetected"
  )
)


write.csv(
  global_cyp,
  file.path(
    results_dir,
    "08_All_CYP_Global_Expression.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. PRINT TOP CYP GENES
# ============================================================

cat("\nTOP CYP GENES BY POSITIVE CELLS:\n")

print(
  head(
    global_cyp,
    30
  )
)


cat("\nCYP DETECTION CLASSES:\n")

print(
  table(
    global_cyp$Detection_Class
  )
)


# ============================================================
# 11. FUNCTION: RAW EXPRESSION BY GROUP
# ============================================================

summarize_gene_by_group <- function(
    counts_matrix,
    genes,
    groups,
    group_name
) {
  
  groups <- as.character(groups)
  
  group_levels <- unique(groups)
  
  result_list <- list()
  
  counter <- 1
  
  
  for (group in group_levels) {
    
    cells <- which(
      groups == group
    )
    
    
    if (length(cells) == 0) {
      
      next
    }
    
    
    for (gene in genes) {
      
      values <- counts_matrix[
        gene,
        cells,
        drop = TRUE
      ]
      
      
      result_list[[counter]] <- data.frame(
        
        Group =
          group,
        
        Gene =
          gene,
        
        Cells =
          length(cells),
        
        Positive_Cells =
          sum(values > 0),
        
        Percent_Positive =
          100 * mean(values > 0),
        
        Mean_Raw_Count =
          mean(values),
        
        Mean_Count_Positive_Cells =
          ifelse(
            sum(values > 0) > 0,
            mean(values[values > 0]),
            0
          ),
        
        Total_UMI =
          sum(values),
        
        stringsAsFactors = FALSE
      )
      
      
      counter <- counter + 1
    }
  }
  
  
  result <- do.call(
    rbind,
    result_list
  )
  
  
  colnames(result)[
    colnames(result) == "Group"
  ] <- group_name
  
  
  rownames(result) <- NULL
  
  
  return(result)
}


# ============================================================
# 12. ALL CYP GENES BY REFINED CELL TYPE
# ============================================================

cyp_by_celltype <- summarize_gene_by_group(
  counts_matrix = rna_counts,
  genes = all_cyp_genes,
  groups = obj$refined_cell_type,
  group_name = "Refined_Cell_Type"
)


write.csv(
  cyp_by_celltype,
  file.path(
    results_dir,
    "08_All_CYP_By_Refined_Cell_Type.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. ALL CYP GENES BY CLUSTER
# ============================================================

cyp_by_cluster <- summarize_gene_by_group(
  counts_matrix = rna_counts,
  genes = all_cyp_genes,
  groups = obj$atlas_cluster,
  group_name = "Atlas_Cluster"
)


write.csv(
  cyp_by_cluster,
  file.path(
    results_dir,
    "08_All_CYP_By_Cluster.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. ALL CYP GENES BY SAMPLE
# ============================================================

cyp_by_sample <- summarize_gene_by_group(
  counts_matrix = rna_counts,
  genes = all_cyp_genes,
  groups = obj$sample_id,
  group_name = "Sample"
)


write.csv(
  cyp_by_sample,
  file.path(
    results_dir,
    "08_All_CYP_By_Sample.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. DEFINE AA / EICOSANOID GENE SETS
# ============================================================
#
# This is deliberately broad.
#
# Genes absent from the dataset will be recorded separately.
# We do NOT silently discard them.
# ============================================================

gene_sets <- list(
  
  CYP_AA_Epoxygenase_Hydroxylase = c(
    "Cyp2c29",
    "Cyp2c37",
    "Cyp2c38",
    "Cyp2c39",
    "Cyp2c40",
    "Cyp2c44",
    "Cyp2c50",
    "Cyp2c54",
    "Cyp2c55",
    "Cyp2c65",
    "Cyp2c66",
    "Cyp2c68",
    "Cyp2c69",
    "Cyp2c70",
    "Cyp2j5",
    "Cyp2j6",
    "Cyp2j8",
    "Cyp2j9",
    "Cyp2j11",
    "Cyp2j12",
    "Cyp4a10",
    "Cyp4a12a",
    "Cyp4a12b",
    "Cyp4a14",
    "Cyp4a31",
    "Cyp4a32",
    "Cyp4f13",
    "Cyp4f14",
    "Cyp4f15",
    "Cyp4f16",
    "Cyp4f17",
    "Cyp4f18",
    "Cyp4f39",
    "Cyp4v3"
  ),
  
  Epoxide_Hydrolase = c(
    "Ephx1",
    "Ephx2"
  ),
  
  AA_Release_PLA2 = c(
    "Pla2g4a",
    "Pla2g4b",
    "Pla2g4c",
    "Pla2g4d",
    "Pla2g4e",
    "Pla2g4f",
    "Pla2g6",
    "Pla2g7",
    "Pla2g2a",
    "Pla2g2d",
    "Pla2g2e",
    "Pla2g5",
    "Pla2g12a"
  ),
  
  COX_Prostaglandin = c(
    "Ptgs1",
    "Ptgs2",
    "Ptges",
    "Ptges2",
    "Ptges3",
    "Ptgis",
    "Tbxas1",
    "Hpgds",
    "Ptgds",
    "Hpgd",
    "Ptgr1"
  ),
  
  LOX_Leukotriene = c(
    "Alox5",
    "Alox5ap",
    "Alox12",
    "Alox12b",
    "Alox15",
    "Alox15b",
    "Aloxe3",
    "Lta4h",
    "Ltc4s",
    "Ggt1",
    "Ggt5"
  ),
  
  Eicosanoid_Receptors = c(
    "Ptger1",
    "Ptger2",
    "Ptger3",
    "Ptger4",
    "Ptgdr",
    "Ptgdr2",
    "Ptgfr",
    "Ptgir",
    "Tbxa2r",
    "Ltb4r1",
    "Ltb4r2",
    "Cysltr1",
    "Cysltr2"
  ),
  
  AA_Lipid_Context = c(
    "Acsl4",
    "Acsl1",
    "Acsl3",
    "Acsl5",
    "Lpcat3",
    "Lpcat1",
    "Lpcat2",
    "Lpcat4",
    "Fabp3",
    "Fabp4",
    "Fabp5",
    "Cd36"
  )
)


# ============================================================
# 16. ADD LAB-RELEVANT CYP GENES
# ============================================================
#
# Cyp1b1 is not restricted to the classical CYP-AA panel,
# but is strongly detectable and biologically relevant to
# the broader project.
# ============================================================

lab_priority_genes <- c(
  "Cyp1a1",
  "Cyp1a2",
  "Cyp1b1",
  "Cyp2j6",
  "Cyp2j9",
  "Cyp4f13",
  "Cyp4f16",
  "Cyp4f17",
  "Cyp4f18",
  "Ephx2"
)


# ============================================================
# 17. BUILD GENE-SET AVAILABILITY TABLE
# ============================================================

availability_list <- list()

counter <- 1


for (set_name in names(gene_sets)) {
  
  genes <- unique(
    gene_sets[[set_name]]
  )
  
  
  for (gene in genes) {
    
    availability_list[[counter]] <- data.frame(
      
      Gene_Set =
        set_name,
      
      Gene =
        gene,
      
      Present_In_Matrix =
        gene %in% all_genes,
      
      stringsAsFactors = FALSE
    )
    
    
    counter <- counter + 1
  }
}


gene_set_availability <- do.call(
  rbind,
  availability_list
)


rownames(
  gene_set_availability
) <- NULL


write.csv(
  gene_set_availability,
  file.path(
    results_dir,
    "08_AA_Eicosanoid_Gene_Set_Availability.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. CREATE MASTER AA/EICOSANOID GENE LIST
# ============================================================

aa_genes_requested <- unique(
  c(
    unlist(
      gene_sets,
      use.names = FALSE
    ),
    lab_priority_genes
  )
)


aa_genes_present <- aa_genes_requested[
  aa_genes_requested %in% all_genes
]


aa_genes_missing <- setdiff(
  aa_genes_requested,
  all_genes
)


write.csv(
  data.frame(
    Gene = aa_genes_present
  ),
  file.path(
    results_dir,
    "08_AA_Eicosanoid_Genes_Present.csv"
  ),
  row.names = FALSE
)


write.csv(
  data.frame(
    Gene = aa_genes_missing
  ),
  file.path(
    results_dir,
    "08_AA_Eicosanoid_Genes_Missing.csv"
  ),
  row.names = FALSE
)


cat(
  "\nAA/eicosanoid genes requested:",
  length(aa_genes_requested),
  "\n"
)


cat(
  "Present:",
  length(aa_genes_present),
  "\n"
)


cat(
  "Missing:",
  length(aa_genes_missing),
  "\n"
)


# ============================================================
# 19. GLOBAL AA/EICOSANOID EXPRESSION
# ============================================================

aa_global_list <- vector(
  "list",
  length(aa_genes_present)
)


for (i in seq_along(aa_genes_present)) {
  
  gene <- aa_genes_present[i]
  
  values <- rna_counts[
    gene,
    ,
    drop = TRUE
  ]
  
  
  aa_global_list[[i]] <- data.frame(
    
    Gene =
      gene,
    
    Positive_Cells =
      sum(values > 0),
    
    Percent_Positive =
      100 * mean(values > 0),
    
    Mean_Raw_Count =
      mean(values),
    
    Total_UMI =
      sum(values),
    
    stringsAsFactors = FALSE
  )
}


aa_global <- do.call(
  rbind,
  aa_global_list
)


aa_global <- aa_global[
  order(
    -aa_global$Positive_Cells,
    -aa_global$Total_UMI
  ),
  ,
  drop = FALSE
]


rownames(aa_global) <- NULL


write.csv(
  aa_global,
  file.path(
    results_dir,
    "08_AA_Eicosanoid_Global_Expression.csv"
  ),
  row.names = FALSE