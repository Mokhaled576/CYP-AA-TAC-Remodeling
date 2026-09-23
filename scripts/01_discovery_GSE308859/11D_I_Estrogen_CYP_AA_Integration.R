# ============================================================
# STEP 11D-I
# ESTROGEN SIGNALING x CYP/ARACHIDONIC-ACID INTEGRATION
#
# PURPOSE
# ------------------------------------------------------------
# Determine whether estrogen-related transcriptional biology
# is detectable in the TAC dataset and, if supported, examine
# its relationship to CYP/AA remodeling.
#
# ANALYTICAL LOGIC
# ------------------------------------------------------------
# PART A: FEASIBILITY GATE
#   1. Esr1 detection
#   2. Esr2 detection
#   3. Gper1 detection
#   4. receptor localization by cell type
#   5. receptor trajectories across TAC
#   6. spatial receptor detectability
#
# PART B: ESTROGEN-RELATED TRANSCRIPTIONAL PROGRAMS
#   1. estrogen-response program
#   2. ER-associated transcriptional program
#
# PART C: CYP/AA INTEGRATION
#   1. priority CYP genes x estrogen programs
#   2. AA metabolic arms x estrogen programs
#   3. cell-type-specific relationships
#   4. spatial relationships
#   5. program-adjusted CYP residual relationships
#
# IMPORTANT
# ------------------------------------------------------------
# This dataset contains NO estrogen intervention.
#
# Therefore this analysis CANNOT demonstrate:
#   - estrogen treatment effects
#   - ER-dependent regulation
#   - estrogen causality
#   - estrogen-mediated CYP regulation
#
# It CAN identify:
#   - estrogen-receptor expression
#   - estrogen-related transcriptional states
#   - spatial/cellular co-variation with CYP/AA biology
#
# No condition-level inferential statistics.
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Matrix)


# ============================================================
# 2. PATHS
# ============================================================

project_dir <- "D:/Master/ScRNA seq/TAC"

scrna_file <- file.path(
  project_dir,
  "CLEAN DATA",
  "11D_C_FINAL_CORRECTED_ANNOTATED_ATLAS.rds"
)

spatial_file <- file.path(
  project_dir,
  "CLEAN DATA",
  "11D_H_FOCUSED_CYP_AA_NETWORK_SPATIAL_OBJECTS.rds"
)

results_dir <- file.path(
  project_dir,
  "RESULTS",
  "11_SPATIAL",
  "11D_I_ESTROGEN_CYP_AA_INTEGRATION"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "11_SPATIAL",
  "11D_I_ESTROGEN_CYP_AA_INTEGRATION"
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
# 3. LOAD DATA
# ============================================================

cat("\n============================================\n")
cat("STEP 11D-I\n")
cat("ESTROGEN SIGNALING x CYP/AA INTEGRATION\n")
cat("============================================\n")

scrna <- readRDS(scrna_file)
spatial_list <- readRDS(spatial_file)

cat("\nscRNA cells:", ncol(scrna), "\n")
cat("scRNA genes:", nrow(scrna), "\n")
cat("Spatial objects:", length(spatial_list), "\n")
cat(
  "Spatial spots:",
  sum(sapply(spatial_list, ncol)),
  "\n"
)


# ============================================================
# 4. INTEGRITY
# ============================================================

stopifnot(
  ncol(scrna) == 25186,
  nrow(scrna) == 32285,
  length(unique(scrna$seurat_clusters)) == 29,
  length(unique(scrna$corrected_broad_cell_type)) == 15,
  length(spatial_list) == 4,
  sum(sapply(spatial_list, ncol)) == 7079
)


# ============================================================
# 5. CONDITIONS
# ============================================================

condition_levels <- c(
  "Sham",
  "2W",
  "4W",
  "6W"
)

sample_to_condition <- c(
  "Sham" = "Sham",
  "TAC_2W" = "2W",
  "TAC_4W" = "4W",
  "TAC_6W" = "6W"
)

scrna$Condition <- unname(
  sample_to_condition[
    as.character(scrna$sample_id)
  ]
)

scrna$Condition <- factor(
  scrna$Condition,
  levels = condition_levels
)

if (any(is.na(scrna$Condition))) {
  stop("Condition mapping failed for scRNA.")
}

for (i in seq_along(spatial_list)) {
  
  spatial_list[[i]]$Condition <- factor(
    as.character(spatial_list[[i]]$Condition),
    levels = condition_levels
  )
}


# ============================================================
# 6. EXPRESSION MATRICES
# ============================================================

DefaultAssay(scrna) <- "RNA"

scrna_data <- GetAssayData(
  scrna,
  assay = "RNA",
  layer = "data"
)

shared_genes <- Reduce(
  intersect,
  c(
    list(rownames(scrna)),
    lapply(spatial_list, rownames)
  )
)


# ============================================================
# 7. ESTROGEN RECEPTOR FEASIBILITY PANEL
# ============================================================

estrogen_receptors <- c(
  "Esr1",
  "Esr2",
  "Gper1"
)

receptor_presence <- data.frame(
  
  Gene = estrogen_receptors,
  
  Present_scRNA =
    estrogen_receptors %in%
    rownames(scrna),
  
  Present_All_Spatial =
    estrogen_receptors %in%
    shared_genes,
  
  stringsAsFactors = FALSE
)

cat("\n============================================\n")
cat("ESTROGEN RECEPTOR AVAILABILITY\n")
cat("============================================\n")

print(
  receptor_presence,
  row.names = FALSE
)

write.csv(
  receptor_presence,
  file.path(
    results_dir,
    "11D_I_Estrogen_Receptor_Availability.csv"
  ),
  row.names = FALSE
)

available_receptors <- intersect(
  estrogen_receptors,
  shared_genes
)


# ============================================================
# 8. HELPER: EXPRESSION SUMMARY
# ============================================================

summarize_expression <- function(
    expression_matrix,
    genes
) {
  
  genes <- intersect(
    genes,
    rownames(expression_matrix)
  )
  
  if (length(genes) == 0) {
    
    return(
      data.frame()
    )
  }
  
  bind_rows(
    lapply(
      genes,
      function(g) {
        
        x <- as.numeric(
          expression_matrix[
            g,
            ,
            drop = TRUE
          ]
        )
        
        data.frame(
          
          Gene = g,
          
          Mean_Expression =
            mean(x),
          
          Median_Expression =
            median(x),
          
          Percent_Detected =
            mean(x > 0) * 100,
          
          Detected =
            sum(x > 0),
          
          Total =
            length(x),
          
          stringsAsFactors = FALSE
        )
      }
    )
  )
}


# ============================================================
# 9. GLOBAL scRNA RECEPTOR DETECTION
# ============================================================

global_receptor_detection <-
  summarize_expression(
    scrna_data,
    estrogen_receptors
  )

cat("\n============================================\n")
cat("GLOBAL scRNA ESTROGEN RECEPTOR DETECTION\n")
cat("============================================\n")

print(
  global_receptor_detection,
  row.names = FALSE
)

write.csv(
  global_receptor_detection,
  file.path(
    results_dir,
    "11D_I_Global_scRNA_Estrogen_Receptor_Detection.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. RECEPTOR DETECTION BY CONDITION
# ============================================================

receptor_condition_list <- list()
counter <- 1

for (condition_now in condition_levels) {
  
  cells_now <- colnames(scrna)[
    scrna$Condition == condition_now
  ]
  
  temp <- summarize_expression(
    scrna_data[
      ,
      cells_now,
      drop = FALSE
    ],
    estrogen_receptors
  )
  
  if (nrow(temp) > 0) {
    
    temp$Condition <- condition_now
    temp$Cells <- length(cells_now)
    
    receptor_condition_list[[counter]] <- temp
    counter <- counter + 1
  }
}

receptor_condition <- bind_rows(
  receptor_condition_list
)

receptor_condition$Condition <- factor(
  receptor_condition$Condition,
  levels = condition_levels
)

cat("\n============================================\n")
cat("ESTROGEN RECEPTORS BY CONDITION\n")
cat("============================================\n")

print(
  as.data.frame(
    receptor_condition
  ),
  row.names = FALSE
)

write.csv(
  receptor_condition,
  file.path(
    results_dir,
    "11D_I_scRNA_Estrogen_Receptors_By_Condition.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. RECEPTOR DETECTION BY CELL TYPE x CONDITION
# ============================================================

cell_types <- sort(
  unique(
    as.character(
      scrna$corrected_broad_cell_type
    )
  )
)

receptor_celltype_list <- list()
counter <- 1

for (condition_now in condition_levels) {
  
  for (celltype_now in cell_types) {
    
    cells_now <- colnames(scrna)[
      scrna$Condition == condition_now &
        scrna$corrected_broad_cell_type == celltype_now
    ]
    
    if (length(cells_now) == 0) {
      next
    }
    
    temp <- summarize_expression(
      scrna_data[
        ,
        cells_now,
        drop = FALSE
      ],
      estrogen_receptors
    )
    
    if (nrow(temp) == 0) {
      next
    }
    
    temp$Condition <- condition_now
    temp$Cell_Type <- celltype_now
    temp$Cells <- length(cells_now)
    
    receptor_celltype_list[[counter]] <- temp
    counter <- counter + 1
  }
}

receptor_celltype <- bind_rows(
  receptor_celltype_list
)

receptor_celltype$Condition <- factor(
  receptor_celltype$Condition,
  levels = condition_levels
)

write.csv(
  receptor_celltype,
  file.path(
    results_dir,
    "11D_I_scRNA_Estrogen_Receptors_By_CellType_Condition.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. TOP RECEPTOR-EXPRESSING CELL TYPES
# ============================================================

top_receptor_celltypes <- receptor_celltype %>%
  
  group_by(
    Condition,
    Gene
  ) %>%
  
  slice_max(
    order_by = Percent_Detected,
    n = 3,
    with_ties = FALSE
  ) %>%
  
  ungroup() %>%
  
  arrange(
    Gene,
    Condition,
    desc(Percent_Detected)
  )

cat("\n============================================\n")
cat("TOP ESTROGEN-RECEPTOR CELL TYPES\n")
cat("============================================\n")

print(
  as.data.frame(
    top_receptor_celltypes
  ),
  row.names = FALSE
)

write.csv(
  top_receptor_celltypes,
  file.path(
    results_dir,
    "11D_I_Top_Estrogen_Receptor_CellTypes.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. SPATIAL RECEPTOR DETECTION
# ============================================================

spatial_receptor_list <- list()
counter <- 1

for (i in seq_along(spatial_list)) {
  
  obj <- spatial_list[[i]]
  
  DefaultAssay(obj) <- "Spatial"
  
  spatial_data <- GetAssayData(
    obj,
    assay = "Spatial",
    layer = "data"
  )
  
  condition_now <- unique(
    as.character(
      obj$Condition
    )
  )
  
  temp <- summarize_expression(
    spatial_data,
    estrogen_receptors
  )
  
  if (nrow(temp) > 0) {
    
    temp$Condition <- condition_now
    temp$Spots <- ncol(obj)
    
    spatial_receptor_list[[counter]] <- temp
    counter <- counter + 1
  }
}

spatial_receptor <- bind_rows(
  spatial_receptor_list
)

spatial_receptor$Condition <- factor(
  spatial_receptor$Condition,
  levels = condition_levels
)

cat("\n============================================\n")
cat("SPATIAL ESTROGEN RECEPTOR DETECTION\n")
cat("============================================\n")

print(
  as.data.frame(
    spatial_receptor
  ),
  row.names = FALSE
)

write.csv(
  spatial_receptor,
  file.path(
    results_dir,
    "11D_I_Spatial_Estrogen_Receptor_Detection.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. FEASIBILITY CLASSIFICATION
#
# This is an analytical quality-control classification.
#
# It is NOT a biological significance test.
# ============================================================

receptor_feasibility <- global_receptor_detection %>%
  
  mutate(
    
    Detection_Class =
      case_when(
        
        Percent_Detected >= 5 ~
          "Robust",
        
        Percent_Detected >= 1 ~
          "Moderate",
        
        Percent_Detected > 0 ~
          "Sparse",
        
        TRUE ~
          "Undetected"
      )
  )

cat("\n============================================\n")
cat("ESTROGEN RECEPTOR FEASIBILITY CLASSIFICATION\n")
cat("============================================\n")

print(
  receptor_feasibility,
  row.names = FALSE
)

write.csv(
  receptor_feasibility,
  file.path(
    results_dir,
    "11D_I_Estrogen_Receptor_Feasibility.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. DEFINE ESTROGEN-RELATED TRANSCRIPTIONAL PROGRAMS
#
# IMPORTANT:
# These are compact predefined transcriptional programs.
#
# CYP genes are deliberately NOT included.
# AA/eicosanoid genes are deliberately NOT included.
#
# This reduces circularity when relating estrogen programs
# to CYP/AA biology.
# ============================================================

estrogen_candidate_sets <- list(
  
  Estrogen_Response_Early = c(
    "Tff1",
    "Pgr",
    "Greb1",
    "Igfbp4",
    "Krt19",
    "Ccnd1",
    "Myc",
    "Fos",
    "Jun",
    "Egr1",
    "Nrip1",
    "Klf4"
  ),
  
  Estrogen_Response_Late = c(
    "Pgr",
    "Greb1",
    "Igfbp4",
    "Ccnd1",
    "Myc",
    "Bcl2",
    "Socs2",
    "Cav1",
    "Klf9",
    "Nrip1",
    "Tgfa",
    "Ret"
  ),
  
  ER_Transcriptional_Coregulation = c(
    "Ncoa1",
    "Ncoa2",
    "Ncoa3",
    "Ncor1",
    "Ncor2",
    "Med1",
    "Crebbp",
    "Ep300",
    "Nrip1",
    "Foxa1",
    "Gata3",
    "Sp1"
  )
)


# ============================================================
# 16. REMOVE CYP / AA GENES FROM ESTROGEN PROGRAMS
# ============================================================

strict_cyp_genes <- grep(
  "^Cyp[0-9]",
  shared_genes,
  value = TRUE
)

aa_exclusion_candidates <- c(
  
  "Pla2g4a",
  "Pla2g4b",
  "Pla2g4c",
  "Pla2g6",
  "Pla2g2a",
  "Pla2g5",
  
  "Ephx1",
  "Ephx2",
  
  "Ptgs1",
  "Ptgs2",
  "Ptges",
  "Ptges2",
  "Ptges3",
  "Ptgis",
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

aa_exclusion_genes <- intersect(
  aa_exclusion_candidates,
  shared_genes
)

estrogen_exclusion <- unique(
  c(
    strict_cyp_genes,
    aa_exclusion_genes
  )
)

estrogen_sets <- lapply(
  estrogen_candidate_sets,
  function(x) {
    
    genes_now <- intersect(
      x,
      shared_genes
    )
    
    setdiff(
      genes_now,
      estrogen_exclusion
    )
  }
)


# ============================================================
# 17. ESTROGEN PROGRAM AVAILABILITY
# ============================================================

estrogen_set_summary <- data.frame(
  
  Estrogen_Program =
    names(estrogen_sets),
  
  Candidate_Genes =
    sapply(
      estrogen_candidate_sets,
      length
    ),
  
  Available_Genes =
    sapply(
      estrogen_sets,
      length
    ),
  
  stringsAsFactors = FALSE
)

cat("\n============================================\n")
cat("ESTROGEN PROGRAM GENE AVAILABILITY\n")
cat("============================================\n")

print(
  estrogen_set_summary,
  row.names = FALSE
)

write.csv(
  estrogen_set_summary,
  file.path(
    results_dir,
    "11D_I_Estrogen_Program_Gene_Availability.csv"
  ),
  row.names = FALSE
)

estrogen_gene_table <- bind_rows(
  lapply(
    names(estrogen_sets),
    function(x) {
      
      data.frame(
        Estrogen_Program = x,
        Gene = estrogen_sets[[x]],
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  estrogen_gene_table,
  file.path(
    results_dir,
    "11D_I_Final_Estrogen_Gene_Sets.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. HELPER: PROGRAM SCORES
# ============================================================

calculate_scores <- function(
    expression_matrix,
    gene_sets,
    prefix
) {
  
  result <- lapply(
    names(gene_sets),
    function(set_name) {
      
      genes_now <- intersect(
        gene_sets[[set_name]],
        rownames(expression_matrix)
      )
      
      if (length(genes_now) == 0) {
        
        return(
          rep(
            NA_real_,
            ncol(expression_matrix)
          )
        )
      }
      
      Matrix::colMeans(
        expression_matrix[
          genes_now,
          ,
          drop = FALSE
        ]
      )
    }
  )
  
  result <- as.data.frame(
    result,
    check.names = FALSE
  )
  
  colnames(result) <- paste0(
    prefix,
    names(gene_sets)
  )
  
  result
}


# ============================================================
# 19. SCORE ESTROGEN PROGRAMS IN scRNA
# ============================================================

scrna_estrogen_scores <- calculate_scores(
  scrna_data,
  estrogen_sets,
  "Estrogen_"
)

rownames(scrna_estrogen_scores) <-
  colnames(scrna)

for (score_name in
     colnames(scrna_estrogen_scores)) {
  
  scrna[[score_name]] <-
    scrna_estrogen_scores[
      colnames(scrna),
      score_name
    ]
}


# ============================================================
# 20. scRNA ESTROGEN PROGRAM TRAJECTORIES
# ============================================================

estrogen_score_cols <- colnames(
  scrna_estrogen_scores
)

scrna_estrogen_long <- scrna@meta.data %>%
  
  select(
    Condition,
    corrected_broad_cell_type,
    all_of(estrogen_score_cols)
  ) %>%
  
  pivot_longer(
    cols = all_of(estrogen_score_cols),
    names_to = "Estrogen_Program",
    values_to = "Score"
  ) %>%
  
  mutate(
    Estrogen_Program =
      sub(
        "^Estrogen_",
        "",
        Estrogen_Program
      )
  )

scrna_estrogen_condition <-
  scrna_estrogen_long %>%
  
  group_by(
    Condition,
    Estrogen_Program
  ) %>%
  
  summarise(
    Cells = n(),
    Mean_Score =
      mean(
        Score,
        na.rm = TRUE
      ),
    Median_Score =
      median(
        Score,
        na.rm = TRUE
      ),
    SD_Score =
      sd(
        Score,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

write.csv(
  scrna_estrogen_condition,
  file.path(
    results_dir,
    "11D_I_scRNA_Estrogen_Program_Trajectories.csv"
  ),
  row.names = FALSE
)

cat("\n============================================\n")
cat("scRNA ESTROGEN-PROGRAM TRAJECTORIES\n")
cat("============================================\n")

print(
  as.data.frame(
    scrna_estrogen_condition
  ),
  row.names = FALSE