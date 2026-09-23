# ============================================================
# STEP 11D-H
# FOCUSED CYP–ARACHIDONIC ACID METABOLIC NETWORK
#
# PURPOSE
# ------------------------------------------------------------
# Reconstruct the CYP/AA metabolic architecture across:
#
#   Sham -> TAC 2W -> TAC 4W -> TAC 6W
#
# using:
#   1. corrected scRNA-seq atlas
#   2. spatial transcriptomics
#
# QUESTIONS
# ------------------------------------------------------------
# 1. How are epoxygenase-related CYPs remodeled?
# 2. How is Ephx2 remodeled relative to epoxygenase genes?
# 3. How is the CYP4F/CYP4A hydroxylase arm remodeled?
# 4. How do CYP branches relate to PLA2, COX and LOX biology?
# 5. Which cell types carry each AA metabolic branch?
# 6. Is there evidence of temporal/spatial "handoff" between
#    CYP/AA branches?
#
# IMPORTANT INTERPRETATION
# ------------------------------------------------------------
# These are transcriptomic expression programs.
#
# They are NOT:
#   - direct EET concentrations
#   - direct HETE concentrations
#   - enzymatic activity measurements
#   - metabolic flux
#   - EET/DHET biochemical ratios
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
  "11D_G_MECHANISTICALLY_SCORED_SPATIAL_OBJECTS.rds"
)

results_dir <- file.path(
  project_dir,
  "RESULTS",
  "11_SPATIAL",
  "11D_H_FOCUSED_CYP_AA_NETWORK"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "11_SPATIAL",
  "11D_H_FOCUSED_CYP_AA_NETWORK"
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
cat("STEP 11D-H\n")
cat("FOCUSED CYP–AA METABOLIC NETWORK\n")
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
    as.character(
      spatial_list[[i]]$Condition
    ),
    levels = condition_levels
  )
}


# ============================================================
# 6. DEFINE AA METABOLIC ARCHITECTURE
#
# Candidate genes are deliberately defined BEFORE looking at
# condition-specific results.
#
# The script will retain only genes actually present.
# ============================================================

aa_candidate_sets <- list(
  
  AA_Liberation = c(
    "Pla2g4a",
    "Pla2g4b",
    "Pla2g4c",
    "Pla2g6",
    "Pla2g2a",
    "Pla2g5"
  ),
  
  CYP_Epoxygenase = c(
    "Cyp2j5",
    "Cyp2j6",
    "Cyp2j8",
    "Cyp2j9",
    "Cyp2j11",
    "Cyp2j12",
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
    "Cyp2c67",
    "Cyp2c68",
    "Cyp2c69",
    "Cyp2c70"
  ),
  
  Epoxide_Hydrolysis = c(
    "Ephx1",
    "Ephx2"
  ),
  
  CYP_Omega_Hydroxylase = c(
    "Cyp4a10",
    "Cyp4a12a",
    "Cyp4a12b",
    "Cyp4a14",
    "Cyp4f13",
    "Cyp4f14",
    "Cyp4f15",
    "Cyp4f16",
    "Cyp4f17",
    "Cyp4f18",
    "Cyp4f39"
  ),
  
  COX_Arm = c(
    "Ptgs1",
    "Ptgs2",
    "Ptges",
    "Ptges2",
    "Ptges3",
    "Ptgis",
    "Tbxas1"
  ),
  
  LOX_Arm = c(
    "Alox5",
    "Alox5ap",
    "Alox12",
    "Alox12b",
    "Alox15",
    "Alox15b"
  ),
  
  Eicosanoid_Receptors = c(
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


# ============================================================
# 7. SHARED GENES
# ============================================================

shared_genes <- Reduce(
  intersect,
  c(
    list(rownames(scrna)),
    lapply(spatial_list, rownames)
  )
)

aa_sets <- lapply(
  aa_candidate_sets,
  function(x) {
    intersect(x, shared_genes)
  }
)

aa_set_summary <- data.frame(
  AA_Arm = names(aa_sets),
  Candidate_Genes = sapply(
    aa_candidate_sets,
    length
  ),
  Available_Genes = sapply(
    aa_sets,
    length
  ),
  stringsAsFactors = FALSE
)

write.csv(
  aa_set_summary,
  file.path(
    results_dir,
    "11D_H_AA_Arm_Gene_Availability.csv"
  ),
  row.names = FALSE
)

aa_gene_table <- bind_rows(
  lapply(
    names(aa_sets),
    function(x) {
      data.frame(
        AA_Arm = x,
        Gene = aa_sets[[x]],
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  aa_gene_table,
  file.path(
    results_dir,
    "11D_H_Final_AA_Gene_Sets.csv"
  ),
  row.names = FALSE
)

cat("\n============================================\n")
cat("AA METABOLIC ARM AVAILABILITY\n")
cat("============================================\n")

print(
  aa_set_summary,
  row.names = FALSE
)


# ============================================================
# 8. PRIORITY INDIVIDUAL GENES
# ============================================================

priority_genes <- c(
  "Cyp1b1",
  "Cyp2j6",
  "Cyp2j9",
  "Cyp4f13",
  "Cyp4f16",
  "Cyp4f17",
  "Cyp4f18",
  "Ephx2"
)

priority_available <- intersect(
  priority_genes,
  shared_genes
)

cat(
  "\nPriority genes available:",
  length(priority_available),
  "/",
  length(priority_genes),
  "\n"
)

print(priority_available)


# ============================================================
# 9. CORE AA GENES FOR GENE-LEVEL ANALYSIS
# ============================================================

core_aa_genes <- unique(
  c(
    priority_available,
    unlist(aa_sets)
  )
)

core_aa_genes <- intersect(
  core_aa_genes,
  shared_genes
)

write.csv(
  data.frame(
    Gene = core_aa_genes
  ),
  file.path(
    results_dir,
    "11D_H_Core_AA_Genes.csv"
  ),
  row.names = FALSE
)

cat(
  "Core AA genes available:",
  length(core_aa_genes),
  "\n"
)


# ============================================================
# 10. HELPER: EXPRESSION SUMMARY
# ============================================================

summarize_gene_expression <- function(
    expression_matrix,
    genes
) {
  
  genes <- intersect(
    genes,
    rownames(expression_matrix)
  )
  
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
          Mean_Expression = mean(x),
          Median_Expression = median(x),
          Percent_Detected = mean(x > 0) * 100,
          stringsAsFactors = FALSE
        )
      }
    )
  )
}


# ============================================================
# 11. scRNA GENE-LEVEL AA LANDSCAPE BY CONDITION
# ============================================================

DefaultAssay(scrna) <- "RNA"

scrna_data <- GetAssayData(
  scrna,
  assay = "RNA",
  layer = "data"
)

scrna_condition_gene_list <- list()
counter <- 1

for (condition_now in condition_levels) {
  
  cells_now <- colnames(scrna)[
    scrna$Condition == condition_now
  ]
  
  expression_now <- scrna_data[
    ,
    cells_now,
    drop = FALSE
  ]
  
  temp <- summarize_gene_expression(
    expression_now,
    core_aa_genes
  )
  
  temp$Condition <- condition_now
  temp$Cells <- length(cells_now)
  
  scrna_condition_gene_list[[counter]] <- temp
  counter <- counter + 1
}

scrna_gene_condition <- bind_rows(
  scrna_condition_gene_list
)

scrna_gene_condition$Condition <- factor(
  scrna_gene_condition$Condition,
  levels = condition_levels
)

write.csv(
  scrna_gene_condition,
  file.path(
    results_dir,
    "11D_H_scRNA_AA_Gene_Trajectories.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. scRNA GENE-LEVEL LANDSCAPE BY CELL TYPE x CONDITION
# ============================================================

cell_types <- sort(
  unique(
    as.character(
      scrna$corrected_broad_cell_type
    )
  )
)

celltype_gene_list <- list()
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
    
    expression_now <- scrna_data[
      ,
      cells_now,
      drop = FALSE
    ]
    
    temp <- summarize_gene_expression(
      expression_now,
      core_aa_genes
    )
    
    temp$Condition <- condition_now
    temp$Cell_Type <- celltype_now
    temp$Cells <- length(cells_now)
    
    celltype_gene_list[[counter]] <- temp
    counter <- counter + 1
  }
}

scrna_celltype_gene <- bind_rows(
  celltype_gene_list
)

scrna_celltype_gene$Condition <- factor(
  scrna_celltype_gene$Condition,
  levels = condition_levels
)

write.csv(
  scrna_celltype_gene,
  file.path(
    results_dir,
    "11D_H_scRNA_AA_Genes_By_CellType_Condition.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. HELPER: AA ARM SCORES
#
# Mean normalized expression.
# ============================================================

calculate_scores <- function(
    expression_matrix,
    gene_sets
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
    "AA_",
    names(gene_sets)
  )
  
  result
}


# ============================================================
# 14. SCORE scRNA AA ARMS
# ============================================================

scrna_aa_scores <- calculate_scores(
  scrna_data,
  aa_sets
)

rownames(scrna_aa_scores) <- colnames(scrna)

for (score_name in colnames(scrna_aa_scores)) {
  
  scrna[[score_name]] <-
    scrna_aa_scores[
      colnames(scrna),
      score_name
    ]
}


# ============================================================
# 15. scRNA AA ARM SUMMARY
# ============================================================

aa_score_columns <- colnames(
  scrna_aa_scores
)

scrna_meta <- scrna@meta.data %>%
  mutate(
    Cell = rownames(scrna@meta.data)
  )

scrna_arm_long <- scrna_meta %>%
  
  select(
    Cell,
    Condition,
    corrected_broad_cell_type,
    all_of(aa_score_columns)
  ) %>%
  
  pivot_longer(
    cols = all_of(aa_score_columns),
    names_to = "AA_Arm",
    values_to = "Score"
  ) %>%
  
  mutate(
    AA_Arm = sub(
      "^AA_",
      "",
      AA_Arm
    )
  )


# ============================================================
# 16. WHOLE-ATLAS AA ARM TRAJECTORIES
# ============================================================

scrna_arm_condition <- scrna_arm_long %>%
  
  group_by(
    Condition,
    AA_Arm
  ) %>%
  
  summarise(
    Cells = n(),
    Mean_Score = mean(
      Score,
      na.rm = TRUE
    ),
    Median_Score = median(
      Score,
      na.rm = TRUE
    ),
    SD_Score = sd(
      Score,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write.csv(
  scrna_arm_condition,
  file.path(
    results_dir,
    "11D_H_scRNA_AA_Arm_Trajectories.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. CELL-TYPE-SPECIFIC AA ARM TRAJECTORIES
# ============================================================

scrna_arm_celltype <- scrna_arm_long %>%
  
  group_by(
    Condition,
    corrected_broad_cell_type,
    AA_Arm
  ) %>%
  
  summarise(
    Cells = n(),
    Mean_Score = mean(
      Score,
      na.rm = TRUE
    ),
    Median_Score = median(
      Score,
      na.rm = TRUE
    ),
    SD_Score = sd(
      Score,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write.csv(
  scrna_arm_celltype,
  file.path(
    results_dir,
    "11D_H_scRNA_AA_Arms_By_CellType.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. 4W AA ARM CHANGE vs SHAM WITHIN CELL TYPE
# ============================================================

sham_arm <- scrna_arm_celltype %>%
  
  filter(
    Condition == "Sham"
  ) %>%
  
  select(
    corrected_broad_cell_type,
    AA_Arm,
    Sham_Mean = Mean_Score
  )

arm_change <- scrna_arm_celltype %>%
  
  left_join(
    sham_arm,
    by = c(
      "corrected_broad_cell_type",
      "AA_Arm"
    )
  ) %>%
  
  mutate(
    Delta_vs_Sham =
      Mean_Score - Sham_Mean,
    
    Fold_vs_Sham =
      (Mean_Score + 1e-8) /
      (Sham_Mean + 1e-8),
    
    Log2FC_vs_Sham =
      log2(Fold_vs_Sham)
  )

arm_4w <- arm_change %>%
  
  filter(
    Condition == "4W",
    Cells >= 30
  ) %>%
  
  arrange(
    desc(
      abs(
        Delta_vs_Sham
      )
    )
  )

write.csv(
  arm_4w,
  file.path(
    results_dir,
    "11D_H_scRNA_4W_AA_Arm_State.csv"
  ),
  row.names = FALSE
)

cat("\n============================================\n")
cat("4W scRNA AA-ARM REMODELING\n")
cat("============================================\n")

print(
  as.data.frame(
    head(
      arm_4w,
      50
    )
  ),
  row.names = FALSE
)


# ============================================================
# 19. FIGURE: 4W AA ARM CELL-TYPE HEATMAP
# ============================================================

p_aa_4w <- ggplot(
  arm_4w,
  aes(
    x = AA_Arm,
    y = corrected_broad_cell_type,
    fill = Delta_vs_Sham
  )
) +
  geom_tile() +
  labs(
    title =
      "Cell-type-specific AA metabolic remodeling at 4W TAC",
    subtitle =
      "Difference in mean transcriptomic AA-arm score from Sham",
    x = NULL,
    y = NULL,
    fill = "Delta score"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  file.path(
    figures_dir,
    "11D_H_scRNA_4W_AA_Arm_Heatmap.png"
  ),
  p_aa_4w,
  width = 11,
  height = 8,
  dpi = 400
)


# ============================================================
# 20. EPOXYGENASE–EPHX2 EXPRESSION BALANCE
#
# IMPORTANT:
# This is NOT an EET/DHET ratio.
#
# It is a transcriptomic expression index comparing the
# epoxygenase-related gene score with Ephx2 expression.
# ============================================================

epoxygenase_col <- "AA_CYP_Epoxygenase"

if (!epoxygenase_col %in%
    colnames(scrna@meta.data)) {
  stop("Epoxygenase score missing.")
}

if (!"Ephx2" %in% rownames(scrna_data)) {
  stop("Ephx2 missing.")
}

scrna$Ephx2_NormalizedExpression <-
  as.numeric(
    scrna_data[
      "Ephx2",
      ,
      drop = TRUE
    ]
  )

scrna$Epoxygenase_Ephx2_Balance <-
  scrna[[epoxygenase_col]][, 1] -
  scrna$Ephx2_NormalizedExpression


# ============================================================
# 21. SUMMARIZE EPOXYGENASE–EPHX2 BALANCE
# ============================================================

balance_summary <- scrna@meta.data %>%
  
  group_by(
    Condition,
    corrected_broad_cell_type
  ) %>%
  
  summarise(
    Cells = n(),
    
    Mean_Epoxygenase =
      mean(
        AA_CYP_Epoxygenase,
        na.rm = TRUE
      ),
    
    Mean_Ephx2 =
      mean(
        Ephx2_NormalizedExpression,
        na.rm = TRUE
      ),
    
    Mean_Expression_Balance =
      mean(
        Epoxygenase_Ephx2_Balance,
        na.rm = TRUE
      ),
    
    Median_Expression_Balance =
      median(
        Epoxygenase_Ephx2_Balance,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

write.csv(
  balance_summary,
  file.path(
    results_dir,
    "11D_H_scRNA_Epoxygenase_Ephx2_Expression_Balance.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 22. PRIORITY CYP FAMILY HANDOFF
#
# Determine which cell types contribute detected expression
# for each priority CYP/AA gene at each condition.
# ============================================================

handoff_list <- list()
counter <- 1

for (condition_now in condition_levels) {
  
  for (gene_now in priority_available) {
    
    cells_condition <- colnames(scrna)[
      scrna$Condition == condition_now
    ]
    
    x <- as.numeric(
      scrna_data[
        gene_now,
        cells_condition,
        drop = TRUE
      ]
    )
    
    celltypes_now <- as.character(
      scrna$corrected_broad_cell_type[
        match(
          cells_condition,