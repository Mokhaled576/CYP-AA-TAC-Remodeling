# ============================================================
# STEP 11D-D
# CORRECTED CYP / ARACHIDONIC ACID / EICOSANOID
# CELL-TYPE REANALYSIS
#
# PURPOSE
# ------------------------------------------------------------
# Recalculate all annotation-dependent CYP/AA results using the
# corrected frozen atlas generated in Step 11D-C.
#
# This script:
#   1. Uses corrected broad/refined annotations
#   2. Does NOT alter the atlas
#   3. Does NOT recluster or filter cells
#   4. Quantifies CYP genes by corrected cell type
#   5. Quantifies priority CYP/AA genes across conditions
#   6. Separates:
#        A. cell-composition changes
#        B. within-cell-type expression/detection changes
#   7. Produces publication-oriented tables and figures
#
# IMPORTANT
# ------------------------------------------------------------
# Condition comparisons are DESCRIPTIVE.
# Cells are NOT treated as biological replicates.
#
# INPUT
# ------------------------------------------------------------
# 11D_C_FINAL_CORRECTED_ANNOTATED_ATLAS.rds
#
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

input_file <- file.path(
  project_dir,
  "CLEAN DATA",
  "11D_C_FINAL_CORRECTED_ANNOTATED_ATLAS.rds"
)

results_dir <- file.path(
  project_dir,
  "RESULTS",
  "11_SPATIAL",
  "11D_D_CORRECTED_CYP_AA_CELLTYPE"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "11_SPATIAL",
  "11D_D_CORRECTED_CYP_AA_CELLTYPE"
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
# 3. LOAD CORRECTED FROZEN ATLAS
# ============================================================

cat("\n============================================\n")
cat("STEP 11D-D\n")
cat("CORRECTED CYP/AA CELL-TYPE REANALYSIS\n")
cat("============================================\n")

scrna <- readRDS(input_file)

DefaultAssay(scrna) <- "RNA"

cat("Input:", input_file, "\n")
cat("Cells:", ncol(scrna), "\n")
cat("Genes:", nrow(scrna), "\n")


# ============================================================
# 4. VERIFY CORRECTED ATLAS
# ============================================================

required_metadata <- c(
  "seurat_clusters",
  "sample_id",
  "corrected_broad_cell_type",
  "corrected_refined_cell_type",
  "corrected_annotation_confidence"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(scrna@meta.data)
)

if (length(missing_metadata) > 0) {
  
  stop(
    paste(
      "Missing required corrected metadata:",
      paste(missing_metadata, collapse = ", ")
    )
  )
}

if (ncol(scrna) != 25186) {
  stop("Unexpected cell count.")
}

cluster_ids <- sort(
  unique(
    as.numeric(
      as.character(
        scrna$seurat_clusters
      )
    )
  )
)

if (!setequal(
  cluster_ids,
  as.numeric(0:28)
)) {
  stop("Frozen cluster structure is not 0-28.")
}

cat("Corrected atlas verified.\n")


# ============================================================
# 5. IDENTIFY CONDITION COLUMN
#
# Prefer explicit condition metadata if present.
# Otherwise derive condition from sample_id.
# ============================================================
# ============================================================
# 5-6. DEFINE AND STANDARDIZE CONDITION
# CORRECTED VERSION
# ============================================================

cat("\n============================================\n")
cat("SAMPLE IDs BEFORE CONDITION MAPPING\n")
cat("============================================\n")

print(
  table(
    meta$sample_id,
    useNA = "ifany"
  )
)

# Convert sample_id explicitly to character
sample_text <- as.character(
  meta$sample_id
)

# Direct mapping from the actual sample IDs
meta$Condition <- case_when(
  
  sample_text == "Sham" ~ "Sham",
  
  sample_text == "TAC_2W" ~ "2W",
  
  sample_text == "TAC_4W" ~ "4W",
  
  sample_text == "TAC_6W" ~ "6W",
  
  TRUE ~ NA_character_
)

# Expected biological order
condition_levels <- c(
  "Sham",
  "2W",
  "4W",
  "6W"
)

meta$Condition <- factor(
  meta$Condition,
  levels = condition_levels
)


# ============================================================
# STRICT VALIDATION
# ============================================================

if (any(is.na(meta$Condition))) {
  
  bad_samples <- unique(
    sample_text[
      is.na(meta$Condition)
    ]
  )
  
  stop(
    paste(
      "Condition assignment failed for sample_id:",
      paste(
        bad_samples,
        collapse = ", "
      )
    )
  )
}

# Check that all four expected conditions exist
observed_conditions <- unique(
  as.character(
    meta$Condition
  )
)

if (!setequal(
  observed_conditions,
  condition_levels
)) {
  
  stop(
    paste(
      "Expected Sham, 2W, 4W, 6W but found:",
      paste(
        observed_conditions,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# ADD CONDITION TO SEURAT OBJECT
# ============================================================

scrna$Condition <- meta$Condition


# ============================================================
# VERIFY MAPPING
# ============================================================

condition_mapping_check <- data.frame(
  sample_id = as.character(
    scrna$sample_id
  ),
  Condition = as.character(
    scrna$Condition
  )
) %>%
  distinct() %>%
  arrange(
    factor(
      Condition,
      levels = condition_levels
    )
  )

cat("\n============================================\n")
cat("CONDITION MAPPING CHECK\n")
cat("============================================\n")

print(
  condition_mapping_check,
  row.names = FALSE
)

cat("\nCondition mapping successfully verified.\n")

# ============================================================
# 7. CONDITION COUNTS
# ============================================================

condition_counts <- data.frame(
  Condition = condition_levels,
  Cells = as.numeric(
    table(
      factor(
        scrna$Condition,
        levels = condition_levels
      )
    )
  )
)

write.csv(
  condition_counts,
  file.path(
    results_dir,
    "11D_D_Condition_Cell_Counts.csv"
  ),
  row.names = FALSE
)

cat("\n============================================\n")
cat("CONDITION CELL COUNTS\n")
cat("============================================\n")

print(
  condition_counts,
  row.names = FALSE
)


# ============================================================
# 8. DEFINE CYP GENE UNIVERSE
#
# Same strict CYP definition used previously.
# ============================================================

cyp_genes <- grep(
  "^Cyp[0-9]",
  rownames(scrna),
  value = TRUE
)

cyp_genes <- sort(
  unique(
    cyp_genes
  )
)

cat(
  "\nStrict CYP genes detected:",
  length(cyp_genes),
  "\n"
)

if (length(cyp_genes) == 0) {
  stop("No CYP genes detected.")
}

write.csv(
  data.frame(
    Gene = cyp_genes
  ),
  file.path(
    results_dir,
    "11D_D_CYP_Gene_Universe.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9. PRIORITY CYP/AA GENES
# ============================================================

priority_genes_requested <- c(
  "Cyp1b1",
  "Cyp2j6",
  "Cyp2j9",
  "Cyp4f13",
  "Cyp4f16",
  "Cyp4f17",
  "Cyp4f18",
  "Ephx2"
)

priority_genes <- intersect(
  priority_genes_requested,
  rownames(scrna)
)

missing_priority <- setdiff(
  priority_genes_requested,
  priority_genes
)

cat(
  "Priority genes present:",
  paste(
    priority_genes,
    collapse = ", "
  ),
  "\n"
)

if (length(missing_priority) > 0) {
  
  cat(
    "Priority genes missing:",
    paste(
      missing_priority,
      collapse = ", "
    ),
    "\n"
  )
}


# ============================================================
# 10. GET NORMALIZED RNA EXPRESSION
# ============================================================

rna_data <- GetAssayData(
  scrna,
  assay = "RNA",
  layer = "data"
)

cat(
  "Normalized RNA matrix:",
  nrow(rna_data),
  "genes x",
  ncol(rna_data),
  "cells\n"
)


# ============================================================
# 11. WHOLE-ATLAS CYP DETECTION
#
# This is repeated only as a consistency check.
# ============================================================

cyp_matrix <- rna_data[
  cyp_genes,
  ,
  drop = FALSE
]

whole_cyp <- data.frame(
  
  Gene = cyp_genes,
  
  Cells_Expressing =
    as.numeric(
      Matrix::rowSums(
        cyp_matrix > 0
      )
    ),
  
  Percent_Expressing =
    as.numeric(
      Matrix::rowMeans(
        cyp_matrix > 0
      ) * 100
    ),
  
  Mean_Expression =
    as.numeric(
      Matrix::rowMeans(
        cyp_matrix
      )
    ),
  
  stringsAsFactors = FALSE
)

whole_cyp <- whole_cyp %>%
  arrange(
    desc(Percent_Expressing),
    desc(Mean_Expression)
  )

write.csv(
  whole_cyp,
  file.path(
    results_dir,
    "11D_D_Whole_Atlas_CYP_Detection.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. CYP DETECTION CLASSIFICATION
#
# Retain the same descriptive classification logic used in the
# previous CYP landscape where possible:
#
# Robust      >= 5%
# Moderate    >= 1% and < 5%
# Sparse      > 0 and < 1%
# Undetected  = 0
# ============================================================

whole_cyp <- whole_cyp %>%
  
  mutate(
    
    Detection_Class = case_when(
      
      Percent_Expressing >= 5 ~
        "Robust",
      
      Percent_Expressing >= 1 ~
        "Moderate",
      
      Percent_Expressing > 0 ~
        "Sparse",
      
      TRUE ~
        "Undetected"
    )
  )

write.csv(
  whole_cyp,
  file.path(
    results_dir,
    "11D_D_Whole_Atlas_CYP_Detection_Classified.csv"
  ),
  row.names = FALSE
)

cat("\n============================================\n")
cat("CYP DETECTION CLASSIFICATION\n")
cat("============================================\n")

print(
  as.data.frame(
    whole_cyp %>%
      count(
        Detection_Class
      )
  ),
  row.names = FALSE
)


# ============================================================
# 13. HELPER FUNCTION:
# GENE STATISTICS BY GROUP
# ============================================================

calculate_gene_group_stats <- function(
    expression_matrix,
    genes,
    groups
) {
  
  group_levels <- unique(
    as.character(
      groups
    )
  )
  
  output_list <- vector(
    "list",
    length(group_levels)
  )
  
  names(output_list) <- group_levels
  
  for (grp in group_levels) {
    
    cells_now <- which(
      as.character(groups) == grp
    )
    
    mat_now <- expression_matrix[
      genes,
      cells_now,
      drop = FALSE
    ]
    
    output_list[[grp]] <- data.frame(
      
      Group = grp,
      
      Gene = genes,
      
      Total_Cells =
        length(cells_now),
      
      Cells_Expressing =
        as.numeric(
          Matrix::rowSums(
            mat_now > 0
          )
        ),
      
      Percent_Expressing =
        as.numeric(
          Matrix::rowMeans(
            mat_now > 0
          ) * 100
        ),
      
      Mean_Expression =
        as.numeric(
          Matrix::rowMeans(
            mat_now
          )
        ),
      
      Mean_Expression_Among_Expressing =
        vapply(
          seq_len(nrow(mat_now)),
          function(i) {
            
            x <- mat_now[i, ]
            
            positive <- x[x > 0]
            
            if (length(positive) == 0) {
              return(0)
            }
            
            mean(positive)
          },
          numeric(1)
        ),
      
      stringsAsFactors = FALSE
    )
  }
  
  bind_rows(
    output_list
  )
}


# ============================================================
# 14. ALL CYP GENES × CORRECTED BROAD CELL TYPE
# ============================================================

broad_stats <- calculate_gene_group_stats(
  expression_matrix = rna_data,
  genes = cyp_genes,
  groups = scrna$corrected_broad_cell_type
)

broad_stats <- broad_stats %>%
  rename(
    Corrected_Broad_Cell_Type = Group
  )

write.csv(
  broad_stats,
  file.path(
    results_dir,
    "11D_D_ALL_CYP_BY_CORRECTED_BROAD_CELLTYPE.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. ALL CYP GENES × CORRECTED REFINED CELL TYPE
# ============================================================

refined_stats <- calculate_gene_group_stats(
  expression_matrix = rna_data,
  genes = cyp_genes,
  groups = scrna$corrected_refined_cell_type
)

refined_stats <- refined_stats %>%
  rename(
    Corrected_Refined_Cell_Type = Group
  )

write.csv(
  refined_stats,
  file.path(
    results_dir,
    "11D_D_ALL_CYP_BY_CORRECTED_REFINED_CELLTYPE.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. PRIORITY GENES × BROAD CELL TYPE
# ============================================================

priority_broad <- broad_stats %>%
  
  filter(
    Gene %in% priority_genes
  ) %>%
  
  arrange(
    Gene,
    desc(Percent_Expressing)
  )

write.csv(
  priority_broad,
  file.path(
    results_dir,
    "11D_D_PRIORITY_GENES_BY_BROAD_CELLTYPE.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. PRIORITY GENES × REFINED CELL TYPE
# ============================================================

priority_refined <- refined_stats %>%
  
  filter(
    Gene %in% priority_genes
  ) %>%
  
  arrange(
    Gene,
    desc(Percent_Expressing)
  )

write.csv(
  priority_refined,
  file.path(
    results_dir,
    "11D_D_PRIORITY_GENES_BY_REFINED_CELLTYPE.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. HIGHEST CORRECTED BROAD LOCALIZATION
# ============================================================

priority_broad_top <- priority_broad %>%
  
  group_by(
    Gene
  ) %>%
  
  slice_max(
    order_by = Percent_Expressing,
    n = 1,
    with_ties = FALSE
  ) %>%
  
  ungroup() %>%
  
  arrange(
    Gene
  )

write.csv(
  priority_broad_top,
  file.path(
    results_dir,
    "11D_D_PRIORITY_TOP_BROAD_LOCALIZATION.csv"
  ),
  row.names = FALSE
)

cat("\n============================================\n")
cat("CORRECTED PRIORITY-GENE BROAD LOCALIZATION\n")
cat("============================================\n")

print(
  as.data.frame(
    priority_broad_top %>%
      select(
        Gene,
        Corrected_Broad_Cell_Type,
        Total_Cells,
        Cells_Expressing,
        Percent_Expressing,
        Mean_Expression
      )
  ),
  row.names = FALSE
)


# ============================================================
# 19. HIGHEST CORRECTED REFINED LOCALIZATION
# ============================================================

priority_refined_top <- priority_refined %>%
  
  group_by(
    Gene
  ) %>%
  
  slice_max(
    order_by = Percent_Expressing,
    n = 1,
    with_ties = FALSE
  ) %>%
  
  ungroup() %>%
  
  arrange(
    Gene
  )

write.csv(
  priority_refined_top,
  file.path(
    results_dir,
    "11D_D_PRIORITY_TOP_REFINED_LOCALIZATION.csv"
  ),
  row.names = FALSE
)

cat("\n============================================\n")
cat("CORRECTED PRIORITY-GENE REFINED LOCALIZATION\n")
cat("============================================\n")

print(
  as.data.frame(
    priority_refined_top %>%
      select(
        Gene,
        Corrected_Refined_Cell_Type,
        Total_Cells,
        Cells_Expressing,
        Percent_Expressing,
        Mean_Expression
      )
  ),
  row.names = FALSE
)


# ============================================================
# 20. CELL-TYPE COMPOSITION BY CONDITION
#
# Descriptive only.
# ============================================================

composition_broad <- scrna@meta.data %>%
  
  mutate(
    Condition =
      as.character(
        scrna$Condition
      )
  ) %>%
  
  count(
    Condition,
    corrected_broad_cell_type,
    name = "Cells"
  ) %>%
  
  group_by(
    Condition
  ) %>%
  
  mutate(
    Total_Cells_Condition =
      sum(Cells),
    
    Percent_of_Condition =
      100 * Cells /
      Total_Cells_Condition
  ) %>%
  
  ungroup()

write.csv(
  composition_broad,
  file.path(
    results_dir,
    "11D_D_BROAD_CELLTYPE_COMPOSITION_BY_CONDITION.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 21. REFINED COMPOSITION BY CONDITION
# ============================================================

composition_refined <- scrna@meta.data %>%
  
  mutate(
    Condition =
      as.character(
        scrna$Condition
      )
  ) %>%
  
  count(
    Condition,
    corrected_refined_cell_type,
    name = "Cells"
  ) %>%
  
  group_by(
    Condition
  ) %>%
  
  mutate(
    Total_Cells_Condition =
      sum(Cells),
    
    Percent_of_Condition =
      100 * Cells /
      Total_Cells_Condition
  ) %>%
  
  ungroup()

write.csv(
  composition_refined,
  file.path(
    results_dir,
    "11D_D_REFINED_CELLTYPE_COMPOSITION_BY_CONDITION.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 22. PRIORITY GENES:
# CONDITION × BROAD CELL TYPE
#
# This is the key within-cell-type analysis.
# ============================================================

priority_matrix <- rna_data[
  priority_genes,
  ,
  drop = FALSE
]

combo_broad <- paste(
  as.character(
    scrna$Condition
  ),
  as.character(
    scrna$corrected_broad_cell_type
  ),
  sep = "|||"
)

priority_condition_broad <-
  calculate_gene_group_stats(
    expression_matrix = priority_matrix,
    genes = priority_genes,
    groups = combo_broad
  )

priority_condition_broad <- priority_condition_broad %>%
  
  separate(
    Group,
    into = c(
      "Condition",
      "Corrected_Broad_Cell_Type"
    ),
    sep = "\\|\\|\\|",
    remove = TRUE
  )

priority_condition_broad$Condition <- factor(
  priority_condition_broad$Condition,
  levels = condition_levels
)

priority_condition_broad <-
  priority_condition_broad %>%
  
  arrange(
    Gene,
    Corrected_Broad_Cell_Type,
    Condition
  )

write.csv(
  priority_condition_broad,
  file.path(
    results_dir,
    "11D_D_PRIORITY_CONDITION_BY_BROAD_CELLTYPE.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 23. PRIORITY GENES:
# CONDITION × REFINED CELL TYPE