# ============================================================
# V08 — PREDEFINED CYP/AA INDEPENDENT VALIDATION
# GSE155882
# ============================================================
#
# INPUT:
#   V07_GSE155882_FINAL_ANNOTATED_ATLAS.rds
#
# PURPOSE:
#   Independently evaluate PREDEFINED CYP/AA findings from
#   GSE308859 in GSE155882.
#
# DESIGN PRINCIPLES
# -----------------
# 1. Validation targets were predefined before V08.
# 2. Biological replicate = sequencing library, NOT cell.
# 3. Sham n = 2 libraries; TAC n = 2 libraries.
# 4. No cell-level inferential Sham-vs-TAC p-values.
# 5. Replicate-level pseudobulk is used for expression summaries.
# 6. Detection is also summarized independently per library.
# 7. Cluster 10 remains separate and cannot support a
#    fibroblast- or myeloid-specific validation claim.
# 8. Cyp4f17 is expected to be poorly evaluable because of
#    extreme sparsity established before V08.
# 9. No genes are added/removed based on favorable results.
# 10. No reclustering, integration, filtering, or annotation.
#
# IMPORTANT:
#   This step evaluates expression/detection/localization.
#   It does NOT measure CYP enzymatic activity, EET/HETE levels,
#   AA flux, or metabolite production.
#
# ============================================================


# ============================================================
# 1. CLEAN WORKSPACE
# ============================================================

rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  max.print = 100000
)

set.seed(20260921)


# ============================================================
# 2. PACKAGES
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})


# ============================================================
# 3. PATHS
# ============================================================

project_root <- "D:/Master/ScRNA seq/TAC/VALIDATION/GSE155882"

input_file <- file.path(
  project_root,
  "CLEAN DATA",
  "V07_GSE155882_FINAL_ANNOTATED_ATLAS.rds"
)

results_dir <- file.path(
  project_root,
  "RESULTS",
  "V08_PREDEFINED_CYP_AA_VALIDATION"
)

figures_dir <- file.path(
  project_root,
  "FIGURES",
  "V08_PREDEFINED_CYP_AA_VALIDATION"
)

png_dir <- file.path(
  figures_dir,
  "PNG"
)

output_file <- file.path(
  project_root,
  "CLEAN DATA",
  "V08_GSE155882_CYP_AA_VALIDATION_ATLAS.rds"
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  png_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 4. LOAD FROZEN V07 ATLAS
# ============================================================

if (!file.exists(input_file)) {
  stop("V07 final annotated atlas not found.")
}

obj <- readRDS(input_file)

input_n <- ncol(obj)
input_cells <- colnames(obj)

cat("\n============================================\n")
cat("V08 — PREDEFINED CYP/AA INDEPENDENT VALIDATION\n")
cat("============================================\n")

cat("\nCells:", input_n, "\n")
cat("Genes:", nrow(obj), "\n")


# ============================================================
# 5. INPUT INTEGRITY
# ============================================================

required_metadata <- c(
  "sample_id",
  "condition",
  "replicate",
  "scDblFinder_class",
  "V07_primary_cluster",
  "V07_broad_cell_type",
  "V07_refined_cell_type",
  "V07_annotation_confidence",
  "V07_annotation_status",
  "V07_annotation_version"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(obj@meta.data)
)

if (length(missing_metadata) > 0) {
  stop(
    paste(
      "Missing required metadata:",
      paste(missing_metadata, collapse = ", ")
    )
  )
}

if (input_n != 28388) {
  stop("Unexpected V07 cell count.")
}

if (!all(as.character(obj$scDblFinder_class) == "singlet")) {
  stop("Non-singlet cells detected.")
}

if (!all(
  obj$V07_annotation_version ==
  "GSE155882_Validation_Atlas_v1_Final"
)) {
  stop("Unexpected V07 annotation version.")
}

if (!"RNA" %in% Assays(obj)) {
  stop("RNA assay missing.")
}


# ============================================================
# 6. JOIN RNA LAYERS
# ============================================================

DefaultAssay(obj) <- "RNA"

obj <- JoinLayers(
  obj,
  assay = "RNA"
)

counts_mat <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "counts"
)

data_mat <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "data"
)


# ============================================================
# 7. PREDEFINED PRIORITY GENES
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

priority_available <- priority_genes[
  priority_genes %in% rownames(obj)
]

priority_missing <- setdiff(
  priority_genes,
  rownames(obj)
)

cat("\nPredefined priority genes:\n")
print(priority_genes)

cat("\nAvailable priority genes:\n")
print(priority_available)

cat("\nMissing priority genes:\n")
print(priority_missing)

if (length(priority_missing) > 0) {
  stop(
    paste(
      "Predefined priority genes missing:",
      paste(priority_missing, collapse = ", ")
    )
  )
}


# ============================================================
# 8. LOCK EXPECTED DISCOVERY LOCALIZATION
# ============================================================
#
# These expectations come from the corrected GSE308859
# discovery atlas and are declared BEFORE V08 results.
#
# They are NOT derived from GSE155882.
#
# ============================================================

discovery_expectations <- data.frame(
  
  Gene = priority_genes,
  
  Discovery_Primary_Context = c(
    "Fibroblast",
    "Fibroblast/Endothelial",
    "Fibroblast/Endothelial",
    "Resident macrophage/Smooth muscle",
    "Monocyte",
    "Endothelial",
    "Monocyte/Neutrophil",
    "Endothelial/Fibroblast"
  ),
  
  Discovery_Key_Observation = c(
    "Strong fibroblast contribution; marked temporal remodeling",
    "Predominantly fibroblast contribution with temporal remodeling",
    "Fibroblast contribution with endothelial redistribution at 4W",
    "Macrophage/smooth-muscle-associated localization and remodeling",
    "Strong monocyte-associated expression",
    "Low-expression endothelial-associated signal",
    "Strong immune redistribution, especially monocyte/neutrophil compartments",
    "Major endothelial-associated signal with marked 4W loss and later fibroblast redistribution"
  ),
  
  stringsAsFactors = FALSE
)

write.csv(
  discovery_expectations,
  file.path(
    results_dir,
    "V08_Predefined_Discovery_Expectations.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9. SAMPLE AUDIT
# ============================================================

sample_audit <- obj@meta.data %>%
  transmute(
    sample_id = as.character(sample_id),
    condition = as.character(condition),
    replicate = as.character(replicate)
  ) %>%
  distinct() %>%
  arrange(
    condition,
    sample_id
  )

sample_cell_counts <- obj@meta.data %>%
  transmute(
    sample_id = as.character(sample_id),
    condition = as.character(condition)
  ) %>%
  count(
    sample_id,
    condition,
    name = "Cells"
  )

sample_audit <- left_join(
  sample_audit,
  sample_cell_counts,
  by = c(
    "sample_id",
    "condition"
  )
)

write.csv(
  sample_audit,
  file.path(
    results_dir,
    "V08_Sample_Audit.csv"
  ),
  row.names = FALSE
)

cat("\nSample audit:\n")
print(sample_audit, row.names = FALSE)


# ============================================================
# 10. WHOLE-ATLAS DETECTION — EACH LIBRARY SEPARATELY
# ============================================================

sample_ids <- unique(
  as.character(obj$sample_id)
)

whole_detection_list <- list()
counter <- 1

for (sid in sample_ids) {
  
  cells_sid <- colnames(obj)[
    as.character(obj$sample_id) == sid
  ]
  
  condition_sid <- unique(
    as.character(
      obj$condition[
        as.character(obj$sample_id) == sid
      ]
    )
  )
  
  mat <- counts_mat[
    priority_genes,
    cells_sid,
    drop = FALSE
  ]
  
  positive <- Matrix::rowSums(
    mat > 0
  )
  
  whole_detection_list[[counter]] <- data.frame(
    Sample = sid,
    Condition = condition_sid,
    Gene = priority_genes,
    Cells = length(cells_sid),
    Positive_Cells = as.numeric(positive),
    Percent_Detected =
      100 * as.numeric(positive) / length(cells_sid),
    stringsAsFactors = FALSE
  )
  
  counter <- counter + 1
}

whole_detection <- bind_rows(
  whole_detection_list
)

write.csv(
  whole_detection,
  file.path(
    results_dir,
    "V08_Whole_Atlas_Priority_Gene_Detection_by_Library.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. WHOLE-ATLAS NORMALIZED EXPRESSION BY LIBRARY
# ============================================================

whole_expression_list <- list()
counter <- 1

for (sid in sample_ids) {
  
  cells_sid <- colnames(obj)[
    as.character(obj$sample_id) == sid
  ]
  
  condition_sid <- unique(
    as.character(
      obj$condition[
        as.character(obj$sample_id) == sid
      ]
    )
  )
  
  mat <- data_mat[
    priority_genes,
    cells_sid,
    drop = FALSE
  ]
  
  avg_expr <- Matrix::rowMeans(
    mat
  )
  
  whole_expression_list[[counter]] <- data.frame(
    Sample = sid,
    Condition = condition_sid,
    Gene = priority_genes,
    Mean_Normalized_Expression =
      as.numeric(avg_expr),
    stringsAsFactors = FALSE
  )
  
  counter <- counter + 1
}

whole_expression <- bind_rows(
  whole_expression_list
)

write.csv(
  whole_expression,
  file.path(
    results_dir,
    "V08_Whole_Atlas_Normalized_Expression_by_Library.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. WHOLE-ATLAS RAW-COUNT PSEUDOBULK
# ============================================================
#
# Counts are summed per biological library.
#
# We calculate CPM for descriptive replicate-aware comparison.
# No cell-level statistical testing.
#
# ============================================================

pseudobulk_counts_list <- list()
counter <- 1

for (sid in sample_ids) {
  
  cells_sid <- colnames(obj)[
    as.character(obj$sample_id) == sid
  ]
  
  condition_sid <- unique(
    as.character(
      obj$condition[
        as.character(obj$sample_id) == sid
      ]
    )
  )
  
  total_counts <- Matrix::rowSums(
    counts_mat[, cells_sid, drop = FALSE]
  )
  
  library_size <- sum(total_counts)
  
  priority_counts <- total_counts[
    priority_genes
  ]
  
  priority_cpm <-
    1e6 * priority_counts / library_size
  
  pseudobulk_counts_list[[counter]] <- data.frame(
    Sample = sid,
    Condition = condition_sid,
    Gene = priority_genes,
    Raw_Pseudobulk_Count =
      as.numeric(priority_counts),
    Library_Size =
      as.numeric(library_size),
    CPM =
      as.numeric(priority_cpm),
    Log2_CPM =
      log2(as.numeric(priority_cpm) + 1),
    stringsAsFactors = FALSE
  )
  
  counter <- counter + 1
}

whole_pseudobulk <- bind_rows(
  pseudobulk_counts_list
)

write.csv(
  whole_pseudobulk,
  file.path(
    results_dir,
    "V08_Whole_Atlas_Pseudobulk_by_Library.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. CONDITION-LEVEL DESCRIPTIVE EFFECTS
# ============================================================
#
# With only 2 libraries per condition:
# report effect size + replicate values,
# NOT inferential p-values.
#
# ============================================================

condition_pb_summary <- whole_pseudobulk %>%
  group_by(
    Gene,
    Condition
  ) %>%
  summarise(
    N_Libraries = n(),
    Mean_CPM = mean(CPM),
    Median_CPM = median(CPM),
    Mean_Log2_CPM = mean(Log2_CPM),
    .groups = "drop"
  )

condition_pb_wide <- condition_pb_summary %>%
  select(
    Gene,
    Condition,
    Mean_CPM
  ) %>%
  pivot_wider(
    names_from = Condition,
    values_from = Mean_CPM
  )

if (!all(c("Sham", "TAC") %in% colnames(condition_pb_wide))) {
  stop(
    paste(
      "Expected condition columns Sham and TAC.",
      "Found:",
      paste(colnames(condition_pb_wide), collapse = ", ")
    )
  )
}

condition_pb_effect <- condition_pb_wide %>%
  mutate(
    TAC_vs_Sham_Log2FC =
      log2(
        (TAC + 0.01) /
          (Sham + 0.01)
      ),
    Direction =
      case_when(
        TAC_vs_Sham_Log2FC > 0 ~ "Higher_in_TAC",
        TAC_vs_Sham_Log2FC < 0 ~ "Lower_in_TAC",
        TRUE ~ "No_change"
      )
  )

write.csv(
  condition_pb_effect,
  file.path(
    results_dir,
    "V08_Whole_Atlas_TAC_vs_Sham_Pseudobulk_Effect.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. REPLICATE CONCORDANCE
# ============================================================
#
# Direction is assessed independently in replicate pairs:
#
#   TAC_Rep1 vs Sham_Rep1
#   TAC_Rep2 vs Sham_Rep2
#
# This pairing is used ONLY as a simple concordance diagnostic.
# It is not assumed to represent experimental matching.
#
# ============================================================

pb_pair <- whole_pseudobulk %>%
  select(
    Sample,
    Gene,
    CPM
  ) %>%
  pivot_wider(
    names_from = Sample,
    values_from = CPM
  )

required_sample_columns <- c(
  "Sham_Rep1",
  "Sham_Rep2",
  "TAC_Rep1",
  "TAC_Rep2"
)

if (!all(required_sample_columns %in% colnames(pb_pair))) {
  stop(
    paste(
      "Missing expected sample columns:",
      paste(
        setdiff(
          required_sample_columns,
          colnames(pb_pair)
        ),
        collapse = ", "
      )
    )
  )
}

replicate_concordance <- pb_pair %>%
  mutate(
    Pair1_Log2FC =
      log2(
        (TAC_Rep1 + 0.01) /
          (Sham_Rep1 + 0.01)
      ),
    
    Pair2_Log2FC =
      log2(
        (TAC_Rep2 + 0.01) /
          (Sham_Rep2 + 0.01)
      ),
    
    Pair1_Direction =
      case_when(
        Pair1_Log2FC > 0 ~ "Higher_in_TAC",
        Pair1_Log2FC < 0 ~ "Lower_in_TAC",
        TRUE ~ "No_change"
      ),
    
    Pair2_Direction =
      case_when(
        Pair2_Log2FC > 0 ~ "Higher_in_TAC",
        Pair2_Log2FC < 0 ~ "Lower_in_TAC",
        TRUE ~ "No_change"
      ),
    
    Direction_Concordant =
      Pair1_Direction == Pair2_Direction
  )

write.csv(
  replicate_concordance,
  file.path(
    results_dir,
    "V08_Whole_Atlas_Replicate_Direction_Concordance.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. DEFINE VALIDATION CELL-TYPE GROUPS
# ============================================================
#
# We retain the frozen V07 broad annotation.
#
# IMPORTANT:
# "Ambiguous/Mixed" remains its own group.
#
# ============================================================

broad_types <- unique(
  as.character(obj$V07_broad_cell_type)
)

broad_types <- sort(
  broad_types
)

cat("\nBroad validation groups:\n")
print(broad_types)


# ============================================================
# 16. CELL COUNTS BY LIBRARY × BROAD CELL TYPE
# ============================================================

celltype_sample_counts <- obj@meta.data %>%
  transmute(
    Sample =
      as.character(sample_id),
    Condition =
      as.character(condition),
    Broad_Cell_Type =
      as.character(V07_broad_cell_type)
  ) %>%
  count(
    Sample,
    Condition,
    Broad_Cell_Type,
    name = "Cells"
  )

write.csv(
  celltype_sample_counts,
  file.path(
    results_dir,
    "V08_Cell_Counts_by_Library_and_Broad_Cell_Type.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. CELL-TYPE × LIBRARY DETECTION
# ============================================================

ct_detection_list <- list()
counter <- 1

for (sid in sample_ids) {
  
  condition_sid <- unique(
    as.character(
      obj$condition[
        as.character(obj$sample_id) == sid
      ]
    )
  )
  
  for (ct in broad_types) {
    
    cells_ct <- colnames(obj)[
      as.character(obj$sample_id) == sid &
        as.character(obj$V07_broad_cell_type) == ct
    ]
    
    n_cells <- length(cells_ct)
    
    if (n_cells == 0) {
      next
    }
    
    mat <- counts_mat[
      priority_genes,
      cells_ct,
      drop = FALSE
    ]
    
    positive <- Matrix::rowSums(
      mat > 0
    )
    
    ct_detection_list[[counter]] <- data.frame(
      Sample = sid,
      Condition = condition_sid,
      Broad_Cell_Type = ct,
      Gene = priority_genes,
      Cells = n_cells,
      Positive_Cells =
        as.numeric(positive),
      Percent_Detected =
        100 * as.numeric(positive) / n_cells,
      stringsAsFactors = FALSE
    )
    
    counter <- counter + 1
  }
}

ct_detection <- bind_rows(
  ct_detection_list
)

write.csv(
  ct_detection,
  file.path(
    results_dir,
    "V08_Priority_Gene_Detection_by_CellType_and_Library.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. CELL-TYPE × LIBRARY NORMALIZED EXPRESSION
# ============================================================

ct_expression_list <- list()
counter <- 1

for (sid in sample_ids) {
  
  condition_sid <- unique(
    as.character(
      obj$condition[
        as.character(obj$sample_id) == sid
      ]
    )
  )
  
  for (ct in broad_types) {
    
    cells_ct <- colnames(obj)[
      as.character(obj$sample_id) == sid &
        as.character(obj$V07_broad_cell_type) == ct
    ]
    
    if (length(cells_ct) == 0) {
      next
    }
    
    mat <- data_mat[
      priority_genes,
      cells_ct,
      drop = FALSE
    ]
    
    avg <- Matrix::rowMeans(
      mat
    )
    
    ct_expression_list[[counter]] <- data.frame(
      Sample = sid,
      Condition = condition_sid,
      Broad_Cell_Type = ct,
      Gene = priority_genes,
      Cells = length(cells_ct),
      Mean_Normalized_Expression =
        as.numeric(avg),
      stringsAsFactors = FALSE
    )
    
    counter <- counter + 1
  }
}

ct_expression <- bind_rows(
  ct_expression_list
)

write.csv(
  ct_expression,
  file.path(
    results_dir,
    "V08_Priority_Gene_Normalized_Expression_by_CellType_and_Library.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. CELL-TYPE PSEUDOBULK
# ============================================================
#
# Sum raw counts independently for each:
#
#   library × broad cell type
#
# CPM denominator = all genes within that cell-type pseudobulk.
#
# ============================================================

ct_pb_list <- list()
counter <- 1

for (sid in sample_ids) {
  
  condition_sid <- unique(
    as.character(
      obj$condition[
        as.character(obj$sample_id) == sid
      ]
    )
  )
  
  for (ct in broad_types) {
    
    cells_ct <- colnames(obj)[
      as.character(obj$sample_id) == sid &
        as.character(obj$V07_broad_cell_type) == ct
    ]
    
    n_cells <- length(cells_ct)
    
    if (n_cells == 0) {
      next
    }
    
    summed_counts <- Matrix::rowSums(
      counts_mat[
        ,
        cells_ct,
        drop = FALSE
      ]
    )
    
    library_size <- sum(
      summed_counts
    )
    
    priority_counts <-
      summed_counts[
        priority_genes
      ]
    
    priority_cpm <-
      1e6 *
      priority_counts /
      library_size
    
    ct_pb_list[[counter]] <- data.frame(
      Sample = sid,
      Condition = condition_sid,
      Broad_Cell_Type = ct,
      Cells = n_cells,
      Gene = priority_genes,
      Raw_Pseudobulk_Count =
        as.numeric(priority_counts),
      Pseudobulk_Library_Size =
        as.numeric(library_size),
      CPM =
        as.numeric(priority_cpm),
      Log2_CPM =
        log2(
          as.numeric(priority_cpm) + 1
        ),
      stringsAsFactors = FALSE
    )
    
    counter <- counter + 1
  }
}

ct_pseudobulk <- bind_rows(
  ct_pb_list
)

write.csv(
  ct_pseudobulk,
  file.path(
    results_dir,
    "V08_CellType_Pseudobulk_by_Library.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 20. PSEUDOBULK EVALUABILITY RULES
# ============================================================
#
# Conservative thresholds for interpretation:
#
#   >= 50 cells in each available library
#   AND
#   >= 10 total positive cells across all four libraries
#
# These rules are for whether we INTERPRET a cell-type/gene