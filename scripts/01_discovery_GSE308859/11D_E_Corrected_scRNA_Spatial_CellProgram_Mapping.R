# ============================================================
# STEP 11D-E
# CORRECTED scRNA -> SPATIAL CELL-PROGRAM MAPPING
#
# PURPOSE
# ------------------------------------------------------------
# Build broad cell-program signatures from the FINAL corrected
# scRNA atlas and project them onto the frozen spatial dataset.
#
# CRITICAL PRINCIPLES
# ------------------------------------------------------------
# 1. Corrected annotations only.
# 2. Original scRNA clusters are unchanged.
# 3. Spatial QC is unchanged.
# 4. CYP / eicosanoid genes are excluded from mapping signatures.
# 5. Signatures are validated in scRNA before spatial projection.
# 6. Spatial scores are interpreted as PROGRAM SCORES, not
#    literal cell fractions.
# 7. No condition-level inferential statistics are performed,
#    because there is one deposited spatial library/condition.
#
# INPUTS
# ------------------------------------------------------------
# scRNA:
# 11D_C_FINAL_CORRECTED_ANNOTATED_ATLAS.rds
#
# Spatial:
# 11C_CYP_AA_ANALYZED_SPATIAL_OBJECTS.rds
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

scrna_file <- file.path(
  project_dir,
  "CLEAN DATA",
  "11D_C_FINAL_CORRECTED_ANNOTATED_ATLAS.rds"
)

spatial_file <- file.path(
  project_dir,
  "CLEAN DATA",
  "11C_CYP_AA_ANALYZED_SPATIAL_OBJECTS.rds"
)

results_dir <- file.path(
  project_dir,
  "RESULTS",
  "11_SPATIAL",
  "11D_E_CORRECTED_CELL_PROGRAM_MAPPING"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "11_SPATIAL",
  "11D_E_CORRECTED_CELL_PROGRAM_MAPPING"
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
# 3. LOAD scRNA ATLAS
# ============================================================

cat("\n============================================\n")
cat("STEP 11D-E\n")
cat("CORRECTED scRNA -> SPATIAL PROGRAM MAPPING\n")
cat("============================================\n")

scrna <- readRDS(scrna_file)

DefaultAssay(scrna) <- "RNA"

cat("\nscRNA atlas loaded.\n")
cat("Cells:", ncol(scrna), "\n")
cat("Genes:", nrow(scrna), "\n")


# ============================================================
# 4. VERIFY CORRECTED scRNA METADATA
# ============================================================

required_scrna_meta <- c(
  "seurat_clusters",
  "sample_id",
  "corrected_broad_cell_type",
  "corrected_refined_cell_type"
)

missing_scrna_meta <- setdiff(
  required_scrna_meta,
  colnames(scrna@meta.data)
)

if (length(missing_scrna_meta) > 0) {
  
  stop(
    paste(
      "Missing scRNA metadata:",
      paste(
        missing_scrna_meta,
        collapse = ", "
      )
    )
  )
}

if (ncol(scrna) != 25186) {
  stop("Unexpected scRNA cell count.")
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
  stop("Expected frozen clusters 0-28.")
}

cat("Corrected scRNA atlas verified.\n")


# ============================================================
# 5. LOAD SPATIAL OBJECT LIST
# ============================================================

spatial_list <- readRDS(
  spatial_file
)

if (!is.list(spatial_list)) {
  stop("Spatial input is expected to be a list of Seurat objects.")
}

cat("\nSpatial objects loaded:", length(spatial_list), "\n")

print(
  names(spatial_list)
)


# ============================================================
# 6. VERIFY SPATIAL OBJECTS
# ============================================================

spatial_counts_before <- sapply(
  spatial_list,
  ncol
)

cat("\n============================================\n")
cat("SPATIAL SPOT COUNTS\n")
cat("============================================\n")

print(
  spatial_counts_before
)

if (sum(spatial_counts_before) != 7079) {
  
  stop(
    paste(
      "Unexpected total spatial spot count:",
      sum(spatial_counts_before),
      "Expected 7079."
    )
  )
}


# ============================================================
# 7. STANDARDIZE SPATIAL CONDITION NAMES
#
# Use list names first.
# ============================================================

standardize_condition <- function(x) {
  
  x <- as.character(x)
  
  case_when(
    
    grepl(
      "sham",
      x,
      ignore.case = TRUE
    ) ~ "Sham",
    
    grepl(
      "2w|2_w|2week|2-week|tac2",
      x,
      ignore.case = TRUE
    ) ~ "2W",
    
    grepl(
      "4w|4_w|4week|4-week|tac4",
      x,
      ignore.case = TRUE
    ) ~ "4W",
    
    grepl(
      "6w|6_w|6week|6-week|tac6",
      x,
      ignore.case = TRUE
    ) ~ "6W",
    
    TRUE ~ NA_character_
  )
}


condition_levels <- c(
  "Sham",
  "2W",
  "4W",
  "6W"
)


# ============================================================
# 8. ASSIGN CONDITION TO EACH SPATIAL OBJECT
# ============================================================

spatial_conditions <- standardize_condition(
  names(spatial_list)
)

if (any(is.na(spatial_conditions))) {
  
  cat(
    "\nCould not infer all conditions from list names.\n"
  )
  
  cat(
    "List names:",
    paste(
      names(spatial_list),
      collapse = ", "
    ),
    "\n"
  )
  
  # Try common metadata columns
  spatial_conditions <- rep(
    NA_character_,
    length(spatial_list)
  )
  
  for (i in seq_along(spatial_list)) {
    
    obj <- spatial_list[[i]]
    
    candidate_cols <- c(
      "sample_id",
      "sample",
      "condition",
      "Condition",
      "orig.ident"
    )
    
    candidate_cols <- candidate_cols[
      candidate_cols %in%
        colnames(obj@meta.data)
    ]
    
    if (length(candidate_cols) > 0) {
      
      for (candidate in candidate_cols) {
        
        vals <- unique(
          as.character(
            obj@meta.data[[candidate]]
          )
        )
        
        vals <- vals[
          !is.na(vals)
        ]
        
        if (length(vals) > 0) {
          
          mapped <- standardize_condition(
            vals
          )
          
          mapped <- unique(
            mapped[
              !is.na(mapped)
            ]
          )
          
          if (length(mapped) == 1) {
            
            spatial_conditions[i] <- mapped
            
            break
          }
        }
      }
    }
  }
}

if (any(is.na(spatial_conditions))) {
  
  stop(
    paste(
      "Could not assign condition to spatial object(s):",
      paste(
        names(spatial_list)[
          is.na(spatial_conditions)
        ],
        collapse = ", "
      )
    )
  )
}

if (!setequal(
  spatial_conditions,
  condition_levels
)) {
  
  stop(
    paste(
      "Expected spatial conditions Sham, 2W, 4W, 6W; found:",
      paste(
        spatial_conditions,
        collapse = ", "
      )
    )
  )
}

cat("\n============================================\n")
cat("SPATIAL CONDITION MAPPING\n")
cat("============================================\n")

condition_map <- data.frame(
  
  Spatial_Object =
    names(spatial_list),
  
  Condition =
    spatial_conditions,
  
  Spots =
    spatial_counts_before,
  
  stringsAsFactors = FALSE
)

print(
  condition_map,
  row.names = FALSE
)

write.csv(
  condition_map,
  file.path(
    results_dir,
    "11D_E_Spatial_Condition_Mapping.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9. ADD CONDITION METADATA TO SPATIAL OBJECTS
# ============================================================

for (i in seq_along(spatial_list)) {
  
  spatial_list[[i]]$Condition <-
    factor(
      spatial_conditions[i],
      levels = condition_levels
    )
}


# ============================================================
# 10. DEFINE EXCLUDED PATHWAY GENES
#
# Signatures must NOT be driven by the pathway being tested.
# ============================================================

all_genes <- rownames(scrna)

cyp_excluded <- grep(
  "^Cyp[0-9]",
  all_genes,
  value = TRUE
)

eicosanoid_excluded_requested <- c(
  
  # Epoxide hydrolases
  "Ephx1",
  "Ephx2",
  
  # PLA2 / AA release
  "Pla2g4a",
  "Pla2g4b",
  "Pla2g4c",
  "Pla2g6",
  "Pla2g2a",
  "Pla2g5",
  
  # COX pathway
  "Ptgs1",
  "Ptgs2",
  "Ptges",
  "Ptges2",
  "Ptges3",
  "Ptgis",
  "Tbxas1",
  
  # LOX pathway
  "Alox5",
  "Alox5ap",
  "Alox12",
  "Alox12b",
  "Alox15",
  "Alox15b",
  
  # Leukotriene / prostaglandin metabolism
  "Lta4h",
  "Ltc4s",
  "Hpgds",
  "Hpgd",
  
  # Selected eicosanoid receptors
  "Ptger1",
  "Ptger2",
  "Ptger3",
  "Ptger4",
  "Ptgdr",
  "Ptgdr2",
  "Tbxa2r",
  "Cysltr1",
  "Cysltr2",
  "Ltb4r1",
  "Ltb4r2"
)

eicosanoid_excluded <- intersect(
  eicosanoid_excluded_requested,
  all_genes
)

excluded_signature_genes <- sort(
  unique(
    c(
      cyp_excluded,
      eicosanoid_excluded
    )
  )
)

cat(
  "\nCYP genes excluded from signatures:",
  length(cyp_excluded),
  "\n"
)

cat(
  "Additional eicosanoid genes excluded:",
  length(eicosanoid_excluded),
  "\n"
)

cat(
  "Total pathway genes excluded:",
  length(excluded_signature_genes),
  "\n"
)

write.csv(
  data.frame(
    Gene = excluded_signature_genes
  ),
  file.path(
    results_dir,
    "11D_E_Excluded_Pathway_Genes.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. SET CORRECTED BROAD IDENTITIES
# ============================================================

Idents(scrna) <- "corrected_broad_cell_type"

broad_types <- sort(
  unique(
    as.character(
      scrna$corrected_broad_cell_type
    )
  )
)

cat(
  "\nCorrected broad cell types:",
  length(broad_types),
  "\n"
)

print(
  broad_types
)

if (length(broad_types) != 15) {
  stop("Expected 15 corrected broad cell types.")
}


# ============================================================
# 12. FIND BROAD CELL-TYPE MARKERS
#
# Positive markers only.
#
# IMPORTANT:
# These p-values are used ONLY for marker ranking within this
# technical signature-building step, NOT for biological
# condition-level inference.
# ============================================================

cat("\n============================================\n")
cat("FINDING CORRECTED BROAD CELL-TYPE MARKERS\n")
cat("============================================\n")

broad_markers <- FindAllMarkers(
  object = scrna,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = TRUE
)

write.csv(
  broad_markers,
  file.path(
    results_dir,
    "11D_E_All_Broad_CellType_Markers_RAW.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. IDENTIFY MARKER FC COLUMN
# ============================================================

fc_candidates <- c(
  "avg_log2FC",
  "avg_logFC"
)

fc_col <- fc_candidates[
  fc_candidates %in%
    colnames(broad_markers)
]

if (length(fc_col) == 0) {
  stop("Could not identify marker fold-change column.")
}

fc_col <- fc_col[1]

cat(
  "Marker FC column:",
  fc_col,
  "\n"
)


# ============================================================
# 14. REMOVE CYP / EICOSANOID GENES FROM MARKERS
# ============================================================

broad_markers_clean <- broad_markers %>%
  
  filter(
    !gene %in% excluded_signature_genes
  )

write.csv(
  broad_markers_clean,
  file.path(
    results_dir,
    "11D_E_All_Broad_CellType_Markers_PathwayExcluded.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. REQUIRE SPATIAL GENE AVAILABILITY
#
# A mapping signature is only useful if the genes are present
# in the spatial assay.
# ============================================================

spatial_gene_sets <- lapply(
  spatial_list,
  rownames
)

shared_spatial_genes <- Reduce(
  intersect,
  spatial_gene_sets
)

shared_mapping_genes <- intersect(
  rownames(scrna),
  shared_spatial_genes
)

cat(
  "\nGenes shared by scRNA and all spatial objects:",
  length(shared_mapping_genes),
  "\n"
)

broad_markers_clean <- broad_markers_clean %>%
  
  filter(
    gene %in% shared_mapping_genes
  )


# ============================================================
# 16. BUILD INITIAL SIGNATURES
#
# Select up to 30 strongest markers/cell type.
# Require positive FC.
# ============================================================

broad_markers_clean <-
  broad_markers_clean %>%
  
  filter(
    .data[[fc_col]] > 0
  )

signature_table <- broad_markers_clean %>%
  
  group_by(
    cluster
  ) %>%
  
  arrange(
    desc(
      .data[[fc_col]]
    ),
    .by_group = TRUE
  ) %>%
  
  slice_head(
    n = 30
  ) %>%
  
  ungroup()

write.csv(
  signature_table,
  file.path(
    results_dir,
    "11D_E_Initial_Top30_Signature_Markers.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. CONVERT TO SIGNATURE LIST
# ============================================================

signature_list <- split(
  signature_table$gene,
  as.character(
    signature_table$cluster
  )
)

signature_list <- lapply(
  signature_list,
  unique
)

signature_sizes <- data.frame(
  
  Cell_Type =
    names(signature_list),
  
  Signature_Genes =
    lengths(signature_list),
  
  stringsAsFactors = FALSE
)

cat("\n============================================\n")
cat("INITIAL SIGNATURE SIZES\n")
cat("============================================\n")

print(
  signature_sizes,
  row.names = FALSE
)

write.csv(
  signature_sizes,
  file.path(
    results_dir,
    "11D_E_Initial_Signature_Sizes.csv"
  ),
  row.names = FALSE
)

if (any(signature_sizes$Signature_Genes < 5)) {
  
  warning(
    paste(
      "Some cell types have fewer than 5 signature genes:",
      paste(
        signature_sizes$Cell_Type[
          signature_sizes$Signature_Genes < 5
        ],
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 18. CHECK CROSS-SIGNATURE DUPLICATION
#
# Genes occurring in many signatures are less discriminative.
# ============================================================

signature_membership <- bind_rows(
  lapply(
    names(signature_list),
    function(ct) {
      
      data.frame(
        Cell_Type = ct,
        Gene = signature_list[[ct]],
        stringsAsFactors = FALSE
      )
    }
  )
)

gene_signature_frequency <- signature_membership %>%
  
  count(
    Gene,
    name = "Number_of_Signatures"
  ) %>%
  
  arrange(
    desc(Number_of_Signatures),
    Gene
  )

write.csv(
  gene_signature_frequency,
  file.path(
    results_dir,
    "11D_E_Signature_Gene_Sharing.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. REMOVE MARKERS SHARED BY >2 SIGNATURES
#
# Conservative specificity filter.
# ============================================================

genes_too_shared <- gene_signature_frequency %>%
  
  filter(
    Number_of_Signatures > 2
  ) %>%
  
  pull(
    Gene
  )

signature_list_specific <- lapply(
  signature_list,
  function(x) {
    
    setdiff(
      x,
      genes_too_shared
    )
  }
)


# ============================================================
# 20. ENSURE SIGNATURE SIZE AFTER SPECIFICITY FILTER
#
# If a population becomes too small, retain its original
# signature rather than creating an unusable program.
# ============================================================

for (ct in names(signature_list_specific)) {
  
  if (length(signature_list_specific[[ct]]) < 5) {
    
    warning(
      paste(
        ct,
        "had <5 genes after specificity filtering;",
        "retaining original pathway-excluded signature."
      )
    )
    
    signature_list_specific[[ct]] <-
      signature_list[[ct]]
  }
}

signature_list <- signature_list_specific

final_signature_sizes <- data.frame(
  
  Cell_Type =
    names(signature_list),
  
  Signature_Genes =
    lengths(signature_list),
  
  stringsAsFactors = FALSE
)

write.csv(
  final_signature_sizes,
  file.path(
    results_dir,
    "11D_E_Final_Signature_Sizes.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 21. SAVE FINAL SIGNATURE MEMBERSHIP
# ============================================================

final_signature_membership <- bind_rows(
  lapply(
    names(signature_list),
    function(ct) {
      
      data.frame(
        
        Cell_Type = ct,
        
        Gene =
          signature_list[[ct]],
        
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  final_signature_membership,
  file.path(
    results_dir,
    "11D_E_FINAL_BROAD_CELL_PROGRAM_SIGNATURES.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 22. PRINT FINAL SIGNATURES
# ============================================================

cat("\n============================================\n")
cat("FINAL BROAD CELL-PROGRAM SIGNATURES\n")
cat("============================================\n")

for (ct in names(signature_list)) {
  
  cat("\n", ct, "\n", sep = "")
  
  cat(
    paste(
      signature_list[[ct]],
      collapse = ", "
    ),
    "\n"
  )
}


# ============================================================
# 23. VALIDATE SIGNATURES IN scRNA
#
# Use mean expression of signature genes rather than
# AddModuleScore to keep scoring directly reproducible across
# scRNA and spatial modalities.
# ============================================================

scrna_data <- GetAssayData(
  scrna,
  assay = "RNA",
  layer = "data"
)

scrna_score_matrix <- matrix(
  NA_real_,
  nrow = ncol(scrna),
  ncol = length(signature_list)
)

rownames(scrna_score_matrix) <- colnames(scrna)
colnames(scrna_score_matrix) <- names(signature_list)

for (ct in names(signature_list)) {
  
  genes_now <- intersect(
    signature_list[[ct]],
    rownames(scrna_data)
  )
  
  scrna_score_matrix[, ct] <-
    as.numeric(
      Matrix::colMeans(
        scrna_data[
          genes_now,
          ,
          drop = FALSE
        ]
      )
    )
}


# ============================================================
# 24. scRNA SIGNATURE VALIDATION SUMMARY
#
# Expected score = score of each cell type's own signature.
# Compare expected population with all other cells.
# ============================================================

validation_list <- list()
