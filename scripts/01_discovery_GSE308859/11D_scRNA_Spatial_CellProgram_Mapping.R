# ============================================================
# STEP 11D
# scRNA -> SPATIAL CELL-PROGRAM MAPPING
# AND CYP-AA SPATIAL NICHE ANALYSIS
#
# PART 1:
#   1. Load frozen scRNA atlas
#   2. Load frozen Step 11C spatial objects
#   3. Verify object integrity
#   4. Verify frozen annotation
#   5. Derive broad cell-type marker signatures
#   6. Remove problematic marker classes
#   7. Rank/select robust markers
#   8. Export signatures for review
#
# IMPORTANT:
# - NO scRNA reclustering
# - NO annotation changes
# - NO spatial re-QC
# - NO spot removal
# - NO CYP-based marker selection
# - NO condition-level inferential statistics
# - Spatial spots are NOT biological replicates
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(12345)

cat("\n============================================\n")
cat("STEP 11D — PART 1\n")
cat("scRNA -> SPATIAL CELL-PROGRAM MAPPING\n")
cat("============================================\n")


# ============================================================
# 2. PATHS
# ============================================================

project_dir <- "D:/Master/ScRNA seq/TAC"

scrna_file <- file.path(
  project_dir,
  "CLEAN DATA",
  "07C_FINAL_FROZEN_ANNOTATED_ATLAS.rds"
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
  "11D_CELL_PROGRAM_MAPPING"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "11_SPATIAL",
  "11D_CELL_PROGRAM_MAPPING"
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

if (!file.exists(scrna_file)) {
  stop(
    paste0(
      "Frozen scRNA atlas not found:\n",
      scrna_file
    )
  )
}

if (!file.exists(spatial_file)) {
  stop(
    paste0(
      "Step 11C spatial object not found:\n",
      spatial_file
    )
  )
}


# ============================================================
# 3. LOAD OBJECTS
# ============================================================

cat("\nLoading frozen scRNA atlas...\n")

scrna <- readRDS(
  scrna_file
)

cat("Loading Step 11C spatial objects...\n")

spatial_objects <- readRDS(
  spatial_file
)

sample_order <- c(
  "Sham",
  "TAC_2W",
  "TAC_4W",
  "TAC_6W"
)

if (!all(sample_order %in% names(spatial_objects))) {
  stop(
    paste0(
      "Spatial object list does not contain all expected samples.\n",
      "Observed: ",
      paste(names(spatial_objects), collapse = ", ")
    )
  )
}

spatial_objects <- spatial_objects[sample_order]


# ============================================================
# 4. VERIFY FROZEN scRNA ATLAS
# ============================================================

cat("\n============================================\n")
cat("VERIFYING FROZEN scRNA ATLAS\n")
cat("============================================\n")

cat(
  "\nscRNA cells:",
  ncol(scrna),
  "\n"
)

cat(
  "scRNA genes:",
  nrow(scrna),
  "\n"
)

if (ncol(scrna) != 25186) {
  stop(
    paste0(
      "Frozen atlas cell count changed. Expected 25186; observed ",
      ncol(scrna),
      "."
    )
  )
}

if (nrow(scrna) != 32285) {
  stop(
    paste0(
      "Frozen atlas gene count changed. Expected 32285; observed ",
      nrow(scrna),
      "."
    )
  )
}

required_metadata <- c(
  "sample_id",
  "broad_cell_type",
  "refined_cell_type",
  "annotation_confidence"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(scrna@meta.data)
)

if (length(missing_metadata) > 0) {
  stop(
    paste0(
      "Missing frozen atlas metadata: ",
      paste(
        missing_metadata,
        collapse = ", "
      )
    )
  )
}

cat("\nFrozen refined cell types:\n")

refined_counts <- table(
  scrna$refined_cell_type
)

print(
  refined_counts
)


# ============================================================
# 5. VERIFY FROZEN SPATIAL OBJECTS
# ============================================================

cat("\n============================================\n")
cat("VERIFYING SPATIAL OBJECTS\n")
cat("============================================\n")

expected_spots <- c(
  Sham = 1657,
  TAC_2W = 2040,
  TAC_4W = 1212,
  TAC_6W = 2170
)

spatial_integrity_list <- list()

for (sid in sample_order) {
  
  sobj <- spatial_objects[[sid]]
  
  assay_names <- names(sobj@assays)
  
  if (!any(assay_names == "Spatial")) {
    stop(
      paste0(
        "Spatial assay missing in ",
        sid
      )
    )
  }
  
  DefaultAssay(sobj) <- "Spatial"
  
  counts_mat <- GetAssayData(
    object = sobj,
    assay = "Spatial",
    layer = "counts"
  )
  
  data_mat <- GetAssayData(
    object = sobj,
    assay = "Spatial",
    layer = "data"
  )
  
  if (ncol(sobj) != expected_spots[[sid]]) {
    stop(
      paste0(
        sid,
        ": expected ",
        expected_spots[[sid]],
        " spots but found ",
        ncol(sobj),
        "."
      )
    )
  }
  
  if (!identical(colnames(counts_mat), colnames(data_mat))) {
    stop(
      paste0(
        sid,
        ": counts/data spot names are not identical."
      )
    )
  }
  
  if (!identical(colnames(sobj), colnames(counts_mat))) {
    stop(
      paste0(
        sid,
        ": object/counts spot names are not identical."
      )
    )
  }
  
  spatial_integrity_list[[sid]] <- data.frame(
    Sample = sid,
    Spots = ncol(sobj),
    Genes = nrow(sobj),
    Counts_Data_Spots_Identical =
      identical(
        colnames(counts_mat),
        colnames(data_mat)
      ),
    Object_Counts_Spots_Identical =
      identical(
        colnames(sobj),
        colnames(counts_mat)
      ),
    stringsAsFactors = FALSE
  )
}

spatial_integrity <- bind_rows(
  spatial_integrity_list
)

write.csv(
  spatial_integrity,
  file.path(
    results_dir,
    "11D_Spatial_Object_Integrity.csv"
  ),
  row.names = FALSE
)

print(
  as.data.frame(spatial_integrity),
  row.names = FALSE
)


# ============================================================
# 6. DEFINE BROAD CELL TYPES FOR SPATIAL MAPPING
#
# We use broad types first because Visium spots contain
# multiple cells and broad signatures are more robust.
#
# Refined immune subtypes will be analyzed later as a
# secondary layer.
# ============================================================

cat("\n============================================\n")
cat("PREPARING BROAD CELL-TYPE SIGNATURES\n")
cat("============================================\n")

broad_counts <- scrna@meta.data %>%
  count(
    broad_cell_type,
    name = "Cells"
  ) %>%
  arrange(
    desc(Cells)
  )

print(
  as.data.frame(broad_counts),
  row.names = FALSE
)

write.csv(
  broad_counts,
  file.path(
    results_dir,
    "11D_Broad_Cell_Type_Counts.csv"
  ),
  row.names = FALSE
)

# Require enough cells to derive a stable reference signature.

min_cells_for_signature <- 100

eligible_cell_types <- broad_counts %>%
  filter(
    Cells >= min_cells_for_signature
  ) %>%
  pull(
    broad_cell_type
  )

cat(
  "\nCell types eligible for broad signature derivation:\n"
)

print(
  eligible_cell_types
)

if (length(eligible_cell_types) < 5) {
  stop(
    "Too few broad cell types passed the signature-size criterion."
  )
}


# ============================================================
# 7. DERIVE scRNA CELL-TYPE MARKERS
#
# IMPORTANT:
# These markers are derived from the frozen scRNA atlas.
# They are NOT selected using CYP/eicosanoid behavior.
#
# FindAllMarkers is used here to identify cell-identity
# signatures, NOT to make TAC-vs-Sham biological claims.
# ============================================================

cat("\n============================================\n")
cat("DERIVING BROAD CELL-TYPE MARKERS\n")
cat("============================================\n")

DefaultAssay(scrna) <- "RNA"

# Use the frozen broad annotation as identity.

Idents(scrna) <- scrna$broad_cell_type

# Restrict marker derivation to eligible broad populations.

cells_for_markers <- WhichCells(
  scrna,
  idents = eligible_cell_types
)

scrna_marker_reference <- subset(
  scrna,
  cells = cells_for_markers
)

Idents(scrna_marker_reference) <-
  scrna_marker_reference$broad_cell_type

cat(
  "\nCells used for marker derivation:",
  ncol(scrna_marker_reference),
  "\n"
)

cat(
  "Cell types used:",
  length(unique(Idents(scrna_marker_reference))),
  "\n"
)

# ------------------------------------------------------------
# Marker detection
#
# Only positive markers.
# min.pct = 0.20 helps avoid extremely sparse genes.
# logfc.threshold = 0.25 avoids weak identity markers.
#
# These p-values are used only as marker-ranking information
# within the reference atlas, not as TAC condition inference.
# ------------------------------------------------------------

broad_markers_all <- FindAllMarkers(
  object = scrna_marker_reference,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.20,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  return.thresh = 0.05,
  verbose = TRUE
)

if (nrow(broad_markers_all) == 0) {
  stop(
    "FindAllMarkers returned zero broad cell-type markers."
  )
}

write.csv(
  broad_markers_all,
  file.path(
    results_dir,
    "11D_All_Broad_CellType_Markers_Unfiltered.csv"
  ),
  row.names = FALSE
)

cat(
  "\nTotal unfiltered markers:",
  nrow(broad_markers_all),
  "\n"
)


# ============================================================
# 8. BUILD ROBUST SPATIAL-MAPPING SIGNATURES
#
# Remove genes likely to reflect technical state rather than
# cell identity.
#
# Also remove CYP/eicosanoid genes from the signatures.
# This is critical:
#
# We later test association between cell programs and CYP-AA.
# Therefore CYP-AA genes must NOT define the cell-type scores.
# ============================================================

cat("\n============================================\n")
cat("BUILDING ROBUST MAPPING SIGNATURES\n")
cat("============================================\n")


# ------------------------------------------------------------
# 8A. Define genes excluded from cell-identity signatures
# ------------------------------------------------------------

all_genes <- rownames(scrna_marker_reference)

mitochondrial_genes <- grep(
  "^mt-",
  all_genes,
  value = TRUE
)

ribosomal_genes <- grep(
  "^Rp[sl]",
  all_genes,
  value = TRUE
)

hemoglobin_genes <- grep(
  "^Hb[ab]",
  all_genes,
  value = TRUE
)

cell_cycle_genes <- unique(
  c(
    Seurat::cc.genes.updated.2019$s.genes,
    Seurat::cc.genes.updated.2019$g2m.genes
  )
)

# Strict CYP family exclusion

cyp_genes <- grep(
  "^Cyp[0-9]",
  all_genes,
  value = TRUE
)

# Eicosanoid-related genes excluded from identity signatures
# to prevent circularity in downstream CYP-AA niche analysis.

eicosanoid_exclusion <- c(
  "Ephx1",
  "Ephx2",
  
  "Pla2g4a",
  "Pla2g4b",
  "Pla2g4c",
  "Pla2g4d",
  "Pla2g4e",
  "Pla2g4f",
  "Pla2g6",
  "Pla2g7",
  "Pla2g2d",
  "Pla2g2e",
  "Pla2g5",
  "Pla2g12a",
  
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
  "Ptgr1",
  
  "Alox5",
  "Alox5ap",
  "Alox12",
  "Alox15",
  "Aloxe3",
  "Lta4h",
  "Ltc4s",
  "Ggt1",
  "Ggt5",
  
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
  "Cysltr2",
  
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

excluded_signature_genes <- unique(
  c(
    mitochondrial_genes,
    ribosomal_genes,
    hemoglobin_genes,
    cell_cycle_genes,
    cyp_genes,
    eicosanoid_exclusion
  )
)

cat(
  "\nGenes excluded from identity signatures:",
  length(excluded_signature_genes),
  "\n"
)


# ------------------------------------------------------------
# 8B. Detect fold-change column automatically
# ------------------------------------------------------------

possible_fc_columns <- c(
  "avg_log2FC",
  "avg_logFC"
)

fc_column <- intersect(
  possible_fc_columns,
  colnames(broad_markers_all)
)

if (length(fc_column) == 0) {
  stop(
    paste0(
      "Could not identify marker fold-change column. Columns are: ",
      paste(
        colnames(broad_markers_all),
        collapse = ", "
      )
    )
  )
}

fc_column <- fc_column[1]

cat(
  "Using marker FC column:",
  fc_column,
  "\n"
)


# ------------------------------------------------------------
# 8C. Filter markers
# ------------------------------------------------------------

broad_markers_filtered <- broad_markers_all %>%
  filter(
    !gene %in% excluded_signature_genes
  ) %>%
  filter(
    pct.1 >= 0.25
  ) %>%
  filter(
    .data[[fc_column]] >= 0.50
  )


# ------------------------------------------------------------
# 8D. Keep genes present in spatial data
# ------------------------------------------------------------

spatial_reference_genes <- rownames(
  spatial_objects[[sample_order[1]]]
)

broad_markers_filtered <- broad_markers_filtered %>%
  filter(
    gene %in% spatial_reference_genes
  )


# ------------------------------------------------------------
# 8E. Rank markers within each cell type
#
# Ranking prioritizes:
# 1. adjusted p-value
# 2. fold change
# 3. pct.1 - pct.2 specificity
# ------------------------------------------------------------

broad_markers_filtered <- broad_markers_filtered %>%
  mutate(
    Detection_Specificity =
      pct.1 - pct.2
  ) %>%
  arrange(
    cluster,
    p_val_adj,
    desc(.data[[fc_column]]),
    desc(Detection_Specificity)
  )


# ------------------------------------------------------------
# 8F. Select top 50 markers per broad cell type
# ------------------------------------------------------------

markers_per_cell_type <- 50

broad_signature_markers <- broad_markers_filtered %>%
  group_by(
    cluster
  ) %>%
  slice_head(
    n = markers_per_cell_type
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 8G. Prevent marker overlap between cell-type signatures
#
# A gene occurring in multiple signatures can make spatial
# scores artificially correlated.
# Keep only genes assigned to ONE broad cell type.
# ------------------------------------------------------------

marker_gene_frequency <- broad_signature_markers %>%
  count(
    gene,
    name = "CellType_Assignments"
  )

unique_identity_genes <- marker_gene_frequency %>%
  filter(
    CellType_Assignments == 1
  ) %>%
  pull(
    gene
  )

broad_signature_markers_unique <- broad_signature_markers %>%
  filter(
    gene %in% unique_identity_genes
  )


# ------------------------------------------------------------
# 8H. Count final signature genes
# ------------------------------------------------------------

signature_counts <- broad_signature_markers_unique %>%
  count(
    cluster,
    name = "Signature_Genes"
  ) %>%
  arrange(
    desc(Signature_Genes)
  )

cat(
  "\nFinal unique marker counts by broad cell type:\n"
)

print(
  as.data.frame(signature_counts),
  row.names = FALSE
)


# ------------------------------------------------------------
# 8I. Safety check
# ------------------------------------------------------------

low_signature_types <- signature_counts %>%
  filter(
    Signature_Genes < 10
  )

if (nrow(low_signature_types) > 0) {
  
  warning(
    paste0(
      "At least one broad cell type has fewer than 10 ",
      "unique signature genes. Review before spatial scoring."
    )
  )
  
  print(
    as.data.frame(low_signature_types),
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# 8J. Build signature list
# ------------------------------------------------------------

broad_signature_list <- split(
  broad_signature_markers_unique$gene,
  broad_signature_markers_unique$cluster
)


# ------------------------------------------------------------
# 8K. Save all signature outputs
# ------------------------------------------------------------

write.csv(
  broad_markers_filtered,
  file.path(
    results_dir,
    "11D_Broad_Markers_Filtered_All.csv"
  ),
  row.names = FALSE
)

write.csv(
  broad_signature_markers_unique,
  file.path(
    results_dir,
    "11D_Final_Broad_Spatial_Signature_Markers.csv"
  ),
  row.names = FALSE
)

write.csv(
  signature_counts,
  file.path(
    results_dir,
    "11D_Final_Broad_Signature_Counts.csv"
  ),
  row.names = FALSE
)

saveRDS(
  broad_signature_list,
  file.path(
    results_dir,
    "11D_Final_Broad_Signature_List.rds"
  )
)


# ------------------------------------------------------------
# 8L. Print actual marker genes
# ------------------------------------------------------------

cat("\n============================================\n")
cat("FINAL BROAD CELL-TYPE SIGNATURES\n")
cat("============================================\n")

for (cell_type in names(broad_signature_list)) {
  
  cat(
    "\n--- ",
    cell_type,
    " ---\n",
    sep = ""
  )
  
  cat(
    paste(
      broad_signature_list[[cell_type]],
      collapse = ", "
    ),
    "\n"
  )
}


# ------------------------------------------------------------
# 8M. Final Part-1 status
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 11D PART 1 COMPLETE\n")
cat("============================================\n")

cat(
  "Frozen scRNA cells:",
  ncol(scrna),
  "\n"
)

cat(
  "Frozen spatial spots:",
  sum(
    vapply(
      spatial_objects,
      ncol,
      numeric(1)
    )
  ),
  "\n"
)

cat(
  "Broad cell types with signatures:",
  length(broad_signature_list),
  "\n"
)

cat(
  "No scRNA clusters were changed.\n"
)

cat(
  "No cell-type annotations were changed.\n"
)

cat(
  "No spatial spots were removed.\n"
)

cat(
  "CYP/eicosanoid genes were excluded from cell-type signature derivation.\n"
)

cat(
  "No condition-level statistical testing was performed.\n"
)

cat("============================================\n")