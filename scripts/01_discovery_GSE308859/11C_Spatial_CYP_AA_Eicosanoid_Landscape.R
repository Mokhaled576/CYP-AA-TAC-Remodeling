# ============================================================
# STEP 11C
# SPATIAL CYP–ARACHIDONIC ACID / EICOSANOID LANDSCAPE
# GSE308859
#
# INPUT:
#   11B_FILTERED_NORMALIZED_SPATIAL_OBJECTS.rds
#
# PURPOSE:
#   1. Quantify all CYP genes across spatial samples
#   2. Characterize CYP-AA / eicosanoid pathway genes
#   3. Map priority CYP-AA genes spatially
#   4. Quantify spot-level detection by sample
#   5. Generate descriptive temporal trajectories
#   6. Generate pathway-level spatial module maps
#   7. Identify spatially informative genes for later
#      scRNA-spatial integration
#
# IMPORTANT:
#   - Spatial QC is FROZEN.
#   - No spots are filtered here.
#   - No cell-type labels are transferred here.
#   - No condition-level inferential tests are performed.
#   - Spots are NOT biological replicates.
#   - All TAC-vs-Sham comparisons are DESCRIPTIVE.
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
  library(scales)
})

set.seed(1234)


# ============================================================
# 2. PROJECT DIRECTORIES
# ============================================================

project_dir <- "D:/Master/ScRNA seq/TAC"

input_file <- file.path(
  project_dir,
  "CLEAN DATA",
  "11B_FILTERED_NORMALIZED_SPATIAL_OBJECTS.rds"
)

results_dir <- file.path(
  project_dir,
  "RESULTS",
  "11_SPATIAL",
  "11C_CYP_AA_EICOSANOID"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "11_SPATIAL",
  "11C_CYP_AA_EICOSANOID"
)

clean_dir <- file.path(
  project_dir,
  "CLEAN DATA"
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

dir.create(
  clean_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. LOAD FROZEN POST-QC SPATIAL OBJECTS
# ============================================================

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Input file does not exist:\n",
      input_file
    )
  )
}

spatial_objects <- readRDS(input_file)

if (!is.list(spatial_objects)) {
  stop(
    "Expected a list of Seurat spatial objects."
  )
}

cat("\nLoaded spatial samples:\n")
print(names(spatial_objects))


# ============================================================
# 4. DEFINE SAMPLE ORDER
# ============================================================

sample_order <- c(
  "Sham",
  "TAC_2W",
  "TAC_4W",
  "TAC_6W"
)

missing_samples <- setdiff(
  sample_order,
  names(spatial_objects)
)

if (length(missing_samples) > 0) {
  stop(
    paste(
      "Missing spatial samples:",
      paste(
        missing_samples,
        collapse = ", "
      )
    )
  )
}

spatial_objects <- spatial_objects[sample_order]


# ============================================================
# 5. VERIFY OBJECT INTEGRITY
# ============================================================
# ============================================================
# 5. VERIFY OBJECT INTEGRITY
# ============================================================

integrity_list <- list()

for (sid in sample_order) {
  
  cat("\nChecking spatial object:", sid, "\n")
  
  sobj <- spatial_objects[[sid]]
  
  # ----------------------------------------------------------
  # Confirm Seurat object
  # ----------------------------------------------------------
  
  if (!inherits(sobj, "Seurat")) {
    stop(
      paste0(
        sid,
        " is not a Seurat object."
      )
    )
  }
  
  # ----------------------------------------------------------
  # Confirm Spatial assay
  #
  # Avoid Assays(sobj) %in% logic because some
  # SeuratObject versions return a non-vector assay container.
  # ----------------------------------------------------------
  
  assay_names <- names(sobj@assays)
  
  if (!is.character(assay_names)) {
    assay_names <- as.character(assay_names)
  }
  
  if (!any(assay_names == "Spatial")) {
    stop(
      paste0(
        sid,
        " does not contain a Spatial assay. Available assays: ",
        paste(
          assay_names,
          collapse = ", "
        )
      )
    )
  }
  
  DefaultAssay(sobj) <- "Spatial"
  
  # ----------------------------------------------------------
  # Retrieve raw counts
  # ----------------------------------------------------------
  
  counts_mat <- GetAssayData(
    object = sobj,
    assay = "Spatial",
    layer = "counts"
  )
  
  # ----------------------------------------------------------
  # Retrieve normalized data
  # ----------------------------------------------------------
  
  data_mat <- GetAssayData(
    object = sobj,
    assay = "Spatial",
    layer = "data"
  )
  
  # ----------------------------------------------------------
  # Integrity checks
  # ----------------------------------------------------------
  
  counts_data_cellnames_identical <- identical(
    colnames(counts_mat),
    colnames(data_mat)
  )
  
  object_counts_cellnames_identical <- identical(
    colnames(sobj),
    colnames(counts_mat)
  )
  
  counts_data_genes_identical <- identical(
    rownames(counts_mat),
    rownames(data_mat)
  )
  
  # ----------------------------------------------------------
  # Store results
  # ----------------------------------------------------------
  
  integrity_list[[sid]] <- data.frame(
    Sample = sid,
    Spots_Object = ncol(sobj),
    Genes_Object = nrow(sobj),
    Spots_Counts = ncol(counts_mat),
    Genes_Counts = nrow(counts_mat),
    Spots_Data = ncol(data_mat),
    Genes_Data = nrow(data_mat),
    Counts_Data_Cellnames_Identical =
      counts_data_cellnames_identical,
    Object_Counts_Cellnames_Identical =
      object_counts_cellnames_identical,
    Counts_Data_Genes_Identical =
      counts_data_genes_identical,
    stringsAsFactors = FALSE
  )
  
  cat(
    "  Spots:",
    ncol(sobj),
    "\n"
  )
  
  cat(
    "  Genes:",
    nrow(sobj),
    "\n"
  )
  
  cat(
    "  Counts/data spot names identical:",
    counts_data_cellnames_identical,
    "\n"
  )
  
  cat(
    "  Object/counts spot names identical:",
    object_counts_cellnames_identical,
    "\n"
  )
  
  cat(
    "  Counts/data gene names identical:",
    counts_data_genes_identical,
    "\n"
  )
}

integrity_table <- dplyr::bind_rows(
  integrity_list
)

# ------------------------------------------------------------
# Strict integrity validation
# ------------------------------------------------------------

if (
  any(
    integrity_table$Spots_Object !=
    integrity_table$Spots_Counts
  )
) {
  stop(
    "Object and counts spot numbers are inconsistent."
  )
}

if (
  any(
    integrity_table$Spots_Counts !=
    integrity_table$Spots_Data
  )
) {
  stop(
    "Counts and normalized-data spot numbers are inconsistent."
  )
}

if (
  any(
    integrity_table$Genes_Counts !=
    integrity_table$Genes_Data
  )
) {
  stop(
    "Counts and normalized-data gene numbers are inconsistent."
  )
}

if (
  !all(
    integrity_table$Counts_Data_Cellnames_Identical
  )
) {
  stop(
    "Counts and normalized-data spot names are not identical."
  )
}

if (
  !all(
    integrity_table$Object_Counts_Cellnames_Identical
  )
) {
  stop(
    "Object and counts spot names are not identical."
  )
}

if (
  !all(
    integrity_table$Counts_Data_Genes_Identical
  )
) {
  stop(
    "Counts and normalized-data gene names are not identical."
  )
}

# ------------------------------------------------------------
# Save integrity table
# ------------------------------------------------------------

write.csv(
  integrity_table,
  file.path(
    results_dir,
    "11C_Spatial_Object_Integrity.csv"
  ),
  row.names = FALSE
)

cat("\n============================================\n")
cat("SPATIAL OBJECT INTEGRITY CHECK\n")
cat("============================================\n")

print(
  as.data.frame(integrity_table),
  row.names = FALSE
)

cat("\nAll four spatial objects passed integrity checks.\n")
cat("============================================\n")

# ============================================================
# 6. DEFINE GENE SETS
# ============================================================

priority_cyp_aa <- c(
  "Cyp1b1",
  "Cyp2j6",
  "Cyp2j9",
  "Cyp4f13",
  "Cyp4f16",
  "Cyp4f17",
  "Cyp4f18",
  "Ephx2"
)

gene_sets <- list(
  
  CYP_AA_Epoxygenase_Hydroxylase = c(
    "Cyp1a1",
    "Cyp1a2",
    "Cyp1b1",
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
    "Cyp2c70",
    "Cyp2j5",
    "Cyp2j6",
    "Cyp2j8",
    "Cyp2j9",
    "Cyp2j11",
    "Cyp2j12",
    "Cyp2j13",
    "Cyp4a10",
    "Cyp4a12a",
    "Cyp4a12b",
    "Cyp4a14",
    "Cyp4f13",
    "Cyp4f14",
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
    "Aloxe3",
    "Lta4h",
    "Ltc4s",
    "Ggt1",
    "Ggt5",
    "Ggt7"
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
# 7. VERIFY GENE AVAILABILITY
# ============================================================

reference_genes <- rownames(
  spatial_objects[[sample_order[1]]]
)

all_cyp_genes <- grep(
  "^Cyp[0-9]",
  reference_genes,
  value = TRUE
)

cat(
  "\nNumber of CYP genes in spatial dataset:",
  length(all_cyp_genes),
  "\n"
)

gene_availability_list <- list()

for (module_name in names(gene_sets)) {
  
  requested <- unique(
    gene_sets[[module_name]]
  )
  
  present <- requested[
    requested %in% reference_genes
  ]
  
  absent <- setdiff(
    requested,
    reference_genes
  )
  
  gene_availability_list[[module_name]] <- data.frame(
    Module = module_name,
    Requested_Genes = length(requested),
    Present_Genes = length(present),
    Missing_Genes = length(absent),
    Present = paste(
      present,
      collapse = ";"
    ),
    Missing = paste(
      absent,
      collapse = ";"
    ),
    stringsAsFactors = FALSE
  )
}

gene_availability <- bind_rows(
  gene_availability_list
)

write.csv(
  gene_availability,
  file.path(
    results_dir,
    "11C_Pathway_Gene_Availability.csv"
  ),
  row.names = FALSE
)

priority_availability <- data.frame(
  Gene = priority_cyp_aa,
  Present = priority_cyp_aa %in% reference_genes,
  stringsAsFactors = FALSE
)

write.csv(
  priority_availability,
  file.path(
    results_dir,
    "11C_Priority_CYP_AA_Gene_Availability.csv"
  ),
  row.names = FALSE
)

cat("\nPriority CYP-AA availability:\n")
print(
  priority_availability,
  row.names = FALSE
)


# ============================================================
# 8. HELPER FUNCTION:
#    GENE EXPRESSION SUMMARY
# ============================================================

summarize_gene_expression <- function(
    sobj,
    genes,
    sample_name
) {
  
  DefaultAssay(sobj) <- "Spatial"
  
  counts_mat <- GetAssayData(
    sobj,
    assay = "Spatial",
    layer = "counts"
  )
  
  data_mat <- GetAssayData(
    sobj,
    assay = "Spatial",
    layer = "data"
  )
  
  genes_present <- intersect(
    genes,
    rownames(counts_mat)
  )
  
  if (length(genes_present) == 0) {
    return(
      data.frame()
    )
  }
  
  counts_sub <- counts_mat[
    genes_present,
    ,
    drop = FALSE
  ]
  
  data_sub <- data_mat[
    genes_present,
    ,
    drop = FALSE
  ]
  
  positive_spots <- Matrix::rowSums(
    counts_sub > 0
  )
  
  total_umi <- Matrix::rowSums(
    counts_sub
  )
  
  mean_raw_count <- Matrix::rowMeans(
    counts_sub
  )
  
  mean_normalized <- Matrix::rowMeans(
    data_sub
  )
  
  total_spots <- ncol(
    counts_sub
  )
  
  data.frame(
    Sample = sample_name,
    Gene = genes_present,
    Total_Spots = total_spots,
    Positive_Spots = as.numeric(
      positive_spots
    ),
    Percent_Positive = as.numeric(
      positive_spots
    ) / total_spots * 100,
    Mean_Raw_Count = as.numeric(
      mean_raw_count
    ),
    Mean_Normalized_Expression = as.numeric(
      mean_normalized
    ),
    Total_UMI = as.numeric(
      total_umi
    ),
    stringsAsFactors = FALSE
  )
}


# ============================================================
# 9. ALL CYP SPATIAL EXPRESSION SUMMARY
# ============================================================

all_cyp_summary_list <- list()

for (sid in sample_order) {
  
  cat(
    "\nSummarizing all CYP genes:",
    sid,
    "\n"
  )
  
  all_cyp_summary_list[[sid]] <-
    summarize_gene_expression(
      spatial_objects[[sid]],
      all_cyp_genes,
      sid
    )
}

all_cyp_summary <- bind_rows(
  all_cyp_summary_list
)

all_cyp_summary$Sample <- factor(
  all_cyp_summary$Sample,
  levels = sample_order
)

write.csv(
  all_cyp_summary,
  file.path(
    results_dir,
    "11C_All_CYP_Spatial_Expression_By_Sample.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. GLOBAL CYP SPATIAL DETECTION
# ============================================================
# ============================================================
# 10. GLOBAL CYP SPATIAL DETECTION
# ============================================================

cat("\n============================================\n")
cat("GLOBAL CYP SPATIAL DETECTION\n")
cat("============================================\n")

# ------------------------------------------------------------
# 10A. Verify required columns
# ------------------------------------------------------------

required_cyp_columns <- c(
  "Sample",
  "Gene",
  "Total_Spots",
  "Positive_Spots",
  "Total_UMI",
  "Mean_Normalized_Expression"
)

missing_cyp_columns <- setdiff(
  required_cyp_columns,
  colnames(all_cyp_summary)
)

if (length(missing_cyp_columns) > 0) {
  stop(
    paste0(
      "Missing required columns in all_cyp_summary: ",
      paste(
        missing_cyp_columns,
        collapse = ", "
      )
    )
  )
}

# ------------------------------------------------------------
# 10B. Create weighted-expression numerator BEFORE summarising
# ------------------------------------------------------------

cyp_global_input <- all_cyp_summary %>%
  mutate(
    Normalized_Expression_Weighted_Numerator =
      Mean_Normalized_Expression * Total_Spots
  )

# ------------------------------------------------------------
# 10C. Collapse across the four spatial samples
# ------------------------------------------------------------

global_cyp_summary <- cyp_global_input %>%
  group_by(Gene) %>%
  summarise(
    Global_Total_Spots =
      sum(
        Total_Spots,
        na.rm = TRUE
      ),
    
    Global_Positive_Spots =
      sum(
        Positive_Spots,
        na.rm = TRUE
      ),
    
    Global_Total_UMI =
      sum(
        Total_UMI,
        na.rm = TRUE
      ),
    
    Weighted_Normalized_Numerator =
      sum(
        Normalized_Expression_Weighted_Numerator,
        na.rm = TRUE
      ),
    
    Samples_Detected =
      sum(
        Positive_Spots > 0,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  mutate(
    Percent_Positive =
      100 *
      Global_Positive_Spots /
      Global_Total_Spots,
    
    Mean_Raw_Count =
      Global_Total_UMI /
      Global_Total_Spots,
    
    Mean_Normalized_Expression =
      Weighted_Normalized_Numerator /
      Global_Total_Spots
  ) %>%
  select(
    Gene,
    Total_Spots = Global_Total_Spots,
    Positive_Spots = Global_Positive_Spots,
    Percent_Positive,
    Total_UMI = Global_Total_UMI,
    Mean_Raw_Count,
    Mean_Normalized_Expression,
    Samples_Detected
  ) %>%
  arrange(
    desc(Percent_Positive),
    desc(Total_UMI)
  )

# ------------------------------------------------------------
# 10D. Sanity checks
# ------------------------------------------------------------

expected_total_spots <- sum(
  vapply(
    spatial_objects,
    ncol,
    numeric(1)
  )
)

cat(
  "\nExpected total retained spatial spots:",
  expected_total_spots,
  "\n"
)

cat(
  "Total spots represented per CYP gene should equal:",
  expected_total_spots,
  "\n"
)

if (
  any(
    global_cyp_summary$Total_Spots !=
    expected_total_spots
  )
) {
  warning(
    paste0(
      "At least one CYP gene does not contain the expected ",
      "total number of spots. Inspect all_cyp_summary before ",
      "biological interpretation."
    )
  )
}

if (
  any(
    global_cyp_summary$Positive_Spots >
    global_cyp_summary$Total_Spots
  )
) {
  stop(
    "Positive_Spots exceeds Total_Spots for at least one CYP gene."
  )
}

if (
  any(
    global_cyp_summary$Percent_Positive < 0 |
    global_cyp_summary$Percent_Positive > 100
  )
) {
  stop(
    "Invalid CYP detection percentage detected."
  )
}

if (
  any(
    global_cyp_summary$Samples_Detected < 0 |
    global_cyp_summary$Samples_Detected > length(sample_order)
  )
) {
  stop(
    "Invalid Samples_Detected value detected."
  )
}

# ------------------------------------------------------------
# 10E. Assign descriptive spatial detection classes
#
# These are descriptive categories only.
# They are NOT statistical significance categories.
# ------------------------------------------------------------

global_cyp_summary <- global_cyp_summary %>%
  mutate(
    Spatial_Detection_Class =
      case_when(
        Percent_Positive >= 10 ~
          "Robust",
        
        Percent_Positive >= 2 ~
          "Moderate",
        
        Percent_Positive > 0 ~
          "Sparse",
        
        TRUE ~
          "Undetected"
      )
  )

# ------------------------------------------------------------
# 10F. Save complete global CYP table
# ------------------------------------------------------------

write.csv(
  global_cyp_summary,
  file.path(
    results_dir,
    "11C_Global_CYP_Spatial_Detection.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10G. Save detected CYP genes only
# ------------------------------------------------------------

detected_cyp_summary <- global_cyp_summary %>%
  filter(
    Positive_Spots > 0
  )

write.csv(
  detected_cyp_summary,
  file.path(
    results_dir,
    "11C_Detected_CYPs_Spatial.csv"
  ),
  row.names = FALSE
)