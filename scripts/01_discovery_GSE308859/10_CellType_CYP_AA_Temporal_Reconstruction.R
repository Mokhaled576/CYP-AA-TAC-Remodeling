# ============================================================
# STEP 10
# CELL-TYPE-SPECIFIC TEMPORAL RECONSTRUCTION OF THE
# CYP–ARACHIDONIC ACID / EICOSANOID NETWORK
#
# GSE308859 TAC CARDIAC scRNA-seq
#
# PRIMARY QUESTION:
# Are CYP/eicosanoid changes caused mainly by changing
# cellular composition, or do individual cardiac cell types
# themselves remodel these pathways during TAC progression?
#
# IMPORTANT:
#   - Frozen cell identities are preserved
#   - No reclustering
#   - No reannotation
#   - No cell-level condition p-values
#   - One scRNA library per condition
#   - All temporal comparisons are DESCRIPTIVE
# ============================================================


# ============================================================
# 1. CLEAN ENVIRONMENT
# ============================================================

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
set.seed(12345)


# ============================================================
# 2. LOAD PACKAGES
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

project_dir <- "D:/Master/ScRNA seq/TAC"

clean_dir <- file.path(
  project_dir,
  "CLEAN DATA"
)

results_dir <- file.path(
  project_dir,
  "RESULTS",
  "10_CELLTYPE_CYP_AA_TEMPORAL"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "10_CELLTYPE_CYP_AA_TEMPORAL"
)

notes_dir <- file.path(
  project_dir,
  "NOTES"
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
  notes_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 4. LOAD STEP 09 OBJECT
# ============================================================

input_file <- file.path(
  clean_dir,
  "09_BIOLOGICAL_PROGRAM_SCORED_ATLAS.rds"
)

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Step 09 object not found:\n",
      input_file
    )
  )
}

cat("\nLoading Step 09 object...\n")

obj <- readRDS(input_file)

cat("Object loaded successfully.\n")
cat("Cells:", ncol(obj), "\n")
cat("Genes:", nrow(obj), "\n")


# ============================================================
# 5. VERIFY FROZEN METADATA
# ============================================================

required_metadata <- c(
  "sample_id",
  "refined_cell_type"
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

cat("\nRefined cell types:\n")
print(table(obj$refined_cell_type))

cat("\nSamples:\n")
print(table(obj$sample_id))


# ============================================================
# 6. STANDARDIZE CONDITION
# ============================================================

condition_map <- c(
  "Sham" = "Sham",
  "TAC_2W" = "TAC_2W",
  "TAC_4W" = "TAC_4W",
  "TAC_6W" = "TAC_6W"
)

obj$condition_10 <- unname(
  condition_map[
    as.character(obj$sample_id)
  ]
)

if (any(is.na(obj$condition_10))) {
  
  bad_samples <- unique(
    as.character(obj$sample_id)[
      is.na(obj$condition_10)
    ]
  )
  
  stop(
    paste(
      "Unrecognized sample IDs:",
      paste(bad_samples, collapse = ", ")
    )
  )
}

condition_order <- c(
  "Sham",
  "TAC_2W",
  "TAC_4W",
  "TAC_6W"
)

obj$condition_10 <- factor(
  obj$condition_10,
  levels = condition_order
)


# ============================================================
# 7. VERIFY RNA ASSAY
# ============================================================

assay_names <- names(obj@assays)

cat("\nAssays present:\n")
print(assay_names)

if (!("RNA" %in% assay_names)) {
  stop("RNA assay is absent.")
}

DefaultAssay(obj) <- "RNA"

cat("\nDefault assay:", DefaultAssay(obj), "\n")


# ============================================================
# 8. GET RAW COUNTS
# ============================================================

counts <- GetAssayData(
  object = obj,
  assay = "RNA",
  layer = "counts"
)

cat("\nRaw count matrix:\n")
cat("Genes:", nrow(counts), "\n")
cat("Cells:", ncol(counts), "\n")


# ============================================================
# 9. GET NORMALIZED RNA DATA
# ============================================================

rna_layers <- Layers(
  obj[["RNA"]]
)

cat("\nRNA layers:\n")
print(rna_layers)

if (!("data" %in% rna_layers)) {
  
  cat(
    "\nNormalized RNA data absent.",
    "\nRunning LogNormalize...\n"
  )
  
  obj <- NormalizeData(
    obj,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
}

norm_data <- GetAssayData(
  object = obj,
  assay = "RNA",
  layer = "data"
)

cat("\nNormalized expression matrix ready.\n")


# ============================================================
# 10. DEFINE CYP–AA / EICOSANOID NETWORK
# ============================================================
#
# Keep pathway arms separate.
# We do NOT collapse the biology into one arbitrary score.
# ============================================================

network_genes <- list(
  
  CYP_EPOXYGENASE = c(
    "Cyp2j5",
    "Cyp2j6",
    "Cyp2j8",
    "Cyp2j9",
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
    "Cyp2c70"
  ),
  
  CYP_OMEGA_HYDROXYLASE = c(
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
    "Cyp4f39",
    "Cyp4v3"
  ),
  
  OTHER_RELEVANT_CYP = c(
    "Cyp1a1",
    "Cyp1a2",
    "Cyp1b1",
    "Cyp2e1",
    "Cyp2u1"
  ),
  
  EPOXIDE_HYDROLASE = c(
    "Ephx1",
    "Ephx2"
  ),
  
  AA_RELEASE = c(
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
    "Pla2g12a"
  ),
  
  COX_PROSTANOID = c(
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
  
  LOX_LEUKOTRIENE = c(
    "Alox5",
    "Alox5ap",
    "Alox12",
    "Alox15",
    "Aloxe3",
    "Lta4h",
    "Ltc4s",
    "Ggt1",
    "Ggt5"
  ),
  
  EICOSANOID_RECEPTORS = c(
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
  
  LIPID_CONTEXT = c(
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
# 11. CHECK GENE AVAILABILITY
# ============================================================

all_genes <- rownames(counts)

availability_list <- lapply(
  names(network_genes),
  function(pathway) {
    
    requested <- unique(
      network_genes[[pathway]]
    )
    
    present <- requested[
      requested %in% all_genes
    ]
    
    missing <- requested[
      !requested %in% all_genes
    ]
    
    data.frame(
      Pathway = pathway,
      Requested = length(requested),
      Present = length(present),
      Missing = length(missing),
      Present_Genes = paste(
        present,
        collapse = ";"
      ),
      Missing_Genes = paste(
        missing,
        collapse = ";"
      ),
      stringsAsFactors = FALSE
    )
  }
)

availability_table <- bind_rows(
  availability_list
)

write.csv(
  availability_table,
  file.path(
    results_dir,
    "10_Network_Gene_Availability.csv"
  ),
  row.names = FALSE
)

cat("\nNetwork gene availability:\n")

print(
  availability_table[
    ,
    c(
      "Pathway",
      "Requested",
      "Present",
      "Missing"
    )
  ]
)


# ============================================================
# 12. CREATE GENE-TO-PATHWAY LOOKUP
# ============================================================

gene_pathway_lookup <- bind_rows(
  lapply(
    names(network_genes),
    function(pathway) {
      
      genes_i <- network_genes[[pathway]]
      
      genes_i <- genes_i[
        genes_i %in% all_genes
      ]
      
      data.frame(
        Gene = genes_i,
        Pathway = pathway,
        stringsAsFactors = FALSE
      )
    }
  )
)

gene_pathway_lookup <- gene_pathway_lookup %>%
  distinct()

write.csv(
  gene_pathway_lookup,
  file.path(
    results_dir,
    "10_Gene_Pathway_Lookup.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. DEFINE ALL PRESENT NETWORK GENES
# ============================================================

network_present <- unique(
  gene_pathway_lookup$Gene
)

cat(
  "\nTotal unique network genes present:",
  length(network_present),
  "\n"
)

if (length(network_present) == 0) {
  stop("No network genes were found.")
}


# ============================================================
# 14. SUBSET MATRICES TO NETWORK GENES ONLY
# ============================================================
#
# This avoids converting the full 32,285 x 25,186 matrix
# into a dense object.
# ============================================================

counts_network <- counts[
  network_present,
  ,
  drop = FALSE
]

norm_network <- norm_data[
  network_present,
  ,
  drop = FALSE
]

cat("\nNetwork matrices created.\n")
cat("Genes:", nrow(counts_network), "\n")
cat("Cells:", ncol(counts_network), "\n")


# ============================================================
# 15. CELL-TYPE ORDER
# ============================================================

preferred_celltype_order <- c(
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

actual_celltypes <- unique(
  as.character(obj$refined_cell_type)
)

celltype_order <- preferred_celltype_order[
  preferred_celltype_order %in% actual_celltypes
]

extra_celltypes <- setdiff(
  actual_celltypes,
  celltype_order
)

celltype_order <- c(
  celltype_order,
  extra_celltypes
)


# ============================================================
# 16. BUILD CELL METADATA TABLE
# ============================================================

cell_metadata <- data.frame(
  Cell = colnames(obj),
  Cell_Type = as.character(
    obj$refined_cell_type
  ),
  Condition = as.character(
    obj$condition_10
  ),
  stringsAsFactors = FALSE
)

rownames(cell_metadata) <- cell_metadata$Cell


# ============================================================
# 17. CELL COUNTS BY CELL TYPE × CONDITION
# ============================================================

cell_counts <- cell_metadata %>%
  count(
    Cell_Type,
    Condition,
    name = "Cells"
  )

cell_counts$Reliable_20Cells <- (
  cell_counts$Cells >= 20
)

cell_counts$Reliable_25Cells <- (
  cell_counts$Cells >= 25
)

write.csv(
  cell_counts,
  file.path(
    results_dir,
    "10_Cell_Counts_By_CellType_Condition.csv"
  ),
  row.names = FALSE
)

cat("\nCell counts by cell type and condition:\n")
print(cell_counts)


# ============================================================
# 18. CELL COMPOSITION
# ============================================================
#
# This represents relative captured cell representation,
# NOT direct biological abundance.
# ============================================================

composition <- cell_counts %>%
  group_by(
    Condition
  ) %>%
  mutate(
    Total_Captured_Cells = sum(Cells),
    Percent_Captured_Cells = 100 * Cells /
      Total_Captured_Cells
  ) %>%
  ungroup()

write.csv(
  composition,
  file.path(
    results_dir,
    "10_Relative_Captured_Cell_Representation.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. FUNCTION FOR GENE-LEVEL SUMMARY
# ============================================================

summarize_gene_group <- function(
    count_matrix,
    normalized_matrix,
    selected_cells,
    genes
) {
  
  if (length(selected_cells) == 0) {
    return(NULL)
  }
  
  count_sub <- count_matrix[
    genes,
    selected_cells,
    drop = FALSE
  ]
  
  norm_sub <- normalized_matrix[
    genes,
    selected_cells,
    drop = FALSE
  ]
  
  n_cells <- length(selected_cells)
  
  positive_cells <- Matrix::rowSums(
    count_sub > 0
  )
  
  percent_positive <- 100 * positive_cells /
    n_cells
  
  total_umi <- Matrix::rowSums(
    count_sub
  )
  
  mean_raw_count <- Matrix::rowMeans(
    count_sub
  )
  
  mean_normalized <- Matrix::rowMeans(
    norm_sub
  )
  
  positive_mean_raw <- rep(
    NA_real_,
    length(genes)
  )
  
  positive_mean_normalized <- rep(
    NA_real_,
    length(genes)
  )
  
  for (g in seq_along(genes)) {
    
    raw_values <- count_sub[
      g,
      ,
      drop = TRUE
    ]
    
    norm_values <- norm_sub[
      g,
      ,
      drop = TRUE
    ]
    
    positive_index <- raw_values > 0
    
    if (any(positive_index)) {
      
      positive_mean_raw[g] <- mean(
        raw_values[positive_index]
      )
      
      positive_mean_normalized[g] <- mean(
        norm_values[positive_index]
      )
    }
  }
  
  data.frame(
    Gene = genes,
    Cells = n_cells,
    Positive_Cells = as.numeric(
      positive_cells
    ),
    Percent_Expressing = as.numeric(
      percent_positive
    ),
    Mean_Raw_Count = as.numeric(
      mean_raw_count
    ),
    Mean_Normalized_Expression = as.numeric(
      mean_normalized
    ),
    Positive_Cell_Mean_Raw = positive_mean_raw,
    Positive_Cell_Mean_Normalized =
      positive_mean_normalized,
    Total_UMI = as.numeric(
      total_umi
    ),
    stringsAsFactors = FALSE
  )
}


# ============================================================
# 20. GENE × CELL TYPE × CONDITION SUMMARY
# ============================================================

summary_list <- list()

counter <- 1

for (celltype in celltype_order) {
  
  for (condition in condition_order) {
    
    selected_cells <- cell_metadata$Cell[
      cell_metadata$Cell_Type == celltype &
        cell_metadata$Condition == condition
    ]
    
    cat(
      "\nProcessing:",
      celltype,
      "|",
      condition,
      "| cells:",
      length(selected_cells),
      "\n"
    )
    
    if (length(selected_cells) == 0) {
      next
    }
    
    temp <- summarize_gene_group(
      count_matrix = counts_network,
      normalized_matrix = norm_network,
      selected_cells = selected_cells,
      genes = network_present
    )
    
    temp$Cell_Type <- celltype
    temp$Condition <- condition
    
    temp$Reliable_20Cells <- (
      temp$Cells >= 20
    )
    
    temp$Reliable_25Cells <- (
      temp$Cells >= 25
    )
    
    summary_list[[counter]] <- temp
    
    counter <- counter + 1
  }
}

gene_summary <- bind_rows(
  summary_list
)

gene_summary <- gene_summary %>%
  left_join(
    gene_pathway_lookup,
    by = "Gene"
  )

gene_summary$Condition <- factor(
  gene_summary$Condition,
  levels = condition_order
)

gene_summary$Cell_Type <- factor(
  gene_summary$Cell_Type,
  levels = celltype_order
)

write.csv(
  gene_summary,
  file.path(
    results_dir,
    "10_All_Eicosanoid_Genes_CellType_Condition.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 21. SHAM BASELINE TABLE
# ============================================================

sham_baseline <- gene_summary %>%
  filter(
    Condition == "Sham"
  ) %>%
  select(
    Gene,
    Pathway,
    Cell_Type,
    Sham_Cells = Cells,
    Sham_Positive_Cells = Positive_Cells,
    Sham_Percent_Expressing = Percent_Expressing,
    Sham_Mean_Raw_Count = Mean_Raw_Count,
    Sham_Mean_Normalized_Expression =
      Mean_Normalized_Expression,
    Sham_Positive_Cell_Mean_Normalized =
      Positive_Cell_Mean_Normalized
  )

write.csv(
  sham_baseline,
  file.path(
    results_dir,
    "10_Sham_Baseline_By_Gene_CellType.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 22. JOIN SHAM BASELINE TO ALL CONDITIONS
# ============================================================

temporal_change <- gene_summary %>%
  left_join(
    sham_baseline,
    by = c(
      "Gene",
      "Pathway",
      "Cell_Type"
    )
  )


# ============================================================
# 23. CALCULATE DESCRIPTIVE CHANGE VS SHAM
# ============================================================

detection_pseudocount <- 0.1
expression_pseudocount <- 0.01

temporal_change <- temporal_change %>%
  mutate(
    
    Delta_Percent_Expressing =
      Percent_Expressing -
      Sham_Percent_Expressing,
    
    Log2_Detection_Ratio_vs_Sham =
      log2(
        (Percent_Expressing +
           detection_pseudocount) /
          (Sham_Percent_Expressing +
             detection_pseudocount)
      ),
    
    Delta_Mean_Normalized =
      Mean_Normalized_Expression -
      Sham_Mean_Normalized_Expression,
    
    Log2_Mean_Normalized_Ratio_vs_Sham =
      log2(
        (Mean_Normalized_Expression +
           expression_pseudocount) /
          (Sham_Mean_Normalized_Expression +
             expression_pseudocount)
      )
  )

write.csv(
  temporal_change,
  file.path(
    results_dir,
    "10_Eicosanoid_Temporal_Change_vs_Sham.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 24. PRIORITY GENES
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

priority_genes <- priority_genes[
  priority_genes %in% network_present
]

cat("\nPriority genes detected:\n")
print(priority_genes)

priority_table <- temporal_change %>%
  filter(
    Gene %in% priority_genes
  )

write.csv(
  priority_table,
  file.path(
    results_dir,
    "10_Priority_CYP_AA_CellType_Temporal.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 25. PRIORITY GENE DETECTION TRAJECTORIES
# ============================================================

priority_plot_data <- gene_summary %>%
  filter(
    Gene %in% priority_genes
  ) %>%
  filter(