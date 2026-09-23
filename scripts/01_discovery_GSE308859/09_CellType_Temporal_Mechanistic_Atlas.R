# ============================================================
# STEP 09
# CELL-TYPE-RESOLVED BIOLOGICAL PROGRAM LANDSCAPE
# GSE308859 TAC CARDIAC scRNA-seq
#
# PURPOSE:
#   Map major biological programs relevant to CYP-mediated
#   arachidonic-acid metabolism and cardiac remodeling across
#   the frozen cardiac cell atlas and across TAC time.
#
# IMPORTANT:
#   - Uses the frozen Step 07C/Step 08 atlas
#   - Does NOT recluster
#   - Does NOT change cell-type annotations
#   - Descriptive analysis only because there is one scRNA
#     library per condition
# ============================================================


# ============================================================
# 1. CLEAR ENVIRONMENT
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
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})


# ============================================================
# 3. DEFINE PATHS
# ============================================================

project_dir <- "D:/Master/ScRNA seq/TAC"

clean_dir <- file.path(
  project_dir,
  "CLEAN DATA"
)

results_dir <- file.path(
  project_dir,
  "RESULTS",
  "09_BIOLOGICAL_PROGRAMS"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "09_BIOLOGICAL_PROGRAMS"
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
# 4. FIND AND LOAD INPUT OBJECT
# ============================================================

candidate_files <- c(
  file.path(
    clean_dir,
    "08_CYP_LANDSCAPE_ATLAS.rds"
  ),
  file.path(
    clean_dir,
    "08_CYP_ANNOTATED_ATLAS.rds"
  ),
  file.path(
    clean_dir,
    "07C_FINAL_FROZEN_ANNOTATED_ATLAS.rds"
  )
)

existing_files <- candidate_files[
  file.exists(candidate_files)
]

if (length(existing_files) == 0) {
  stop(
    paste0(
      "Could not find the Step 08 or Step 07C atlas in:\n",
      clean_dir
    )
  )
}

input_file <- existing_files[1]

cat("\nLoading:\n")
cat(input_file, "\n\n")

obj <- readRDS(input_file)

cat("Object loaded successfully.\n")
cat("Cells:", ncol(obj), "\n")
cat("Genes:", nrow(obj), "\n")


# ============================================================
# 5. VERIFY REQUIRED METADATA
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

cat("\nRequired metadata verified.\n")

cat("\nRefined cell-type counts:\n")
print(
  table(obj$refined_cell_type)
)

cat("\nSample counts:\n")
print(
  table(obj$sample_id)
)


# ============================================================
# 6. STANDARDIZE CONDITION METADATA
# ============================================================

sample_values <- unique(
  as.character(obj$sample_id)
)

cat("\nSample IDs detected:\n")
print(sample_values)

sample_map <- c(
  "Sham" = "Sham",
  "TAC_2W" = "TAC_2W",
  "TAC_4W" = "TAC_4W",
  "TAC_6W" = "TAC_6W"
)

obj$condition_09 <- unname(
  sample_map[
    as.character(obj$sample_id)
  ]
)

if (any(is.na(obj$condition_09))) {
  
  stop(
    paste0(
      "Some sample_id values could not be mapped. ",
      "Detected sample IDs: ",
      paste(sample_values, collapse = ", ")
    )
  )
}

condition_order <- c(
  "Sham",
  "TAC_2W",
  "TAC_4W",
  "TAC_6W"
)

obj$condition_09 <- factor(
  obj$condition_09,
  levels = condition_order
)

cat("\nCondition counts:\n")
print(
  table(obj$condition_09)
)


# ============================================================
# ============================================================
# 7. SET AND VERIFY RNA ASSAY
# ============================================================

assay_names <- names(obj@assays)

cat("\nAssays present in object:\n")
print(assay_names)

if (!("RNA" %in% assay_names)) {
  stop(
    paste0(
      "RNA assay is not present in the object. ",
      "Available assays: ",
      paste(assay_names, collapse = ", ")
    )
  )
}

DefaultAssay(obj) <- "RNA"

cat("\nDefault assay successfully set to:\n")
print(DefaultAssay(obj))

cat("\nRNA assay dimensions:\n")
cat("Genes:", nrow(obj[["RNA"]]), "\n")
cat("Cells:", ncol(obj[["RNA"]]), "\n")
# ============================================================
# 8. VERIFY NORMALIZED RNA DATA
# ============================================================

rna_layers <- Layers(
  obj[["RNA"]]
)

cat("\nRNA assay layers:\n")
print(rna_layers)

if (!"data" %in% rna_layers) {
  
  cat(
    "\nNormalized RNA data layer not found.",
    "\nRunning NormalizeData()...\n"
  )
  
  obj <- NormalizeData(
    obj,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
}

cat("\nRNA data ready.\n")


# ============================================================
# 9. DEFINE BIOLOGICAL PROGRAM GENE SETS
# ============================================================
#
# These are targeted mechanistic programs relevant to:
#
#   CYP-AA metabolism
#   eicosanoid biology
#   inflammation
#   NF-kB
#   oxidative stress
#   Nrf2
#   ferroptosis
#   lipid peroxidation
#   fibrosis
#   TGF-beta
#   hypoxia
#   apoptosis
#   mitochondrial stress
#   PPAR/lipid metabolism
#   estrogen signaling
#
# We will automatically retain only genes actually present
# in this mouse dataset.
# ============================================================


program_genes <- list(
  
  CYP_AA_EPOXYGENASE = c(
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
  
  CYP_AA_OMEGA_HYDROXYLASE = c(
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
  
  EET_METABOLISM = c(
    "Ephx1",
    "Ephx2"
  ),
  
  AA_RELEASE = c(
    "Pla2g4a",
    "Pla2g4b",
    "Pla2g4c",
    "Pla2g6",
    "Pla2g7"
  ),
  
  COX_PROSTANOID = c(
    "Ptgs1",
    "Ptgs2",
    "Ptges",
    "Ptges2",
    "Ptges3",
    "Ptgis",
    "Tbxas1"
  ),
  
  LOX_LEUKOTRIENE = c(
    "Alox5",
    "Alox5ap",
    "Alox12",
    "Alox12b",
    "Alox15",
    "Alox15b",
    "Lta4h",
    "Ltc4s"
  ),
  
  NF_KB_SIGNALING = c(
    "Nfkb1",
    "Nfkb2",
    "Rela",
    "Relb",
    "Rel",
    "Nfkbia",
    "Nfkbib",
    "Nfkbie",
    "Ikbkb",
    "Ikbkg",
    "Tnfaip3"
  ),
  
  PRO_INFLAMMATORY = c(
    "Tnf",
    "Il1a",
    "Il1b",
    "Il6",
    "Ccl2",
    "Ccl3",
    "Ccl4",
    "Ccl5",
    "Ccl7",
    "Cxcl1",
    "Cxcl2",
    "Cxcl9",
    "Cxcl10",
    "Cxcl12",
    "Nlrp3"
  ),
  
  INFLAMMASOME = c(
    "Nlrp3",
    "Pycard",
    "Casp1",
    "Il1b",
    "Il18",
    "Gsdmd"
  ),
  
  NRF2_ANTIOXIDANT = c(
    "Nfe2l2",
    "Keap1",
    "Hmox1",
    "Nqo1",
    "Gclc",
    "Gclm",
    "Gsr",
    "Txnrd1",
    "Srxn1",
    "Sqstm1"
  ),
  
  ROS_DETOXIFICATION = c(
    "Sod1",
    "Sod2",
    "Sod3",
    "Cat",
    "Gpx1",
    "Gpx3",
    "Gpx4",
    "Prdx1",
    "Prdx2",
    "Prdx3",
    "Prdx4",
    "Prdx5",
    "Prdx6",
    "Txn1",
    "Txn2"
  ),
  
  ROS_GENERATION = c(
    "Cyba",
    "Cybb",
    "Ncf1",
    "Ncf2",
    "Ncf4",
    "Nox1",
    "Nox2",
    "Nox4",
    "Duox1",
    "Duox2"
  ),
  
  FERROPTOSIS_DEFENSE = c(
    "Gpx4",
    "Slc7a11",
    "Slc3a2",
    "Aifm2",
    "Gch1",
    "Fth1",
    "Ftl1",
    "Nfe2l2",
    "Hmox1"
  ),
  
  LIPID_PEROXIDATION = c(
    "Acsl4",
    "Lpcat3",
    "Alox5",
    "Alox12",
    "Alox15",
    "Por",
    "Sat1"
  ),
  
  TGF_BETA_FIBROSIS = c(
    "Tgfb1",
    "Tgfb2",
    "Tgfb3",
    "Tgfbr1",
    "Tgfbr2",
    "Smad2",
    "Smad3",
    "Smad4",
    "Col1a1",
    "Col1a2",
    "Col3a1",
    "Fn1",
    "Postn",
    "Ctgf",
    "Acta2"
  ),
  
  ECM_REMODELING = c(
    "Col1a1",
    "Col1a2",
    "Col3a1",
    "Col5a1",
    "Col5a2",
    "Fn1",
    "Postn",
    "Timp1",
    "Timp2",
    "Mmp2",
    "Mmp3",
    "Mmp9",
    "Mmp14"
  ),
  
  HYPOXIA_HIF1 = c(
    "Hif1a",
    "Epas1",
    "Vegfa",
    "Slc2a1",
    "Ldha",
    "Pdk1",
    "Bnip3",
    "Ndrg1"
  ),
  
  APOPTOSIS_PRO = c(
    "Bax",
    "Bak1",
    "Bbc3",
    "Pmaip1",
    "Casp3",
    "Casp7",
    "Casp8",
    "Casp9",
    "Fas"
  ),
  
  APOPTOSIS_SURVIVAL = c(
    "Bcl2",
    "Bcl2l1",
    "Mcl1",
    "Birc2",
    "Birc3",
    "Xiap"
  ),
  
  MITOCHONDRIAL_BIOGENESIS = c(
    "Ppargc1a",
    "Ppargc1b",
    "Tfam",
    "Nrf1",
    "Esrra"
  ),
  
  FATTY_ACID_OXIDATION = c(
    "Ppara",
    "Ppard",
    "Cpt1a",
    "Cpt1b",
    "Cpt2",
    "Acadm",
    "Acadvl",
    "Hadha",
    "Hadhb",
    "Echs1"
  ),
  
  PPAR_SIGNALING = c(
    "Ppara",
    "Ppard",
    "Pparg",
    "Rxra",
    "Rxrb",
    "Ppargc1a",
    "Ppargc1b",
    "Fabp3",
    "Fabp4",
    "Cd36",
    "Cpt1a",
    "Cpt1b",
    "Acox1"
  ),
  
  ESTROGEN_SIGNALING = c(
    "Esr1",
    "Esr2",
    "Gper1",
    "Ncoa1",
    "Ncoa2",
    "Ncoa3",
    "Nrip1",
    "GreB1"
  ),
  
  MAPK_STRESS = c(
    "Mapk1",
    "Mapk3",
    "Mapk8",
    "Mapk9",
    "Mapk14",
    "Jun",
    "Fos",
    "Atf2"
  ),
  
  JAK_STAT = c(
    "Jak1",
    "Jak2",
    "Stat1",
    "Stat3",
    "Stat5a",
    "Stat5b",
    "Socs1",
    "Socs3"
  )
)


# ============================================================
# 10. CHECK WHICH GENES ARE PRESENT
# ============================================================

all_genes <- rownames(obj)

gene_set_qc <- lapply(
  names(program_genes),
  function(program_name) {
    
    requested <- unique(
      program_genes[[program_name]]
    )
    
    present <- requested[
      requested %in% all_genes
    ]
    
    missing <- requested[
      !requested %in% all_genes
    ]
    
    data.frame(
      Program = program_name,
      Requested_Genes = length(requested),
      Present_Genes = length(present),
      Missing_Genes = length(missing),
      Percent_Present = round(
        100 * length(present) / length(requested),
        2
      ),
      Present_Gene_List = paste(
        present,
        collapse = ";"
      ),
      Missing_Gene_List = paste(
        missing,
        collapse = ";"
      ),
      stringsAsFactors = FALSE
    )
  }
)

gene_set_qc <- bind_rows(
  gene_set_qc
)

write.csv(
  gene_set_qc,
  file.path(
    results_dir,
    "09_Gene_Set_QC.csv"
  ),
  row.names = FALSE
)

cat("\nGene-set coverage:\n")
print(
  gene_set_qc[
    ,
    c(
      "Program",
      "Requested_Genes",
      "Present_Genes",
      "Percent_Present"
    )
  ]
)


# ============================================================
# 11. CREATE FILTERED GENE SETS
# ============================================================

program_genes_present <- lapply(
  program_genes,
  function(x) {
    unique(
      x[x %in% all_genes]
    )
  }
)

program_genes_present <- program_genes_present[
  lengths(program_genes_present) >= 3
]

if (length(program_genes_present) == 0) {
  stop(
    "No biological program retained at least 3 genes."
  )
}

cat(
  "\nPrograms retained for scoring:",
  length(program_genes_present),
  "\n"
)

print(
  names(program_genes_present)
)


# ============================================================
# ============================================================
# 12. SCORE BIOLOGICAL PROGRAMS
# ============================================================
#
# AddModuleScore creates one score per program.
# Each program is scored independently.
# ============================================================

score_columns <- character(0)

for (i in seq_along(program_genes_present)) {
  
  program_name <- names(program_genes_present)[i]
  
  genes_i <- program_genes_present[[program_name]]
  
  cat(
    "\nScoring:",
    program_name,
    "| genes:",
    length(genes_i),
    "\n"
  )
  
  temporary_prefix <- paste0("TMP09_", i, "_")
  
  obj <- AddModuleScore(
    object = obj,
    features = list(genes_i),
    assay = "RNA",
    name = temporary_prefix,
    seed = 12345 + i
  )
  
  temporary_column <- paste0(temporary_prefix, "1")
  final_column <- paste0("Program09_", program_name)
  
  if (!(temporary_column %in% colnames(obj@meta.data))) {
    stop(
      paste0(
        "Expected AddModuleScore column was not created: ",
        temporary_column
      )
    )
  }
  
  obj@meta.data[[final_column]] <- obj@meta.data[[temporary_column]]
  
  obj@meta.data[[temporary_column]] <- NULL
  
  score_columns <- c(
    score_columns,
    final_column
  )
}

cat("\n========================================\n")
cat("PROGRAM SCORING COMPLETE\n")
cat("========================================\n")

cat("\nTotal programs scored:", length(score_columns), "\n")

cat("\nScore columns created:\n")
print(score_columns)

cat("\nChecking for NA values:\n")
print(
  colSums(
    is.na(
      obj@meta.data[, score_columns, drop = FALSE]
    )
  )
)
# ============================================================
# 13. VERIFY SCORE COLUMNS
# ============================================================

missing_scores <- setdiff(
  paste0(
    "Program09_",
    names(program_genes_present)
  ),
  colnames(obj@meta.data)
)

if (length(missing_scores) > 0) {
  stop(
    paste(
      "Missing score columns:",
      paste(missing_scores, collapse = ", ")
    )
  )
}

cat("\nAll program scores verified.\n")


# ============================================================
# 14. CREATE CELL-LEVEL SCORE TABLE
# ============================================================

cell_score_table <- obj@meta.data %>%
  tibble::rownames_to_column(
    "Cell"
  ) %>%
  select(
    Cell,
    sample_id,
    condition_09,
    refined_cell_type,
    all_of(score_columns)
  )

write.csv(
  cell_score_table,
  file.path(
    results_dir,
    "09_Cell_Level_Program_Scores.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. SUMMARIZE SCORES BY CELL TYPE
# ============================================================

celltype_summary <- obj@meta.data %>%
  group_by(
    refined_cell_type
  ) %>%
  summarise(
    Cells = n(),
    across(
      all_of(score_columns),
      list(
        Mean = ~ mean(
          .x,
          na.rm = TRUE
        ),
        Median = ~ median(
          .x,
          na.rm = TRUE
        )
      )
    ),
    .groups = "drop"
  )

write.csv(
  celltype_summary,
  file.path(
    results_dir,
    "09_Program_Scores_By_Cell_Type.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. SUMMARIZE SCORES BY CONDITION
# ============================================================

condition_summary <- obj@meta.data %>%
  group_by(
    condition_09
  ) %>%
  summarise(
    Cells = n(),
    across(
      all_of(score_columns),
      list(
        Mean = ~ mean(
          .x,
          na.rm = TRUE
        ),
        Median = ~ median(
          .x,
          na.rm = TRUE
        )
      )
    ),
    .groups = "drop"
  )

write.csv(
  condition_summary,
  file.path(
    results_dir,
    "09_Program_Scores_By_Condition.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. SUMMARIZE BY CELL TYPE AND CONDITION
# ============================================================

celltype_condition_summary <- obj@meta.data %>%
  group_by(
    refined_cell_type,
    condition_09
  ) %>%
  summarise(
    Cells = n(),
    across(
      all_of(score_columns),
      list(
        Mean = ~ mean(
          .x,
          na.rm = TRUE
        ),
        Median = ~ median(
          .x,
          na.rm = TRUE
        )
      )
    ),
    .groups = "drop"
  )

write.csv(
  celltype_condition_summary,
  file.path(
    results_dir,
    "09_Program_Scores_By_Cell_Type_Condition.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. CREATE LONG-FORM MEAN SCORE TABLE
# ============================================================

mean_columns <- paste0(
  score_columns,
  "_Mean"
)

celltype_condition_long <- celltype_condition_summary %>%
  select(
    refined_cell_type,
    condition_09,
    Cells,
    all_of(mean_columns)
  ) %>%
  pivot_longer(
    cols = all_of(mean_columns),
    names_to = "Program",
    values_to = "Mean_Score"
  )

celltype_condition_long$Program <- gsub(
  "^Program09_",
  "",
  celltype_condition_long$Program
)

celltype_condition_long$Program <- gsub(
  "_Mean$",
  "",
  celltype_condition_long$Program
)

write.csv(
  celltype_condition_long,
  file.path(