# ============================================================
# STEP 11D-G
# MECHANISTIC PATHWAY INTEGRATION OF CYP/AA REMODELING
#
# PURPOSE
# ------------------------------------------------------------
# Integrate the corrected scRNA-seq atlas and spatial
# transcriptomics to determine which predefined biological
# programs accompany CYP/AA remodeling during TAC progression.
#
# PRIMARY QUESTIONS
# ------------------------------------------------------------
# 1. Which mechanistic programs change across Sham -> 2W -> 4W
#    -> 6W?
#
# 2. Which corrected cell types carry these programs?
#
# 3. Which pathways accompany the 4W CYP/AA spatial state?
#
# 4. Which CYP/AA genes spatially co-vary with these pathways?
#
# IMPORTANT
# ------------------------------------------------------------
# - Pathways are predefined before analysis.
# - CYP/AA genes are excluded from pathway scores.
# - No cells/spots are removed.
# - No reclustering.
# - No annotation changes.
# - No condition-level inferential statistics.
# - Correlations are descriptive within sections.
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
  "11D_F_CELL_PROGRAM_ADJUSTED_CYP_AA_SPATIAL_OBJECTS.rds"
)

results_dir <- file.path(
  project_dir,
  "RESULTS",
  "11_SPATIAL",
  "11D_G_MECHANISTIC_PATHWAY_INTEGRATION"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "11_SPATIAL",
  "11D_G_MECHANISTIC_PATHWAY_INTEGRATION"
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
cat("STEP 11D-G\n")
cat("MECHANISTIC PATHWAY INTEGRATION\n")
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
# 4. BASIC INTEGRITY
# ============================================================

if (ncol(scrna) != 25186) {
  stop("Unexpected scRNA cell count.")
}

if (nrow(scrna) != 32285) {
  stop("Unexpected scRNA gene count.")
}

if (length(unique(scrna$seurat_clusters)) != 29) {
  stop("Expected 29 frozen clusters.")
}

if (length(unique(scrna$corrected_broad_cell_type)) != 15) {
  stop("Expected 15 corrected broad cell types.")
}

if (length(spatial_list) != 4) {
  stop("Expected four spatial objects.")
}

if (sum(sapply(spatial_list, ncol)) != 7079) {
  stop("Expected 7079 spatial spots.")
}


# ============================================================
# 5. CONDITION STANDARDIZATION
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

if (!"sample_id" %in% colnames(scrna@meta.data)) {
  stop("sample_id missing from scRNA metadata.")
}

scrna$Condition <- unname(
  sample_to_condition[
    as.character(scrna$sample_id)
  ]
)

if (any(is.na(scrna$Condition))) {
  stop("Could not map all scRNA sample IDs.")
}

scrna$Condition <- factor(
  scrna$Condition,
  levels = condition_levels
)

for (i in seq_along(spatial_list)) {
  
  if (!"Condition" %in%
      colnames(spatial_list[[i]]@meta.data)) {
    stop("Condition missing from spatial metadata.")
  }
  
  spatial_list[[i]]$Condition <- factor(
    as.character(
      spatial_list[[i]]$Condition
    ),
    levels = condition_levels
  )
}


# ============================================================
# 6. PRIORITY CYP/AA GENES
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

if (!all(priority_genes %in% rownames(scrna))) {
  stop("One or more priority genes missing from scRNA.")
}

if (!all(
  priority_genes %in%
  Reduce(
    intersect,
    lapply(spatial_list, rownames)
  )
)) {
  stop("One or more priority genes missing from spatial data.")
}


# ============================================================
# 7. PREDEFINED MECHANISTIC GENE SETS
#
# These are deliberately compact and interpretable.
#
# They are NOT intended to represent every gene in each
# biological pathway.
# ============================================================

mechanistic_sets <- list(
  
  Fibrosis_ECM = c(
    "Col1a1",
    "Col1a2",
    "Col3a1",
    "Postn",
    "Fn1",
    "Ccn2",
    "Lox",
    "Timp1",
    "Thbs1",
    "Sparc"
  ),
  
  TGFb_Response = c(
    "Tgfb1",
    "Tgfbr1",
    "Tgfbr2",
    "Smad2",
    "Smad3",
    "Smad7",
    "Serpine1",
    "Ccn2",
    "Thbs1",
    "Timp1"
  ),
  
  NFkB_Inflammation = c(
    "Nfkb1",
    "Rela",
    "Nfkbia",
    "Tnf",
    "Il1b",
    "Il6",
    "Ccl2",
    "Cxcl2",
    "Ptgs2",
    "Icam1"
  ),
  
  Nrf2_Antioxidant = c(
    "Nfe2l2",
    "Keap1",
    "Hmox1",
    "Nqo1",
    "Gclc",
    "Gclm",
    "Txnrd1",
    "Srxn1",
    "Gsta1",
    "Gstm1"
  ),
  
  Hypoxia_HIF = c(
    "Hif1a",
    "Egln1",
    "Egln3",
    "Vegfa",
    "Slc2a1",
    "Pdk1",
    "Ldha",
    "Bnip3",
    "Ndrg1",
    "Adm"
  ),
  
  Mitochondrial_OXPHOS = c(
    "Ndufs1",
    "Ndufs2",
    "Sdha",
    "Sdhb",
    "Uqcrc1",
    "Uqcrc2",
    "Cox4i1",
    "Cox5a",
    "Atp5f1a",
    "Atp5f1b"
  ),
  
  Fatty_Acid_Oxidation = c(
    "Ppara",
    "Ppargc1a",
    "Cpt1a",
    "Cpt1b",
    "Cpt2",
    "Acadm",
    "Acadvl",
    "Hadha",
    "Hadhb",
    "Cd36"
  ),
  
  Lipid_Peroxidation_Ferroptosis = c(
    "Acsl4",
    "Lpcat3",
    "Alox15",
    "Alox5",
    "Slc7a11",
    "Gpx4",
    "Fth1",
    "Ftl1",
    "Hmox1",
    "Ncoa4"
  )
)


# ============================================================
# 8. REMOVE CYP/AA TARGET GENES FROM MECHANISTIC SETS
#
# Prevent circularity.
# ============================================================

strict_cyp_genes <- grep(
  "^Cyp[0-9]",
  rownames(scrna),
  value = TRUE
)

aa_exclusion <- unique(
  c(
    strict_cyp_genes,
    priority_genes,
    "Ephx1",
    "Ephx2",
    "Pla2g4a",
    "Pla2g6",
    "Ptgs1",
    "Ptgs2",
    "Alox5",
    "Alox12",
    "Alox15",
    "Alox15b"
  )
)

mechanistic_sets <- lapply(
  mechanistic_sets,
  function(x) {
    setdiff(
      x,
      aa_exclusion
    )
  }
)


# ============================================================
# 9. REQUIRE GENES PRESENT IN BOTH MODALITIES
# ============================================================

shared_all_genes <- Reduce(
  intersect,
  c(
    list(rownames(scrna)),
    lapply(
      spatial_list,
      rownames
    )
  )
)

mechanistic_sets_present <- lapply(
  mechanistic_sets,
  function(x) {
    intersect(
      x,
      shared_all_genes
    )
  }
)

gene_set_table <- bind_rows(
  lapply(
    names(mechanistic_sets_present),
    function(set_name) {
      
      data.frame(
        Pathway = set_name,
        Gene = mechanistic_sets_present[[set_name]],
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  gene_set_table,
  file.path(
    results_dir,
    "11D_G_Final_Mechanistic_Gene_Sets.csv"
  ),
  row.names = FALSE
)

gene_set_summary <- data.frame(
  Pathway = names(mechanistic_sets_present),
  Genes_Available = sapply(
    mechanistic_sets_present,
    length
  ),
  stringsAsFactors = FALSE
)

cat("\n============================================\n")
cat("FINAL MECHANISTIC GENE SETS\n")
cat("============================================\n")

print(
  gene_set_summary,
  row.names = FALSE
)

if (any(gene_set_summary$Genes_Available < 5)) {
  warning(
    "At least one mechanistic gene set contains fewer than five genes."
  )
}


# ============================================================
# 10. HELPER:
# MEAN NORMALIZED EXPRESSION SCORE
#
# We use transparent mean normalized expression rather than
# AddModuleScore so the same definition can be applied to both
# scRNA and spatial data.
# ============================================================

calculate_pathway_scores <- function(
    expression_matrix,
    gene_sets
) {
  
  score_list <- lapply(
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
  
  score_df <- as.data.frame(
    score_list,
    check.names = FALSE
  )
  
  colnames(score_df) <- paste0(
    "Mechanism_",
    names(gene_sets)
  )
  
  score_df
}


# ============================================================
# 11. SCORE scRNA
# ============================================================

DefaultAssay(scrna) <- "RNA"

scrna_data <- GetAssayData(
  scrna,
  assay = "RNA",
  layer = "data"
)

scrna_mechanism_scores <-
  calculate_pathway_scores(
    expression_matrix = scrna_data,
    gene_sets = mechanistic_sets_present
  )

rownames(scrna_mechanism_scores) <-
  colnames(scrna)

for (score_name in
     colnames(scrna_mechanism_scores)) {
  
  scrna[[score_name]] <-
    scrna_mechanism_scores[
      colnames(scrna),
      score_name
    ]
}


# ============================================================
# 12. scRNA PATHWAY SUMMARY:
# CONDITION x BROAD CELL TYPE
# ============================================================

mechanism_columns <- colnames(
  scrna_mechanism_scores
)

scrna_meta <- scrna@meta.data %>%
  mutate(
    Cell = rownames(scrna@meta.data)
  )

scrna_long <- scrna_meta %>%
  
  select(
    Cell,
    Condition,
    corrected_broad_cell_type,
    all_of(mechanism_columns)
  ) %>%
  
  pivot_longer(
    cols = all_of(mechanism_columns),
    names_to = "Mechanism",
    values_to = "Score"
  ) %>%
  
  mutate(
    Mechanism =
      sub(
        "^Mechanism_",
        "",
        Mechanism
      )
  )

scrna_summary <- scrna_long %>%
  
  group_by(
    Condition,
    corrected_broad_cell_type,
    Mechanism
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
    
    Q25 =
      quantile(
        Score,
        0.25,
        na.rm = TRUE
      ),
    
    Q75 =
      quantile(
        Score,
        0.75,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

write.csv(
  scrna_summary,
  file.path(
    results_dir,
    "11D_G_scRNA_Mechanism_By_Condition_CellType.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. scRNA WHOLE-ATLAS CONDITION SUMMARY
# ============================================================

scrna_condition_summary <- scrna_long %>%
  
  group_by(
    Condition,
    Mechanism
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
  scrna_condition_summary,
  file.path(
    results_dir,
    "11D_G_scRNA_WholeAtlas_Mechanism_Trajectories.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. scRNA CHANGE vs SHAM WITHIN EACH CELL TYPE
# ============================================================

scrna_sham <- scrna_summary %>%
  
  filter(
    Condition == "Sham"
  ) %>%
  
  select(
    corrected_broad_cell_type,
    Mechanism,
    Sham_Mean_Score = Mean_Score
  )

scrna_change <- scrna_summary %>%
  
  left_join(
    scrna_sham,
    by = c(
      "corrected_broad_cell_type",
      "Mechanism"
    )
  ) %>%
  
  mutate(
    
    Delta_vs_Sham =
      Mean_Score -
      Sham_Mean_Score,
    
    Fold_vs_Sham =
      (Mean_Score + 1e-8) /
      (Sham_Mean_Score + 1e-8),
    
    Log2_Fold_vs_Sham =
      log2(
        Fold_vs_Sham
      )
  )

write.csv(
  scrna_change,
  file.path(
    results_dir,
    "11D_G_scRNA_Mechanism_Change_vs_Sham.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. PRINT STRONGEST 4W scRNA CHANGES
#
# Only populations with >=30 cells in that condition.
# Ranking is descriptive.
# ============================================================

scrna_4w_changes <- scrna_change %>%
  
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
  scrna_4w_changes,
  file.path(
    results_dir,
    "11D_G_scRNA_4W_Mechanistic_State.csv"
  ),
  row.names = FALSE
)

cat("\n============================================\n")
cat("TOP 4W scRNA MECHANISTIC CHANGES\n")
cat("============================================\n")

print(
  as.data.frame(
    head(
      scrna_4w_changes,
      40
    )
  ),
  row.names = FALSE
)


# ============================================================
# 16. FIGURE:
# 4W scRNA CELL-TYPE x MECHANISM HEATMAP
#
# Delta from Sham.
# ============================================================

plot_scrna_4w <- scrna_4w_changes %>%
  
  filter(
    Cells >= 30
  )

p_scrna_4w <- ggplot(
  plot_scrna_4w,
  aes(
    x = Mechanism,
    y = corrected_broad_cell_type,
    fill = Delta_vs_Sham
  )
) +
  geom_tile() +
  labs(
    title =
      "Cell-type-specific mechanistic remodeling at 4W TAC",
    subtitle =
      "Mean pathway score difference from Sham; descriptive only",
    x = NULL,
    y = NULL,
    fill = "Δ score"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )

ggsave(
  file.path(
    figures_dir,
    "11D_G_scRNA_4W_CellType_Mechanism_Heatmap.png"
  ),
  p_scrna_4w,
  width = 11,
  height = 8,
  dpi = 400
)


# ============================================================
# 17. SCORE SPATIAL DATA
# ============================================================

spatial_scored <- spatial_list

spatial_summary_list <- list()

spatial_counter <- 1

for (i in seq_along(spatial_scored)) {
  
  obj <- spatial_scored[[i]]
  
  DefaultAssay(obj) <- "Spatial"
  
  spatial_data <- GetAssayData(
    obj,
    assay = "Spatial",
    layer = "data"
  )
  
  scores_now <- calculate_pathway_scores(
    expression_matrix = spatial_data,
    gene_sets = mechanistic_sets_present
  )
  
  rownames(scores_now) <- colnames(obj)
  
  for (score_name in colnames(scores_now)) {
    
    obj[[score_name]] <-
      scores_now[
        colnames(obj),
        score_name
      ]
  }
  
  condition_now <- unique(
    as.character(
      obj$Condition
    )
  )
  
  if (length(condition_now) != 1) {
    stop("Non-unique spatial condition.")
  }
  
  summary_now <- data.frame(
    Condition = condition_now,
    Mechanism = sub(
      "^Mechanism_",
      "",
      colnames(scores_now)
    ),
    Mean_Score = sapply(
      scores_now,
      mean,
      na.rm = TRUE
    ),
    Median_Score = sapply(
      scores_now,
      median,
      na.rm = TRUE
    ),
    SD_Score = sapply(
      scores_now,
      sd,
      na.rm = TRUE
    ),
    Q25 = sapply(
      scores_now,
      quantile,
      probs = 0.25,
      na.rm = TRUE
    ),
    Q75 = sapply(
      scores_now,
      quantile,
      probs = 0.75,
      na.rm = TRUE
    ),
    Spots = ncol(obj),
    stringsAsFactors = FALSE
  )
  
  spatial_summary_list[[spatial_counter]] <-
    summary_now
  
  spatial_counter <- spatial_counter + 1
  
  spatial_scored[[i]] <- obj
}


# ============================================================
# 18. SPATIAL MECHANISM SUMMARY
# ============================================================

spatial_summary <- bind_rows(
  spatial_summary_list
)

spatial_summary$Condition <- factor(
  spatial_summary$Condition,
  levels = condition_levels
)

spatial_summary <- spatial_summary %>%
  arrange(
    Mechanism,
    Condition
  )

write.csv(
  spatial_summary,
  file.path(
    results_dir,
    "11D_G_Spatial_Mechanism_Trajectories.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. SPATIAL CHANGE vs SHAM
# ============================================================

spatial_sham <- spatial_summary %>%
  
  filter(
    Condition == "Sham"
  ) %>%
  
  select(
    Mechanism,
    Sham_Mean_Score = Mean_Score
  )

spatial_change <- spatial_summary %>%
  
  left_join(
    spatial_sham,
    by = "Mechanism"
  ) %>%
  
  mutate(
    
    Delta_vs_Sham =
      Mean_Score -
      Sham_Mean_Score,
    
    Fold_vs_Sham =
      (Mean_Score + 1e-8) /
      (Sham_Mean_Score + 1e-8),
    
    Log2_Fold_vs_Sham =
      log2(
        Fold_vs_Sham
      )
  )

write.csv(
  spatial_change,
  file.path(
    results_dir,
    "11D_G_Spatial_Mechanism_Change_vs_Sham.csv"
  ),