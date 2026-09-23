# ============================================================
# V04 — GSE155882
# NORMALIZATION + CONSENSUS HVGs + PCA +
# UNINTEGRATED ATLAS STRUCTURE
#
# INPUT:
#   V03_GSE155882_SINGLET_OBJECTS.rds
#
# PURPOSE:
#   1. Load frozen V03 singlets
#   2. Preserve four biological libraries separately
#   3. Log-normalize each library independently
#   4. Identify 2,000 HVGs independently per library
#   5. Construct a replicate-aware consensus HVG set
#   6. Merge normalized libraries
#   7. Scale consensus HVGs
#   8. Run PCA
#   9. Audit PC variance / choose exploratory PC range
#  10. Run UNINTEGRATED UMAP
#  11. Examine sample/condition structure BEFORE integration
#  12. Save all important figures as 300-dpi PNG
#
# IMPORTANT:
#   - NO cells removed
#   - NO doublet detection rerun
#   - NO CYP/AA-driven feature selection
#   - NO batch integration
#   - NO clustering
#   - NO cell-type annotation
#   - NO differential expression
#   - NO condition-level statistical testing
#
# V04 is a technical atlas-construction step.
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(Seurat)
library(Matrix)
library(ggplot2)


# ============================================================
# 2. PROJECT PATHS
# ============================================================

project_dir <- "D:/Master/ScRNA seq/TAC/VALIDATION/GSE155882"

clean_dir <- file.path(
  project_dir,
  "CLEAN DATA"
)

results_dir <- file.path(
  project_dir,
  "RESULTS",
  "V04_NORMALIZATION_PCA_UNINTEGRATED"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "V04_NORMALIZATION_PCA_UNINTEGRATED"
)

png_dir <- file.path(
  figures_dir,
  "PNG"
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


input_file <- file.path(
  clean_dir,
  "V03_GSE155882_SINGLET_OBJECTS.rds"
)

output_file <- file.path(
  clean_dir,
  "V04_GSE155882_NORMALIZED_PCA_UNINTEGRATED_ATLAS.rds"
)


# ============================================================
# 3. LOAD FROZEN V03 SINGLET OBJECTS
# ============================================================

if (!file.exists(input_file)) {
  stop(
    paste0(
      "V03 singlet object not found:\n",
      input_file
    )
  )
}

singlet_objects <- readRDS(
  input_file
)

expected_samples <- c(
  "Sham_Rep1",
  "Sham_Rep2",
  "TAC_Rep1",
  "TAC_Rep2"
)

if (!identical(
  sort(names(singlet_objects)),
  sort(expected_samples)
)) {
  stop(
    "Unexpected sample structure in V03 object."
  )
}


cat("\n============================================\n")
cat("V04 — NORMALIZATION + PCA + UNINTEGRATED ATLAS\n")
cat("============================================\n")


input_counts <- sapply(
  singlet_objects,
  ncol
)

cat(
  "Cells entering V04:",
  sum(input_counts),
  "\n"
)

print(
  data.frame(
    Sample_ID = names(input_counts),
    Cells = as.integer(input_counts),
    row.names = NULL
  )
)


# ============================================================
# 4. PRE-V04 INTEGRITY AUDIT
# ============================================================

pre_audit_list <- list()

for (sample_name in names(singlet_objects)) {
  
  obj <- singlet_objects[[sample_name]]
  
  pre_audit_list[[sample_name]] <- data.frame(
    Sample_ID = sample_name,
    Condition = unique(obj$condition),
    Replicate = unique(obj$replicate),
    Cells = ncol(obj),
    Genes = nrow(obj),
    Median_nFeature_RNA = median(
      obj$nFeature_RNA,
      na.rm = TRUE
    ),
    Median_nCount_RNA = median(
      obj$nCount_RNA,
      na.rm = TRUE
    ),
    Median_percent_mt = median(
      obj$percent.mt,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
}

pre_audit <- do.call(
  rbind,
  pre_audit_list
)

rownames(pre_audit) <- NULL

write.csv(
  pre_audit,
  file.path(
    results_dir,
    "V04_Pre_Normalization_Audit.csv"
  ),
  row.names = FALSE
)

cat("\nPRE-NORMALIZATION AUDIT\n")

print(
  pre_audit,
  row.names = FALSE
)


# ============================================================
# 5. VERIFY ALL V03 CELLS ARE SINGLET
# ============================================================

all_singlets <- all(
  sapply(
    singlet_objects,
    function(obj) {
      all(
        obj$scDblFinder_class ==
          "singlet"
      )
    }
  )
)

if (!all_singlets) {
  stop(
    "V04 input contains cells not classified as singlets."
  )
}


# ============================================================
# 6. NORMALIZE EACH LIBRARY INDEPENDENTLY
#
# Standard log-normalization:
#   counts / total counts * 10,000
#   followed by log1p transformation.
#
# No regression is performed here.
#
# We deliberately do NOT regress:
#   condition
#   sample
#   nCount_RNA
#   percent.mt
#
# because doing so could remove genuine biology.
# ============================================================

normalized_objects <- list()

for (sample_name in names(singlet_objects)) {
  
  cat(
    "\nNormalizing:",
    sample_name,
    "\n"
  )
  
  obj <- singlet_objects[[sample_name]]
  
  DefaultAssay(obj) <- "RNA"
  
  obj <- NormalizeData(
    object = obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  
  obj <- FindVariableFeatures(
    object = obj,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  
  normalized_objects[[sample_name]] <- obj
}


# ============================================================
# 7. AUDIT HVGs BY LIBRARY
# ============================================================

hvg_lists <- lapply(
  normalized_objects,
  VariableFeatures
)

hvg_counts <- sapply(
  hvg_lists,
  length
)

cat("\n============================================\n")
cat("HVG COUNTS BY LIBRARY\n")
cat("============================================\n")

print(hvg_counts)


if (!all(hvg_counts == 2000)) {
  warning(
    "At least one library did not return exactly 2,000 HVGs."
  )
}


# ============================================================
# 8. BUILD REPLICATE-AWARE CONSENSUS HVG TABLE
#
# We record how many of the four libraries independently
# selected each gene as variable.
#
# Primary consensus:
#   genes selected in >=2 libraries.
#
# If this yields fewer than 2,000 genes, we supplement using
# Seurat SelectIntegrationFeatures ranking.
#
# This does NOT perform integration.
# ============================================================

all_hvgs <- unique(
  unlist(
    hvg_lists,
    use.names = FALSE
  )
)

hvg_frequency <- sapply(
  all_hvgs,
  function(gene_now) {
    sum(
      sapply(
        hvg_lists,
        function(x) {
          gene_now %in% x
        }
      )
    )
  }
)

hvg_frequency_table <- data.frame(
  Gene = names(hvg_frequency),
  Libraries_Selected = as.integer(
    hvg_frequency
  ),
  stringsAsFactors = FALSE
)

hvg_frequency_table <- hvg_frequency_table[
  order(
    -hvg_frequency_table$Libraries_Selected,
    hvg_frequency_table$Gene
  ),
]

rownames(
  hvg_frequency_table
) <- NULL


write.csv(
  hvg_frequency_table,
  file.path(
    results_dir,
    "V04_HVG_Replicate_Frequency.csv"
  ),
  row.names = FALSE
)


cat("\n============================================\n")
cat("HVG REPLICATE FREQUENCY\n")
cat("============================================\n")

print(
  table(
    hvg_frequency_table$Libraries_Selected
  )
)


consensus_hvgs <- hvg_frequency_table$Gene[
  hvg_frequency_table$Libraries_Selected >= 2
]


cat(
  "\nGenes selected as HVG in >=2 libraries:",
  length(consensus_hvgs),
  "\n"
)


# ============================================================
# 9. SELECT FINAL 2,000 CONSENSUS FEATURES
#
# Seurat's SelectIntegrationFeatures is used ONLY as a
# replicate-aware feature-ranking utility.
#
# NO integration is performed in V04.
# ============================================================

ranked_features <- SelectIntegrationFeatures(
  object.list = normalized_objects,
  nfeatures = 2000
)


if (length(consensus_hvgs) >= 2000) {
  
  final_hvgs <- ranked_features[
    ranked_features %in%
      consensus_hvgs
  ]
  
  if (length(final_hvgs) < 2000) {
    
    supplement <- setdiff(
      ranked_features,
      final_hvgs
    )
    
    final_hvgs <- c(
      final_hvgs,
      head(
        supplement,
        2000 - length(final_hvgs)
      )
    )
  }
  
  final_hvgs <- head(
    final_hvgs,
    2000
  )
  
} else {
  
  final_hvgs <- ranked_features
}


final_hvgs <- unique(
  final_hvgs
)


if (length(final_hvgs) != 2000) {
  stop(
    paste(
      "Final HVG set contains",
      length(final_hvgs),
      "genes instead of 2000."
    )
  )
}


write.csv(
  data.frame(
    Gene = final_hvgs,
    stringsAsFactors = FALSE
  ),
  file.path(
    results_dir,
    "V04_Final_2000_Consensus_HVGs.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. CYP/AA HVG AUDIT — DESCRIPTIVE ONLY
#
# IMPORTANT:
# We are NOT adding/removing CYP genes from HVGs.
# This merely records whether predefined genes naturally
# entered the unbiased HVG set.
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

priority_hvg_audit <- data.frame(
  Gene = priority_genes,
  In_Final_HVG_Set =
    priority_genes %in%
    final_hvgs,
  stringsAsFactors = FALSE
)

write.csv(
  priority_hvg_audit,
  file.path(
    results_dir,
    "V04_Priority_Gene_HVG_Audit.csv"
  ),
  row.names = FALSE
)

cat("\nPriority-gene HVG audit — descriptive only:\n")

print(
  priority_hvg_audit,
  row.names = FALSE
)


# ============================================================
# 11. MERGE NORMALIZED LIBRARIES
#
# Replicate identity is preserved in metadata.
#
# IMPORTANT:
# No integration occurs here.
# ============================================================

merged_obj <- merge(
  x = normalized_objects[["Sham_Rep1"]],
  y = list(
    normalized_objects[["Sham_Rep2"]],
    normalized_objects[["TAC_Rep1"]],
    normalized_objects[["TAC_Rep2"]]
  ),
  add.cell.ids = c(
    "ShamR1",
    "ShamR2",
    "TACR1",
    "TACR2"
  ),
  project = "GSE155882_Validation"
)


cat(
  "\nMerged cells:",
  ncol(merged_obj),
  "\n"
)

cat(
  "Merged genes:",
  nrow(merged_obj),
  "\n"
)


# ============================================================
# 12. JOIN SEURAT V5 RNA LAYERS
#
# merge() may create sample-specific layers.
# JoinLayers creates unified RNA layers for downstream
# normalization-based analyses.
# ============================================================

merged_obj <- JoinLayers(
  merged_obj,
  assay = "RNA"
)

DefaultAssay(
  merged_obj
) <- "RNA"


# ============================================================
# 13. VERIFY METADATA AFTER MERGE
# ============================================================

required_metadata <- c(
  "sample_id",
  "condition",
  "replicate",
  "nFeature_RNA",
  "nCount_RNA",
  "percent.mt",
  "scDblFinder_score",
  "scDblFinder_class"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(
    merged_obj@meta.data
  )
)

if (length(missing_metadata) > 0) {
  
  stop(
    paste(
      "Missing required metadata after merge:",
      paste(
        missing_metadata,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 14. SET CONSENSUS VARIABLE FEATURES
# ============================================================

VariableFeatures(
  merged_obj
) <- final_hvgs


# ============================================================
# 15. SCALE CONSENSUS HVGs
#
# No covariates are regressed.
# ============================================================

cat("\nScaling consensus HVGs...\n")

merged_obj <- ScaleData(
  object = merged_obj,
  features = final_hvgs,
  verbose = FALSE
)


# ============================================================
# 16. PCA
# ============================================================

cat("\nRunning PCA...\n")

set.seed(20260920)

merged_obj <- RunPCA(
  object = merged_obj,
  features = final_hvgs,
  npcs = 50,
  verbose = FALSE
)


# ============================================================
# 17. PCA VARIANCE AUDIT
# ============================================================

pca_stdev <- Stdev(
  merged_obj,
  reduction = "pca"
)

variance <- pca_stdev^2

percent_variance <- (
  variance /
    sum(variance)
) * 100

cumulative_variance <- cumsum(
  percent_variance
)

pc_table <- data.frame(
  PC = seq_along(
    percent_variance
  ),
  Percent_Variance = percent_variance,
  Cumulative_Variance = cumulative_variance,
  stringsAsFactors = FALSE
)

write.csv(
  pc_table,
  file.path(
    results_dir,
    "V04_PCA_Variance_Explained.csv"
  ),
  row.names = FALSE
)


cat("\n============================================\n")
cat("PCA VARIANCE — FIRST 30 PCs\n")
cat("============================================\n")

print(
  head(
    pc_table,
    30
  ),
  row.names = FALSE
)


# ============================================================
# 18. DATA-SUPPORTED EXPLORATORY PC SELECTION
#
# We use two descriptive criteria:
#
# Criterion A:
#   first PC where cumulative variance >= 90%
#
# Criterion B:
#   last PC contributing >= 0.5% variance
#
# For exploratory UMAP we use the larger of these,
# constrained to 10–40 PCs.
#
# This is NOT a biological hypothesis-driven choice.
# ============================================================

pc_90_candidates <- which(
  cumulative_variance >= 90
)

if (length(pc_90_candidates) > 0) {
  pc_90 <- min(
    pc_90_candidates
  )
} else {
  pc_90 <- length(
    cumulative_variance
  )
}


pc_half_percent_candidates <- which(
  percent_variance >= 0.5
)

if (length(
  pc_half_percent_candidates
) > 0) {
  
  pc_half_percent <- max(
    pc_half_percent_candidates
  )
  
} else {
  
  pc_half_percent <- 10
}


selected_npcs <- max(
  pc_90,
  pc_half_percent
)

selected_npcs <- max(
  10,
  selected_npcs
)

selected_npcs <- min(
  40,
  selected_npcs
)


cat(
  "\nPC reaching >=90% cumulative variance:",
  pc_90,
  "\n"
)

cat(
  "Last PC contributing >=0.5% variance:",
  pc_half_percent,
  "\n"
)

cat(
  "Selected exploratory PC count:",
  selected_npcs,
  "\n"
)


pc_selection <- data.frame(
  Criterion = c(
    "PC_at_90_percent_cumulative_variance",
    "Last_PC_with_at_least_0.5_percent_variance",
    "Final_exploratory_PC_count"
  ),
  Value = c(
    pc_90,
    pc_half_percent,
    selected_npcs
  ),
  stringsAsFactors = FALSE
)

write.csv(
  pc_selection,
  file.path(
    results_dir,
    "V04_PC_Selection_Audit.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. PCA ELBOW PNG
# ============================================================

p_elbow <- ElbowPlot(
  merged_obj,
  ndims = 50
) +
  ggtitle(
    "GSE155882 — PCA Elbow Plot"
  ) +
  theme_classic(
    base_size = 12
  )


ggsave(
  filename = file.path(
    png_dir,
    "V04_PCA_ElbowPlot.png"
  ),
  plot = p_elbow,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)


# ============================================================
# 20. PCA VARIANCE PNG
# ============================================================

p_variance <- ggplot(
  pc_table,
  aes(
    x = PC,
    y = Percent_Variance
  )
) +
  geom_line() +
  geom_point(
    size = 1.5
  ) +
  geom_vline(
    xintercept = selected_npcs,
    linetype = "dashed"
  ) +
  labs(
    title = "Variance Explained by Principal Components",
    x = "Principal component",
    y = "Variance explained (%)"
  ) +
  theme_classic(
    base_size = 12
  )


ggsave(
  filename = file.path(
    png_dir,
    "V04_PCA_Percent_Variance.png"
  ),
  plot = p_variance,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)


# ============================================================
# 21. PCA CUMULATIVE VARIANCE PNG
# ============================================================

p_cumulative <- ggplot(
  pc_table,
  aes(
    x = PC,
    y = Cumulative_Variance
  )
) +
  geom_line() +
  geom_point(
    size = 1.5
  ) +
  geom_hline(
    yintercept = 90,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = selected_npcs,
    linetype = "dashed"
  ) +
  labs(
    title = "Cumulative PCA Variance",
    x = "Principal component",
    y = "Cumulative variance (%)"
  ) +
  theme_classic(
    base_size = 12
  )


ggsave(
  filename = file.path(
    png_dir,
    "V04_PCA_Cumulative_Variance.png"
  ),
  plot = p_cumulative,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)


# ============================================================
# 22. PCA SAMPLE STRUCTURE
# ============================================================

p_pca_sample <- DimPlot(
  merged_obj,
  reduction = "pca",
  group.by = "sample_id",
  pt.size = 0.25
) +
  ggtitle(
    "PCA — Sample Identity"
  ) +
  theme_classic(
    base_size = 12
  )


ggsave(
  filename = file.path(
    png_dir,
    "V04_PCA_by_Sample.png"
  ),
  plot = p_pca_sample,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)


# ============================================================
# 23. PCA CONDITION STRUCTURE
# ============================================================

p_pca_condition <- DimPlot(
  merged_obj,
  reduction = "pca",
  group.by = "condition",
  pt.size = 0.25
) +
  ggtitle(
    "PCA — Condition"
  ) +
  theme_classic(
    base_size = 12
  )


ggsave(
  filename = file.path(
    png_dir,
    "V04_PCA_by_Condition.png"
  ),
  plot = p_pca_condition,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"