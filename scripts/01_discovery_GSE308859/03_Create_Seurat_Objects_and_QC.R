# ============================================================
# GSE308859 TAC scRNA-seq PROJECT
# STEP 03: CREATE SEURAT OBJECTS + PRE-QC ASSESSMENT
#
# Samples:
#   Sham
#   TAC 2 weeks
#   TAC 4 weeks
#   TAC 6 weeks
#
# IMPORTANT:
#   NO CELLS ARE FILTERED IN THIS SCRIPT.
#   This script is diagnostic only.
# ============================================================


# ============================================================
# 0. CLEAN ENVIRONMENT
# ============================================================

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)


# ============================================================
# 1. PROJECT DIRECTORIES
# ============================================================

project_dir <- "D:/Master/ScRNA seq/TAC"

raw_dir <- file.path(
  project_dir,
  "RAW DATA",
  "GSE308859",
  "GSE308859_RAW"
)

clean_dir <- file.path(
  project_dir,
  "CLEAN DATA"
)

results_dir <- file.path(
  project_dir,
  "RESULTS",
  "03_PRE_QC"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "03_PRE_QC"
)

notes_dir <- file.path(
  project_dir,
  "NOTES"
)


# Create folders if they do not exist

dirs_to_create <- c(
  clean_dir,
  results_dir,
  figures_dir,
  notes_dir
)

for (d in dirs_to_create) {
  
  if (!dir.exists(d)) {
    
    dir.create(
      d,
      recursive = TRUE
    )
  }
}


cat("\nProject directory:\n")
cat(project_dir, "\n")

cat("\nRaw data directory exists:\n")
print(dir.exists(raw_dir))


# ============================================================
# 2. INSTALL / LOAD REQUIRED PACKAGES
# ============================================================

required_packages <- c(
  "Seurat",
  "Matrix",
  "ggplot2",
  "patchwork"
)


for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    install.packages(pkg)
  }
}


library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)


cat("\n========================================\n")
cat("PACKAGE VERSIONS\n")
cat("========================================\n")

cat("R version:", R.version.string, "\n")
cat("Seurat:", as.character(packageVersion("Seurat")), "\n")
cat("Matrix:", as.character(packageVersion("Matrix")), "\n")
cat("ggplot2:", as.character(packageVersion("ggplot2")), "\n")
cat("patchwork:", as.character(packageVersion("patchwork")), "\n")


# ============================================================
# 3. DEFINE THE FOUR SAMPLES
# ============================================================

samples <- data.frame(
  
  sample_id = c(
    "Sham",
    "TAC_2W",
    "TAC_4W",
    "TAC_6W"
  ),
  
  prefix = c(
    "GSM9254695_Sham_scRNA",
    "GSM9254692_TAC2w_scRNA",
    "GSM9254693_TAC4w_scRNA",
    "GSM9254694_TAC6w_scRNA"
  ),
  
  timepoint = c(
    "Sham",
    "2W",
    "4W",
    "6W"
  ),
  
  condition = c(
    "Sham",
    "TAC",
    "TAC",
    "TAC"
  ),
  
  stringsAsFactors = FALSE
)


print(samples)


# ============================================================
# 4. FUNCTION TO READ ONE SAMPLE
# ============================================================

read_10x_sample <- function(prefix, sample_id) {
  
  cat("\n\n========================================\n")
  cat("READING SAMPLE:", sample_id, "\n")
  cat("========================================\n")
  
  matrix_file <- file.path(
    raw_dir,
    paste0(prefix, "_matrix.mtx.gz")
  )
  
  feature_file <- file.path(
    raw_dir,
    paste0(prefix, "_features.tsv.gz")
  )
  
  barcode_file <- file.path(
    raw_dir,
    paste0(prefix, "_barcodes.tsv.gz")
  )
  
  
  # Check all files exist
  
  if (!file.exists(matrix_file)) {
    stop("Matrix file not found: ", matrix_file)
  }
  
  if (!file.exists(feature_file)) {
    stop("Feature file not found: ", feature_file)
  }
  
  if (!file.exists(barcode_file)) {
    stop("Barcode file not found: ", barcode_file)
  }
  
  
  # Read sparse count matrix
  
  counts <- Matrix::readMM(matrix_file)
  
  
  # Read feature information
  
  features <- read.delim(
    feature_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  
  # Read cell barcodes
  
  barcodes <- read.delim(
    barcode_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  
  # Verify dimensions
  
  if (nrow(counts) != nrow(features)) {
    
    stop(
      "Feature number does not match matrix rows for ",
      sample_id
    )
  }
  
  
  if (ncol(counts) != nrow(barcodes)) {
    
    stop(
      "Barcode number does not match matrix columns for ",
      sample_id
    )
  }
  
  
  # Use gene symbols as row names
  
  gene_symbols <- features[[2]]
  
  rownames(counts) <- make.unique(gene_symbols)
  
  
  # Make cell names sample-specific immediately
  
  cell_names <- paste0(
    sample_id,
    "_",
    barcodes[[1]]
  )
  
  colnames(counts) <- cell_names
  
  
  cat("Genes:", nrow(counts), "\n")
  cat("Cells:", ncol(counts), "\n")
  cat(
    "Non-zero matrix entries:",
    length(counts@x),
    "\n"
  )
  
  
  return(counts)
}


# ============================================================
# 5. CREATE SEURAT OBJECTS
# ============================================================

seurat_list <- list()


for (i in seq_len(nrow(samples))) {
  
  sample_name <- samples$sample_id[i]
  
  counts <- read_10x_sample(
    prefix = samples$prefix[i],
    sample_id = sample_name
  )
  
  
  # IMPORTANT:
  # min.cells = 0
  # min.features = 0
  #
  # We deliberately do NOT filter cells here.
  
  obj <- CreateSeuratObject(
    counts = counts,
    project = sample_name,
    min.cells = 0,
    min.features = 0
  )
  
  
  # Add metadata
  
  obj$sample_id <- sample_name
  
  obj$condition <- samples$condition[i]
  
  obj$timepoint <- samples$timepoint[i]
  
  
  seurat_list[[sample_name]] <- obj
  
  
  rm(counts, obj)
  
  gc()
}


# ============================================================
# 6. CHECK CELL NUMBERS
# ============================================================

cat("\n\n========================================\n")
cat("RAW CELL NUMBERS\n")
cat("========================================\n")


raw_cell_numbers <- data.frame(
  
  Sample = names(seurat_list),
  
  Cells = sapply(
    seurat_list,
    ncol
  ),
  
  stringsAsFactors = FALSE
)


print(raw_cell_numbers)


write.csv(
  raw_cell_numbers,
  file.path(
    results_dir,
    "03_Raw_Cell_Numbers.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 7. CALCULATE MITOCHONDRIAL PERCENTAGE
# ============================================================

# Mouse mitochondrial genes normally begin with "mt-".
#
# We calculate this separately for each sample.

for (sample_name in names(seurat_list)) {
  
  seurat_list[[sample_name]][["percent.mt"]] <-
    PercentageFeatureSet(
      seurat_list[[sample_name]],
      pattern = "^mt-",
      assay = "RNA"
    )
}


# ============================================================
# 8. CALCULATE RIBOSOMAL PERCENTAGE
# ============================================================

for (sample_name in names(seurat_list)) {
  
  seurat_list[[sample_name]][["percent.ribo"]] <-
    PercentageFeatureSet(
      seurat_list[[sample_name]],
      pattern = "^Rp[sl]",
      assay = "RNA"
    )
}


# ============================================================
# 9. CALCULATE HEMOGLOBIN PERCENTAGE
# ============================================================

# Useful diagnostic for blood/erythrocyte-associated signal.

for (sample_name in names(seurat_list)) {
  
  seurat_list[[sample_name]][["percent.hb"]] <-
    PercentageFeatureSet(
      seurat_list[[sample_name]],
      pattern = "^Hb[ab]",
      assay = "RNA"
    )
}


# ============================================================
# 10. VERIFY MITOCHONDRIAL GENES
# ============================================================

mito_genes <- grep(
  "^mt-",
  rownames(seurat_list[[1]]),
  value = TRUE
)


cat("\n\n========================================\n")
cat("MITOCHONDRIAL GENE CHECK\n")
cat("========================================\n")

cat(
  "Number of genes matching ^mt-:",
  length(mito_genes),
  "\n"
)

print(mito_genes)


# ============================================================
# 11. CREATE QC SUMMARY TABLE
# ============================================================

qc_summary_list <- list()


for (sample_name in names(seurat_list)) {
  
  md <- seurat_list[[sample_name]]@meta.data
  
  
  qc_summary_list[[sample_name]] <- data.frame(
    
    Sample = sample_name,
    
    Cells = nrow(md),
    
    Median_nFeature_RNA =
      median(md$nFeature_RNA),
    
    Mean_nFeature_RNA =
      mean(md$nFeature_RNA),
    
    Q1_nFeature_RNA =
      quantile(md$nFeature_RNA, 0.25),
    
    Q3_nFeature_RNA =
      quantile(md$nFeature_RNA, 0.75),
    
    Median_nCount_RNA =
      median(md$nCount_RNA),
    
    Mean_nCount_RNA =
      mean(md$nCount_RNA),
    
    Q1_nCount_RNA =
      quantile(md$nCount_RNA, 0.25),
    
    Q3_nCount_RNA =
      quantile(md$nCount_RNA, 0.75),
    
    Median_percent_mt =
      median(md$percent.mt),
    
    Mean_percent_mt =
      mean(md$percent.mt),
    
    Q1_percent_mt =
      quantile(md$percent.mt, 0.25),
    
    Q3_percent_mt =
      quantile(md$percent.mt, 0.75),
    
    Median_percent_ribo =
      median(md$percent.ribo),
    
    Median_percent_hb =
      median(md$percent.hb),
    
    stringsAsFactors = FALSE
  )
}


qc_summary <- do.call(
  rbind,
  qc_summary_list
)

rownames(qc_summary) <- NULL


cat("\n\n========================================\n")
cat("PRE-QC SUMMARY\n")
cat("========================================\n")

print(qc_summary)


write.csv(
  qc_summary,
  file.path(
    results_dir,
    "03_PreQC_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. MORE DETAILED QC QUANTILES
# ============================================================

quantile_probs <- c(
  0,
  0.01,
  0.05,
  0.10,
  0.25,
  0.50,
  0.75,
  0.90,
  0.95,
  0.99,
  1
)


detailed_qc <- list()


for (sample_name in names(seurat_list)) {
  
  md <- seurat_list[[sample_name]]@meta.data
  
  
  temp <- data.frame(
    
    Sample = sample_name,
    
    Quantile = quantile_probs,
    
    nFeature_RNA =
      as.numeric(
        quantile(
          md$nFeature_RNA,
          probs = quantile_probs
        )
      ),
    
    nCount_RNA =
      as.numeric(
        quantile(
          md$nCount_RNA,
          probs = quantile_probs
        )
      ),
    
    percent_mt =
      as.numeric(
        quantile(
          md$percent.mt,
          probs = quantile_probs
        )
      ),
    
    percent_ribo =
      as.numeric(
        quantile(
          md$percent.ribo,
          probs = quantile_probs
        )
      ),
    
    percent_hb =
      as.numeric(
        quantile(
          md$percent.hb,
          probs = quantile_probs
        )
      )
  )
  
  
  detailed_qc[[sample_name]] <- temp
}


detailed_qc <- do.call(
  rbind,
  detailed_qc
)

rownames(detailed_qc) <- NULL


write.csv(
  detailed_qc,
  file.path(
    results_dir,
    "03_PreQC_Detailed_Quantiles.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. COMBINE METADATA FOR PLOTTING
# ============================================================

metadata_all <- do.call(
  
  rbind,
  
  lapply(
    
    names(seurat_list),
    
    function(sample_name) {
      
      x <- seurat_list[[sample_name]]@meta.data
      
      x$Cell <- rownames(x)
      
      x
      
    }
  )
)


metadata_all$sample_id <- factor(
  metadata_all$sample_id,
  levels = c(
    "Sham",
    "TAC_2W",
    "TAC_4W",
    "TAC_6W"
  )
)


# Save full pre-QC cell metadata

write.csv(
  metadata_all,
  file.path(
    results_dir,
    "03_All_Cells_PreQC_Metadata.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. VIOLIN PLOT — nFeature_RNA
# ============================================================

p_feature <- ggplot(
  metadata_all,
  aes(
    x = sample_id,
    y = nFeature_RNA
  )
) +
  geom_violin(
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA
  ) +
  labs(
    title = "Detected genes per cell before QC",
    x = NULL,
    y = "nFeature_RNA"
  ) +
  theme_classic(base_size = 13)


ggsave(
  filename = file.path(
    figures_dir,
    "03_PreQC_nFeature_RNA_Violin.png"
  ),
  plot = p_feature,
  width = 8,
  height = 6,
  dpi = 600
)


# ============================================================
# 15. VIOLIN PLOT — nCount_RNA
# ============================================================

p_count <- ggplot(
  metadata_all,
  aes(
    x = sample_id,
    y = nCount_RNA
  )
) +
  geom_violin(
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA
  ) +
  labs(
    title = "UMI counts per cell before QC",
    x = NULL,
    y = "nCount_RNA"
  ) +
  theme_classic(base_size = 13)


ggsave(
  filename = file.path(
    figures_dir,
    "03_PreQC_nCount_RNA_Violin.png"
  ),
  plot = p_count,
  width = 8,
  height = 6,
  dpi = 600
)


# ============================================================
# 16. VIOLIN PLOT — MITOCHONDRIAL %
# ============================================================

p_mt <- ggplot(
  metadata_all,
  aes(
    x = sample_id,
    y = percent.mt
  )
) +
  geom_violin(
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA
  ) +
  labs(
    title = "Mitochondrial RNA percentage before QC",
    x = NULL,
    y = "Mitochondrial RNA (%)"
  ) +
  theme_classic(base_size = 13)


ggsave(
  filename = file.path(
    figures_dir,
    "03_PreQC_Percent_Mitochondrial_Violin.png"
  ),
  plot = p_mt,
  width = 8,
  height = 6,
  dpi = 600
)


# ============================================================
# 17. VIOLIN PLOT — RIBOSOMAL %
# ============================================================

p_ribo <- ggplot(
  metadata_all,
  aes(
    x = sample_id,
    y = percent.ribo
  )
) +
  geom_violin(
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA
  ) +
  labs(
    title = "Ribosomal RNA-associated percentage before QC",
    x = NULL,
    y = "Ribosomal genes (%)"
  ) +
  theme_classic(base_size = 13)


ggsave(
  filename = file.path(
    figures_dir,
    "03_PreQC_Percent_Ribosomal_Violin.png"
  ),
  plot = p_ribo,
  width = 8,
  height = 6,
  dpi = 600
)


# ============================================================
# 18. VIOLIN PLOT — HEMOGLOBIN %
# ============================================================

p_hb <- ggplot(
  metadata_all,
  aes(
    x = sample_id,
    y = percent.hb
  )
) +
  geom_violin(
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA
  ) +
  labs(
    title = "Hemoglobin-associated RNA before QC",
    x = NULL,
    y = "Hemoglobin genes (%)"
  ) +
  theme_classic(base_size = 13)


ggsave(
  filename = file.path(
    figures_dir,
    "03_PreQC_Percent_Hemoglobin_Violin.png"
  ),
  plot = p_hb,
  width = 8,
  height = 6,
  dpi = 600
)


# ============================================================
# 19. SCATTER: COUNTS vs FEATURES
# ============================================================

p_count_feature <- ggplot(
  metadata_all,
  aes(
    x = nCount_RNA,
    y = nFeature_RNA
  )
) +
  geom_point(
    alpha = 0.25,
    size = 0.4
  ) +
  facet_wrap(
    ~ sample_id,
    scales = "free"
  ) +
  labs(
    title = "UMI counts versus detected genes before QC",
    x = "nCount_RNA",
    y = "nFeature_RNA"
  ) +
  theme_classic(base_size = 12)


ggsave(
  filename = file.path(
    figures_dir,
    "03_PreQC_Counts_vs_Features.png"
  ),
  plot = p_count_feature,
  width = 10,
  height = 8,
  dpi = 600
)


# ============================================================
# 20. SCATTER: COUNTS vs MITOCHONDRIAL %
# ============================================================

p_count_mt <- ggplot(
  metadata_all,
  aes(
    x = nCount_RNA,
    y = percent.mt
  )
) +
  geom_point(
    alpha = 0.25,
    size = 0.4
  ) +
  facet_wrap(
    ~ sample_id,
    scales = "free_x"
  ) +
  labs(
    title = "UMI counts versus mitochondrial RNA before QC",
    x = "nCount_RNA",
    y = "Mitochondrial RNA (%)"
  ) +
  theme_classic(base_size = 12)


ggsave(
  filename = file.path(
    figures_dir,
    "03_PreQC_Counts_vs_Mitochondrial.png"
  ),
  plot = p_count_mt,
  width = 10,
  height = 8,
  dpi = 600
)


# ============================================================
# 21. HISTOGRAMS FOR EACH QC METRIC
# ============================================================

p_feature_hist <- ggplot(
  metadata_all,
  aes(x = nFeature_RNA)
) +
  geom_histogram(
    bins = 100
  ) +
  facet_wrap(
    ~ sample_id,
    scales = "free_y"
  ) +
  labs(
    title = "Distribution of detected genes before QC",
    x = "nFeature_RNA",
    y = "Number of cells"
  ) +
  theme_classic(base_size = 12)


ggsave(
  file.path(
    figures_dir,
    "03_PreQC_nFeature_Histogram.png"
  ),
  p_feature_hist,
  width = 10,
  height = 7,
  dpi = 600
)


p_count_hist <- ggplot(
  metadata_all,
  aes(x = nCount_RNA)
) +
  geom_histogram(
    bins = 100
  ) +
  facet_wrap(
    ~ sample_id,
    scales = "free_y"
  ) +
  labs(
    title = "Distribution of UMI counts before QC",
    x = "nCount_RNA",
    y = "Number of cells"
  ) +
  theme_classic(base_size = 12)

