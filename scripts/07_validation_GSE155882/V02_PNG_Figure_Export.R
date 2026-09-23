# ============================================================
# GSE155882
# V02A + V02B PNG FIGURE EXPORT
#
# PURPOSE:
#   Recreate diagnostic figures from already-frozen
#   V02A and V02B objects and save them as individual
#   high-resolution PNG files.
#
# IMPORTANT:
#   - DOES NOT rerun EmptyDrops
#   - DOES NOT rerun QC
#   - DOES NOT change any saved object
#   - DOES NOT filter any additional cells
#   - Figure export only
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(Seurat)
library(Matrix)


# ============================================================
# 2. PATHS
# ============================================================

project_dir <- "D:/Master/ScRNA seq/TAC/VALIDATION/GSE155882"

clean_dir <- file.path(
  project_dir,
  "CLEAN DATA"
)

figure_root <- file.path(
  project_dir,
  "FIGURES"
)

v02a_png_dir <- file.path(
  figure_root,
  "V02A_CELL_BARCODE_RECOVERY",
  "PNG"
)

v02b_png_dir <- file.path(
  figure_root,
  "V02B_FINAL_QC",
  "PNG"
)

dir.create(
  v02a_png_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  v02b_png_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. INPUT FILES
# ============================================================

v01_file <- file.path(
  clean_dir,
  "V01_GSE155882_RAW_SEURAT_OBJECTS.rds"
)

v02a_file <- file.path(
  clean_dir,
  "V02A_GSE155882_CANDIDATE_CELL_OBJECTS.rds"
)

v02b_file <- file.path(
  clean_dir,
  "V02B_GSE155882_FINAL_QC_OBJECTS.rds"
)


required_files <- c(
  v01_file,
  v02a_file,
  v02b_file
)


if (!all(file.exists(required_files))) {
  
  missing_files <- required_files[
    !file.exists(required_files)
  ]
  
  stop(
    paste(
      "Missing required file(s):",
      paste(
        missing_files,
        collapse = "\n"
      ),
      sep = "\n"
    )
  )
}


# ============================================================
# 4. LOAD FROZEN OBJECTS
# ============================================================

raw_objects <- readRDS(
  v01_file
)

candidate_objects <- readRDS(
  v02a_file
)

qc_objects <- readRDS(
  v02b_file
)


expected_samples <- c(
  "Sham_Rep1",
  "Sham_Rep2",
  "TAC_Rep1",
  "TAC_Rep2"
)


stopifnot(
  identical(
    sort(names(raw_objects)),
    sort(expected_samples)
  )
)

stopifnot(
  identical(
    sort(names(candidate_objects)),
    sort(expected_samples)
  )
)

stopifnot(
  identical(
    sort(names(qc_objects)),
    sort(expected_samples)
  )
)


cat("\n============================================\n")
cat("V02 PNG FIGURE EXPORT\n")
cat("============================================\n")

cat("Objects loaded successfully.\n")
cat("No analysis will be rerun.\n\n")


# ============================================================
# 5. HELPER — RAW COUNTS
# ============================================================

get_counts <- function(obj) {
  
  GetAssayData(
    obj,
    assay = "RNA",
    layer = "counts"
  )
}


# ============================================================
# 6. V02A — BARCODE RANK PNGs
# ============================================================

cat("Creating V02A barcode-rank PNGs...\n")


for (sample_name in expected_samples) {
  
  counts <- get_counts(
    raw_objects[[sample_name]]
  )
  
  total_umi <- Matrix::colSums(
    counts
  )
  
  total_umi <- sort(
    total_umi,
    decreasing = TRUE
  )
  
  
  output_png <- file.path(
    v02a_png_dir,
    paste0(
      "V02A_",
      sample_name,
      "_Barcode_Rank.png"
    )
  )
  
  
  png(
    filename = output_png,
    width = 7,
    height = 6,
    units = "in",
    res = 300,
    bg = "white"
  )
  
  
  plot(
    seq_along(total_umi),
    total_umi + 1,
    log = "xy",
    pch = 16,
    cex = 0.25,
    xlab = "Barcode rank",
    ylab = "Total UMI + 1",
    main = paste(
      sample_name,
      "Barcode Rank"
    )
  )
  
  
  abline(
    h = c(
      10,
      50,
      100,
      200,
      500,
      1000
    ),
    lty = 3
  )
  
  
  dev.off()
}


# ============================================================
# 7. V02A — FEATURE vs UMI PNGs
# ============================================================

cat("Creating V02A feature-vs-UMI PNGs...\n")


for (sample_name in expected_samples) {
  
  counts <- get_counts(
    raw_objects[[sample_name]]
  )
  
  
  total_umi <- Matrix::colSums(
    counts
  )
  
  
  detected_genes <- Matrix::colSums(
    counts > 0
  )
  
  
  positive <- total_umi > 0
  
  
  output_png <- file.path(
    v02a_png_dir,
    paste0(
      "V02A_",
      sample_name,
      "_Features_vs_UMI.png"
    )
  )
  
  
  png(
    filename = output_png,
    width = 7,
    height = 6,
    units = "in",
    res = 300,
    bg = "white"
  )
  
  
  plot(
    total_umi[positive],
    detected_genes[positive],
    log = "xy",
    pch = 16,
    cex = 0.25,
    xlab = "Total UMI",
    ylab = "Detected genes",
    main = paste(
      sample_name,
      "Features vs UMI"
    )
  )
  
  
  dev.off()
}


# ============================================================
# 8. V02A — CANDIDATE CELL QC PNGs
#
# One PNG per sample with three panels.
# ============================================================

cat("Creating V02A candidate-cell QC PNGs...\n")


for (sample_name in expected_samples) {
  
  obj <- candidate_objects[[sample_name]]
  
  
  output_png <- file.path(
    v02a_png_dir,
    paste0(
      "V02A_",
      sample_name,
      "_Candidate_Cell_QC.png"
    )
  )
  
  
  png(
    filename = output_png,
    width = 12,
    height = 4,
    units = "in",
    res = 300,
    bg = "white"
  )
  
  
  old_par <- par(
    mfrow = c(
      1,
      3
    )
  )
  
  
  hist(
    obj$nFeature_RNA,
    breaks = 80,
    main = paste(
      sample_name,
      "nFeature_RNA"
    ),
    xlab = "Detected genes"
  )
  
  
  hist(
    obj$nCount_RNA,
    breaks = 80,
    main = paste(
      sample_name,
      "nCount_RNA"
    ),
    xlab = "UMIs"
  )
  
  
  hist(
    obj$percent.mt,
    breaks = 80,
    main = paste(
      sample_name,
      "percent.mt"
    ),
    xlab = "Mitochondrial %"
  )
  
  
  par(
    old_par
  )
  
  
  dev.off()
}


# ============================================================
# 9. V02B — PRE vs POST QC PNGs
#
# One PNG per sample.
# Top row = before final QC
# Bottom row = after final QC
# ============================================================

cat("Creating V02B pre-vs-post QC PNGs...\n")


for (sample_name in expected_samples) {
  
  pre_obj <- candidate_objects[[sample_name]]
  
  post_obj <- qc_objects[[sample_name]]
  
  
  output_png <- file.path(
    v02b_png_dir,
    paste0(
      "V02B_",
      sample_name,
      "_Pre_vs_Post_QC.png"
    )
  )
  
  
  png(
    filename = output_png,
    width = 12,
    height = 8,
    units = "in",
    res = 300,
    bg = "white"
  )
  
  
  old_par <- par(
    mfrow = c(
      2,
      3
    )
  )
  
  
  # --------------------------
  # PRE-QC
  # --------------------------
  
  hist(
    pre_obj$nFeature_RNA,
    breaks = 80,
    main = paste(
      sample_name,
      "PRE - Features"
    ),
    xlab = "nFeature_RNA"
  )
  
  
  hist(
    pre_obj$nCount_RNA,
    breaks = 80,
    main = paste(
      sample_name,
      "PRE - UMIs"
    ),
    xlab = "nCount_RNA"
  )
  
  
  hist(
    pre_obj$percent.mt,
    breaks = 80,
    main = paste(
      sample_name,
      "PRE - mt%"
    ),
    xlab = "percent.mt"
  )
  
  
  # --------------------------
  # POST-QC
  # --------------------------
  
  hist(
    post_obj$nFeature_RNA,
    breaks = 80,
    main = paste(
      sample_name,
      "POST - Features"
    ),
    xlab = "nFeature_RNA"
  )
  
  
  hist(
    post_obj$nCount_RNA,
    breaks = 80,
    main = paste(
      sample_name,
      "POST - UMIs"
    ),
    xlab = "nCount_RNA"
  )
  
  
  hist(
    post_obj$percent.mt,
    breaks = 80,
    main = paste(
      sample_name,
      "POST - mt%"
    ),
    xlab = "percent.mt"
  )
  
  
  par(
    old_par
  )
  
  
  dev.off()
}


# ============================================================
# 10. V02B — POST-QC FEATURE vs COUNT PNGs
# ============================================================

cat("Creating V02B feature-vs-count PNGs...\n")


for (sample_name in expected_samples) {
  
  obj <- qc_objects[[sample_name]]
  
  
  output_png <- file.path(
    v02b_png_dir,
    paste0(
      "V02B_",
      sample_name,
      "_Post_QC_Feature_vs_Count.png"
    )
  )
  
  
  png(
    filename = output_png,
    width = 7,
    height = 6,
    units = "in",
    res = 300,
    bg = "white"
  )
  
  
  plot(
    obj$nCount_RNA,
    obj$nFeature_RNA,
    log = "xy",
    pch = 16,
    cex = 0.35,
    xlab = "nCount_RNA",
    ylab = "nFeature_RNA",
    main = paste(
      sample_name,
      "Post-QC"
    )
  )
  
  
  dev.off()
}


# ============================================================
# 11. V02B — CELL RETENTION SUMMARY PNG
# ============================================================

before_counts <- sapply(
  candidate_objects[
    expected_samples
  ],
  ncol
)


after_counts <- sapply(
  qc_objects[
    expected_samples
  ],
  ncol
)


retention_percent <- (
  100 *
    after_counts /
    before_counts
)


output_png <- file.path(
  v02b_png_dir,
  "V02B_Cell_Retention_Summary.png"
)


png(
  filename = output_png,
  width = 8,
  height = 6,
  units = "in",
  res = 300,
  bg = "white"
)


barplot(
  retention_percent,
  names.arg = expected_samples,
  ylim = c(
    0,
    100
  ),
  ylab = "Cells retained (%)",
  xlab = "Library",
  main = "V02B Final QC Cell Retention",
  las = 2
)


abline(
  h = 100,
  lty = 3
)


dev.off()


# ============================================================
# 12. V02B — FINAL CELL COUNTS PNG
# ============================================================

output_png <- file.path(
  v02b_png_dir,
  "V02B_Final_Cell_Counts.png"
)


png(
  filename = output_png,
  width = 8,
  height = 6,
  units = "in",
  res = 300,
  bg = "white"
)


barplot(
  after_counts,
  names.arg = expected_samples,
  ylab = "Final QC cells",
  xlab = "Library",
  main = "GSE155882 Final QC Cell Counts",
  las = 2
)


dev.off()


# ============================================================
# 13. VERIFY PNG FILES
# ============================================================

v02a_png_files <- list.files(
  v02a_png_dir,
  pattern = "\\.png$",
  full.names = TRUE
)


v02b_png_files <- list.files(
  v02b_png_dir,
  pattern = "\\.png$",
  full.names = TRUE
)


png_audit <- data.frame(
  
  Stage = c(
    rep(
      "V02A",
      length(v02a_png_files)
    ),
    rep(
      "V02B",
      length(v02b_png_files)
    )
  ),
  
  File = c(
    basename(
      v02a_png_files
    ),
    basename(
      v02b_png_files
    )
  ),
  
  Exists = file.exists(
    c(
      v02a_png_files,
      v02b_png_files
    )
  ),
  
  Size_MB = round(
    file.info(
      c(
        v02a_png_files,
        v02b_png_files
      )
    )$size /
      1024^2,
    3
  ),
  
  stringsAsFactors = FALSE
)


cat("\n============================================\n")
cat("PNG EXPORT AUDIT\n")
cat("============================================\n")

print(
  png_audit,
  row.names = FALSE
)


# ============================================================
# 14. FINAL STATUS
# ============================================================

cat("\n\n============================================\n")
cat("PNG EXPORT COMPLETE\n")
cat("============================================\n")

cat(
  "V02A PNG files:",
  length(
    v02a_png_files
  ),
  "\n"
)

cat(
  "V02B PNG files:",
  length(
    v02b_png_files
  ),
  "\n"
)

cat(
  "Total PNG files:",
  length(
    v02a_png_files
  ) +
    length(
      v02b_png_files
    ),
  "\n"
)

cat(
  "\nV02A PNG folder:\n",
  v02a_png_dir,
  "\n"
)

cat(
  "\nV02B PNG folder:\n",
  v02b_png_dir,
  "\n"
)

cat("\nAnalysis objects modified: NO\n")
cat("EmptyDrops rerun: NO\n")
cat("QC rerun: NO\n")
cat("Cells added/removed: NO\n")

cat("============================================\n")