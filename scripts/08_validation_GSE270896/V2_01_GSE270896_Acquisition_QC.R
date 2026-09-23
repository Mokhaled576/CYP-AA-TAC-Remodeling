# ============================================================
# VALIDATION #2 — GSE270896
# V2-01: ACQUISITION + SAMPLE AUDIT + RAPID QC
# ============================================================
#
# Study:
#   Whole-heart snRNA-seq
#   Male C57BL/6J mice
#   10 weeks after TAC or Sham
#
# INCLUDED:
#   GSM8354334 = WT Sham replicate 1
#   GSM8354335 = WT Sham replicate 2
#   GSM8354336 = WT TAC replicate 1
#   GSM8354337 = WT TAC replicate 2
#
# EXCLUDED:
#   ob/ob Sham/TAC samples
#
# IMPORTANT:
#   - Biological replicate = mouse/library
#   - No CYP-driven filtering
#   - No integration
#   - No clustering
#   - No annotation
#   - No DE
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
  timeout = 3600
)

set.seed(20260921)


# ============================================================
# 2. PACKAGES
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
})


# ============================================================
# 3. PROJECT DIRECTORIES
# ============================================================

project_root <-
  "D:/Master/ScRNA seq/TAC/VALIDATION/GSE270896"

raw_dir <-
  file.path(
    project_root,
    "RAW DATA"
  )

clean_dir <-
  file.path(
    project_root,
    "CLEAN DATA"
  )

scripts_dir <-
  file.path(
    project_root,
    "SCRIPTS"
  )

results_dir <-
  file.path(
    project_root,
    "RESULTS",
    "V2_01_ACQUISITION_QC"
  )

figures_dir <-
  file.path(
    project_root,
    "FIGURES",
    "V2_01_ACQUISITION_QC",
    "PNG"
  )

dirs <- c(
  raw_dir,
  clean_dir,
  scripts_dir,
  results_dir,
  figures_dir
)

for (x in dirs) {
  dir.create(
    x,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ============================================================
# 4. LOCK SAMPLE DESIGN
# ============================================================

sample_table <- data.frame(
  
  sample_id = c(
    "Sham_Rep1",
    "Sham_Rep2",
    "TAC_Rep1",
    "TAC_Rep2"
  ),
  
  gsm = c(
    "GSM8354334",
    "GSM8354335",
    "GSM8354336",
    "GSM8354337"
  ),
  
  geo_prefix = c(
    "WT_S1",
    "WT_S2",
    "WT_T1",
    "WT_T2"
  ),
  
  condition = c(
    "Sham",
    "Sham",
    "TAC",
    "TAC"
  ),
  
  replicate = c(
    1,
    2,
    1,
    2
  ),
  
  genotype = rep(
    "C57BL/6J_WT",
    4
  ),
  
  tissue = rep(
    "Whole_heart",
    4
  ),
  
  modality = rep(
    "snRNA-seq",
    4
  ),
  
  timepoint = rep(
    "10_weeks_post_surgery",
    4
  ),
  
  stringsAsFactors = FALSE
)

write.csv(
  sample_table,
  file.path(
    results_dir,
    "V2_01_Locked_Sample_Design.csv"
  ),
  row.names = FALSE
)

cat("\nLOCKED SAMPLE DESIGN\n")
print(
  sample_table,
  row.names = FALSE
)


# ============================================================
# 5. LOCK PREDEFINED VALIDATION GENES
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

write.csv(
  data.frame(
    Gene = priority_genes
  ),
  file.path(
    results_dir,
    "V2_01_Predefined_Priority_Genes.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 6. GEO DOWNLOAD INFORMATION
# ============================================================

base_url <-
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8354nnn"

file_manifest <- list()

counter <- 1

for (i in seq_len(nrow(sample_table))) {
  
  gsm <-
    sample_table$gsm[i]
  
  prefix <-
    sample_table$geo_prefix[i]
  
  gsm_folder <-
    paste0(
      substr(gsm, 1, 7),
      "nnn"
    )
  
  sample_url <-
    paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/samples/",
      gsm_folder,
      "/",
      gsm,
      "/suppl/"
    )
  
  files <- c(
    paste0(
      gsm,
      "_",
      prefix,
      "_barcodes.tsv.gz"
    ),
    paste0(
      gsm,
      "_",
      prefix,
      "_features.tsv.gz"
    ),
    paste0(
      gsm,
      "_",
      prefix,
      "_matrix.mtx.gz"
    )
  )
  
  types <- c(
    "barcodes",
    "features",
    "matrix"
  )
  
  for (j in seq_along(files)) {
    
    file_manifest[[counter]] <-
      data.frame(
        sample_id =
          sample_table$sample_id[i],
        
        gsm = gsm,
        
        geo_prefix = prefix,
        
        file_type =
          types[j],
        
        filename =
          files[j],
        
        url =
          paste0(
            sample_url,
            files[j]
          ),
        
        local_path =
          file.path(
            raw_dir,
            files[j]
          ),
        
        stringsAsFactors = FALSE
      )
    
    counter <- counter + 1
  }
}

file_manifest <-
  bind_rows(
    file_manifest
  )

write.csv(
  file_manifest,
  file.path(
    results_dir,
    "V2_01_File_Manifest.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 7. DOWNLOAD FILES
# ============================================================

cat("\n============================================\n")
cat("DOWNLOADING GEO MATRICES\n")
cat("============================================\n")

for (i in seq_len(nrow(file_manifest))) {
  
  dest <-
    file_manifest$local_path[i]
  
  url <-
    file_manifest$url[i]
  
  if (
    file.exists(dest) &&
    file.info(dest)$size > 0
  ) {
    
    cat(
      "\nAlready present:",
      basename(dest),
      "\n"
    )
    
  } else {
    
    cat(
      "\nDownloading:",
      basename(dest),
      "\n"
    )
    
    download.file(
      url = url,
      destfile = dest,
      mode = "wb",
      quiet = FALSE
    )
  }
}


# ============================================================
# 8. VERIFY DOWNLOADS
# ============================================================

file_manifest$Exists <-
  file.exists(
    file_manifest$local_path
  )

file_manifest$Size_Bytes <-
  ifelse(
    file_manifest$Exists,
    file.info(
      file_manifest$local_path
    )$size,
    NA
  )

file_manifest$Nonempty <-
  file_manifest$Exists &
  file_manifest$Size_Bytes > 0

write.csv(
  file_manifest,
  file.path(
    results_dir,
    "V2_01_Download_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(file_manifest$Nonempty)) {
  
  print(
    file_manifest[
      !file_manifest$Nonempty,
      ,
      drop = FALSE
    ]
  )
  
  stop(
    "One or more GEO files failed to download."
  )
}


# ============================================================
# 9. READ FOUR MATRICES
# ============================================================

objects_raw <- list()

for (i in seq_len(nrow(sample_table))) {
  
  sid <-
    sample_table$sample_id[i]
  
  gsm <-
    sample_table$gsm[i]
  
  prefix <-
    sample_table$geo_prefix[i]
  
  matrix_file <-
    file.path(
      raw_dir,
      paste0(
        gsm,
        "_",
        prefix,
        "_matrix.mtx.gz"
      )
    )
  
  feature_file <-
    file.path(
      raw_dir,
      paste0(
        gsm,
        "_",
        prefix,
        "_features.tsv.gz"
      )
    )
  
  barcode_file <-
    file.path(
      raw_dir,
      paste0(
        gsm,
        "_",
        prefix,
        "_barcodes.tsv.gz"
      )
    )
  
  cat(
    "\nReading:",
    sid,
    "\n"
  )
  
  mat <-
    ReadMtx(
      mtx = matrix_file,
      features = feature_file,
      cells = barcode_file,
      feature.column = 2,
      cell.column = 1,
      unique.features = TRUE
    )
  
  obj <-
    CreateSeuratObject(
      counts = mat,
      project = sid,
      min.cells = 0,
      min.features = 0
    )
  
  colnames(obj) <-
    paste0(
      sid,
      "_",
      colnames(obj)
    )
  
  obj$sample_id <-
    sid
  
  obj$gsm <-
    gsm
  
  obj$condition <-
    sample_table$condition[i]
  
  obj$replicate <-
    sample_table$replicate[i]
  
  obj$genotype <-
    "C57BL/6J_WT"
  
  obj$tissue <-
    "Whole_heart"
  
  obj$modality <-
    "snRNA-seq"
  
  obj$timepoint <-
    "10_weeks_post_surgery"
  
  objects_raw[[sid]] <-
    obj
  
  rm(
    mat,
    obj
  )
  
  gc()
}


# ============================================================
# 10. RAW MATRIX AUDIT
# ============================================================

raw_audit <- bind_rows(
  lapply(
    names(objects_raw),
    function(sid) {
      
      x <-
        objects_raw[[sid]]
      
      data.frame(
        
        Sample = sid,
        
        Condition =
          unique(
            as.character(
              x$condition
            )
          ),
        
        Genes =
          nrow(x),
        
        Barcodes =
          ncol(x),
        
        Median_nFeature =
          median(
            x$nFeature_RNA
          ),
        
        Median_nCount =
          median(
            x$nCount_RNA
          ),
        
        Barcodes_nCount_ge_100 =
          sum(
            x$nCount_RNA >= 100
          ),
        
        Barcodes_nCount_ge_500 =
          sum(
            x$nCount_RNA >= 500
          ),
        
        Barcodes_nCount_ge_1000 =
          sum(
            x$nCount_RNA >= 1000
          ),
        
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  raw_audit,
  file.path(
    results_dir,
    "V2_01_Raw_Matrix_Audit.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. SHARED GENE UNIVERSE
# ============================================================

gene_sets <-
  lapply(
    objects_raw,
    rownames
  )

shared_genes <-
  Reduce(
    intersect,
    gene_sets
  )

gene_universe_audit <-
  data.frame(
    Sample =
      names(objects_raw),
    
    Genes =
      sapply(
        objects_raw,
        nrow
      ),
    
    Shared_Genes =
      length(shared_genes),
    
    Identical_to_First =
      sapply(
        objects_raw,
        function(x) {
          identical(
            rownames(x),
            rownames(
              objects_raw[[1]]
            )
          )
        }
      ),
    
    stringsAsFactors = FALSE
  )

write.csv(
  gene_universe_audit,
  file.path(
    results_dir,
    "V2_01_Gene_Universe_Audit.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. PRIORITY GENE AVAILABILITY
# ============================================================

priority_availability <-
  data.frame(
    Gene =
      priority_genes,
    
    Present =
      priority_genes %in%
      shared_genes,
    
    stringsAsFactors = FALSE
  )

write.csv(
  priority_availability,
  file.path(
    results_dir,
    "V2_01_Priority_Gene_Availability.csv"
  ),
  row.names = FALSE
)

if (!all(priority_availability$Present)) {
  
  warning(
    "One or more predefined genes are absent."
  )
}


# ============================================================
# 13. STRICT CYP FAMILY AUDIT
# ============================================================

strict_cyp_genes <-
  grep(
    "^Cyp[0-9]",
    shared_genes,
    value = TRUE
  )

write.csv(
  data.frame(
    Gene = strict_cyp_genes
  ),
  file.path(
    results_dir,
    "V2_01_Strict_CYP_Genes.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. MITOCHONDRIAL NAMING AUDIT
# ============================================================

mito_patterns <- c(
  "^mt-",
  "^Mt-",
  "^MT-"
)

mito_audit <- data.frame()

for (pat in mito_patterns) {
  
  hits <-
    grep(
      pat,
      shared_genes,
      value = TRUE
    )
  
  mito_audit <-
    bind_rows(
      mito_audit,
      data.frame(
        Pattern = pat,
        N_Genes =
          length(hits),
        Example =
          paste(
            head(
              hits,
              10
            ),
            collapse = "; "
          ),
        stringsAsFactors = FALSE
      )
    )
}

write.csv(
  mito_audit,
  file.path(
    results_dir,
    "V2_01_Mitochondrial_Naming_Audit.csv"
  ),
  row.names = FALSE
)

best_mito_pattern <-
  mito_audit$Pattern[
    which.max(
      mito_audit$N_Genes
    )
  ]

cat(
  "\nSelected mitochondrial pattern:",
  best_mito_pattern,
  "\n"
)


# ============================================================
# 15. CALCULATE QC METRICS
# ============================================================

for (sid in names(objects_raw)) {
  
  objects_raw[[sid]][["percent.mt"]] <-
    PercentageFeatureSet(
      objects_raw[[sid]],
      pattern =
        best_mito_pattern
    )
}


# ============================================================
# 16. RAW QC QUANTILES
# ============================================================

qc_quantile_list <- list()

counter <- 1

for (sid in names(objects_raw)) {
  
  x <-
    objects_raw[[sid]]
  
  for (
    metric in c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt"
    )
  ) {
    
    q <-
      quantile(
        x[[metric]][, 1],
        probs = c(
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
        ),
        na.rm = TRUE
      )
    
    qc_quantile_list[[counter]] <-
      data.frame(
        Sample = sid,
        Metric = metric,
        Quantile =
          names(q),
        Value =
          as.numeric(q),
        stringsAsFactors = FALSE
      )
    
    counter <- counter + 1
  }
}

qc_quantiles <-
  bind_rows(
    qc_quantile_list
  )

write.csv(
  qc_quantiles,
  file.path(
    results_dir,
    "V2_01_Raw_QC_Quantiles.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. DETERMINE WHETHER MATRICES ARE ALREADY CELL-FILTERED
# ============================================================
#
# GEO/Cell Ranger supplementary matrices are often filtered
# feature-barcode matrices, but we verify from the data rather
# than assuming.
#
# If the median library has extremely low RNA complexity,
# STOP before applying biological QC.
#
# ============================================================

matrix_status <-
  raw_audit %>%
  mutate(
    Likely_Cell_Filtered =
      Median_nFeature >= 100 &
      Median_nCount >= 200
  )

write.csv(
  matrix_status,
  file.path(
    results_dir,
    "V2_01_Matrix_Filtering_Status.csv"
  ),
  row.names = FALSE
)

all_filtered <-
  all(
    matrix_status$Likely_Cell_Filtered
  )

cat(
  "\nMatrices appear already cell-filtered:",
  all_filtered,
  "\n"
)

if (!all_filtered) {
  
  saveRDS(
    objects_raw,
    file.path(
      clean_dir,
      "V2_01_GSE270896_RAW_OBJECTS_NEED_BARCODE_RECOVERY.rds"
    ),
    compress = FALSE
  )
  
  stop(
    paste(
      "At least one matrix does not look cell-filtered.",
      "Do NOT apply the QC section below.",
      "Send me the console output and we will perform",
      "barcode recovery/EmptyDrops first."
    )
  )
}


# ============================================================
# 18. CONSERVATIVE snRNA QC THRESHOLDS
# ============================================================
#
# Lower:
#   nFeature >= 200
#   nCount   >= 300
#
# Upper:
#   median + 5*MAD independently per library
#
# Mito:
#   nuclei generally have low mitochondrial RNA.
#   We use sample-specific median + 3*MAD,
#   constrained between 5% and 15%.
#
# No CYP/AA genes influence these thresholds.
#
# ============================================================

threshold_list <- list()

counter <- 1

for (sid in names(objects_raw)) {
  
  x <-
    objects_raw[[sid]]@meta.data
  
  feature_upper <-
    median(
      x$nFeature_RNA
    ) +
    5 *
    mad(
      x$nFeature_RNA
    )
  
  count_upper <-
    median(
      x$nCount_RNA
    ) +
    5 *
    mad(
      x$nCount_RNA
    )
  
  mt_raw <-
    median(
      x$percent.mt
    ) +
    3 *
    mad(
      x$percent.mt
    )
  
  mt_final <-
    max(
      5,
      min(
        15,
        mt_raw
      )