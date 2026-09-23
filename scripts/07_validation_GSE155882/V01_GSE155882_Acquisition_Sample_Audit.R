# ============================================================
# V01 — GSE155882
# INDEPENDENT TAC VALIDATION
# DATA ACQUISITION + SAMPLE AUDIT + RAW OBJECT CONSTRUCTION
#
# PRIMARY VALIDATION:
#   Sham rep1
#   Sham rep2
#   TAC rep1
#   TAC rep2
#
# IMPORTANT:
#   - No QC filtering in V01
#   - No clustering
#   - No integration
#   - No differential expression
#   - No CYP-based filtering
#   - Replicates remain separate
# ============================================================


# ============================================================
# 1. PACKAGE
# ============================================================

library(Seurat)
library(Matrix)


# ============================================================
# 2. PROJECT PATHS
# ============================================================

project_dir <-
  "D:/Master/ScRNA seq/TAC/VALIDATION/GSE155882"

raw_dir <-
  file.path(project_dir, "RAW DATA")

clean_dir <-
  file.path(project_dir, "CLEAN DATA")

scripts_dir <-
  file.path(project_dir, "SCRIPTS")

results_dir <-
  file.path(project_dir, "RESULTS")

figures_dir <-
  file.path(project_dir, "FIGURES")


dirs <- c(
  project_dir,
  raw_dir,
  clean_dir,
  scripts_dir,
  results_dir,
  figures_dir
)

for (d in dirs) {
  
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ============================================================
# 3. LOCK PRIMARY VALIDATION SAMPLES
# ============================================================

sample_manifest <- data.frame(
  
  GSM = c(
    "GSM4715045",
    "GSM4715046",
    "GSM4715047",
    "GSM4715048"
  ),
  
  Sample_ID = c(
    "Sham_Rep1",
    "Sham_Rep2",
    "TAC_Rep1",
    "TAC_Rep2"
  ),
  
  Condition = c(
    "Sham",
    "Sham",
    "TAC",
    "TAC"
  ),
  
  Replicate = c(
    1,
    2,
    1,
    2
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  sample_manifest,
  file.path(
    results_dir,
    "V01_Primary_Validation_Sample_Manifest.csv"
  ),
  row.names = FALSE
)


cat("\n============================================\n")
cat("GSE155882 PRIMARY VALIDATION SAMPLES\n")
cat("============================================\n")

print(
  sample_manifest,
  row.names = FALSE
)


# ============================================================
# 4. GEO FILE NAMES
# ============================================================

file_manifest <- data.frame(
  
  GSM = rep(
    sample_manifest$GSM,
    each = 3
  ),
  
  Sample_ID = rep(
    sample_manifest$Sample_ID,
    each = 3
  ),
  
  File_Type = rep(
    c(
      "barcodes",
      "features",
      "matrix"
    ),
    times = 4
  ),
  
  stringsAsFactors = FALSE
)


prefix_map <- c(
  
  GSM4715045 =
    "Sham_scRNA_rep1",
  
  GSM4715046 =
    "Sham_scRNA_rep2",
  
  GSM4715047 =
    "TAC_scRNA_rep1",
  
  GSM4715048 =
    "TAC_scRNA_rep2"
)


file_manifest$Filename <- NA_character_


for (i in seq_len(nrow(file_manifest))) {
  
  gsm_now <-
    file_manifest$GSM[i]
  
  prefix_now <-
    prefix_map[[gsm_now]]
  
  type_now <-
    file_manifest$File_Type[i]
  
  extension_now <-
    ifelse(
      type_now == "matrix",
      "matrix.mtx.gz",
      paste0(type_now, ".tsv.gz")
    )
  
  file_manifest$Filename[i] <-
    paste0(
      gsm_now,
      "_",
      prefix_now,
      "_",
      extension_now
    )
}


# ============================================================
# 5. GEO DOWNLOAD URLs
# ============================================================

file_manifest$URL <- paste0(
  
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/",
  
  substr(
    file_manifest$GSM,
    1,
    7
  ),
  
  "nnn/",
  
  file_manifest$GSM,
  
  "/suppl/",
  
  file_manifest$Filename
)


file_manifest$Local_Path <- file.path(
  raw_dir,
  file_manifest$Filename
)


write.csv(
  file_manifest,
  file.path(
    results_dir,
    "V01_GEO_File_Manifest.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 6. DOWNLOAD FILES
# ============================================================

cat("\n============================================\n")
cat("DOWNLOADING GSE155882 scRNA FILES\n")
cat("============================================\n")


options(timeout = 3600)


for (i in seq_len(nrow(file_manifest))) {
  
  destination <-
    file_manifest$Local_Path[i]
  
  url_now <-
    file_manifest$URL[i]
  
  if (
    file.exists(destination) &&
    file.info(destination)$size > 0
  ) {
    
    cat(
      "\nAlready exists:",
      basename(destination),
      "\n"
    )
    
    next
  }
  
  
  cat(
    "\nDownloading:",
    basename(destination),
    "\n"
  )
  
  
  download.file(
    url = url_now,
    destfile = destination,
    mode = "wb",
    quiet = FALSE
  )
}


# ============================================================
# 7. VERIFY DOWNLOADS
# ============================================================

file_manifest$Exists <-
  file.exists(
    file_manifest$Local_Path
  )


file_manifest$Size_MB <-
  round(
    file.info(
      file_manifest$Local_Path
    )$size /
      1024^2,
    3
  )


cat("\n============================================\n")
cat("DOWNLOAD AUDIT\n")
cat("============================================\n")

print(
  file_manifest[
    ,
    c(
      "GSM",
      "Sample_ID",
      "File_Type",
      "Exists",
      "Size_MB"
    )
  ],
  row.names = FALSE
)


write.csv(
  file_manifest,
  file.path(
    results_dir,
    "V01_Download_Audit.csv"
  ),
  row.names = FALSE
)


if (!all(file_manifest$Exists)) {
  
  stop(
    "One or more GEO files were not downloaded."
  )
}


if (any(file_manifest$Size_MB <= 0)) {
  
  stop(
    "One or more downloaded files are empty."
  )
}


# ============================================================
# 8. HELPER TO READ ONE SAMPLE
# ============================================================

read_geo_sample <- function(
    gsm,
    sample_id,
    condition,
    replicate_number
) {
  
  sample_files <-
    file_manifest[
      file_manifest$GSM == gsm,
      ,
      drop = FALSE
    ]
  
  
  barcode_file <-
    sample_files$Local_Path[
      sample_files$File_Type ==
        "barcodes"
    ]
  
  
  feature_file <-
    sample_files$Local_Path[
      sample_files$File_Type ==
        "features"
    ]
  
  
  matrix_file <-
    sample_files$Local_Path[
      sample_files$File_Type ==
        "matrix"
    ]
  
  
  cat(
    "\nReading:",
    sample_id,
    "\n"
  )
  
  
  counts <- ReadMtx(
    
    mtx = matrix_file,
    
    cells = barcode_file,
    
    features = feature_file,
    
    feature.column = 2,
    
    unique.features = TRUE
  )
  
  
  obj <- CreateSeuratObject(
    
    counts = counts,
    
    project = "GSE155882",
    
    min.cells = 0,
    
    min.features = 0
  )
  
  
  obj$GSE <-
    "GSE155882"
  
  obj$GSM <-
    gsm
  
  obj$sample_id <-
    sample_id
  
  obj$condition <-
    condition
  
  obj$replicate <-
    replicate_number
  
  obj$validation_role <-
    "Independent_validation"
  
  
  obj[["percent.mt"]] <-
    PercentageFeatureSet(
      obj,
      pattern = "^mt-"
    )
  
  
  return(obj)
}


# ============================================================
# 9. CONSTRUCT FOUR RAW OBJECTS
# ============================================================
# ============================================================
# 9. CONSTRUCT FOUR RAW OBJECTS
# ============================================================

raw_objects <- list()

for (i in seq_len(nrow(sample_manifest))) {
  
  sample_now <- sample_manifest[i, ]
  
  raw_objects[[sample_now$Sample_ID]] <- read_geo_sample(
    gsm = sample_now$GSM,
    sample_id = sample_now$Sample_ID,
    condition = sample_now$Condition,
    replicate_number = sample_now$Replicate
  )
}
# ============================================================
# 10. RAW SAMPLE DIMENSIONS
# ============================================================

raw_dimensions <- data.frame(
  
  Sample_ID =
    names(raw_objects),
  
  Condition =
    sapply(
      raw_objects,
      function(x) {
        unique(x$condition)
      }
    ),
  
  Replicate =
    sapply(
      raw_objects,
      function(x) {
        unique(x$replicate)
      }
    ),
  
  Genes =
    sapply(
      raw_objects,
      nrow
    ),
  
  Raw_Cells =
    sapply(
      raw_objects,
      ncol
    ),
  
  stringsAsFactors = FALSE
)


cat("\n============================================\n")
cat("RAW OBJECT DIMENSIONS\n")
cat("============================================\n")

print(
  raw_dimensions,
  row.names = FALSE
)


write.csv(
  raw_dimensions,
  file.path(
    results_dir,
    "V01_Raw_Object_Dimensions.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. RAW QC SUMMARY — NO FILTERING
# ============================================================

qc_summary_list <- list()


for (sample_name in names(raw_objects)) {
  
  obj <-
    raw_objects[[sample_name]]
  
  
  qc_summary_list[[sample_name]] <-
    data.frame(
      
      Sample_ID =
        sample_name,
      
      Condition =
        unique(obj$condition),
      
      Replicate =
        unique(obj$replicate),
      
      Cells =
        ncol(obj),
      
      Median_nFeature_RNA =
        median(
          obj$nFeature_RNA
        ),
      
      Mean_nFeature_RNA =
        mean(
          obj$nFeature_RNA
        ),
      
      Median_nCount_RNA =
        median(
          obj$nCount_RNA
        ),
      
      Mean_nCount_RNA =
        mean(
          obj$nCount_RNA
        ),
      
      Median_percent_mt =
        median(
          obj$percent.mt
        ),
      
      Mean_percent_mt =
        mean(
          obj$percent.mt
        ),
      
      stringsAsFactors = FALSE
    )
}


qc_summary <- do.call(
  rbind,
  qc_summary_list
)


rownames(qc_summary) <- NULL


cat("\n============================================\n")
cat("RAW QC SUMMARY — NO FILTERING\n")
cat("============================================\n")

print(
  qc_summary,
  row.names = FALSE
)


write.csv(
  qc_summary,
  file.path(
    results_dir,
    "V01_Raw_QC_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. CHECK GENE UNIVERSE CONSISTENCY
# ============================================================

gene_lists <-
  lapply(
    raw_objects,
    rownames
  )


shared_genes <-
  Reduce(
    intersect,
    gene_lists
  )


gene_universe_summary <- data.frame(
  
  Sample_ID =
    c(
      names(raw_objects),
      "Shared_All_4"
    ),
  
  Genes =
    c(
      sapply(
        raw_objects,
        nrow
      ),
      length(shared_genes)
    ),
  
  stringsAsFactors = FALSE
)


cat("\n============================================\n")
cat("GENE UNIVERSE AUDIT\n")
cat("============================================\n")

print(
  gene_universe_summary,
  row.names = FALSE
)


write.csv(
  gene_universe_summary,
  file.path(
    results_dir,
    "V01_Gene_Universe_Audit.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. PREDEFINED VALIDATION TARGETS
#
# LOCKED BEFORE EXPRESSION ANALYSIS
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


priority_availability <- data.frame(
  
  Gene =
    priority_genes,
  
  Present_All_4 =
    priority_genes %in%
    shared_genes,
  
  stringsAsFactors = FALSE
)


cat("\n============================================\n")
cat("PREDEFINED CYP/AA TARGET AVAILABILITY\n")
cat("============================================\n")

print(
  priority_availability,
  row.names = FALSE
)


write.csv(
  priority_availability,
  file.path(
    results_dir,
    "V01_Priority_CYP_AA_Target_Availability.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. STRICT CYP FAMILY AUDIT
# ============================================================

strict_cyp <-
  grep(
    "^Cyp[0-9]",
    shared_genes,
    value = TRUE
  )


strict_cyp <- sort(
  strict_cyp
)


cat(
  "\nStrict CYP genes shared across all samples:",
  length(strict_cyp),
  "\n"
)


write.csv(
  data.frame(
    Gene = strict_cyp
  ),
  file.path(
    results_dir,
    "V01_Strict_CYP_Genes_Shared.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. AA NETWORK GENE AVAILABILITY
# ============================================================

aa_candidate_sets <- list(
  
  AA_Liberation = c(
    "Pla2g4a",
    "Pla2g4b",
    "Pla2g4c",
    "Pla2g6",
    "Pla2g2a",
    "Pla2g5"
  ),
  
  CYP_Epoxygenase = c(
    "Cyp2j5",
    "Cyp2j6",
    "Cyp2j8",
    "Cyp2j9",
    "Cyp2j11",
    "Cyp2j12",
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
    "Cyp2c70"
  ),
  
  Epoxide_Hydrolysis = c(
    "Ephx1",
    "Ephx2"
  ),
  
  CYP_Omega_Hydroxylase = c(
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
  
  COX_Arm = c(
    "Ptgs1",
    "Ptgs2",
    "Ptges",
    "Ptges2",
    "Ptges3",
    "Ptgis",
    "Tbxas1"
  ),
  
  LOX_Arm = c(
    "Alox5",
    "Alox5ap",
    "Alox12",
    "Alox12b",
    "Alox15",
    "Alox15b"
  ),
  
  Eicosanoid_Receptors = c(
    "Ptger1",
    "Ptger2",
    "Ptger3",
    "Ptger4",
    "Ptgdr",
    "Ptgdr2",
    "Tbxa2r",
    "Ltb4r1",
    "Ltb4r2",
    "Cysltr1",
    "Cysltr2"
  )
)


aa_availability_list <- list()


for (set_name in names(aa_candidate_sets)) {
  
  genes_now <-
    aa_candidate_sets[[set_name]]
  
  
  aa_availability_list[[set_name]] <-
    data.frame(
      
      AA_Arm =
        set_name,
      
      Gene =
        genes_now,
      
      Present_All_4 =
        genes_now %in%
        shared_genes,
      
      stringsAsFactors = FALSE
    )
}


aa_availability <- do.call(
  rbind,
  aa_availability_list
)


rownames(aa_availability) <- NULL


write.csv(
  aa_availability,
  file.path(
    results_dir,
    "V01_AA_Network_Gene_Availability.csv"
  ),
  row.names = FALSE
)


aa_arm_summary <-
  aggregate(
    
    Present_All_4 ~ AA_Arm,
    
    data = aa_availability,
    
    FUN = function(x) {
      c(
        Available = sum(x),
        Candidate = length(x)
      )
    }
  )


cat("\n============================================\n")
cat("AA NETWORK AVAILABILITY\n")
cat("============================================\n")


for (set_name in names(aa_candidate_sets)) {
  
  candidate_n <-
    length(
      aa_candidate_sets[[set_name]]
    )
  
  available_n <-
    sum(
      aa_candidate_sets[[set_name]]
      %in%
        shared_genes
    )
  
  cat(
    set_name,
    ":",
    available_n,
    "/",
    candidate_n,
    "\n"
  )
}


# ============================================================
# 16. RAW PRIORITY GENE DETECTION
#
# Descriptive feasibility only.
# NO biological interpretation yet.
# ============================================================

priority_detection_list <- list()


for (sample_name in names(raw_objects)) {
  
  obj <-
    raw_objects[[sample_name]]
  
  
  counts <-
    GetAssayData(
      obj,
      assay = "RNA",
      layer = "counts"
    )
  