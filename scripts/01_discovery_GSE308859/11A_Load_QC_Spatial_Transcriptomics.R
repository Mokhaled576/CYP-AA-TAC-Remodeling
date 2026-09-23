# ============================================================
# STEP 11A
# LOAD, RECONSTRUCT, AND QC THE FOUR VISIUM SAMPLES
# GSE308859
#
# PURPOSE:
#   1. Discover all spatial files
#   2. Reconstruct each Visium sample
#   3. Attach tissue coordinates
#   4. Verify image/scalefactor files
#   5. Calculate spot-level QC
#   6. Generate diagnostic QC plots
#   7. Save raw reconstructed spatial objects
#
# IMPORTANT:
#   - NO CYP interpretation
#   - NO spatial deconvolution
#   - NO biological inference
#   - NO spot filtering
#   - QC thresholds are NOT selected in this step
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
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(jsonlite)
  library(png)
})


# ============================================================
# 3. DEFINE PROJECT PATHS
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
  "11_SPATIAL",
  "11A_QC"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "11_SPATIAL",
  "11A_QC"
)

notes_dir <- file.path(
  project_dir,
  "NOTES"
)

dir.create(
  clean_dir,
  recursive = TRUE,
  showWarnings = FALSE
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

if (!dir.exists(raw_dir)) {
  stop(
    paste(
      "Raw spatial-data directory does not exist:",
      raw_dir
    )
  )
}

cat("\n========================================\n")
cat("PROJECT PATHS VERIFIED\n")
cat("========================================\n")

cat("\nRaw directory:\n")
cat(raw_dir, "\n")

cat("\nResults directory:\n")
cat(results_dir, "\n")

cat("\nFigures directory:\n")
cat(figures_dir, "\n")


# ============================================================
# 4. DEFINE THE FOUR SPATIAL SAMPLES
# ============================================================

sample_info <- data.frame(
  sample_id = c(
    "Sham",
    "TAC_2W",
    "TAC_4W",
    "TAC_6W"
  ),
  geo_id = c(
    "GSM9254699",
    "GSM9254696",
    "GSM9254697",
    "GSM9254698"
  ),
  geo_label = c(
    "Sham_ST",
    "TAC2w_ST",
    "TAC4w_ST",
    "TAC6w_ST"
  ),
  stringsAsFactors = FALSE
)

cat("\n========================================\n")
cat("SPATIAL SAMPLES\n")
cat("========================================\n")

print(sample_info)


# ============================================================
# 5. INVENTORY ALL FILES
# ============================================================

all_files <- list.files(
  raw_dir,
  full.names = TRUE,
  recursive = TRUE
)

if (length(all_files) == 0) {
  stop("No files were found in the raw-data directory.")
}

spatial_files <- all_files[
  grepl(
    "_ST",
    basename(all_files),
    ignore.case = TRUE
  )
]

spatial_inventory <- data.frame(
  File = spatial_files,
  Filename = basename(spatial_files),
  Size_MB = round(
    file.info(spatial_files)$size / 1024^2,
    4
  ),
  stringsAsFactors = FALSE
)

write.csv(
  spatial_inventory,
  file.path(
    results_dir,
    "11A_Spatial_File_Inventory.csv"
  ),
  row.names = FALSE
)

cat("\n========================================\n")
cat("SPATIAL FILE INVENTORY\n")
cat("========================================\n")

print(spatial_inventory$Filename)


# ============================================================
# 6. HELPER FUNCTION TO FIND SAMPLE FILE
# ============================================================

find_sample_file <- function(
    geo_id,
    required_pattern
) {
  
  candidates <- all_files[
    grepl(
      geo_id,
      basename(all_files),
      fixed = TRUE
    ) &
      grepl(
        required_pattern,
        basename(all_files),
        ignore.case = TRUE
      )
  ]
  
  if (length(candidates) == 0) {
    return(NA_character_)
  }
  
  candidates[1]
}


# ============================================================
# 7. IDENTIFY MATRIX / BARCODE / FEATURE FILES
# ============================================================

sample_files <- list()

for (i in seq_len(nrow(sample_info))) {
  
  sid <- sample_info$sample_id[i]
  gid <- sample_info$geo_id[i]
  
  cat("\n----------------------------------------\n")
  cat("Discovering expression files:", sid, "\n")
  cat("----------------------------------------\n")
  
  matrix_file <- find_sample_file(
    gid,
    "matrix.mtx"
  )
  
  barcode_file <- find_sample_file(
    gid,
    "barcodes.tsv"
  )
  
  feature_file <- find_sample_file(
    gid,
    "features.tsv"
  )
  
  if (is.na(matrix_file)) {
    stop(
      paste(
        "Matrix file not found for",
        sid
      )
    )
  }
  
  if (is.na(barcode_file)) {
    stop(
      paste(
        "Barcode file not found for",
        sid
      )
    )
  }
  
  if (is.na(feature_file)) {
    stop(
      paste(
        "Feature file not found for",
        sid
      )
    )
  }
  
  sample_files[[sid]] <- list(
    matrix = matrix_file,
    barcodes = barcode_file,
    features = feature_file
  )
  
  cat("Matrix:\n")
  cat(matrix_file, "\n")
  
  cat("Barcodes:\n")
  cat(barcode_file, "\n")
  
  cat("Features:\n")
  cat(feature_file, "\n")
}


# ============================================================
# 8. IDENTIFY POSITION / SCALEFACTOR / IMAGE FILES
# ============================================================

spatial_support_files <- list()

for (i in seq_len(nrow(sample_info))) {
  
  sid <- sample_info$sample_id[i]
  gid <- sample_info$geo_id[i]
  
  sample_candidates <- all_files[
    grepl(
      gid,
      basename(all_files),
      fixed = TRUE
    )
  ]
  
  position_candidates <- sample_candidates[
    grepl(
      "position",
      basename(sample_candidates),
      ignore.case = TRUE
    )
  ]
  
  scalefactor_candidates <- sample_candidates[
    grepl(
      "scalefactor",
      basename(sample_candidates),
      ignore.case = TRUE
    )
  ]
  
  image_candidates <- sample_candidates[
    grepl(
      "\\.(png|jpg|jpeg)$",
      basename(sample_candidates),
      ignore.case = TRUE
    )
  ]
  
  spatial_support_files[[sid]] <- list(
    positions = position_candidates,
    scalefactors = scalefactor_candidates,
    images = image_candidates
  )
  
  cat("\n----------------------------------------\n")
  cat("Spatial support files:", sid, "\n")
  cat("----------------------------------------\n")
  
  cat("\nPosition candidates:\n")
  print(position_candidates)
  
  cat("\nScalefactor candidates:\n")
  print(scalefactor_candidates)
  
  cat("\nImage candidates:\n")
  print(image_candidates)
}


# ============================================================
# 9. SAVE DISCOVERED FILE MAP
# ============================================================

file_map_list <- list()

for (sid in sample_info$sample_id) {
  
  file_map_list[[sid]] <- data.frame(
    Sample = sid,
    Matrix = sample_files[[sid]]$matrix,
    Barcodes = sample_files[[sid]]$barcodes,
    Features = sample_files[[sid]]$features,
    Position_Files = paste(
      spatial_support_files[[sid]]$positions,
      collapse = ";"
    ),
    Scalefactor_Files = paste(
      spatial_support_files[[sid]]$scalefactors,
      collapse = ";"
    ),
    Image_Files = paste(
      spatial_support_files[[sid]]$images,
      collapse = ";"
    ),
    stringsAsFactors = FALSE
  )
}

file_map <- bind_rows(
  file_map_list
)

write.csv(
  file_map,
  file.path(
    results_dir,
    "11A_Spatial_File_Map.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. FUNCTION TO READ GEO SPATIAL MATRIX
# ============================================================

read_geo_spatial_matrix <- function(
    matrix_file,
    barcode_file,
    feature_file,
    sample_id
) {
  
  cat("\n----------------------------------------\n")
  cat("Reading expression matrix:", sample_id, "\n")
  cat("----------------------------------------\n")
  
  mat <- Matrix::readMM(
    matrix_file
  )
  
  barcodes <- read.delim(
    barcode_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  features <- read.delim(
    feature_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  cat(
    "Matrix dimensions:",
    nrow(mat),
    "x",
    ncol(mat),
    "\n"
  )
  
  cat(
    "Feature rows:",
    nrow(features),
    "\n"
  )
  
  cat(
    "Barcode rows:",
    nrow(barcodes),
    "\n"
  )
  
  if (nrow(mat) != nrow(features)) {
    stop(
      paste(
        sample_id,
        "- feature count does not match matrix rows."
      )
    )
  }
  
  if (ncol(mat) != nrow(barcodes)) {
    stop(
      paste(
        sample_id,
        "- barcode count does not match matrix columns."
      )
    )
  }
  
  if (ncol(features) >= 2) {
    
    gene_names <- as.character(
      features[, 2]
    )
    
  } else {
    
    gene_names <- as.character(
      features[, 1]
    )
  }
  
  missing_gene_names <- is.na(gene_names) |
    gene_names == ""
  
  if (any(missing_gene_names)) {
    
    gene_names[missing_gene_names] <- as.character(
      features[
        missing_gene_names,
        1
      ]
    )
  }
  
  gene_names <- make.unique(
    gene_names
  )
  
  spot_names <- paste0(
    sample_id,
    "_",
    as.character(barcodes[, 1])
  )
  
  rownames(mat) <- gene_names
  colnames(mat) <- spot_names
  
  return(mat)
}


# ============================================================
# 11. CREATE RAW SEURAT OBJECTS
# ============================================================

spatial_objects <- list()

for (sid in sample_info$sample_id) {
  
  mat <- read_geo_spatial_matrix(
    matrix_file = sample_files[[sid]]$matrix,
    barcode_file = sample_files[[sid]]$barcodes,
    feature_file = sample_files[[sid]]$features,
    sample_id = sid
  )
  
  sobj <- CreateSeuratObject(
    counts = mat,
    assay = "Spatial",
    project = paste0(
      "GSE308859_",
      sid,
      "_Spatial"
    ),
    min.cells = 0,
    min.features = 0
  )
  
  sobj$sample_id <- sid
  sobj$condition <- sid
  sobj$spatial_analysis_step <- "11A_RAW"
  
  spatial_objects[[sid]] <- sobj
  
  cat(
    "\nCreated:",
    sid,
    "| spots:",
    ncol(sobj),
    "| genes:",
    nrow(sobj),
    "\n"
  )
}

rm(mat)
gc()


# ============================================================
# 12. FUNCTION TO READ POSITION FILE
# ============================================================

read_position_file <- function(
    position_files,
    sample_id
) {
  
  if (length(position_files) == 0) {
    
    warning(
      paste(
        "No position file found for",
        sample_id
      )
    )
    
    return(NULL)
  }
  
  position_file <- position_files[1]
  
  cat(
    "\nReading position file for",
    sample_id,
    ":\n"
  )
  
  cat(
    position_file,
    "\n"
  )
  
  first_line <- readLines(
    position_file,
    n = 1,
    warn = FALSE
  )
  
  has_header <- grepl(
    "barcode|in_tissue|array_row|array_col",
    first_line,
    ignore.case = TRUE
  )
  
  pos <- read.csv(
    position_file,
    header = has_header,
    stringsAsFactors = FALSE
  )
  
  if (!has_header) {
    
    if (ncol(pos) < 6) {
      stop(
        paste(
          "Position file has fewer than 6 columns:",
          sample_id
        )
      )
    }
    
    colnames(pos)[1:6] <- c(
      "barcode",
      "in_tissue",
      "array_row",
      "array_col",
      "pxl_row_in_fullres",
      "pxl_col_in_fullres"
    )
    
  } else {
    
    names(pos) <- tolower(
      names(pos)
    )
    
    if (!"barcode" %in% names(pos)) {
      
      barcode_candidate <- grep(
        "barcode",
        names(pos),
        value = TRUE
      )
      
      if (length(barcode_candidate) >= 1) {
        names(pos)[names(pos) == barcode_candidate[1]] <- "barcode"
      }
    }
    
    if (!"in_tissue" %in% names(pos)) {
      
      tissue_candidate <- grep(
        "in.*tissue",
        names(pos),
        value = TRUE
      )
      
      if (length(tissue_candidate) >= 1) {
        names(pos)[names(pos) == tissue_candidate[1]] <- "in_tissue"
      }
    }
    
    if (!"array_row" %in% names(pos)) {
      
      row_candidate <- grep(
        "array.*row",
        names(pos),
        value = TRUE
      )
      
      if (length(row_candidate) >= 1) {
        names(pos)[names(pos) == row_candidate[1]] <- "array_row"
      }
    }
    
    if (!"array_col" %in% names(pos)) {
      
      col_candidate <- grep(
        "array.*col",
        names(pos),
        value = TRUE
      )
      
      if (length(col_candidate) >= 1) {
        names(pos)[names(pos) == col_candidate[1]] <- "array_col"
      }
    }
    
    if (!"pxl_row_in_fullres" %in% names(pos)) {
      
      pixel_row_candidate <- grep(
        "pxl.*row|pixel.*row",
        names(pos),
        value = TRUE
      )
      
      if (length(pixel_row_candidate) >= 1) {
        names(pos)[names(pos) == pixel_row_candidate[1]] <- "pxl_row_in_fullres"
      }
    }
    
    if (!"pxl_col_in_fullres" %in% names(pos)) {
      
      pixel_col_candidate <- grep(
        "pxl.*col|pixel.*col",
        names(pos),
        value = TRUE
      )
      
      if (length(pixel_col_candidate) >= 1) {
        names(pos)[names(pos) == pixel_col_candidate[1]] <- "pxl_col_in_fullres"
      }
    }
  }
  
  return(pos)
}


# ============================================================
# 13. ATTACH POSITION METADATA
# ============================================================

position_tables <- list()

for (sid in sample_info$sample_id) {
  
  pos <- read_position_file(
    position_files = spatial_support_files[[sid]]$positions,
    sample_id = sid
  )
  
  if (is.null(pos)) {
    next
  }
  
  required_position_columns <- c(
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres"
  )
  
  missing_position_columns <- setdiff(
    required_position_columns,
    colnames(pos)
  )
  
  if (length(missing_position_columns) > 0) {
    
    stop(
      paste(
        sid,
        "- missing position columns:",
        paste(
          missing_position_columns,
          collapse = ", "
        )
      )
    )
  }
  
  pos$Cell <- paste0(
    sid,
    "_",
    pos$barcode
  )
  
  rownames(pos) <- pos$Cell
  
  object_cells <- colnames(
    spatial_objects[[sid]]
  )
  
  matched_cells <- intersect(
    object_cells,
    rownames(pos)
  )
  
  cat(
    "\n",
    sid,
    " - position matches: ",
    length(matched_cells),
    "/",
    length(object_cells),
    "\n",
    sep = ""
  )
  
  if (length(matched_cells) == 0) {
    
    stop(
      paste(
        "No spatial barcodes matched for",
        sid
      )
    )
  }
  
  missing_cells <- setdiff(
    object_cells,
    rownames(pos)
  )
  
  if (length(missing_cells) > 0) {
    
    warning(
      paste(
        sid,
        "has",
        length(missing_cells),
        "expression barcodes without position metadata."
      )
    )
  }
  
  metadata_to_add <- pos[
    object_cells,
    required_position_columns,
    drop = FALSE
  ]
  
  spatial_objects[[sid]] <- AddMetaData(
    object = spatial_objects[[sid]],
    metadata = metadata_to_add
  )
  
  position_tables[[sid]] <- pos
  
  write.csv(
    pos,
    file.path(
      results_dir,
      paste0(
        "11A_",
        sid,
        "_Spatial_Positions.csv"
      )
    ),
    row.names = FALSE
  )
}


# ============================================================
# 14. CALCULATE SPOT QC
# ============================================================

for (sid in sample_info$sample_id) {
  
  sobj <- spatial_objects[[sid]]
  
  assay_names <- names(sobj@assays)
  
  if (!"Spatial" %in% assay_names) {
    
    stop(
      paste(
        "Spatial assay is not present for",
        sid
      )
    )
  }
  
  DefaultAssay(sobj) <- "Spatial"
  
  sobj[["percent.mt"]] <- PercentageFeatureSet(
    object = sobj,
    pattern = "^mt-"
  )
  
  sobj[["percent.ribo"]] <- PercentageFeatureSet(
    object = sobj,
    pattern = "^Rp[sl]"
  )
  
  sobj[["percent.hb"]] <- PercentageFeatureSet(
    object = sobj,
    pattern = "^Hb[ab]"
  )
  
  spatial_objects[[sid]] <- sobj
  
  cat(
    "\nQC metrics calculated for:",
    sid,
    "\n"
  )
}

cat("\n========================================\n")
cat("SPOT QC METRICS CALCULATED\n")
cat("========================================\n")

for (sid in sample_info$sample_id) {
  
  cat(
    "\n",
    sid,
    "\n",
    sep = ""
  )
  
  print(
    summary(
      spatial_objects[[sid]]@meta.data[
        ,
        c(
          "nFeature_Spatial",
          "nCount_Spatial",
          "percent.mt",
          "percent.ribo",
          "percent.hb"
        ),
        drop = FALSE
      ]
    )
  )
}


# ============================================================
# 15. CREATE SPOT-LEVEL QC TABLE
# ============================================================

qc_tables <- list()

for (sid in sample_info$sample_id) {
  
  sobj <- spatial_objects[[sid]]
  
  md <- sobj@meta.data
  
  md$Spot <- rownames(md)
  md$Sample <- sid
  
  qc_tables[[sid]] <- md
}

all_qc <- bind_rows(
  qc_tables
)

write.csv(
  all_qc,
  file.path(
    results_dir,
    "11A_All_Spatial_Spot_QC.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. QC SUMMARY
# ============================================================

qc_summary <- all_qc %>%
  group_by(Sample) %>%
  summarise(
    Spots = n(),
    
    In_Tissue_Spots = sum(
      in_tissue == 1,
      na.rm = TRUE
    ),
    
    Percent_In_Tissue = 100 * mean(
      in_tissue == 1,