# ============================================================
# GSE308859 TAC scRNA-seq PROJECT
# STEP 2: CYP EXPRESSION FEASIBILITY
# Sham vs TAC 2W vs TAC 4W vs TAC 6W
# ============================================================

rm(list = ls())
gc()

# ------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------

if (!requireNamespace("Matrix", quietly = TRUE)) {
  install.packages("Matrix")
}

library(Matrix)


# ------------------------------------------------------------
# 2. DATA DIRECTORY
# ------------------------------------------------------------

data_dir <- "D:/Master/ScRNA seq/TAC/RAW DATA/GSE308859/GSE308859_RAW"

stopifnot(dir.exists(data_dir))


# ------------------------------------------------------------
# 3. DEFINE SAMPLES
# ------------------------------------------------------------

samples <- data.frame(
  condition = c("Sham", "TAC_2W", "TAC_4W", "TAC_6W"),
  prefix = c(
    "GSM9254695_Sham_scRNA",
    "GSM9254692_TAC2w_scRNA",
    "GSM9254693_TAC4w_scRNA",
    "GSM9254694_TAC6w_scRNA"
  ),
  stringsAsFactors = FALSE
)

print(samples)


# ------------------------------------------------------------
# 4. FUNCTION TO READ ONE 10X SAMPLE
# ------------------------------------------------------------

read_sample <- function(prefix, condition) {
  
  cat("\n========================================\n")
  cat("READING:", condition, "\n")
  cat("========================================\n")
  
  matrix_file <- file.path(
    data_dir,
    paste0(prefix, "_matrix.mtx.gz")
  )
  
  feature_file <- file.path(
    data_dir,
    paste0(prefix, "_features.tsv.gz")
  )
  
  barcode_file <- file.path(
    data_dir,
    paste0(prefix, "_barcodes.tsv.gz")
  )
  
  # Check files
  stopifnot(
    file.exists(matrix_file),
    file.exists(feature_file),
    file.exists(barcode_file)
  )
  
  # Read matrix
  mat <- readMM(matrix_file)
  
  # Read features
  features <- read.delim(
    feature_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Read barcodes
  barcodes <- read.delim(
    barcode_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Add names
  rownames(mat) <- make.unique(features[[2]])
  colnames(mat) <- barcodes[[1]]
  
  cat("Genes:", nrow(mat), "\n")
  cat("Cells:", ncol(mat), "\n")
  cat("Non-zero counts:", length(mat@x), "\n")
  
  return(mat)
}


# ------------------------------------------------------------
# 5. LOAD ALL FOUR MATRICES
# ------------------------------------------------------------

matrices <- list()

for (i in seq_len(nrow(samples))) {
  
  matrices[[samples$condition[i]]] <- read_sample(
    prefix = samples$prefix[i],
    condition = samples$condition[i]
  )
}


# ------------------------------------------------------------
# 6. IDENTIFY ALL REAL CYP GENES
# ------------------------------------------------------------

all_genes <- unique(
  unlist(
    lapply(matrices, rownames)
  )
)

# Require CYP followed by a number.
# This excludes Cypt1, Cypt2, etc.
all_cyps <- sort(
  grep(
    "^Cyp[0-9]",
    all_genes,
    value = TRUE,
    ignore.case = TRUE
  )
)

cat("\n========================================\n")
cat("CYP GENES TO TEST\n")
cat("========================================\n")

print(all_cyps)

cat("\nNumber of CYP genes tested:",
    length(all_cyps), "\n")


# ------------------------------------------------------------
# 7. FUNCTION TO CALCULATE CYP EXPRESSION
# ------------------------------------------------------------

calculate_cyp_expression <- function(mat, condition) {
  
  genes_present <- intersect(all_cyps, rownames(mat))
  
  submat <- mat[genes_present, , drop = FALSE]
  
  # Number of cells with at least 1 UMI
  expressing_cells <- Matrix::rowSums(submat > 0)
  
  # Percent of cells expressing gene
  percent_cells <- expressing_cells / ncol(mat) * 100
  
  # Mean raw UMI count across ALL cells
  mean_count <- Matrix::rowMeans(submat)
  
  # Total UMI count
  total_umi <- Matrix::rowSums(submat)
  
  result <- data.frame(
    Gene = genes_present,
    Condition = condition,
    Total_Cells = ncol(mat),
    Expressing_Cells = expressing_cells,
    Percent_Expressing = percent_cells,
    Mean_Raw_Count = mean_count,
    Total_UMI = total_umi,
    stringsAsFactors = FALSE
  )
  
  rownames(result) <- NULL
  
  return(result)
}


# ------------------------------------------------------------
# 8. CALCULATE FOR ALL CONDITIONS
# ------------------------------------------------------------

results_list <- lapply(
  names(matrices),
  function(condition) {
    
    calculate_cyp_expression(
      matrices[[condition]],
      condition
    )
  }
)

cyp_results <- do.call(
  rbind,
  results_list
)


# ------------------------------------------------------------
# 9. ADD EPHX2
# ------------------------------------------------------------

calculate_gene_expression <- function(mat, gene, condition) {
  
  if (!(gene %in% rownames(mat))) {
    return(NULL)
  }
  
  x <- mat[gene, ]
  
  data.frame(
    Gene = gene,
    Condition = condition,
    Total_Cells = ncol(mat),
    Expressing_Cells = sum(x > 0),
    Percent_Expressing = sum(x > 0) / ncol(mat) * 100,
    Mean_Raw_Count = mean(x),
    Total_UMI = sum(x),
    stringsAsFactors = FALSE
  )
}

ephx2_results <- do.call(
  rbind,
  lapply(
    names(matrices),
    function(condition) {
      calculate_gene_expression(
        matrices[[condition]],
        "Ephx2",
        condition
      )
    }
  )
)


# ------------------------------------------------------------
# 10. ONLY CYPs WITH ACTUAL EXPRESSION
# ------------------------------------------------------------

expressed_cyps <- cyp_results[
  cyp_results$Total_UMI > 0,
]

expressed_cyps <- expressed_cyps[
  order(
    expressed_cyps$Condition,
    -expressed_cyps$Percent_Expressing
  ),
]


# ------------------------------------------------------------
# 11. PRINT SUMMARY
# ------------------------------------------------------------

cat("\n\n========================================\n")
cat("NUMBER OF EXPRESSED CYPs BY CONDITION\n")
cat("========================================\n")

print(
  aggregate(
    Gene ~ Condition,
    data = expressed_cyps,
    FUN = function(x) length(unique(x))
  )
)


cat("\n\n========================================\n")
cat("TOP EXPRESSED CYPs IN EACH CONDITION\n")
cat("========================================\n")

for (condition in names(matrices)) {
  
  cat("\n\n########", condition, "########\n")
  
  temp <- expressed_cyps[
    expressed_cyps$Condition == condition,
  ]
  
  temp <- temp[
    order(-temp$Percent_Expressing),
  ]
  
  print(
    head(
      temp[
        ,
        c(
          "Gene",
          "Expressing_Cells",
          "Percent_Expressing",
          "Mean_Raw_Count",
          "Total_UMI"
        )
      ],
      25
    ),
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# 12. OUR AA-RELATED CYP FAMILIES
# ------------------------------------------------------------

aa_cyps <- cyp_results[
  grepl(
    "^Cyp2c|^Cyp2j|^Cyp4a|^Cyp4f",
    cyp_results$Gene,
    ignore.case = TRUE
  ),
]

aa_cyps_expressed <- aa_cyps[
  aa_cyps$Total_UMI > 0,
]

cat("\n\n========================================\n")
cat("ARACHIDONIC-ACID-RELATED CYP FAMILIES\n")
cat("========================================\n")

print(
  aa_cyps_expressed[
    order(
      aa_cyps_expressed$Gene,
      aa_cyps_expressed$Condition
    ),
  ],
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. EPHX2 RESULTS
# ------------------------------------------------------------

cat("\n\n========================================\n")
cat("EPHX2 EXPRESSION\n")
cat("========================================\n")

print(
  ephx2_results,
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. SAVE RESULTS
# ------------------------------------------------------------

results_dir <- "D:/Master/ScRNA seq/TAC/RESULTS"

if (!dir.exists(results_dir)) {
  dir.create(
    results_dir,
    recursive = TRUE
  )
}

write.csv(
  cyp_results,
  file.path(
    results_dir,
    "02_All_CYP_Expression_Raw.csv"
  ),
  row.names = FALSE
)

write.csv(
  expressed_cyps,
  file.path(
    results_dir,
    "02_Expressed_CYPs_Only.csv"
  ),
  row.names = FALSE
)

write.csv(
  aa_cyps_expressed,
  file.path(
    results_dir,
    "02_AA_Related_CYPs_Expressed.csv"
  ),
  row.names = FALSE
)

write.csv(
  ephx2_results,
  file.path(
    results_dir,
    "02_EPHX2_Expression.csv"
  ),
  row.names = FALSE
)


cat("\n\n========================================\n")
cat("STEP 2 COMPLETE\n")
cat("========================================\n")

cat("\nResults saved to:\n")
cat(results_dir, "\n")