# ============================================================
# GSE308859 TAC scRNA-seq PROJECT
# STEP 1: CHECK ALL CYP GENES IN THE DATASET
# ============================================================

# Clear environment
rm(list = ls())
gc()

# ------------------------------------------------------------
# 1. SET DATA DIRECTORY
# ------------------------------------------------------------

data_dir <- "D:/Master/ScRNA seq/TAC/RAW DATA/GSE308859/GSE308859_RAW"

# Confirm that R can see the folder
dir.exists(data_dir)

# Show files
files <- list.files(data_dir)
print(files)


# ------------------------------------------------------------
# 2. IDENTIFY THE FOUR scRNA-seq FEATURE FILES
# ------------------------------------------------------------

feature_files <- list.files(
  data_dir,
  pattern = "scRNA_features\\.tsv(\\.gz)?$",
  full.names = TRUE
)

print(feature_files)

cat("\nNumber of scRNA feature files found:",
    length(feature_files), "\n")


# ------------------------------------------------------------
# 3. READ ALL FEATURE FILES
# ------------------------------------------------------------

feature_lists <- lapply(feature_files, function(f) {
  
  x <- read.delim(
    f,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  cat("\n---------------------------------\n")
  cat("FILE:", basename(f), "\n")
  cat("Dimensions:", nrow(x), "x", ncol(x), "\n")
  
  # Display first few rows so we can confirm format
  print(head(x))
  
  return(x)
})


# ------------------------------------------------------------
# 4. COLLECT ALL GENE SYMBOLS
# ------------------------------------------------------------

# 10x features files normally have:
# V1 = Ensembl ID
# V2 = gene symbol
# V3 = feature type

all_genes <- unique(
  unlist(
    lapply(feature_lists, function(x) x[[2]])
  )
)

cat("\nTotal unique genes/features:",
    length(all_genes), "\n")


# ------------------------------------------------------------
# 5. FIND EVERY CYP GENE
# ------------------------------------------------------------

all_cyps <- sort(
  unique(
    grep(
      "^Cyp",
      all_genes,
      value = TRUE,
      ignore.case = TRUE
    )
  )
)

cat("\n========================================\n")
cat("ALL CYP GENES PRESENT IN GSE308859\n")
cat("========================================\n\n")

print(all_cyps)

cat("\nTOTAL NUMBER OF CYP GENES:",
    length(all_cyps), "\n")


# ------------------------------------------------------------
# 6. ALSO CHECK EPHX2
# ------------------------------------------------------------

cat("\n========================================\n")
cat("EPHX2 CHECK\n")
cat("========================================\n")

print(
  grep(
    "^Ephx2$",
    all_genes,
    value = TRUE,
    ignore.case = TRUE
  )
)


# ------------------------------------------------------------
# 7. CYP COUNTS BY FAMILY
# ------------------------------------------------------------

cyp_family <- sub(
  "^([Cc][Yy][Pp][0-9]+).*",
  "\\1",
  all_cyps
)

family_table <- sort(
  table(cyp_family),
  decreasing = TRUE
)

cat("\n========================================\n")
cat("CYP FAMILY SUMMARY\n")
cat("========================================\n")

print(family_table)


# ------------------------------------------------------------
# END STEP 1
# ------------------------------------------------------------