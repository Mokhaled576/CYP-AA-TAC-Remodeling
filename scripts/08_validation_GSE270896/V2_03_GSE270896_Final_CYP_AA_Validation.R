# ============================================================
# VALIDATION #2 — GSE270896
# V2-03: FINAL CYP/AA VALIDATION
# ============================================================
#
# INPUT:
#   V2_02_GSE270896_RAPID_CLUSTERED_ATLAS.rds
#
# PURPOSE:
#   1. Freeze V2-02 biological annotations
#   2. Validate annotations against unbiased cluster evidence
#   3. Quantify the 8 PREDEFINED CYP/AA targets
#   4. Whole-heart replicate-level pseudobulk
#   5. Replicate-level detection
#   6. Cell-type x replicate pseudobulk
#   7. Cell-type detection/localization
#   8. Replicate concordance
#   9. Late-stage validation classification
#  10. Publication-ready validation figures
#
# IMPORTANT:
#   Biological replicate = mouse/library
#   Nuclei are NOT biological replicates.
#
#   There are:
#       2 Sham biological replicates
#       2 TAC biological replicates
#
#   Therefore:
#       - no nucleus-level condition p-values
#       - no pseudoreplication
#       - no target-driven annotation
#       - no reclustering
#       - no re-QC
#       - no batch integration
#
#   Validation is based primarily on:
#       * pseudobulk CPM
#       * biological-replicate concordance
#       * detection/evaluability
#       * cell-type localization
#       * magnitude/direction of remodeling
#
#   This is a 10-week TAC dataset.
#   Interpret as LATE-STAGE INDEPENDENT SUPPORT,
#   not exact replication of the 2/4/6-week discovery trajectory.
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
  max.print = 100000
)

set.seed(20260921)


# ============================================================
# 2. PACKAGES
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})


# ============================================================
# 3. PATHS
# ============================================================

project_root <-
  "D:/Master/ScRNA seq/TAC/VALIDATION/GSE270896"

input_file <-
  file.path(
    project_root,
    "CLEAN DATA",
    "V2_02_GSE270896_RAPID_CLUSTERED_ATLAS.rds"
  )

clean_dir <-
  file.path(
    project_root,
    "CLEAN DATA"
  )

results_dir <-
  file.path(
    project_root,
    "RESULTS",
    "V2_03_FINAL_CYP_AA_VALIDATION"
  )

png_dir <-
  file.path(
    project_root,
    "FIGURES",
    "V2_03_FINAL_CYP_AA_VALIDATION",
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


# ============================================================
# 4. LOAD FROZEN V2-02 ATLAS
# ============================================================

if (!file.exists(input_file)) {
  stop("V2-02 atlas not found.")
}

atlas <-
  readRDS(
    input_file
  )

DefaultAssay(atlas) <- "RNA"

atlas <-
  JoinLayers(
    atlas,
    assay = "RNA"
  )

cat("\n============================================\n")
cat("V2-03 FINAL CYP/AA VALIDATION\n")
cat("============================================\n")

cat(
  "\nInput nuclei:",
  ncol(atlas),
  "\n"
)

cat(
  "Input genes:",
  nrow(atlas),
  "\n"
)


# ============================================================
# 5. VERIFY REQUIRED METADATA
# ============================================================

required_metadata <- c(
  "sample_id",
  "condition",
  "V2_02_primary_cluster"
)

missing_metadata <-
  setdiff(
    required_metadata,
    colnames(atlas@meta.data)
  )

if (length(missing_metadata) > 0) {
  
  stop(
    paste(
      "Missing metadata:",
      paste(
        missing_metadata,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 6. STANDARDIZE SAMPLE + CONDITION LABELS
# ============================================================

atlas$sample_id <-
  as.character(
    unlist(
      atlas@meta.data[["sample_id"]],
      use.names = FALSE
    )
  )

atlas$condition <-
  as.character(
    unlist(
      atlas@meta.data[["condition"]],
      use.names = FALSE
    )
  )

atlas$V2_02_primary_cluster <-
  as.character(
    unlist(
      atlas@meta.data[["V2_02_primary_cluster"]],
      use.names = FALSE
    )
  )

expected_samples <- c(
  "Sham_Rep1",
  "Sham_Rep2",
  "TAC_Rep1",
  "TAC_Rep2"
)

if (!all(expected_samples %in% unique(atlas$sample_id))) {
  stop("Expected four biological libraries not found.")
}

sample_condition_map <- c(
  Sham_Rep1 = "Sham",
  Sham_Rep2 = "Sham",
  TAC_Rep1 = "TAC",
  TAC_Rep2 = "TAC"
)

for (sid in expected_samples) {
  
  observed_condition <-
    unique(
      atlas$condition[
        atlas$sample_id == sid
      ]
    )
  
  if (
    length(observed_condition) != 1 ||
    observed_condition != sample_condition_map[[sid]]
  ) {
    stop(
      paste(
        "Unexpected condition assignment for",
        sid
      )
    )
  }
}


# ============================================================
# 7. LOCK THE 8 PREDEFINED TARGETS
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

missing_priority <-
  setdiff(
    priority_genes,
    rownames(atlas)
  )

if (length(missing_priority) > 0) {
  
  stop(
    paste(
      "Missing priority genes:",
      paste(
        missing_priority,
        collapse = ", "
      )
    )
  )
}

cat(
  "\nAll 8 predefined genes available: TRUE\n"
)


# ============================================================
# 8. FREEZE V2-02 CLUSTER ANNOTATIONS
# ============================================================
#
# Based on:
#   - unbiased markers
#   - canonical cardiac lineage evidence
#   - CYP/AA genes excluded from annotation evidence
#
# No CYP/AA target was used to define these labels.
#
# ============================================================

annotation_map <-
  data.frame(
    
    Cluster =
      as.character(
        0:23
      ),
    
    Broad_Cell_Type = c(
      "Fibroblasts",          # 0
      "Cardiomyocytes",       # 1
      "Endothelial",          # 2
      "Macrophages",          # 3
      "Endothelial",          # 4
      "Fibroblasts",          # 5
      "Pericytes",            # 6
      "Endothelial",          # 7
      "Endothelial",          # 8
      "Cardiomyocytes",       # 9
      "Cardiomyocytes",       # 10
      "Cardiomyocytes",       # 11
      "Endothelial",          # 12
      "Cardiomyocytes",       # 13
      "Endothelial",          # 14
      "Endothelial",          # 15
      "Mesothelial/Epicardial", # 16
      "Smooth muscle",        # 17
      "Cycling",              # 18
      "T cells",              # 19
      "Dendritic cells",      # 20
      "B cells",              # 21
      "Schwann/Neural",       # 22
      "Adipocytes"            # 23
    ),
    
    Refined_Cell_Type = c(
      "Fibroblasts",                 # 0
      "Cardiomyocytes",              # 1
      "Capillary endothelial",       # 2
      "Resident macrophages",        # 3
      "Endothelial",                 # 4
      "Fibroblasts",                 # 5
      "Pericytes",                   # 6
      "Endothelial",                 # 7
      "Remodeling endothelial",      # 8
      "Cardiomyocytes",              # 9
      "Cardiomyocytes",              # 10
      "Cardiomyocytes",              # 11
      "Capillary endothelial",       # 12
      "Cardiomyocytes",              # 13
      "Lymphatic endothelial",       # 14
      "Arterial endothelial",        # 15
      "Mesothelial cells",           # 16
      "Vascular smooth muscle",      # 17
      "Cycling cells",               # 18
      "Gamma-delta T cells",         # 19
      "Xcr1+ cDC1",                  # 20
      "B cells",                     # 21
      "Schwann cells",               # 22
      "Adipocytes"                   # 23
    ),
    
    Annotation_Confidence = c(
      "High",          # 0
      "High",          # 1
      "High",          # 2
      "Very high",     # 3
      "High",          # 4
      "High",          # 5
      "Very high",     # 6
      "High",          # 7
      "Moderate",      # 8
      "High",          # 9
      "Moderate",      # 10
      "High",          # 11
      "Very high",     # 12
      "High",          # 13
      "Very high",     # 14
      "High",          # 15
      "Very high",     # 16
      "Very high",     # 17
      "Very high",     # 18
      "High",          # 19
      "High",          # 20
      "Very high",     # 21
      "Very high",     # 22
      "Very high"      # 23
    ),
    
    stringsAsFactors = FALSE
  )


# ============================================================
# 9. CHECK CLUSTER MAP
# ============================================================
# ============================================================
# 9. CHECK CLUSTER MAP — CORRECTED
# ============================================================

observed_clusters <- sort(
  unique(
    as.integer(
      as.character(
        atlas$V2_02_primary_cluster
      )
    )
  )
)

expected_clusters <- 0:23

cat("\nObserved clusters:\n")
print(observed_clusters)

cat("\nExpected clusters:\n")
print(expected_clusters)

cat(
  "\nNumber of observed clusters:",
  length(observed_clusters),
  "\n"
)

if (
  length(observed_clusters) != 24 ||
  !setequal(
    observed_clusters,
    expected_clusters
  )
) {
  
  stop(
    paste(
      "Observed clusters do not match expected clusters 0:23.",
      "Do not continue."
    )
  )
}

cat(
  "\nCluster-map integrity check passed: TRUE\n"
)

cat(
  "All 24 expected clusters (0-23) are present.\n"
)
# ============================================================
# 10. APPLY FINAL ANNOTATIONS
# ============================================================

broad_lookup <-
  setNames(
    annotation_map$Broad_Cell_Type,
    annotation_map$Cluster
  )

refined_lookup <-
  setNames(
    annotation_map$Refined_Cell_Type,
    annotation_map$Cluster
  )

confidence_lookup <-
  setNames(
    annotation_map$Annotation_Confidence,
    annotation_map$Cluster
  )

atlas$V2_final_broad_cell_type <-
  unname(
    broad_lookup[
      atlas$V2_02_primary_cluster
    ]
  )

atlas$V2_final_refined_cell_type <-
  unname(
    refined_lookup[
      atlas$V2_02_primary_cluster
    ]
  )

atlas$V2_annotation_confidence <-
  unname(
    confidence_lookup[
      atlas$V2_02_primary_cluster
    ]
  )

atlas$V2_annotation_version <-
  "GSE270896_Validation2_v1"

atlas$V2_annotation_step <-
  "V2-03"

atlas$V2_annotation_date <-
  "2026-09-21"

if (
  any(
    is.na(
      atlas$V2_final_broad_cell_type
    )
  )
) {
  stop("NA broad annotations detected.")
}

if (
  any(
    is.na(
      atlas$V2_final_refined_cell_type
    )
  )
) {
  stop("NA refined annotations detected.")
}

write.csv(
  annotation_map,
  file.path(
    results_dir,
    "V2_03_Final_Cluster_Annotation_Map.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. ANNOTATION COUNTS
# ============================================================

broad_counts <-
  as.data.frame(
    table(
      atlas$V2_final_broad_cell_type
    ),
    stringsAsFactors = FALSE
  )

colnames(broad_counts) <-
  c(
    "Broad_Cell_Type",
    "Nuclei"
  )

broad_counts$Percent <-
  100 *
  broad_counts$Nuclei /
  ncol(atlas)

broad_counts <-
  broad_counts[
    order(
      -broad_counts$Nuclei
    ),
    ,
    drop = FALSE
  ]

write.csv(
  broad_counts,
  file.path(
    results_dir,
    "V2_03_Broad_Cell_Type_Counts.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. BROAD CELL TYPE × SAMPLE COUNTS
# ============================================================

broad_sample_counts <-
  as.data.frame(
    table(
      Broad_Cell_Type =
        atlas$V2_final_broad_cell_type,
      Sample =
        atlas$sample_id
    ),
    stringsAsFactors = FALSE
  )

colnames(broad_sample_counts) <-
  c(
    "Broad_Cell_Type",
    "Sample",
    "Nuclei"
  )

write.csv(
  broad_sample_counts,
  file.path(
    results_dir,
    "V2_03_Broad_CellType_by_Sample_Counts.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. CONDITION COMPOSITION
# ============================================================

composition_counts <-
  as.data.frame(
    table(
      Condition =
        atlas$condition,
      Broad_Cell_Type =
        atlas$V2_final_broad_cell_type
    ),
    stringsAsFactors = FALSE
  )

colnames(composition_counts) <-
  c(
    "Condition",
    "Broad_Cell_Type",
    "Nuclei"
  )

condition_totals <-
  tapply(
    composition_counts$Nuclei,
    composition_counts$Condition,
    sum
  )

composition_counts$Percent <-
  100 *
  composition_counts$Nuclei /
  condition_totals[
    composition_counts$Condition
  ]

write.csv(
  composition_counts,
  file.path(
    results_dir,
    "V2_03_Broad_CellType_Condition_Composition.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. EXTRACT RAW COUNTS
# ============================================================

counts <-
  GetAssayData(
    atlas,
    assay = "RNA",
    layer = "counts"
  )

priority_counts <-
  counts[
    priority_genes,
    ,
    drop = FALSE
  ]


# ============================================================
# 15. WHOLE-ATLAS DETECTION BY BIOLOGICAL REPLICATE
# ============================================================

whole_detection_list <- list()

counter <- 1

for (sid in expected_samples) {
  
  cells <-
    colnames(atlas)[
      atlas$sample_id == sid
    ]
  
  mat <-
    priority_counts[
      ,
      cells,
      drop = FALSE
    ]
  
  positive <-
    Matrix::rowSums(
      mat > 0
    )
  
  whole_detection_list[[counter]] <-
    data.frame(
      Sample = sid,
      Condition =
        sample_condition_map[[sid]],
      Gene =
        priority_genes,
      N_Nuclei =
        length(cells),
      Positive_Nuclei =
        as.numeric(
          positive[
            priority_genes
          ]
        ),
      Detection_Percent =
        100 *
        as.numeric(
          positive[
            priority_genes
          ]
        ) /
        length(cells),
      stringsAsFactors = FALSE
    )
  
  counter <- counter + 1
}

whole_detection <-
  do.call(
    rbind,
    whole_detection_list
  )

rownames(whole_detection) <- NULL

write.csv(
  whole_detection,
  file.path(
    results_dir,
    "V2_03_WholeAtlas_Detection_by_Replicate.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. GLOBAL EVALUABILITY
# ============================================================

global_positive <-
  Matrix::rowSums(
    priority_counts > 0
  )

library_positive <-
  sapply(
    expected_samples,
    function(sid) {
      
      cells <-
        colnames(atlas)[
          atlas$sample_id == sid
        ]
      
      Matrix::rowSums(
        priority_counts[
          ,
          cells,
          drop = FALSE
        ] > 0
      )
    }
  )

global_evaluability <-
  data.frame(
    Gene =
      priority_genes,
    
    Total_Positive_Nuclei =
      as.numeric(
        global_positive[
          priority_genes
        ]
      ),
    
    Libraries_with_Positive_Nuclei =
      rowSums(
        library_positive[
          priority_genes,
          ,
          drop = FALSE
        ] > 0
      ),
    
    stringsAsFactors = FALSE
  )

global_evaluability$Evaluability <-
  ifelse(
    global_evaluability$Total_Positive_Nuclei >= 50 &
      global_evaluability$Libraries_with_Positive_Nuclei == 4,
    "Robustly evaluable",
    ifelse(
      global_evaluability$Total_Positive_Nuclei >= 20 &
        global_evaluability$Libraries_with_Positive_Nuclei >= 3,
      "Evaluable with caution",
      "Not robustly evaluable"
    )
  )

write.csv(
  global_evaluability,
  file.path(
    results_dir,
    "V2_03_Global_Gene_Evaluability.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. WHOLE-HEART PSEUDOBULK
# ============================================================
#
# Aggregate RAW COUNTS by biological library.
#
# ============================================================

sample_total_counts <-
  numeric(
    length(expected_samples)
  )

names(sample_total_counts) <-
  expected_samples

whole_pb_counts <-
  matrix(
    0,
    nrow =
      length(priority_genes),
    ncol =
      length(expected_samples),
    dimnames = list(
      priority_genes,
      expected_samples
    )
  )

for (sid in expected_samples) {
  
  cells <-
    colnames(atlas)[
      atlas$sample_id == sid
    ]
  
  whole_pb_counts[
    ,
    sid
  ] <-
    Matrix::rowSums(
      priority_counts[
        ,
        cells,
        drop = FALSE
      ]
    )
  
  sample_total_counts[sid] <-
    sum(
      counts[
        ,
        cells,
        drop = FALSE
      ]
    )
}


# ============================================================
# 18. WHOLE-HEART CPM
# ============================================================

whole_pb_cpm <-
  sweep(
    whole_pb_counts,
    2,
    sample_total_counts,
    "/"
  ) *
  1e6

whole_pb_long_list <- list()
counter <- 1

for (sid in expected_samples) {
  
  whole_pb_long_list[[counter]] <-
    data.frame(
      Gene =
        priority_genes,
      Sample =
        sid,
      Condition =
        sample_condition_map[[sid]],
      Pseudobulk_Count =
        as.numeric(
          whole_pb_counts[
            priority_genes,
            sid
          ]
        ),
      Library_Total_Counts =
        sample_total_counts[
          sid
        ],
      CPM =
        as.numeric(
          whole_pb_cpm[
            priority_genes,
            sid
          ]
        ),
      stringsAsFactors = FALSE
    )
  
  counter <- counter + 1
}

whole_pb_long <-
  do.call(
    rbind,
    whole_pb_long_list
  )

rownames(whole_pb_long) <- NULL

write.csv(
  whole_pb_long,
  file.path(
    results_dir,
    "V2_03_WholeHeart_Pseudobulk_CPM.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. WHOLE-HEART CONDITION SUMMARY
# ============================================================

pseudocount_cpm <- 0.1

whole_summary_list <- list()

counter <- 1

for (g in priority_genes) {
  
  z <-
    whole_pb_long[
      whole_pb_long$Gene == g,
      ,
      drop = FALSE
    ]
  
  sham_cpm <-
    z$CPM[
      match(
        c(
          "Sham_Rep1",
          "Sham_Rep2"
        ),
        z$Sample
      )
    ]
  
  tac_cpm <-
    z$CPM[
      match(
        c(
          "TAC_Rep1",
          "TAC_Rep2"
        ),
        z$Sample
      )
    ]
  