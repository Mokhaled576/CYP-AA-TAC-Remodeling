# ============================================================
# VALIDATION #2 — GSE270896
# V2-02: RAPID SINGLET ATLAS + ANNOTATION EVIDENCE AUDIT
# ============================================================
#
# INPUT:
#   V2_01_GSE270896_QC_OBJECTS.rds
#
# PURPOSE:
#   1. Detect doublets independently in each biological library
#   2. Retain singlets only
#   3. Normalize each library independently
#   4. Build a merged UNINTEGRATED atlas
#   5. Select replicate-supported HVGs
#   6. PCA + UMAP
#   7. Unbiased clustering
#   8. Discover unbiased cluster markers
#   9. Generate canonical cardiac lineage evidence
#  10. Produce cluster/sample/condition audit
#
# IMPORTANT:
#   - CYP/AA genes are EXCLUDED from annotation evidence
#   - CYP/AA genes are NOT used to select clusters
#   - No Sham-vs-TAC differential expression
#   - No final cell-type labels are frozen here
#   - No cells are removed after singlet filtering
#   - Biological replicate remains the library/mouse
#   - No batch integration is performed
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
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scDblFinder)
  library(SingleCellExperiment)
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
    "V2_01_GSE270896_QC_OBJECTS.rds"
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
    "V2_02_RAPID_ATLAS_ANNOTATION_AUDIT"
  )

png_dir <-
  file.path(
    project_root,
    "FIGURES",
    "V2_02_RAPID_ATLAS_ANNOTATION_AUDIT",
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
# 4. LOAD V2-01 OBJECTS
# ============================================================

if (!file.exists(input_file)) {
  stop("V2-01 QC object not found.")
}

objects_qc <-
  readRDS(
    input_file
  )

expected_samples <- c(
  "Sham_Rep1",
  "Sham_Rep2",
  "TAC_Rep1",
  "TAC_Rep2"
)

if (!all(expected_samples %in% names(objects_qc))) {
  stop("Expected four V2-01 libraries not found.")
}

objects_qc <-
  objects_qc[
    expected_samples
  ]

input_cells <-
  sum(
    sapply(
      objects_qc,
      ncol
    )
  )

cat("\n============================================\n")
cat("V2-02 RAPID ATLAS + ANNOTATION AUDIT\n")
cat("============================================\n")

cat(
  "\nV2-01 QC nuclei entering:",
  input_cells,
  "\n"
)


# ============================================================
# 5. LOCK PREDEFINED CYP/AA EXCLUSION SET
# ============================================================
#
# These genes are explicitly excluded from annotation marker
# evidence so the validation remains independent.
#
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

strict_cyp_genes <-
  grep(
    "^Cyp[0-9]",
    rownames(objects_qc[[1]]),
    value = TRUE
  )

aa_context_genes <- c(
  "Ephx1",
  "Ephx2",
  "Pla2g4a",
  "Pla2g4b",
  "Pla2g4c",
  "Pla2g6",
  "Ptgs1",
  "Ptgs2",
  "Alox5",
  "Alox5ap",
  "Alox12",
  "Alox15",
  "Lta4h",
  "Ltb4r1",
  "Ltb4r2",
  "Ptger1",
  "Ptger2",
  "Ptger3",
  "Ptger4",
  "Ptgdr",
  "Ptgfr",
  "Tbxa2r",
  "Edn1",
  "Ednra",
  "Ednrb"
)

annotation_exclusion_genes <-
  unique(
    c(
      strict_cyp_genes,
      aa_context_genes,
      priority_genes
    )
  )

annotation_exclusion_genes <-
  annotation_exclusion_genes[
    annotation_exclusion_genes %in%
      rownames(objects_qc[[1]])
  ]

cat(
  "CYP/AA genes excluded from annotation evidence:",
  length(annotation_exclusion_genes),
  "\n"
)


# ============================================================
# 6. DOUBLETS — EACH LIBRARY INDEPENDENTLY
# ============================================================

singlet_objects <- list()
doublet_audit <- list()

counter <- 1

for (sid in expected_samples) {
  
  cat(
    "\nRunning scDblFinder:",
    sid,
    "\n"
  )
  
  x <- objects_qc[[sid]]
  
  counts <-
    GetAssayData(
      x,
      assay = "RNA",
      layer = "counts"
    )
  
  sce <-
    SingleCellExperiment(
      assays = list(
        counts = counts
      )
    )
  
  set.seed(20260921)
  
  sce <-
    scDblFinder(
      sce
    )
  
  dbl_class <-
    as.character(
      colData(sce)$scDblFinder.class
    )
  
  dbl_score <-
    as.numeric(
      colData(sce)$scDblFinder.score
    )
  
  names(dbl_class) <-
    colnames(sce)
  
  names(dbl_score) <-
    colnames(sce)
  
  x$scDblFinder_class <-
    dbl_class[
      colnames(x)
    ]
  
  x$scDblFinder_score <-
    dbl_score[
      colnames(x)
    ]
  
  n_start <- ncol(x)
  
  n_doublet <-
    sum(
      x$scDblFinder_class == "doublet"
    )
  
  n_singlet <-
    sum(
      x$scDblFinder_class == "singlet"
    )
  
  doublet_audit[[counter]] <-
    data.frame(
      Sample = sid,
      Start_Nuclei = n_start,
      Singlets = n_singlet,
      Doublets = n_doublet,
      Doublet_Percent =
        100 * n_doublet / n_start,
      stringsAsFactors = FALSE
    )
  
  singlet_cells <-
    colnames(x)[
      x$scDblFinder_class == "singlet"
    ]
  
  singlet_objects[[sid]] <-
    subset(
      x,
      cells = singlet_cells
    )
  
  rm(
    counts,
    sce,
    x
  )
  
  gc()
  
  counter <- counter + 1
}

doublet_audit <-
  bind_rows(
    doublet_audit
  )

write.csv(
  doublet_audit,
  file.path(
    results_dir,
    "V2_02_Doublet_Audit.csv"
  ),
  row.names = FALSE
)

singlet_n <-
  sum(
    sapply(
      singlet_objects,
      ncol
    )
  )

cat(
  "\nTotal singlets:",
  singlet_n,
  "\n"
)


# ============================================================
# 7. NORMALIZE EACH LIBRARY INDEPENDENTLY
# ============================================================

for (sid in expected_samples) {
  
  cat(
    "\nNormalizing:",
    sid,
    "\n"
  )
  
  x <- singlet_objects[[sid]]
  
  x <-
    NormalizeData(
      x,
      normalization.method =
        "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )
  
  x <-
    FindVariableFeatures(
      x,
      selection.method = "vst",
      nfeatures = 3000,
      verbose = FALSE
    )
  
  singlet_objects[[sid]] <- x
}


# ============================================================
# 8. REPLICATE-SUPPORTED HVGs
# ============================================================

hvg_list <-
  lapply(
    singlet_objects,
    VariableFeatures
  )

all_hvgs <-
  unique(
    unlist(
      hvg_list
    )
  )

hvg_frequency <-
  data.frame(
    Gene = all_hvgs,
    N_Libraries =
      sapply(
        all_hvgs,
        function(g) {
          sum(
            sapply(
              hvg_list,
              function(v) {
                g %in% v
              }
            )
          )
        }
      ),
    stringsAsFactors = FALSE
  )

hvg_frequency <-
  hvg_frequency %>%
  arrange(
    desc(N_Libraries),
    Gene
  )

write.csv(
  hvg_frequency,
  file.path(
    results_dir,
    "V2_02_HVG_Replicate_Frequency.csv"
  ),
  row.names = FALSE
)

consensus_hvgs <-
  hvg_frequency %>%
  filter(
    N_Libraries >= 2
  ) %>%
  pull(
    Gene
  )

# Exclude predefined CYP/AA genes from clustering features.
consensus_hvgs <-
  setdiff(
    consensus_hvgs,
    annotation_exclusion_genes
  )

# Keep at most 2500 replicate-supported HVGs.
if (length(consensus_hvgs) > 2500) {
  
  consensus_hvgs <-
    consensus_hvgs[
      1:2500
    ]
}

if (length(consensus_hvgs) < 1000) {
  
  warning(
    paste(
      "Fewer than 1000 replicate-supported HVGs.",
      "Proceeding with available features."
    )
  )
}

cat(
  "\nConsensus clustering HVGs:",
  length(consensus_hvgs),
  "\n"
)


# ============================================================
# 9. MERGE SINGLET OBJECTS
# ============================================================

atlas <-
  merge(
    x = singlet_objects[[1]],
    y = singlet_objects[2:4],
    merge.data = TRUE
  )

DefaultAssay(atlas) <- "RNA"

atlas <-
  JoinLayers(
    atlas,
    assay = "RNA"
  )

VariableFeatures(atlas) <-
  consensus_hvgs

cat(
  "\nMerged singlet atlas:",
  ncol(atlas),
  "nuclei x",
  nrow(atlas),
  "genes\n"
)


# ============================================================
# 10. SCALE CONSENSUS HVGs
# ============================================================

atlas <-
  ScaleData(
    atlas,
    features = consensus_hvgs,
    verbose = FALSE
  )


# ============================================================
# 11. PCA
# ============================================================

atlas <-
  RunPCA(
    atlas,
    features = consensus_hvgs,
    npcs = 50,
    verbose = FALSE
  )

pca_stdev <-
  Stdev(atlas, reduction = "pca")

variance_percent <-
  100 *
  pca_stdev^2 /
  sum(
    pca_stdev^2
  )

cumulative_variance <-
  cumsum(
    variance_percent
  )

pca_variance <-
  data.frame(
    PC = seq_along(
      variance_percent
    ),
    Variance_Percent =
      variance_percent,
    Cumulative_Percent =
      cumulative_variance
  )

write.csv(
  pca_variance,
  file.path(
    results_dir,
    "V2_02_PCA_Variance.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. SELECT PCS
# ============================================================

pc90 <-
  which(
    cumulative_variance >= 90
  )[1]

last_half_percent <-
  max(
    which(
      variance_percent >= 0.5
    )
  )

selected_pcs <-
  max(
    pc90,
    last_half_percent
  )

selected_pcs <-
  min(
    selected_pcs,
    40
  )

selected_pcs <-
  max(
    selected_pcs,
    20
  )

cat(
  "\nPC reaching 90%:",
  pc90,
  "\n"
)

cat(
  "Last PC >=0.5% variance:",
  last_half_percent,
  "\n"
)

cat(
  "Selected PCs:",
  selected_pcs,
  "\n"
)


# ============================================================
# 13. UMAP — UNINTEGRATED
# ============================================================

atlas <-
  RunUMAP(
    atlas,
    reduction = "pca",
    dims = 1:selected_pcs,
    reduction.name =
      "umap.unintegrated",
    reduction.key =
      "UMAPUNINT_",
    seed.use = 20260921,
    verbose = FALSE
  )


# ============================================================
# 14. GRAPH
# ============================================================

atlas <-
  FindNeighbors(
    atlas,
    reduction = "pca",
    dims = 1:selected_pcs,
    verbose = FALSE
  )


# ============================================================
# 15. RESOLUTION SWEEP
# ============================================================

resolutions <- c(
  0.2,
  0.4,
  0.6,
  0.8
)

for (res in resolutions) {
  
  atlas <-
    FindClusters(
      atlas,
      resolution = res,
      algorithm = 1,
      random.seed = 20260921,
      verbose = FALSE
    )
}


# ============================================================
# 16. RESOLUTION AUDIT
# ============================================================

resolution_audit_list <- list()
counter <- 1

for (res in resolutions) {
  
  col_name <-
    paste0(
      "RNA_snn_res.",
      res
    )
  
  cl <-
    as.character(
      atlas[[col_name]][, 1]
    )
  
  cluster_sizes <-
    table(cl)
  
  sample_cluster <-
    table(
      cl,
      atlas$sample_id
    )
  
  dominant_fraction <-
    apply(
      sample_cluster,
      1,
      function(z) {
        max(z) / sum(z)
      }
    )
  
  sample_restricted <-
    sum(
      dominant_fraction >= 0.90
    )
  
  resolution_audit_list[[counter]] <-
    data.frame(
      Resolution = res,
      Clusters =
        length(cluster_sizes),
      Minimum_Cluster_Size =
        min(cluster_sizes),
      Median_Cluster_Size =
        median(cluster_sizes),
      Maximum_Cluster_Size =
        max(cluster_sizes),
      Clusters_Under_50 =
        sum(cluster_sizes < 50),
      Clusters_Under_100 =
        sum(cluster_sizes < 100),
      Sample_Restricted_90pct =
        sample_restricted,
      Median_Dominant_Sample_Fraction =
        median(dominant_fraction),
      stringsAsFactors = FALSE
    )
  
  counter <- counter + 1
}

resolution_audit <-
  bind_rows(
    resolution_audit_list
  )

write.csv(
  resolution_audit,
  file.path(
    results_dir,
    "V2_02_Resolution_Audit.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. PRIMARY RESOLUTION
# ============================================================
#
# For rapid validation we use 0.6 as the primary framework.
# We still retain all sweep resolutions for review.
#
# This is NOT target-gene driven.
#
# ============================================================

primary_resolution <- 0.6

primary_column <-
  "RNA_snn_res.0.6"

atlas$V2_02_primary_cluster <-
  as.character(
    atlas[[primary_column]][, 1]
  )

Idents(atlas) <-
  atlas$V2_02_primary_cluster

primary_clusters <-
  sort(
    unique(
      atlas$V2_02_primary_cluster
    )
  )

cat(
  "\nPrimary resolution:",
  primary_resolution,
  "\n"
)

cat(
  "Primary clusters:",
  length(primary_clusters),
  "\n"
)


# ============================================================
# 18. PRIMARY CLUSTER COUNTS
# ============================================================

cluster_counts <-
  as.data.frame(
    table(
      Cluster =
        atlas$V2_02_primary_cluster
    )
  )

colnames(cluster_counts) <-
  c(
    "Cluster",
    "Cells"
  )

cluster_counts$Percent_of_Atlas <-
  100 *
  cluster_counts$Cells /
  ncol(atlas)

write.csv(
  cluster_counts,
  file.path(
    results_dir,
    "V2_02_Primary_Cluster_Counts.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. CLUSTER × SAMPLE
# ============================================================

cluster_sample <-
  as.data.frame.matrix(
    table(
      atlas$V2_02_primary_cluster,
      atlas$sample_id
    )
  )

cluster_sample$Cluster <-
  rownames(
    cluster_sample
  )

cluster_sample <-
  cluster_sample %>%
  select(
    Cluster,
    everything()
  )

rownames(cluster_sample) <- NULL

write.csv(
  cluster_sample,
  file.path(
    results_dir,
    "V2_02_Cluster_by_Sample.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 20. CLUSTER × CONDITION
# ============================================================

cluster_condition <-
  as.data.frame.matrix(
    table(
      atlas$V2_02_primary_cluster,
      atlas$condition
    )
  )

cluster_condition$Cluster <-
  rownames(
    cluster_condition
  )

cluster_condition <-
  cluster_condition %>%
  select(
    Cluster,
    everything()
  )

rownames(cluster_condition) <- NULL

write.csv(
  cluster_condition,
  file.path(
    results_dir,
    "V2_02_Cluster_by_Condition.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 21. SAMPLE DOMINANCE AUDIT
# ============================================================

sample_dominance_list <- list()
counter <- 1

for (cl in primary_clusters) {
  
  cells <-
    colnames(atlas)[
      atlas$V2_02_primary_cluster == cl
    ]
  
  sample_tab <-
    sort(
      table(
        atlas$sample_id[
          match(
            cells,
            colnames(atlas)
          )
        ]
      ),
      decreasing = TRUE
    )
  
  condition_tab <-
    sort(
      table(
        atlas$condition[
          match(
            cells,
            colnames(atlas)
          )
        ]
      ),
      decreasing = TRUE
    )
  
  sample_dominance_list[[counter]] <-
    data.frame(
      Cluster = cl,
      Cells = length(cells),
      
      Dominant_Sample =
        names(sample_tab)[1],
      
      Dominant_Sample_Percent =
        100 *
        as.numeric(sample_tab[1]) /
        length(cells),
      
      Dominant_Condition =
        names(condition_tab)[1],
      
      Dominant_Condition_Percent =
        100 *
        as.numeric(condition_tab[1]) /
        length(cells),
      
      stringsAsFactors = FALSE
    )