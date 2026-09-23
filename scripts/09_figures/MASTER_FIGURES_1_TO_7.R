# ============================================================
# MASTER MANUSCRIPT FIGURE PIPELINE
# FIGURES 1-7
#
# Project:
# Cellular and spatial remodeling of CYP-mediated
# arachidonic-acid metabolism during pressure-overload
# cardiac remodeling
#
# IMPORTANT PRINCIPLES
# ------------------------------------------------------------
# 1. Frozen analyses are preserved.
# 2. No clustering/integration/UMAP recomputation.
# 3. Discovery cells/spots are NOT biological replicates.
# 4. No discovery condition-level inferential statistics.
# 5. Spatial cell-program scores are expression proxies.
# 6. AA pathway scores are transcriptomic programs.
# 7. No claims about metabolite abundance or metabolic flux.
# 8. No causal pathway claims.
# 9. 4W = prominent intermediate remodeling state,
#    NOT "peak disease".
# 10. Validation uses biological-library pseudobulk outputs.
#
# OUTPUT:
# D:/Master/ScRNA seq/TAC/FIGURES/MANUSCRIPT_FINAL/
# ============================================================


# ============================================================
# 0. CLEAN WORKSPACE
# ============================================================

rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  scipen = 999
)


# ============================================================
# 1. PACKAGES
# ============================================================

required_packages <- c(
  "Seurat",
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "stringr",
  "forcats",
  "patchwork",
  "scales",
  "grid",
  "png"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(forcats)
library(patchwork)
library(scales)
library(grid)
library(png)


# ============================================================
# 2. PROJECT PATHS
# ============================================================

ROOT <- "D:/Master/ScRNA seq/TAC"

CLEAN <- file.path(ROOT, "CLEAN DATA")
RESULTS <- file.path(ROOT, "RESULTS")
FIG_ROOT <- file.path(ROOT, "FIGURES", "MANUSCRIPT_FINAL")

VAL1_ROOT <- file.path(
  ROOT,
  "Validation",
  "GSE155882"
)

VAL2_ROOT <- file.path(
  ROOT,
  "Validation",
  "GSE270896"
)

for (i in 1:7) {
  
  x <- file.path(
    FIG_ROOT,
    sprintf("FIGURE_%02d", i)
  )
  
  dir.create(
    file.path(x, "PNG"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    file.path(x, "PDF"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    file.path(x, "TIFF"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    file.path(x, "SOURCE_DATA"),
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ============================================================
# 3. CONSTANTS
# ============================================================

conditions <- c(
  "Sham",
  "2W",
  "4W",
  "6W"
)

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

selected_genes <- c(
  "Cyp1b1",
  "Cyp4f18",
  "Ephx2"
)

broad_order <- c(
  "Cardiomyocytes",
  "Fibroblasts",
  "Endothelial cells",
  "Pericytes",
  "Smooth muscle cells",
  "Macrophages",
  "Monocytes",
  "Myeloid cells",
  "Dendritic cells",
  "Neutrophils",
  "T cells",
  "NK cells",
  "B cells",
  "Platelets",
  "Erythroid cells"
)

refined_order <- c(
  "Cardiomyocytes",
  "Fibroblasts",
  "Activated fibroblasts",
  "Capillary endothelial cells",
  "Arterial endothelial cells",
  "Endothelial cells",
  "Pericytes",
  "Vascular smooth muscle cells",
  "Macrophages",
  "Resident macrophages",
  "Monocytes",
  "Cycling myeloid cells",
  "Dendritic cells",
  "Neutrophils",
  "Inflammatory neutrophils",
  "T cells",
  "Gamma-delta T cells",
  "Immature T cells",
  "NK/cytotoxic lymphocytes",
  "B cells",
  "Immature B cells",
  "Platelets/megakaryocytes",
  "Erythroid cells"
)


# ============================================================
# 4. PALETTES
# ============================================================

broad_palette <- c(
  "Cardiomyocytes" = "#D73027",
  "Fibroblasts" = "#E78AC3",
  "Endothelial cells" = "#1B9E77",
  "Pericytes" = "#A6761D",
  "Smooth muscle cells" = "#E6AB02",
  "Macrophages" = "#7570B3",
  "Monocytes" = "#666666",
  "Myeloid cells" = "#B15928",
  "Dendritic cells" = "#A6CEE3",
  "Neutrophils" = "#E7298A",
  "T cells" = "#1F78B4",
  "NK cells" = "#33A02C",
  "B cells" = "#6A3D9A",
  "Platelets" = "#FB9A99",
  "Erythroid cells" = "#B2182B"
)

condition_palette <- c(
  "Sham" = "#4D4D4D",
  "2W" = "#4C78A8",
  "4W" = "#E45756",
  "6W" = "#72B7B2"
)


# ============================================================
# 5. GLOBAL THEME
# ============================================================

theme_manuscript <- theme_classic(base_size = 10) +
  theme(
    axis.title = element_text(
      face = "bold",
      colour = "black"
    ),
    axis.text = element_text(
      colour = "black"
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      colour = "black"
    ),
    legend.title = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold",
      size = 11
    ),
    plot.tag = element_text(
      face = "bold",
      size = 16
    ),
    plot.background = element_rect(
      fill = "white",
      colour = NA
    )
  )


# ============================================================
# 6. HELPER FUNCTIONS
# ============================================================

find_file <- function(
    root,
    pattern,
    required = TRUE
) {
  
  files <- list.files(
    root,
    recursive = TRUE,
    full.names = TRUE
  )
  
  hit <- files[
    grepl(
      pattern,
      basename(files),
      ignore.case = TRUE
    )
  ]
  
  if (length(hit) == 0) {
    
    if (required) {
      stop(
        "Could not find file matching: ",
        pattern
      )
    }
    
    return(NA_character_)
  }
  
  hit[1]
}


read_csv_safe <- function(path) {
  
  if (
    length(path) == 0 ||
    is.na(path) ||
    !file.exists(path)
  ) {
    return(NULL)
  }
  
  readr::read_csv(
    path,
    show_col_types = FALSE
  )
}


find_col <- function(
    df,
    candidates,
    required = TRUE
) {
  
  if (is.null(df)) {
    if (required) stop("Data frame is NULL.")
    return(NULL)
  }
  
  nms <- colnames(df)
  
  low <- tolower(nms)
  
  cand_low <- tolower(candidates)
  
  exact <- which(
    low %in% cand_low
  )
  
  if (length(exact) > 0) {
    return(nms[exact[1]])
  }
  
  for (x in cand_low) {
    
    hit <- which(
      grepl(
        x,
        low,
        fixed = TRUE
      )
    )
    
    if (length(hit) > 0) {
      return(nms[hit[1]])
    }
  }
  
  if (required) {
    
    stop(
      "Could not identify required column.\nCandidates: ",
      paste(candidates, collapse = ", "),
      "\nAvailable:\n",
      paste(nms, collapse = ", ")
    )
  }
  
  NULL
}


standardize_condition <- function(x) {
  
  x <- as.character(x)
  
  x <- trimws(x)
  
  x <- gsub(
    "TAC[_ -]*2W",
    "2W",
    x,
    ignore.case = TRUE
  )
  
  x <- gsub(
    "TAC[_ -]*4W",
    "4W",
    x,
    ignore.case = TRUE
  )
  
  x <- gsub(
    "TAC[_ -]*6W",
    "6W",
    x,
    ignore.case = TRUE
  )
  
  x <- gsub(
    "^2[ _-]*weeks?$",
    "2W",
    x,
    ignore.case = TRUE
  )
  
  x <- gsub(
    "^4[ _-]*weeks?$",
    "4W",
    x,
    ignore.case = TRUE
  )
  
  x <- gsub(
    "^6[ _-]*weeks?$",
    "6W",
    x,
    ignore.case = TRUE
  )
  
  x[
    grepl(
      "sham",
      x,
      ignore.case = TRUE
    )
  ] <- "Sham"
  
  factor(
    x,
    levels = conditions
  )
}


save_panel <- function(
    plot,
    figure_number,
    panel_name,
    width,
    height
) {
  
  folder <- file.path(
    FIG_ROOT,
    sprintf(
      "FIGURE_%02d",
      figure_number
    )
  )
  
  ggsave(
    file.path(
      folder,
      "PNG",
      paste0(
        panel_name,
        ".png"
      )
    ),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    file.path(
      folder,
      "PDF",
      paste0(
        panel_name,
        ".pdf"
      )
    ),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )
}


save_final_figure <- function(
    plot,
    figure_number,
    width,
    height
) {
  
  folder <- file.path(
    FIG_ROOT,
    sprintf(
      "FIGURE_%02d",
      figure_number
    )
  )
  
  prefix <- paste0(
    "Figure_",
    figure_number,
    "_FINAL"
  )
  
  ggsave(
    file.path(
      folder,
      "PNG",
      paste0(
        prefix,
        "_300dpi.png"
      )
    ),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    file.path(
      folder,
      "TIFF",
      paste0(
        prefix,
        "_600dpi.tiff"
      )
    ),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    file.path(
      folder,
      "PDF",
      paste0(
        prefix,
        ".pdf"
      )
    ),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )
}


save_source <- function(
    df,
    figure_number,
    name
) {
  
  if (is.null(df)) return(invisible(NULL))
  
  write.csv(
    df,
    file.path(
      FIG_ROOT,
      sprintf(
        "FIGURE_%02d",
        figure_number
      ),
      "SOURCE_DATA",
      paste0(
        name,
        ".csv"
      )
    ),
    row.names = FALSE
  )
}


blank_panel <- function(text) {
  
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = text,
      size = 4,
      fontface = "bold"
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void()
}


png_to_plot <- function(filename) {
  
  if (!file.exists(filename)) {
    return(
      blank_panel(
        paste(
          "Missing:",
          basename(filename)
        )
      )
    )
  }
  
  img <- png::readPNG(filename)
  
  grob <- rasterGrob(
    img,
    interpolate = TRUE
  )
  
  ggplot() +
    annotation_custom(
      grob,
      -Inf,
      Inf,
      -Inf,
      Inf
    ) +
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE
    ) +
    theme_void()
}


# ============================================================
# 7. FROZEN OBJECT CHECK
# ============================================================

discovery_atlas_file <- file.path(
  CLEAN,
  "11D_C_FINAL_CORRECTED_ANNOTATED_ATLAS.rds"
)

spatial_cyp_file <- file.path(
  CLEAN,
  "11C_CYP_AA_ANALYZED_SPATIAL_OBJECTS.rds"
)

spatial_program_file <- file.path(
  CLEAN,
  "11D_E_CORRECTED_CELL_PROGRAM_SCORED_SPATIAL_OBJECTS.rds"
)

spatial_adjusted_file <- file.path(
  CLEAN,
  "11D_F_CELL_PROGRAM_ADJUSTED_CYP_AA_SPATIAL_OBJECTS.rds"
)

mechanistic_file <- file.path(
  CLEAN,
  "11D_G_MECHANISTICALLY_SCORED_SPATIAL_OBJECTS.rds"
)

aa_network_file <- file.path(
  CLEAN,
  "11D_H_FOCUSED_CYP_AA_NETWORK_SPATIAL_OBJECTS.rds"
)

estrogen_file <- file.path(
  CLEAN,
  "11D_I_ESTROGEN_CYP_AA_SCORED_SPATIAL_OBJECTS.rds"
)

required_objects <- c(
  discovery_atlas_file,
  spatial_cyp_file,
  spatial_program_file,
  spatial_adjusted_file,
  mechanistic_file,
  aa_network_file
)

if (!all(file.exists(required_objects))) {
  
  stop(
    "One or more required frozen discovery objects are missing."
  )
}

cat("\nFrozen discovery objects: PASS\n")


# ============================================================
# 8. LOAD DISCOVERY ATLAS
# ============================================================

obj <- readRDS(
  discovery_atlas_file
)

stopifnot(
  ncol(obj) == 25186,
  nrow(obj) == 32285
)

cat(
  "Discovery atlas:",
  ncol(obj),
  "cells;",
  nrow(obj),
  "genes\n"
)


# ============================================================
# 9. DISCOVERY METADATA
# ============================================================

stopifnot(
  "timepoint" %in% colnames(obj@meta.data),
  "corrected_broad_cell_type" %in% colnames(obj@meta.data),
  "corrected_refined_cell_type" %in% colnames(obj@meta.data)
)

obj$manuscript_condition <- standardize_condition(
  obj$timepoint
)

obj$manuscript_broad_cell_type <- factor(
  obj$corrected_broad_cell_type,
  levels = broad_order
)

obj$manuscript_refined_cell_type <- factor(
  obj$corrected_refined_cell_type,
  levels = refined_order
)

stopifnot(
  identical(
    as.integer(
      table(
        obj$manuscript_condition
      )
    ),
    c(
      7578L,
      6051L,
      6374L,
      5183L
    )
  )
)


# ============================================================
# 10. UMAP REDUCTION
# ============================================================

available_reductions <- Reductions(obj)

cat(
  "\nAvailable reductions:",
  paste(
    available_reductions,
    collapse = ", "
  ),
  "\n"
)

if (
  "umap.unintegrated" %in%
  available_reductions
) {
  
  manuscript_umap <- "umap.unintegrated"
  
} else if (
  "umap" %in%
  available_reductions
) {
  
  manuscript_umap <- "umap"
  
} else {
  
  stop(
    "No existing UMAP reduction found. ",
    "UMAP will NOT be recomputed."
  )
}

cat(
  "Using frozen UMAP:",
  manuscript_umap,
  "\n"
)


# ============================================================
# ============================================================
# FIGURE 1
# TEMPORAL CARDIAC CELLULAR LANDSCAPE
# ============================================================
# ============================================================


# ============================================================
# F1A STUDY DESIGN
# ============================================================

timeline <- data.frame(
  x = 1:4,
  Condition = conditions
)

F1A <- ggplot(
  timeline,
  aes(
    x = x,
    y = 1
  )
) +
  
  geom_segment(
    aes(
      x = 1,
      xend = 4,
      y = 1,
      yend = 1
    ),
    linewidth = 0.8,
    arrow = arrow(
      length = unit(
        0.15,
        "inches"
      )
    )
  ) +
  
  geom_point(
    shape = 21,
    size = 10,
    stroke = 1,
    fill = "#F4A6A6"
  ) +
  
  geom_text(
    aes(
      y = 1.35,
      label = Condition
    ),
    fontface = "bold",
    size = 4
  ) +
  
  annotate(
    "text",
    x = 2.5,
    y = 0.55,
    label =
      "scRNA-seq + Visium spatial transcriptomics",
    fontface = "bold",
    size = 4
  ) +
  
  annotate(
    "text",
    x = 2.5,
    y = 0.25,
    label =
      "Public deposition: one scRNA-seq and one spatial library per condition",
    size = 3.2
  ) +
  
  annotate(
    "text",
    x = 2.5,
    y = 0.02,
    label =
      "Deposited metadata do not establish whether libraries represent individual animals or pooled material",
    size = 2.8
  ) +
  
  xlim(
    0.5,
    4.5
  ) +
  
  ylim(
    -0.15,
    1.6
  ) +
  
  labs(
    title =
      "Pressure-overload cardiac remodeling"
  ) +
  
  theme_void() +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )


# ============================================================
# F1B BROAD UMAP
# ============================================================

F1B <- DimPlot(
  obj,
  reduction = manuscript_umap,
  group.by = "manuscript_broad_cell_type",
  cols = broad_palette,
  pt.size = 0.25,
  shuffle = TRUE,
  seed = 12345,
  raster = FALSE
) +
  
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  
  guides(
    colour = guide_legend(
      override.aes = list(
        size = 3
      )
    )
  ) +
  
  theme_manuscript


# ============================================================
# F1C REFINED UMAP
# ============================================================

refined_levels_present <- levels(
  droplevels(
    obj$manuscript_refined_cell_type
  )
)

refined_palette <- setNames(
  grDevices::hcl.colors(
    length(refined_levels_present),
    palette = "Dynamic"
  ),
  refined_levels_present
)

F1C <- DimPlot(
  obj,
  reduction = manuscript_umap,
  group.by = "manuscript_refined_cell_type",
  cols = refined_palette,
  pt.size = 0.25,
  shuffle = TRUE,
  seed = 12345,
  raster = FALSE
) +
  
  labs(