# ============================================================
# FINAL MANUSCRIPT FIGURE RE-ASSEMBLY
# Figures 1–7
#
# PURPOSE:
# - Keep existing analyses unchanged
# - Standardize A/B/C/... labels
# - Remove excessive whitespace
# - Improve panel distribution
# - Standardize page dimensions
# - Rebuild Figure 7 with better spacing
#
# OUTPUT:
# D:/Master/ScRNA seq/TAC/FIGURES/MANUSCRIPT_REASSEMBLED
# ============================================================


library(ggplot2)
library(patchwork)
library(grid)
library(png)


# ============================================================
# 1. PATHS
# ============================================================

ROOT <- "D:/Master/ScRNA seq/TAC"

OLD_FINAL <- file.path(
  ROOT,
  "FIGURES",
  "MANUSCRIPT_FINAL"
)

OUT <- file.path(
  ROOT,
  "FIGURES",
  "MANUSCRIPT_REASSEMBLED"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. BASIC HELPERS
# ============================================================

require_file <- function(x) {
  
  if (!file.exists(x)) {
    stop(
      paste0(
        "\nFILE NOT FOUND:\n",
        x,
        "\n"
      )
    )
  }
  
  x
}


# ------------------------------------------------------------
# Read PNG and remove transparent/excess outer margins
# ------------------------------------------------------------

read_png_panel <- function(
    path,
    label = NULL,
    label_size = 18
) {
  
  require_file(path)
  
  img <- png::readPNG(path)
  
  p <- ggplot() +
    
    annotation_custom(
      rasterGrob(
        img,
        interpolate = TRUE
      ),
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf
    ) +
    
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE
    ) +
    
    theme_void() +
    
    theme(
      plot.margin = margin(
        0,
        0,
        0,
        0
      )
    )
  
  if (!is.null(label)) {
    
    p <- p +
      
      annotate(
        "text",
        x = 0.005,
        y = 0.995,
        label = label,
        hjust = 0,
        vjust = 1,
        fontface = "bold",
        size = label_size / ggplot2::.pt
      )
  }
  
  p
}


# ------------------------------------------------------------
# Add panel label to an EXISTING ggplot
# ------------------------------------------------------------

label_plot <- function(
    p,
    label,
    size = 18
) {
  
  p +
    labs(
      tag = label
    ) +
    theme(
      plot.tag = element_text(
        face = "bold",
        size = size
      ),
      plot.tag.position = c(
        0,
        1
      ),
      plot.margin = margin(
        8,
        8,
        8,
        8
      )
    )
}


# ------------------------------------------------------------
# Save final page
# ------------------------------------------------------------

save_final <- function(
    p,
    number,
    width,
    height
) {
  
  prefix <- sprintf(
    "Figure_%02d_FINAL_REASSEMBLED",
    number
  )
  
  ggsave(
    filename = file.path(
      OUT,
      paste0(
        prefix,
        ".png"
      )
    ),
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(
      OUT,
      paste0(
        prefix,
        ".pdf"
      )
    ),
    plot = p,
    width = width,
    height = height,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
  
  cat(
    "\nSAVED: ",
    prefix,
    "\n",
    sep = ""
  )
}


# ============================================================
# 3. FIGURE 1
#
# Existing Figure 1 is scientifically fine.
# Problem: final panel arrangement.
#
# Desired:
# A = study design
# B = broad UMAP
# C = refined UMAP
# D = refined atlas by condition
# E = recovered-cell composition
#
# Because the uploaded final composite already contains all
# panels, preserve it for now rather than rerunning analysis.
#
# IMPORTANT:
# If individual F1A-E objects remain in your current session,
# the block below will make the preferred arrangement.
# ============================================================

if (
  exists("F1A") &&
  exists("F1B") &&
  exists("F1C") &&
  exists("F1D") &&
  exists("F1E")
) {
  
  Figure1_FINAL <- (
    label_plot(F1A, "A") /
      (
        label_plot(F1B, "B") |
          label_plot(F1C, "C")
      ) /
      label_plot(F1D, "D") /
      label_plot(F1E, "E")
  ) +
    plot_layout(
      heights = c(
        0.70,
        1.15,
        1.05,
        0.90
      )
    )
  
  save_final(
    Figure1_FINAL,
    1,
    width = 13,
    height = 15
  )
  
} else {
  
  cat(
    "\nFIGURE 1:\n",
    "F1A-E are not all present in the current R session.\n",
    "Existing Figure 1 is retained; do not rerun analysis.\n",
    sep = ""
  )
}


# ============================================================
# 4. FIGURE 2
#
# Use existing objects if available.
#
# A = whole-atlas trajectories
# B = broad-cell heatmap
# C = within-cell changes
# D = cellular expression contribution
# E = selected trajectories
# ============================================================

if (
  exists("F2A") &&
  exists("F2B") &&
  exists("F2C") &&
  exists("F2D") &&
  exists("F2E")
) {
  
  Figure2_FINAL <- (
    label_plot(F2A, "A") /
      label_plot(F2B, "B") /
      label_plot(F2C, "C") /
      (
        label_plot(F2D, "D") |
          label_plot(F2E, "E")
      )
  ) +
    plot_layout(
      heights = c(
        0.85,
        1.20,
        1.15,
        1.05
      )
    )
  
  save_final(
    Figure2_FINAL,
    2,
    width = 14,
    height = 17
  )
  
} else {
  
  cat(
    "\nFIGURE 2:\n",
    "F2A-E are not all present in the current R session.\n",
    "Existing Figure 2 is retained.\n",
    sep = ""
  )
}


# ============================================================
# 5. FIGURE 3
#
# Reassemble using existing analytical plots F3A-C and
# VERIFIED spatial PNGs D-F.
#
# Major change:
# D/E/F are placed side-by-side rather than vertically.
# ============================================================

SPATIAL_DIR <- file.path(
  ROOT,
  "FIGURES",
  "11_SPATIAL",
  "11C_CYP_AA_EICOSANOID"
)

F3D_file <- file.path(
  SPATIAL_DIR,
  "11C_Spatial_Cyp1b1.png"
)

F3E_file <- file.path(
  SPATIAL_DIR,
  "11C_Spatial_Cyp4f18.png"
)

F3F_file <- file.path(
  SPATIAL_DIR,
  "11C_Spatial_Ephx2.png"
)


if (
  exists("F3A") &&
  exists("F3B") &&
  exists("F3C") &&
  file.exists(F3D_file) &&
  file.exists(F3E_file) &&
  file.exists(F3F_file)
) {
  
  P3A <- label_plot(
    F3A,
    "A"
  )
  
  P3B <- label_plot(
    F3B,
    "B"
  )
  
  P3C <- label_plot(
    F3C,
    "C"
  )
  
  P3D <- read_png_panel(
    F3D_file,
    "D"
  )
  
  P3E <- read_png_panel(
    F3E_file,
    "E"
  )
  
  P3F <- read_png_panel(
    F3F_file,
    "F"
  )
  
  
  Figure3_FINAL <- (
    P3A /
      P3B /
      P3C /
      (
        P3D |
          P3E |
          P3F
      )
  ) +
    plot_layout(
      heights = c(
        0.95,
        1.00,
        0.85,
        1.05
      )
    )
  
  save_final(
    Figure3_FINAL,
    3,
    width = 15,
    height = 16
  )
  
} else {
  
  cat(
    "\nWARNING: Figure 3 inputs incomplete.\n"
  )
}


# ============================================================
# 6. FIGURE 4
#
# Existing Figure 4 has correct content but poor page usage.
#
# Better layout:
#
#        A
#        B
#      C | D
#
# Labels are applied manually and consistently.
# ============================================================

F4_SOURCE_E <- file.path(
  ROOT,
  "FIGURES",
  "11_SPATIAL",
  "11D_E_CORRECTED_CELL_PROGRAM_MAPPING"
)

F4_SOURCE_F <- file.path(
  ROOT,
  "FIGURES",
  "11_SPATIAL",
  "11D_F_CELL_PROGRAM_ADJUSTED_CYP_AA"
)


F4A_file <- file.path(
  F4_SOURCE_E,
  "11D_E_Spatial_CellProgram_Trajectories.png"
)

F4B_file <- file.path(
  F4_SOURCE_E,
  "11D_E_4W_Program_CYP_AA_Association_Heatmap.png"
)

F4C_file <- file.path(
  F4_SOURCE_F,
  "11D_F_CellProgram_Model_R2_Heatmap.png"
)

F4D_file <- file.path(
  F4_SOURCE_F,
  "11D_F_CellProgram_Model_R2_Trajectories.png"
)


P4A <- read_png_panel(
  F4A_file,
  "A"
)

P4B <- read_png_panel(
  F4B_file,
  "B"
)

P4C <- read_png_panel(
  F4C_file,
  "C"
)

P4D <- read_png_panel(
  F4D_file,
  "D"
)


Figure4_FINAL <- (
  P4A /
    P4B /
    (
      P4C |
        P4D
    )
) +
  plot_layout(
    heights = c(
      1.15,
      0.90,
      1.00
    )
  )


save_final(
  Figure4_FINAL,
  4,
  width = 14,
  height = 15
)


# ============================================================
# 7. FIGURE 5
#
# This fixes one of the biggest layout problems.
#
# OLD:
# B
# C
# D
# E
#
# NEW:
#
#       B
#       C
#     D | E
#
# This dramatically reduces empty space.
# ============================================================

F5_SOURCE_H <- file.path(
  ROOT,
  "FIGURES",
  "11_SPATIAL",
  "11D_H_FOCUSED_CYP_AA_NETWORK"
)

F5_SOURCE_G <- file.path(
  ROOT,
  "FIGURES",
  "11_SPATIAL",
  "11D_G_MECHANISTIC_PATHWAY_INTEGRATION"
)


F5B_file <- file.path(
  F5_SOURCE_H,
  "11D_H_Spatial_AA_Arm_Trajectories.png"
)

F5C_file <- file.path(
  F5_SOURCE_H,
  "11D_H_scRNA_4W_AA_Arm_Heatmap.png"
)

F5D_file <- file.path(
  F5_SOURCE_H,
  "11D_H_4W_CYP_AA_Arm_Correlation_Heatmap.png"
)

F5E_file <- file.path(
  F5_SOURCE_G,
  "11D_G_4W_CYP_Mechanism_Correlation_Heatmap.png"
)


P5B <- read_png_panel(
  F5B_file,
  "B"
)

P5C <- read_png_panel(
  F5C_file,
  "C"
)

P5D <- read_png_panel(
  F5D_file,
  "D"
)

P5E <- read_png_panel(
  F5E_file,
  "E"
)


Figure5_FINAL <- (
  P5B /
    P5C /
    (
      P5D |
        P5E
    )
) +
  plot_layout(
    heights = c(
      0.90,
      1.05,
      1.05
    )
  )


save_final(
  Figure5_FINAL,
  5,
  width = 14,
  height = 14
)


# ============================================================
# 8. FIGURE 6
#
# Existing F6B-E ggplot objects are preferable because they
# remain vector-quality.
#
# Better arrangement:
#
#       B | C
#         D
#         E
#
# ============================================================

if (
  exists("F6B") &&
  exists("F6C") &&
  exists("F6D") &&
  exists("F6E")
) {
  
  P6B <- label_plot(
    F6B,
    "B"
  )
  
  P6C <- label_plot(
    F6C,
    "C"
  )
  
  P6D <- label_plot(
    F6D,
    "D"
  )
  
  P6E <- label_plot(
    F6E,
    "E"
  )
  
  
  Figure6_FINAL <- (
    (
      P6B |
        P6C
    ) /
      P6D /
      P6E
  ) +
    plot_layout(
      heights = c(
        1.00,
        0.85,
        0.85
      )
    )
  
  
  save_final(
    Figure6_FINAL,
    6,
    width = 14,
    height = 13
  )
  
} else {
  
  cat(
    "\nWARNING: F6B-E are not all available.\n"
  )
}


# ============================================================
# 9. FIGURE 7 — COMPLETELY REFORMATTED
#
# Same biological interpretation.
# Better page use and hierarchy.
#
# IMPORTANT:
# Lines represent conceptual organization of observed results,
# not demonstrated causal relationships.
# ============================================================


F7 <- ggplot() +
  
  # ----------------------------------------------------------
# TOP ROW — FOUR MAIN COMPONENTS
# ----------------------------------------------------------

annotate(
  "rect",
  xmin = 0.5,
  xmax = 2.8,
  ymin = 7.1,
  ymax = 8.6,
  fill = "grey96",
  colour = "black",
  linewidth = 0.6
) +
  
  annotate(
    "text",
    x = 1.65,
    y = 7.85,
    label = "Pressure overload\n(TAC)",
    fontface = "bold",
    size = 5.0
  ) +
  
  
  annotate(
    "rect",
    xmin = 3.4,
    xmax = 6.1,
    ymin = 7.1,
    ymax = 8.6,
    fill = "grey96",
    colour = "black",
    linewidth = 0.6
  ) +
  
  annotate(
    "text",
    x = 4.75,
    y = 7.85,
    label = "Temporal cardiac\ncellular remodeling",
    fontface = "bold",
    size = 4.8
  ) +
  
  
  annotate(
    "rect",
    xmin = 6.7,
    xmax = 10.1,
    ymin = 7.1,
    ymax = 8.6,
    fill = "grey96",
    colour = "black",
    linewidth = 0.6
  ) +
  
  annotate(
    "text",
    x = 8.4,
    y = 7.85,
    label = paste0(
      "Cell-type-specific CYP/AA\n",
      "transcriptional reorganization"
    ),
    fontface = "bold",
    size = 4.5
  ) +
  
  
  annotate(
    "rect",
    xmin = 10.7,
    xmax = 13.5,
    ymin = 7.1,
    ymax = 8.6,
    fill = "grey96",
    colour = "black",
    linewidth = 0.6
  ) +
  
  annotate(
    "text",
    x = 12.1,
    y = 7.85,
    label = "Spatial CYP/AA\nremodeling",
    fontface = "bold",
    size = 4.8
  ) +
  
  
  # ----------------------------------------------------------
# CONNECTORS — NO ARROWHEADS
# ----------------------------------------------------------

annotate(
  "segment",
  x = 2.8,
  xend = 3.4,
  y = 7.85,
  yend = 7.85,
  linewidth = 0.7
) +
  
  annotate(
    "segment",
    x = 6.1,
    xend = 6.7,
    y = 7.85,
    yend = 7.85,
    linewidth = 0.7
  ) +
  
  annotate(
    "segment",
    x = 10.1,
    xend = 10.7,
    y = 7.85,
    yend = 7.85,
    linewidth = 0.7
  ) +
  
  
  # ----------------------------------------------------------
# CENTRAL 4W STATE
# ----------------------------------------------------------

annotate(
  "rect",
  xmin = 2.2,
  xmax = 11.8,
  ymin = 5.0,
  ymax = 6.25,
  fill = "grey92",
  colour = "black",
  linewidth = 0.7
) +
  
  annotate(
    "text",
    x = 7.0,
    y = 5.63,
    label = paste0(
      "4 weeks: prominent intermediate ",
      "transcriptional/spatial remodeling state"
    ),
    fontface = "bold",
    size = 5.0
  ) +
  
  
  # ----------------------------------------------------------
# REPRESENTATIVE CELL/GENE REMODELING
# ----------------------------------------------------------

annotate(
  "text",
  x = 7.0,
  y = 4.35,
  label = "Representative cell-type-specific remodeling",
  fontface = "bold",
  size = 4.4
) +
  
  annotate(
    "text",
    x = 7.0,
    y = 3.85,
    label = paste0(
      "Endothelial: Cyp1b1, Ephx2     •     ",
      "Macrophages: Cyp4f18     •     ",
      "Monocytes: Cyp4f16     •     ",
      "Fibroblasts: Cyp1b1"
    ),
    size = 3.9
  ) +
  
  
  # ----------------------------------------------------------
# BOTTOM LEFT
# ----------------------------------------------------------

annotate(
  "rect",
  xmin = 0.8,
  xmax = 6.7,
  ymin = 1.6,
  ymax = 3.0,
  fill = "grey98",
  colour = "black",
  linewidth = 0.6
) +
  
  annotate(
    "text",
    x = 3.75,
    y = 2.48,
    label = "Spatial cellular-program context",
    fontface = "bold",
    size = 4.3
  ) +
  
  annotate(
    "text",
    x = 3.75,
    y = 2.03,
    label = paste0(
      "Broad mapped cellular programs captured only a small\n",
      "fraction of spatial CYP/AA heterogeneity (R² ≈ 0.006–0.051)"
    ),
    size = 3.8
  ) +
  
  
  # ----------------------------------------------------------
# BOTTOM RIGHT
# ----------------------------------------------------------

annotate(
  "rect",
  xmin = 7.3,
  xmax = 13.2,
  ymin = 1.6,
  ymax = 3.0,
  fill = "grey98",
  colour = "black",
  linewidth = 0.6
) +
  
  annotate(
    "text",
    x = 10.25,
    y = 2.48,
    label = "Independent TAC datasets",
    fontface = "bold",
    size = 4.3
  ) +
  
  annotate(
    "text",
    x = 10.25,
    y = 2.03,
    label = paste0(
      "Selected components showed independent support:\n",
      "Cyp4f13 ↓     •     Cyp4f18 ↑ in GSE155882     •     Ephx2 ↓"
    ),
    size = 3.8
  ) +
  
  
  # ----------------------------------------------------------
# FOOTNOTE
# ----------------------------------------------------------

annotate(
  "text",
  x = 7.0,
  y = 0.72,
  label = paste0(
    "AA branches represent transcriptomic programs; ",
    "metabolite concentrations, CYP activity, and metabolic flux were not measured."
  ),
  fontface = "italic",