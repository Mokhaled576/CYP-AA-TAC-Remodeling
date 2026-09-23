# ============================================================
# STEP 11D-F
# CELL-PROGRAM-ADJUSTED CYP/AA SPATIAL REMODELING
#
# PURPOSE
# ------------------------------------------------------------
# Determine how strongly spatial CYP/AA expression is captured
# by the corrected 15-program cellular landscape, and quantify
# the remaining cell-program-adjusted spatial signal.
#
# IMPORTANT INTERPRETATION
# ------------------------------------------------------------
# Residual signal does NOT prove cell-intrinsic regulation.
#
# Visium spots contain mixtures of cells and the program scores
# are expression-derived proxies, not measured cell fractions.
#
# Therefore use:
#   "cell-program-adjusted spatial signal"
# or
#   "signal not captured by mapped cellular programs"
#
# Do NOT use:
#   "cell-intrinsic expression"
#   "composition-independent expression"
#
# NO condition-level inferential statistics are performed.
# There is one deposited spatial library per condition.
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Matrix)


# ============================================================
# 2. PATHS
# ============================================================

project_dir <- "D:/Master/ScRNA seq/TAC"

input_file <- file.path(
  project_dir,
  "CLEAN DATA",
  "11D_E_CORRECTED_CELL_PROGRAM_SCORED_SPATIAL_OBJECTS.rds"
)

results_dir <- file.path(
  project_dir,
  "RESULTS",
  "11_SPATIAL",
  "11D_F_CELL_PROGRAM_ADJUSTED_CYP_AA"
)

figures_dir <- file.path(
  project_dir,
  "FIGURES",
  "11_SPATIAL",
  "11D_F_CELL_PROGRAM_ADJUSTED_CYP_AA"
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


# ============================================================
# 3. LOAD 11D-E SPATIAL OBJECTS
# ============================================================

cat("\n============================================\n")
cat("STEP 11D-F\n")
cat("CELL-PROGRAM-ADJUSTED CYP/AA SPATIAL ANALYSIS\n")
cat("============================================\n")

spatial_list <- readRDS(input_file)

if (!is.list(spatial_list)) {
  stop("Expected a list of spatial Seurat objects.")
}

cat("\nSpatial objects:", length(spatial_list), "\n")
print(names(spatial_list))


# ============================================================
# 4. VERIFY SPOT COUNTS
# ============================================================

spot_counts <- sapply(
  spatial_list,
  ncol
)

cat("\nSpatial spot counts:\n")
print(spot_counts)

if (sum(spot_counts) != 7079) {
  stop(
    paste(
      "Unexpected total spatial spot count:",
      sum(spot_counts)
    )
  )
}


# ============================================================
# 5. CONDITION LEVELS
# ============================================================

condition_levels <- c(
  "Sham",
  "2W",
  "4W",
  "6W"
)

for (i in seq_along(spatial_list)) {
  
  if (!"Condition" %in%
      colnames(spatial_list[[i]]@meta.data)) {
    
    stop(
      paste(
        "Condition metadata missing from",
        names(spatial_list)[i]
      )
    )
  }
}

observed_conditions <- unique(
  unlist(
    lapply(
      spatial_list,
      function(obj) {
        as.character(
          unique(obj$Condition)
        )
      }
    )
  )
)

if (!setequal(
  observed_conditions,
  condition_levels
)) {
  stop("Expected Sham, 2W, 4W and 6W.")
}


# ============================================================
# 6. IDENTIFY PROGRAM COLUMNS
# ============================================================

program_sets <- lapply(
  spatial_list,
  function(obj) {
    grep(
      "^Program_",
      colnames(obj@meta.data),
      value = TRUE
    )
  }
)

program_columns <- Reduce(
  intersect,
  program_sets
)

# Exclude derived metadata if any names ever match Program_
program_columns <- setdiff(
  program_columns,
  c(
    "Program",
    "Program_Score"
  )
)

cat(
  "\nShared spatial program columns:",
  length(program_columns),
  "\n"
)

print(program_columns)

if (length(program_columns) != 15) {
  stop(
    paste(
      "Expected 15 shared program columns; found",
      length(program_columns)
    )
  )
}


# ============================================================
# 7. HUMAN-READABLE PROGRAM NAMES
# ============================================================

clean_program_name <- function(x) {
  
  x <- sub(
    "^Program_",
    "",
    x
  )
  
  gsub(
    "\\.",
    " ",
    x
  )
}


program_name_map <- data.frame(
  Program_Column = program_columns,
  Program = clean_program_name(program_columns),
  stringsAsFactors = FALSE
)

write.csv(
  program_name_map,
  file.path(
    results_dir,
    "11D_F_Program_Name_Map.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 8. PRIORITY CYP/AA GENES
# ============================================================

priority_requested <- c(
  "Cyp1b1",
  "Cyp2j6",
  "Cyp2j9",
  "Cyp4f13",
  "Cyp4f16",
  "Cyp4f17",
  "Cyp4f18",
  "Ephx2"
)

shared_genes <- Reduce(
  intersect,
  lapply(
    spatial_list,
    rownames
  )
)

priority_genes <- intersect(
  priority_requested,
  shared_genes
)

cat(
  "\nPriority genes available:",
  length(priority_genes),
  "\n"
)

print(priority_genes)

if (!setequal(
  priority_genes,
  priority_requested
)) {
  
  stop(
    paste(
      "Missing priority genes:",
      paste(
        setdiff(
          priority_requested,
          priority_genes
        ),
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 9. ANALYSIS STRATEGY
#
# For each spatial library and each priority gene:
#
# Y = normalized gene expression across spots
#
# X = 15 corrected spatial program scores
#
# Because these programs are correlated, ordinary regression
# coefficients are NOT interpreted biologically.
#
# We use the model only to estimate:
#
# 1. R2:
#    fraction of spatial expression variation captured by the
#    mapped program landscape.
#
# 2. Adjusted R2.
#
# 3. Residual:
#    expression variation not captured by the 15-program model.
#
# No condition comparison p-values.
# ============================================================


# ============================================================
# 10. HELPER: SAFE MODEL
# ============================================================

fit_program_model <- function(
    gene_expression,
    program_df
) {
  
  model_df <- data.frame(
    GeneExpression =
      as.numeric(gene_expression),
    program_df,
    check.names = FALSE
  )
  
  # Remove zero-variance predictors
  predictor_variance <- sapply(
    model_df[
      ,
      -1,
      drop = FALSE
    ],
    function(x) {
      var(
        x,
        na.rm = TRUE
      )
    }
  )
  
  keep_predictors <- names(
    predictor_variance[
      is.finite(predictor_variance) &
        predictor_variance > 0
    ]
  )
  
  if (length(keep_predictors) == 0) {
    
    return(
      list(
        model = NULL,
        r2 = NA_real_,
        adjusted_r2 = NA_real_,
        residuals = rep(
          NA_real_,
          nrow(model_df)
        ),
        fitted = rep(
          NA_real_,
          nrow(model_df)
        ),
        predictors_used = 0
      )
    )
  }
  
  model_df_use <- model_df[
    ,
    c(
      "GeneExpression",
      keep_predictors
    ),
    drop = FALSE
  ]
  
  model_formula <- reformulate(
    termlabels =
      keep_predictors,
    response =
      "GeneExpression"
  )
  
  fit <- lm(
    model_formula,
    data = model_df_use
  )
  
  fit_summary <- summary(fit)
  
  list(
    model = fit,
    r2 =
      unname(
        fit_summary$r.squared
      ),
    adjusted_r2 =
      unname(
        fit_summary$adj.r.squared
      ),
    residuals =
      residuals(fit),
    fitted =
      fitted(fit),
    predictors_used =
      length(keep_predictors)
  )
}


# ============================================================
# 11. RUN MODELS
# ============================================================

model_summary_list <- list()
residual_long_list <- list()

summary_counter <- 1
residual_counter <- 1

cat("\n============================================\n")
cat("FITTING CELL-PROGRAM MODELS\n")
cat("============================================\n")

for (i in seq_along(spatial_list)) {
  
  obj <- spatial_list[[i]]
  
  DefaultAssay(obj) <- "Spatial"
  
  condition_now <- unique(
    as.character(
      obj$Condition
    )
  )
  
  if (length(condition_now) != 1) {
    stop(
      paste(
        names(spatial_list)[i],
        "has non-unique Condition."
      )
    )
  }
  
  spatial_data <- GetAssayData(
    obj,
    assay = "Spatial",
    layer = "data"
  )
  
  program_df <- obj@meta.data[
    ,
    program_columns,
    drop = FALSE
  ]
  
  # Make syntactically safe names for lm()
  colnames(program_df) <- make.names(
    colnames(program_df),
    unique = TRUE
  )
  
  for (gene_now in priority_genes) {
    
    gene_expression <- as.numeric(
      spatial_data[
        gene_now,
        ,
        drop = TRUE
      ]
    )
    
    detection_percent <-
      100 *
      mean(
        gene_expression > 0
      )
    
    mean_expression <-
      mean(
        gene_expression
      )
    
    median_expression <-
      median(
        gene_expression
      )
    
    expression_variance <-
      var(
        gene_expression
      )
    
    fit_now <- fit_program_model(
      gene_expression =
        gene_expression,
      program_df =
        program_df
    )
    
    model_summary_list[[summary_counter]] <-
      data.frame(
        
        Spatial_Object =
          names(spatial_list)[i],
        
        Condition =
          condition_now,
        
        Gene =
          gene_now,
        
        Spots =
          ncol(obj),
        
        Percent_Expressing =
          detection_percent,
        
        Mean_Expression =
          mean_expression,
        
        Median_Expression =
          median_expression,
        
        Expression_Variance =
          expression_variance,
        
        Predictors_Used =
          fit_now$predictors_used,
        
        R2 =
          fit_now$r2,
        
        Adjusted_R2 =
          fit_now$adjusted_r2,
        
        Unexplained_Fraction =
          ifelse(
            is.na(fit_now$r2),
            NA_real_,
            1 - fit_now$r2
          ),
        
        stringsAsFactors = FALSE
      )
    
    summary_counter <- summary_counter + 1
    
    residual_sd <- sd(
      fit_now$residuals,
      na.rm = TRUE
    )
    
    if (
      is.finite(residual_sd) &&
      residual_sd > 0
    ) {
      
      standardized_residual <-
        fit_now$residuals /
        residual_sd
      
    } else {
      
      standardized_residual <-
        rep(
          NA_real_,
          length(fit_now$residuals)
        )
    }
    
    residual_long_list[[residual_counter]] <-
      data.frame(
        
        Spot =
          colnames(obj),
        
        Spatial_Object =
          names(spatial_list)[i],
        
        Condition =
          condition_now,
        
        Gene =
          gene_now,
        
        Observed_Expression =
          gene_expression,
        
        Fitted_Program_Component =
          fit_now$fitted,
        
        Adjusted_Residual =
          fit_now$residuals,
        
        Standardized_Residual =
          standardized_residual,
        
        stringsAsFactors = FALSE
      )
    
    residual_counter <- residual_counter + 1
  }
  
  cat(
    names(spatial_list)[i],
    "- completed\n"
  )
}


# ============================================================
# 12. COMBINE MODEL RESULTS
# ============================================================

model_summary <- bind_rows(
  model_summary_list
)

model_summary$Condition <- factor(
  model_summary$Condition,
  levels = condition_levels
)

model_summary <- model_summary %>%
  arrange(
    Gene,
    Condition
  )

residual_long <- bind_rows(
  residual_long_list
)

residual_long$Condition <- factor(
  residual_long$Condition,
  levels = condition_levels
)

write.csv(
  model_summary,
  file.path(
    results_dir,
    "11D_F_CellProgram_Model_Summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  residual_long,
  file.path(
    results_dir,
    "11D_F_All_PriorityGene_Adjusted_Residuals.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. PRINT MODEL EXPLANATORY POWER
# ============================================================

cat("\n============================================\n")
cat("CELL-PROGRAM MODEL EXPLANATORY POWER\n")
cat("============================================\n")

print(
  as.data.frame(
    model_summary %>%
      select(
        Gene,
        Condition,
        Spots,
        Percent_Expressing,
        Mean_Expression,
        R2,
        Adjusted_R2,
        Unexplained_Fraction
      )
  ),
  row.names = FALSE
)


# ============================================================
# 14. MODEL R2 HEATMAP
# ============================================================

p_r2 <- ggplot(
  model_summary,
  aes(
    x = Condition,
    y = Gene,
    fill = R2
  )
) +
  geom_tile() +
  labs(
    title =
      "Spatial CYP/AA variation captured by corrected cell programs",
    subtitle =
      "R² from within-section 15-program models",
    x = NULL,
    y = NULL,
    fill = "R²"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )

ggsave(
  file.path(
    figures_dir,
    "11D_F_CellProgram_Model_R2_Heatmap.png"
  ),
  p_r2,
  width = 7,
  height = 6,
  dpi = 400
)


# ============================================================
# 15. MODEL R2 TRAJECTORIES
# ============================================================

p_r2_traj <- ggplot(
  model_summary,
  aes(
    x = Condition,
    y = R2,
    group = Gene
  )
) +
  geom_line() +
  geom_point() +
  facet_wrap(
    ~ Gene,
    scales = "free_y"
  ) +
  labs(
    title =
      "Cell-program explanatory power across TAC progression",
    subtitle =
      "Descriptive within-section R²",
    x = NULL,
    y = "R²"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor =
      element_blank()
  )

ggsave(
  file.path(
    figures_dir,
    "11D_F_CellProgram_Model_R2_Trajectories.png"
  ),
  p_r2_traj,
  width = 12,
  height = 8,
  dpi = 400
)


# ============================================================
# 16. RESIDUAL SUMMARY
#
# NOTE:
# Mean ordinary least-squares residual within each fitted
# section is expected to be approximately zero by construction.
#
# Therefore biologically useful residual summaries are:
# - mean absolute residual
# - residual SD
# - upper/lower residual-tail fractions
#
# We do NOT compare mean residual across conditions.
# ============================================================

residual_summary <- residual_long %>%
  
  group_by(
    Condition,
    Gene
  ) %>%
  
  summarise(
    
    Spots =
      n(),
    
    Mean_Absolute_Residual =
      mean(
        abs(Adjusted_Residual),
        na.rm = TRUE
      ),
    
    Residual_SD =
      sd(
        Adjusted_Residual,
        na.rm = TRUE
      ),
    
    Median_Absolute_Residual =
      median(
        abs(Adjusted_Residual),
        na.rm = TRUE
      ),
    
    Percent_High_Positive_Residual =
      100 *
      mean(
        Standardized_Residual >= 2,
        na.rm = TRUE
      ),
    
    Percent_High_Negative_Residual =
      100 *
      mean(
        Standardized_Residual <= -2,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

residual_summary$Condition <- factor(
  residual_summary$Condition,
  levels = condition_levels
)

residual_summary <- residual_summary %>%
  arrange(
    Gene,
    Condition
  )

write.csv(
  residual_summary,
  file.path(
    results_dir,
    "11D_F_Adjusted_Residual_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. PRINT RESIDUAL SUMMARY
# ============================================================

cat("\n============================================\n")
cat("CELL-PROGRAM-ADJUSTED RESIDUAL SUMMARY\n")
cat("============================================\n")

print(
  as.data.frame(
    residual_summary
  ),
  row.names = FALSE
)


# ============================================================
# 18. RESIDUAL DISPERSION RELATIVE TO SHAM
#
# Descriptive only.
#
# This asks whether spatial heterogeneity not captured by the
# program model becomes larger or smaller than in Sham.
# ============================================================

sham_residual <- residual_summary %>%
  
  filter(
    Condition == "Sham"
  ) %>%
  
  select(
    Gene,
    Sham_Mean_Absolute_Residual =
      Mean_Absolute_Residual,
    Sham_Residual_SD =
      Residual_SD
  )

residual_change <- residual_summary %>%
  
  left_join(
    sham_residual,
    by = "Gene"
  ) %>%
  
  mutate(
    
    AbsoluteResidual_Fold_vs_Sham =
      (Mean_Absolute_Residual + 1e-8) /
      (Sham_Mean_Absolute_Residual + 1e-8),
    
    Log2_AbsoluteResidual_Fold_vs_Sham =
      log2(
        AbsoluteResidual_Fold_vs_Sham
      ),
    
    ResidualSD_Fold_vs_Sham =
      (Residual_SD + 1e-8) /
      (Sham_Residual_SD + 1e-8),
    
    Log2_ResidualSD_Fold_vs_Sham =
      log2(
        ResidualSD_Fold_vs_Sham
      )
  )

write.csv(
  residual_change,
  file.path(
    results_dir,
    "11D_F_Adjusted_Residual_Change_vs_Sham.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. 4W ADJUSTED SPATIAL STATE
# ============================================================

state_4w <- model_summary %>%
  
  filter(
    Condition == "4W"
  ) %>%
  
  select(
    Gene,
    Spots,
    Percent_Expressing,
    Mean_Expression,
    R2,
    Adjusted_R2,
    Unexplained_Fraction
  ) %>%
  
  left_join(
    residual_change %>%
      filter(
        Condition == "4W"
      ) %>%
      select(
        Gene,
        Mean_Absolute_Residual,
        Residual_SD,
        Percent_High_Positive_Residual,
        Percent_High_Negative_Residual,
        Log2_AbsoluteResidual_Fold_vs_Sham,
        Log2_ResidualSD_Fold_vs_Sham
      ),
    by = "Gene"
  ) %>%
  
  arrange(
    desc(
      Unexplained_Fraction
    )
  )

write.csv(
  state_4w,
  file.path(
    results_dir,