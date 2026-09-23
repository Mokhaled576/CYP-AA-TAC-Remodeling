# Script execution order

## Discovery — GSE308859
Run the numbered scripts in `scripts/01_discovery_GSE308859/` in filename order. The workflow progresses through CYP feasibility, Seurat construction/QC, normalization, clustering, annotation, CYP/AA analyses, spatial analyses, annotation audit/correction, corrected cell-type analyses, cell-program mapping, spatial explanatory models, mechanistic/pathway integration, focused AA-network analyses, and exploratory estrogen-related analyses.

The `11D_*` scripts document the audit/correction chain. Where both an earlier and a corrected downstream analysis exist, the corrected/frozen-annotation outputs are the manuscript-relevant versions.

## Validation — GSE155882
Run V01 → V02A → V02B → V03 → V04 → V05 → V06 → V06B → V07 → V08. `V02_PNG_Figure_Export.R` is a figure-export utility.

## Validation — GSE270896
Run V2_01 → V2_02 → V2_03.

## Manuscript figures
Run scripts under `scripts/09_figures/` after their required analysis outputs exist.

## Paths
Submitted scripts preserve original Windows project paths for provenance. Change project-root/path variables only; do not change analysis thresholds or logic unless intentionally performing a new analysis.
