# Expanded analyses version 2.0.0

This folder contains the scientific code and all supporting tables for Cellular and spatial cytochrome P450 and epoxide hydrolase 2 expression during cardiac pressure overload across independent cohorts. It supports either journal-formatted manuscript version and contains no journal submission documents or private reviewer correspondence.

## Download and extract

Download Analysis_Code_v2.zip and every Supporting_Tables_Part ZIP in this folder. Extract all into the SAME empty directory, retaining folder names. The table parts are independent ZIP files, not split binary volumes. Merge the Supporting_Analyses and TAC_Public_Expansion directories when prompted. No table records have been removed or changed by this repackaging. Source-table duplicates shared by the code and table bundles are identical.

Supporting_Analyses contains the scientific discovery diagnostics and earlier external comparisons. TAC_Public_Expansion contains the six additional studies' R code, result tables, input extracts and acquisition manifest. TABLE_INDEX.csv within the table archives inventories the results. SHA256SUMS.csv in this folder verifies the download archives.

## Reproduce

Set TAC_SOURCE_DIR to the original processed Seurat project tree for Supporting_Analyses/scripts. The original preparation scripts remain in the repository's scripts directory; consult docs/RUN_ORDER.md and their local path definitions. For the expansion open TAC_Public_Expansion/Expansion.Rproj and run scripts/run_analyses.R after acquiring the required raw inputs from the acquisition manifest. Some original scripts retain Windows paths and library preferences. Large raw downloads and processed Seurat objects are excluded. A complete clean-room reconstruction has not been verified. All 23 expanded R scripts passed syntax checks; no new biological model was run for this deposit.

## Evidence and provenance

The earlier three datasets are GSE308859, GSE155882 and GSE270896. Added cohorts are GSE296709, GSE120064, GSE95140, GSE262690, PXD061784 and Mendeley Data bfry768hk7 version 1. The supporting tables report concordant, discordant and inconclusive results. Discovery cells and spots are not biological replicates. The independent endothelial study does not confirm the original decreases. Human confounding and sparse detection remain limitations. Whole-heart proteins and three-week lipid abundance do not establish lineage-specific enzyme activity or targeted EET/DHET/HETE concentrations. No estrogen mechanism is retained in the current paper. Earlier historical scripts remain as provenance, not as evidence for current hormone claims.

Code retains the existing MIT license. Third-party datasets retain their original source rights and terms; the code license does not relicense those datasets. Source information and acquisition records accompany the data. AI assistance with code and manuscript preparation is disclosed in the manuscript; human scientific verification remains required.

The earlier archive https://doi.org/10.5281/zenodo.22912931 contains version v1.0.1 and does not contain this expansion. A new exact-version DOI must be verified after the v2.0.0 archive is created.
