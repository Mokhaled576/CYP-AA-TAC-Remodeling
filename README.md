# CYP-AA-TAC-Remodeling

Reproducibility repository for the manuscript:

**Cellular and Spatial Reorganization of the Cytochrome P450–Arachidonic Acid Axis During Pressure-Overload Cardiac Remodeling**

Authors: Mostafa K. Abd El-Aziz and Ayman O. S. El-Kadi

## Overview

This repository contains analysis code, gene-set definitions, annotation resources, figure source-data tables, and reproducibility documentation for a computational study of cytochrome P450 (CYP)-mediated arachidonic acid (AA) transcriptional remodeling during transverse aortic constriction (TAC).

The study integrates:
- single-cell RNA sequencing
- Visium spatial transcriptomics
- longitudinal pressure-overload remodeling
- independent validation datasets

## Public datasets

The analyses use publicly available Gene Expression Omnibus datasets:
- **GSE308859** — discovery dataset with longitudinal scRNA-seq and spatial transcriptomics across Sham, 2-week, 4-week, and 6-week TAC
- **GSE155882** — independent TAC scRNA-seq validation dataset
- **GSE270896** — independent TAC snRNA-seq validation dataset

Raw public data are not duplicated in this repository. Scripts and documentation will specify the required inputs and processing steps.

## Priority CYP/AA panel

The focused CYP/AA panel contains:
- Cyp1b1
- Cyp2j6
- Cyp2j9
- Cyp4f13
- Cyp4f16
- Cyp4f17
- Cyp4f18
- Ephx2

This repository also includes broader AA-metabolic and remodeling-associated transcriptional gene sets used for pathway-level analyses.

## Repository structure

```
.
├── README.md
├── CITATION.cff
├── LICENSE
├── .gitignore
├── annotations/
├── docs/
├── environment/
├── figures/
│   ├── main/
│   └── supplementary/
├── gene_sets/
├── scripts/
│   ├── 01_scRNA_processing/
│   ├── 02_cell_annotation/
│   ├── 03_CYP_AA_analysis/
│   ├── 04_spatial_analysis/
│   ├── 05_AA_programs/
│   ├── 06_estrogen_exploratory/
│   ├── 07_validation_GSE155882/
│   ├── 08_validation_GSE270896/
│   └── 09_figures/
└── source_data/
    ├── main_figures/
    └── supplementary_figures/
```

## Interpretation boundaries

These analyses quantify RNA abundance and transcriptomic program scores. They do **not** directly measure:
- CYP enzyme activity
- EET, DHET, or HETE concentrations
- metabolic flux
- causal regulation
- estrogen dependence

Spatial associations are interpreted as transcriptomic co-variation, not causal or cell-intrinsic effects.

## Reproducibility

The final repository will contain:
- analysis scripts
- frozen cell-type annotation maps
- gene-set definitions
- figure source-data tables
- software/package versions
- computational workflow documentation

## Status

This repository is currently under preparation for manuscript submission.
