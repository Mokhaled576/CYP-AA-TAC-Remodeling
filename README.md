# CYP-AA-TAC-Remodeling

Reproducibility repository for:

**Cellular and Spatial Reorganization of the Cytochrome P450–Arachidonic Acid Axis During Pressure-Overload Cardiac Remodeling**

**Authors:** Mostafa K. Abd El-Aziz and Ayman O. S. El-Kadi

## Study overview
This computational study defines cellular, temporal, and spatial transcriptional remodeling of the cardiac cytochrome P450 (CYP)–arachidonic acid (AA) axis during transverse aortic constriction (TAC).

### Datasets
- **GSE308859** — discovery: longitudinal mouse scRNA-seq + Visium spatial transcriptomics; Sham, 2W, 4W, 6W.
- **GSE155882** — independent scRNA-seq validation; untreated Sham vs TAC biological libraries.
- **GSE270896** — independent snRNA-seq validation; lean C57BL/6J Sham vs TAC, 10 weeks after surgery.

Raw GEO data and large Seurat/RDS objects are intentionally not duplicated here.

## Predefined CYP/AA panel
`Cyp1b1`, `Cyp2j6`, `Cyp2j9`, `Cyp4f13`, `Cyp4f16`, `Cyp4f17`, `Cyp4f18`, and `Ephx2`.

## Repository contents
- `scripts/01_discovery_GSE308859/` — exact submitted discovery scripts, including QC, annotation, CYP/AA, spatial, corrected-annotation, pathway and exploratory estrogen analyses.
- `scripts/07_validation_GSE155882/` — exact submitted GSE155882 validation workflow.
- `scripts/08_validation_GSE270896/` — exact submitted GSE270896 validation workflow.
- `scripts/09_figures/` — final manuscript figure assembly/cleanup scripts.
- `gene_sets/` — priority panel and gene-set documentation.
- `source_data/` — frozen headline source-data tables currently available.
- `docs/` — dataset, workflow, statistical-framework and result documentation.
- `environment/` — commands/instructions for recording the R environment.

## Statistical guardrails
The GSE308859 public deposition contains one scRNA library and one spatial library per condition. Discovery analyses are therefore descriptive; cells and spots are not treated as biological replicates.

For validation, biological replicate is the mouse/library. The GSE155882 and GSE270896 comparisons each use 2 Sham and 2 TAC biological libraries. Interpretation emphasizes effect direction, magnitude, replicate concordance, detection/evaluability and cell-type localization.

## Interpretation boundaries
These analyses quantify RNA abundance and transcriptomic program scores. They do **not** directly establish CYP enzyme activity, EET/DHET/HETE concentrations, metabolic flux, causal regulation, or estrogen dependence. Spatial associations are interpreted as transcriptomic co-variation.

## Reproduction
Scripts retain the original analysis logic and may contain local Windows project paths. To reproduce, clone this repository, download the public GEO inputs, and edit only the project-root/path definitions to match the local filesystem. Run scripts in the order described in `docs/RUN_ORDER.md`.

## Citation
Please cite the associated manuscript when available. Repository citation metadata are provided in `CITATION.cff`.

## License
Code is provided under the MIT License. Public datasets remain subject to their source repositories' terms.
