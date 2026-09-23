# Computational environment

The repository contains the exact analysis scripts, but the exact final R session/package snapshot must be generated on the computer that ran the final analysis.

From the repository root, run:

```r
source("environment/REPRODUCIBILITY_COMMANDS.R")
```

This creates:
- `environment/R_sessionInfo.txt`
- `environment/package_versions.csv`

Commit those two generated files before creating the archival release. They are intentionally not fabricated from script imports because installed package versions must reflect the actual final analysis environment.
