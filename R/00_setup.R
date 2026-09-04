## =============================================================================
## 00_setup.R
## Loads all packages used across the pipeline. Source this once per session,
## before sourcing config.R and running scripts in scripts/.
## =============================================================================

required_packages <- c(
  "data.table", "dplyr", "tidyr", "stringr", "purrr", "tibble", "forcats",
  "scales", "lubridate", "ggplot2", "patchwork",
  "igraph", "ggraph", "networkD3", "TraMineR",
  "lme4", "lmerTest", "emmeans", "broom", "broom.mixed",
  "survival", "survminer", "arsenal", "ggpubr"
)

missing_pkgs <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs)
}
invisible(lapply(required_packages, library, character.only = TRUE))

## tSPM+ and MLHO are not on CRAN.
## install.packages("devtools")
## devtools::install_github("<org>/tSPMPlus")
## devtools::install_github("hestiri/mlho")
suppressWarnings({
  library(tSPMPlus)
  library(mlho)
})

set.seed(42)
