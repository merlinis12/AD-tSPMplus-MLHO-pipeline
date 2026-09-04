# Temporal Sequence Mining and Machine-Learning Pipeline for Predicting Alzheimer's Disease from Longitudinal Clinical Data

This repository contains the analysis pipeline used in the associated manuscript. It is provided
as **pseudo-code / a generalized template** describing the analytical steps, not as a
plug-and-run tool: all references to local file paths, institutional cohort files, and
proprietary/restricted-access datasets have been removed or replaced with generic
placeholders. Users applying this pipeline to their own data must supply their own
input files and adapt column names to match their data dictionary (see `config.R`).

## Overview

The pipeline takes a longitudinal "event table" (one row per patient per clinical event:
diagnosis, medication, lab value, etc.) and an outcome label (e.g., incident diagnosis of
Alzheimer's disease, "AD") and:

1. Mines pairwise temporal event sequences ("A occurs before B") using the
   [tSPM+](https://github.com/JonasHuegel/tSPMPlus_R) algorithm.
2. Performs iterative, cross-validated feature selection over these sequences with
   [MLHO](https://github.com/hestiri/MLHO) to identify the transitions most predictive
   of the outcome.
3. Builds a directed graph of the selected transitions and enumerates candidate
   multi-step pathways, testing their prevalence and association with the outcome.
4. Reconstructs per-patient clinical trajectories (ordered sequences of dysregulated
   domains) and summarizes/visualizes them (Sankey diagrams, sequence-state plots,
   stratified frequency plots) by demographic/genetic subgroup.
5. Models longitudinal cognitive trajectories (linear mixed-effects models) as a
   function of trajectory group, comorbidity burden, sex, and genetic risk.
6. Performs internal/external validation by testing whether identified trajectory
   groups differ in independent biomarker panels (e.g., amyloid, glial injury,
   inflammatory markers, pTau).

## Repository structure

```
.
├── config.R                                   # user-editable paths & parameters (placeholders only)
├── R/                                          # reusable functions, no execution side effects
│   ├── 00_setup.R                              # package loading
│   ├── utils_sequence_mining.R                 # tSPM+/MLHO wrappers, sequence <-> label decoding
│   ├── utils_pathway_analysis.R                # graph construction, path enumeration, path statistics
│   ├── utils_trajectory_analysis.R             # per-patient trajectory reconstruction & summaries
│   ├── utils_stats.R                           # generic LME / group-comparison helpers
│   └── utils_plotting.R                        # generic, group-agnostic plotting helpers
└── scripts/                                    # ordered, runnable pipeline steps
    ├── 01_temporal_sequence_mining.R
    ├── 02_iterative_feature_selection.R
    ├── 03_pathway_network_analysis.R
    ├── 04_patient_trajectory_analysis.R
    ├── 05_longitudinal_cognitive_modeling.R
    └── 06_biomarker_validation.R
```

## Requirements

R (>= 4.2) and the following packages:

`data.table, dplyr, tidyr, ggplot2, stringr, purrr, tibble, forcats, scales, lubridate,
igraph, ggraph, networkD3, TraMineR, patchwork, lme4, lmerTest, emmeans, broom,
broom.mixed, survival, survminer, arsenal, tSPMPlus, mlho`

`tSPMPlus` and `mlho` are not on CRAN; install from their respective GitHub repositories
with `devtools::install_github(...)`.

## How to use this repository

1. Edit `config.R` to point to your own event table, patient metadata/labels, and
   output directory, and to set analysis parameters (sparsity threshold, number of
   MLHO iterations, top-N features, etc.).
2. Source `R/00_setup.R` once per session.
3. Run the scripts in `scripts/` in numeric order. Each script reads the outputs of the
   previous step from `config$output_dir` and writes its own outputs (tables + figures)
   back to the same location.

## Data availability

No data are distributed with this repository. The pipeline was developed and applied to
a restricted-access longitudinal cohort; researchers with appropriate approvals for
comparable data sources can adapt these scripts by editing `config.R` and matching the
expected column names documented at the top of each script.

## Citation

If you use this pipeline, please cite the associated manuscript (details to be added
upon publication) as well as the upstream methods:
tSPM+ (temporal sequence pattern mining) and MLHO (minimal-learning-of-high-order
sequences).
