## =============================================================================
## config.R
## Central place for paths and parameters. Nothing in this repository should
## contain hard-coded local paths or dataset-specific file names outside of
## this file -- edit the values below to point to your own data.
## =============================================================================

config <- list(

  ## ---- input data -----------------------------------------------------------
  ## Long-format event table with (at minimum) the columns:
  ##   patient_num   : unique patient identifier
  ##   phenx         : event/phenotype code (diagnosis, medication, lab, etc.)
  ##   DESCRIPTION   : human-readable label for `phenx`
  ##   start_date    : event date (Date, or numeric days since an arbitrary origin)
  event_table_path      = "data/event_table.RData",

  ## Per-patient metadata: demographics, genetic risk factors, and the
  ## reference/index date used to time-align events (e.g., outcome exam date).
  ## Expected columns include (at minimum):
  ##   patient_num, outcome (0/1), outcome_exam_date, sex, apoe/genetic_risk
  patient_metadata_path = "data/patient_metadata.RData",

  ## Optional lookup describing event categories, used to keep only certain
  ## classes of events (e.g., comorbidities + medications) if desired.
  event_category_lookup_path = "data/event_category_lookup.RData",

  ## ---- outcome / column names -------------------------------------------
  outcome_col        = "alzheimers",   # binary outcome column in patient metadata
  outcome_date_col   = "ad_exam_date", # index/reference date column
  id_col             = "patient_num",

  ## ---- sequence mining (tSPM+) parameters --------------------------------
  sparsity            = 0.05,   # minimum support for a candidate sequence
  n_threads           = NULL,   # NULL -> autodetect available cores

  ## ---- feature selection (MLHO) parameters -------------------------------
  test_fraction       = 0.2,    # held-out fraction for train/test split
  top_n_features      = 200,    # candidate features passed to MLHO per iteration
  n_iterations        = 50,     # number of repeated MLHO/feature-selection runs
  classifier          = "gbm",  # any classifier supported by mlho::mlearn
  cv_folds            = 5,
  min_selection_freq  = 40,     # keep features selected in > this many iterations

  ## ---- pathway analysis parameters ---------------------------------------
  max_path_length      = 7,

  ## ---- longitudinal cognitive modeling -----------------------------------
  ## Repeated-measures cognitive score table: id_col, cog_score, exam_date
  cognitive_measure_path = "data/cognitive_measures.RData",

  ## ---- biomarker validation panels ----------------------------------------
  ## One named entry per independent biomarker source. Each entry needs:
  ##   path       : file path (id + sample_date + biomarker columns)
  ##   id_col     : name of the patient id column in that file
  ##   biomarkers : character vector of biomarker column names to test
  biomarker_panels = list(
    amyloid = list(
      path = "data/biomarkers_amyloid.RData",
      id_col = "id",
      biomarkers = c("amyloid40", "amyloid42")
    ),
    glial_injury = list(
      path = "data/biomarkers_glial.RData",
      id_col = "id",
      biomarkers = c("gfap", "nfl", "tau", "uchl1")
    ),
    inflammation = list(
      path = "data/biomarkers_inflammation.RData",
      id_col = "id",
      biomarkers = c("crp")
    ),
    ptau = list(
      path = "data/biomarkers_ptau.RData",
      id_col = "id",
      biomarkers = c("ptau_181")
    )
  ),

  ## ---- output -------------------------------------------------------------
  output_dir           = "outputs"
)

if (!dir.exists(config$output_dir)) dir.create(config$output_dir, recursive = TRUE)
