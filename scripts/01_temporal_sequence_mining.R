## =============================================================================
## 01_temporal_sequence_mining.R
##
## Step 1 of the pipeline: mine pairwise temporal event transitions from a
## long-format clinical event table using tSPM+.
##
## Expected inputs (see config.R):
##   - event table:      id_col, phenx, DESCRIPTION, start_date
##   - patient metadata: id_col, outcome_col, outcome_date_col, covariates
##
## Outputs (written to config$output_dir):
##   - numeric_event_table.rds   (dbMart list: numeric events + lookups)
##   - candidate_sequences.csv   (all non-sparse pairwise transitions found)
## =============================================================================

source("R/00_setup.R")
source("config.R")
source("R/utils_sequence_mining.R")

## ---- 1. load inputs --------------------------------------------------------
events   <- get(load(config$event_table_path))
patients <- get(load(config$patient_metadata_path))

## Optional: restrict to specific event categories (e.g., comorbidities and
## medications only), using a category lookup table with columns
## `phenx`, `Category`, `DESCRIPTION`.
if (!is.null(config$event_category_lookup_path) &&
    file.exists(config$event_category_lookup_path)) {
  category_lookup <- get(load(config$event_category_lookup_path))
  excluded_categories <- c("Cognitive", "Trend")  # edit as appropriate
  excluded_codes <- category_lookup$phenx[category_lookup$Category %in% excluded_categories]
  events <- dplyr::filter(events, !phenx %in% excluded_codes)
}

## Drop patients with missing values in required covariates (example: a
## genetic-risk covariate used later for stratified analyses).
complete_ids <- patients[[config$id_col]][stats::complete.cases(patients)]
events   <- dplyr::filter(events, .data[[config$id_col]] %in% complete_ids)
patients <- dplyr::filter(patients, .data[[config$id_col]] %in% complete_ids)

## ---- 2. transform to tSPM+'s numeric representation ------------------------
dbmart_numeric <- prepare_numeric_event_table(events, id_col = config$id_col)

## Re-key patient metadata and the event table to the numeric patient IDs
## used internally by tSPM+.
patients <- patients %>%
  dplyr::left_join(dbmart_numeric$patientLookUp, by = config$id_col) %>%
  dplyr::select(-dplyr::all_of(config$id_col)) %>%
  dplyr::rename(!!config$id_col := num_pat_num)

events <- events %>%
  dplyr::left_join(dbmart_numeric$patientLookUp, by = config$id_col) %>%
  dplyr::select(-dplyr::all_of(config$id_col)) %>%
  dplyr::rename(!!config$id_col := num_pat_num)

## ---- 3. mine candidate transitions -----------------------------------------
candidate_sequences <- mine_sequences(
  numeric_events = dbmart_numeric$dbMart,
  sparsity       = config$sparsity,
  output_dir     = config$output_dir,
  n_threads      = config$n_threads
)

## ---- 4. save --------------------------------------------------------------
saveRDS(dbmart_numeric, file.path(config$output_dir, "numeric_event_table.rds"))
saveRDS(patients,       file.path(config$output_dir, "patient_metadata_numeric_id.rds"))
saveRDS(events,         file.path(config$output_dir, "event_table_numeric_id.rds"))
write.csv(candidate_sequences,
         file.path(config$output_dir, "candidate_sequences.csv"), row.names = FALSE)
