## =============================================================================
## utils_sequence_mining.R
## Helper functions for temporal sequence mining (tSPM+) and translating the
## numeric transition codes it returns back into human-readable labels.
## =============================================================================

#' Prepare a long event table for tSPM+ / MLHO
#'
#' @param events data.frame with columns id_col, phenx, start_date
#' @param id_col name of the patient identifier column
#' @return list as produced by tSPMPlus::transformDbMartToNumeric(), containing
#'   a numeric-ID version of the event table plus lookup tables mapping
#'   numeric IDs back to the original patient and event codes.
prepare_numeric_event_table <- function(events, id_col = "patient_num") {
  stopifnot(id_col %in% names(events))
  tSPMPlus::transformDbMartToNumeric(events)
}

#' Extract non-sparse pairwise transitions with tSPM+
#'
#' @param numeric_events data.frame of numeric-coded events (dbMart$dbMart)
#' @param sparsity minimum support threshold
#' @param output_dir directory for any intermediate files written by tSPM+
#' @param n_threads number of parallel threads (defaults to all available cores)
#' @return data.frame of candidate transitions with duration information
mine_sequences <- function(numeric_events,
                           sparsity = 0.05,
                           output_dir = tempdir(),
                           n_threads = NULL,
                           store_intermediate = FALSE) {
  if (is.null(n_threads)) n_threads <- parallel::detectCores()

  sequences <- tSPMPlus::extractNonSparseSequences(
    as.data.frame(numeric_events),
    store_intermediate,
    output_dir,
    "sequence_mining",
    sparsity,
    n_threads
  )

  # tSPM+ encodes each transition's time span as "startDate|endDate"; split it
  # out into usable Date columns and a duration in days.
  sequences[c("dateA", "dateB")] <- do.call(
    rbind, strsplit(sequences$duration_str, "\\|")
  )
  sequences$duration <- as.numeric(as.Date(sequences$dateB) - as.Date(sequences$dateA))
  sequences
}

#' Translate a tSPM+ numeric transition code (or an already-decoded string)
#' back into a human-readable "EventA -> EventB" label.
#'
#' @param sequence_code character scalar; either a pure-numeric transition
#'   code produced by tSPM+, or an already human-readable string (returned
#'   unchanged in that case).
#' @param phenx_lookup lookup mapping numeric phenotype IDs to original codes
#'   (dbMart$phenxLookUp from prepare_numeric_event_table())
#' @param concept_lookup data.frame with columns phenx, DESCRIPTION mapping
#'   original event codes to human-readable descriptions
#' @param phenx_code_length fixed digit-width used by tSPM+ to encode each
#'   individual event code within a concatenated transition code
decode_transition <- function(sequence_code, phenx_lookup, concept_lookup,
                               phenx_code_length = 7) {
  if (grepl("\\D", sequence_code)) {
    return(sequence_code)  # already a decoded / non-numeric label
  }

  split_at <- nchar(sequence_code) - phenx_code_length
  code_a <- as.integer(substring(sequence_code, 1, split_at))
  code_b <- as.integer(substring(sequence_code, split_at + 1))

  codes <- data.frame(num_Phenx = c(code_a, code_b))
  codes <- dplyr::left_join(codes, phenx_lookup, by = "num_Phenx")
  codes <- dplyr::left_join(codes, concept_lookup, by = "phenx")

  paste(codes$DESCRIPTION[1], codes$DESCRIPTION[2], sep = "->")
}

#' Vectorized wrapper around decode_transition()
decode_transitions <- function(sequence_codes, phenx_lookup, concept_lookup,
                                phenx_code_length = 7) {
  vapply(
    sequence_codes,
    decode_transition,
    FUN.VALUE = character(1),
    phenx_lookup = phenx_lookup,
    concept_lookup = concept_lookup,
    phenx_code_length = phenx_code_length
  )
}

#' Run one iteration of MLHO feature selection (MSMR.lite + mlearn) on a
#' train/test split, and return the ranked feature-importance table.
#'
#' @param train_dat,test_dat long-format sequence tables (patient_num, phenx)
#'   with `row` and `value.var` columns already added
#' @param labels data.frame with id_col + a binary outcome column
#' @param covariates optional data.frame of additional predictors (e.g., sex,
#'   genetic risk) to include in the model, keyed by id_col
#' @param outcome_col name of the outcome column expected by mlho::mlearn
#'   (the "aoi" argument)
#' @param top_n number of candidate features passed to MSMR.lite
#' @param classifier classifier name passed to mlho::mlearn
#' @param cv_folds number of cross-validation folds
run_mlho_iteration <- function(train_dat, test_dat, labels,
                               covariates = NULL,
                               outcome_col = "outcome",
                               top_n = 200,
                               classifier = "gbm",
                               cv_folds = 5,
                               note = "mlho_iteration") {

  uniq_train_ids <- as.character(unique(train_dat$patient_num))
  train_sel <- MSMSR.lite(
    MLHO.dat = train_dat, patients = uniq_train_ids,
    sparsity = NA, labels = labels, topn = top_n, multicore = FALSE
  )

  test_sub <- subset(test_dat, test_dat$phenx %in% colnames(train_sel))
  data.table::setDT(test_sub)
  test_sub[, row := .I]
  test_sub$value.var <- 1
  uniq_test_ids <- as.character(unique(test_sub$patient_num))

  test_sel <- MSMSR.lite(
    MLHO.dat = test_sub, patients = uniq_test_ids,
    sparsity = NA, jmi = FALSE, labels = labels
  )

  model <- mlearn(
    train_sel, test_sel,
    dems = covariates,
    classifier = classifier,
    note = note,
    cv = "cv",
    nfold = cv_folds,
    aoi = outcome_col,
    multicore = FALSE
  )

  model
}
