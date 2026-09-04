## =============================================================================
## utils_trajectory_analysis.R
## Reconstruct ordered, per-patient clinical trajectories from a long event
## table and an outcome/index date, and summarize them.
## =============================================================================

#' Collapse a long event table (one row per event) into an ordered trajectory
#' per patient, appending a terminal "outcome / outcome-free" state at the
#' index date.
#'
#' @param events data.frame with id_col, DESCRIPTION, event_date
#' @param patient_index data.frame with id_col, index_date, outcome (0/1),
#'   plus any grouping covariates to carry through (e.g., sex, genetic_risk)
#' @param outcome_labels length-2 character vector: c(positive_label, negative_label)
build_patient_trajectories <- function(events, patient_index,
                                       id_col = "patient_num",
                                       event_date_col = "start_date",
                                       index_date_col = "index_date",
                                       outcome_col = "outcome",
                                       outcome_labels = c("outcome_positive", "outcome_negative"),
                                       covariate_cols = character()) {

  events <- events %>%
    dplyr::left_join(patient_index, by = id_col) %>%
    dplyr::filter(.data[[event_date_col]] < .data[[index_date_col]])

  terminal_events <- patient_index %>%
    dplyr::mutate(
      DESCRIPTION = ifelse(.data[[outcome_col]] == 1, outcome_labels[1], outcome_labels[2]),
      !!event_date_col := .data[[index_date_col]]
    )

  keep_cols <- c(id_col, "DESCRIPTION", event_date_col, outcome_col, covariate_cols)
  df_all <- dplyr::bind_rows(events[keep_cols], terminal_events[keep_cols]) %>%
    dplyr::arrange(.data[[id_col]], .data[[event_date_col]])

  df_all %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::arrange(.data[[event_date_col]], .by_group = TRUE) %>%
    dplyr::mutate(step = dplyr::row_number()) %>%
    dplyr::ungroup()
}

#' Collapse an ordered event table into one row per patient: the sequence of
#' unique event categories (in order of first occurrence), outcome, and
#' number of distinct non-terminal events ("comorbidity burden").
#'
#' @param df_ordered output of build_patient_trajectories()
summarize_patient_sequences <- function(df_ordered,
                                        id_col = "patient_num",
                                        event_date_col = "start_date",
                                        covariate_cols = character()) {

  df_clean <- df_ordered %>%
    dplyr::arrange(.data[[id_col]], .data[[event_date_col]]) %>%
    dplyr::group_by(.data[[id_col]], DESCRIPTION) %>%
    dplyr::slice_head(n = 1) %>%   # keep earliest occurrence of each category
    dplyr::ungroup()

  df_clean %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::arrange(.data[[event_date_col]], .by_group = TRUE) %>%
    dplyr::summarise(
      sequence   = paste(DESCRIPTION, collapse = " -> "),
      outcome    = dplyr::last(DESCRIPTION),
      n_events   = dplyr::n() - 1,
      first_date = min(.data[[event_date_col]]),
      last_date  = max(.data[[event_date_col]]),
      dplyr::across(dplyr::all_of(covariate_cols), dplyr::first),
      .groups = "drop"
    )
}

#' Build a Sankey diagram's link/node tables from an ordered, step-indexed
#' event table.
build_sankey_tables <- function(df_ordered, id_col = "patient_num") {
  df_transitions <- df_ordered %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::arrange(step, .by_group = TRUE) %>%
    dplyr::mutate(
      next_event = dplyr::lead(DESCRIPTION),
      next_step  = dplyr::lead(step)
    ) %>%
    dplyr::filter(!is.na(next_event)) %>%
    dplyr::ungroup()

  links <- df_transitions %>%
    dplyr::mutate(
      source = paste0("Step ", step, ": ", DESCRIPTION),
      target = paste0("Step ", next_step, ": ", next_event)
    ) %>%
    dplyr::count(source, target, name = "value")

  nodes <- data.frame(name = unique(c(links$source, links$target)))
  links$IDsource <- match(links$source, nodes$name) - 1
  links$IDtarget <- match(links$target, nodes$name) - 1

  list(links = links, nodes = nodes)
}
