## =============================================================================
## 06_biomarker_validation.R
##
## Step 6 of the pipeline: internal/external validation. For each of several
## independent biomarker panels (e.g., amyloid, glial injury, inflammatory
## markers, pTau), test whether levels differ across the trajectory-derived
## grouping variables defined in 05_longitudinal_cognitive_modeling.R, in
## outcome-positive patients whose biomarker sample postdates their last
## captured trajectory event. Only statistically significant comparisons are
## plotted.
##
## Inputs:
##   - patient_sequences.csv, plus grouping columns computed as in script 05
##   - one or more biomarker panel tables, each with:
##       id_col, sample_date, <biomarker columns>
##     (paths supplied via config$biomarker_panels; not distributed here)
##
## Outputs, per (group, biomarker) pair:
##   - printed descriptive summary + test result
##   - a boxplot saved to disk, only if the group difference is significant
## =============================================================================

source("R/00_setup.R")
source("config.R")
source("R/utils_stats.R")
source("R/utils_plotting.R")

patient_sequences <- read.csv(file.path(config$output_dir, "patient_sequences.csv"))
grouping_vars <- c("sex", "genetic_risk", "sex_genetic_risk", "seq_group", "n_events_group")

group_colors <- list(
  n_events_group = c("Short (<=1)" = "#1f78b4", "Medium (2-4)" = "#b2df8a", "Long (>4)" = "#33a02c"),
  sex            = c("male" = "#4E589F", "female" = "#F5A4A7"),
  genetic_risk   = c("E4" = "#F9A022", "nonE4" = "#4193CF")
)

#' Load one biomarker panel, restrict it to outcome-positive patients whose
#' sample was collected on/after their last captured trajectory event, attach
#' the grouping variables, and run the group-comparison + conditional-plot
#' pipeline for each biomarker column in that panel.
#'
#' @param panel_path path to the biomarker panel file
#' @param biomarker_cols character vector of biomarker column names to test
#' @param id_col_panel name of the id column in the panel file (often
#'   different from the harmonized cohort id_col)
validate_biomarker_panel <- function(panel_path, biomarker_cols, id_col_panel,
                                     patient_sequences, last_event_dates,
                                     grouping_vars, group_colors, output_dir) {

  panel <- get(load(panel_path))
  panel <- dplyr::filter(panel, .data[[id_col_panel]] %in% patient_sequences[[config$id_col]])

  panel <- panel %>%
    dplyr::inner_join(patient_sequences, by = stats::setNames(config$id_col, id_col_panel)) %>%
    dplyr::inner_join(last_event_dates, by = stats::setNames(config$id_col, id_col_panel)) %>%
    dplyr::filter(outcome == "Outcome", sample_date >= last_event_date)

  for (gv in grouping_vars) {
    for (bio in biomarker_cols) {
      cat("\n== Group:", gv, "| Biomarker:", bio, "==\n")
      print(summarize_by_group(panel, gv, bio))
      print(test_by_group(panel, gv, bio))

      save_path <- file.path(output_dir, paste0("biomarker_", gv, "_", bio, ".svg"))
      plot_if_significant(
        panel, gv, bio, title = paste(bio, "by", gv),
        save_path = save_path, colors = group_colors[[gv]]
      )
    }
  }
}

## Last trajectory-event date per patient, used to ensure the biomarker
## sample was collected after the trajectory being characterized.
last_event_dates <- patient_sequences %>%
  dplyr::select(dplyr::all_of(config$id_col), last_event_date = last_date)

## ---- run validation for each configured biomarker panel --------------------
## config$biomarker_panels is a named list, one entry per panel:
##   list(path = "...", id_col = "...", biomarkers = c("...", "..."))
for (panel_name in names(config$biomarker_panels)) {
  panel_cfg <- config$biomarker_panels[[panel_name]]
  message("Validating biomarker panel: ", panel_name)

  validate_biomarker_panel(
    panel_path        = panel_cfg$path,
    biomarker_cols    = panel_cfg$biomarkers,
    id_col_panel      = panel_cfg$id_col,
    patient_sequences = patient_sequences,
    last_event_dates  = last_event_dates,
    grouping_vars     = grouping_vars,
    group_colors      = group_colors,
    output_dir        = config$output_dir
  )
}
