## =============================================================================
## 05_longitudinal_cognitive_modeling.R
##
## Step 5 of the pipeline: model longitudinal decline in a repeated cognitive
## measure (e.g., MMSE) as a function of time-to-outcome interacted with
## trajectory-derived grouping variables (sex, genetic risk, sex x genetic
## risk, comorbidity-burden group, top trajectory group), using linear
## mixed-effects models with a per-patient random slope.
##
## Inputs:
##   - patient_sequences.csv (from 04_patient_trajectory_analysis.R)
##   - a repeated-measures cognitive-score table with columns:
##       id_col, cog_score, exam_date
##     (path set via config$cognitive_measure_path; not distributed here)
##
## Outputs:
##   - slope_table_<group>.csv / contrast_table_<group>.csv for each grouping
##     variable analyzed
##   - trajectory_plots.svg (combined mean +/- SE trend panels)
## =============================================================================

source("R/00_setup.R")
source("config.R")
source("R/utils_stats.R")
source("R/utils_plotting.R")

patient_sequences <- read.csv(file.path(config$output_dir, "patient_sequences.csv"))

## Repeated cognitive-measure table; replace with your own loader.
## Expected columns: id_col, cog_score, exam_date
cognitive_long <- get(load(config$cognitive_measure_path))

## ---- merge and derive time-to-outcome / grouping variables -----------------
df <- cognitive_long %>%
  dplyr::inner_join(patient_sequences, by = config$id_col) %>%
  dplyr::mutate(
    exam_date  = as.Date(exam_date),
    last_date  = as.Date(last_date),
    time_to_outcome = as.numeric(difftime(exam_date, last_date, units = "days")) / 365.25,
    outcome_group = ifelse(outcome == "Outcome", "Outcome", "Outcome-free")
  )

## Top-N most frequent trajectories among outcome-positive patients; all
## others grouped as "Other".
top_seqs <- df %>%
  dplyr::filter(outcome == "Outcome") %>%
  dplyr::distinct(.data[[config$id_col]], sequence) %>%
  dplyr::count(sequence, name = "n_pat") %>%
  dplyr::arrange(dplyr::desc(n_pat)) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::pull(sequence)

df <- df %>%
  dplyr::mutate(
    seq_group = ifelse(sequence %in% top_seqs, sequence, "Other"),
    n_events_group = dplyr::case_when(
      n_events <= 1 ~ "Short (<=1)",
      n_events <= 4 ~ "Medium (2-4)",
      n_events >  4 ~ "Long (>4)",
      TRUE ~ NA_character_
    ),
    n_events_group = factor(n_events_group, levels = c("Short (<=1)", "Medium (2-4)", "Long (>4)")),
    sex_genetic_risk = paste(sex, genetic_risk)
  )

## ---- fit mixed-effects models for each grouping variable of interest ------
grouping_vars <- c("outcome_group", "seq_group", "n_events_group",
                   "sex", "genetic_risk", "sex_genetic_risk")

results <- list()
for (gv in grouping_vars) {
  message("Fitting LME trend model for group: ", gv)
  fit <- fit_lme_group_trends(
    df, outcome_var = "cog_score", group_var = gv,
    time_var = "time_to_outcome", id_var = config$id_col
  )
  results[[gv]] <- fit

  write.csv(format_slope_table(fit$slopes, gv),
           file.path(config$output_dir, paste0("slope_table_", gv, ".csv")), row.names = FALSE)
  write.csv(format_contrast_table(fit$slope_pairs),
           file.path(config$output_dir, paste0("contrast_table_", gv, ".csv")), row.names = FALSE)
}

## ---- combined trend plots (outcome-positive patients only, as an example) -
group_colors <- list(
  seq_group      = NULL,   # supply a named color vector per level if desired
  n_events_group = c("Short (<=1)" = "#1f78b4", "Medium (2-4)" = "#b2df8a", "Long (>4)" = "#33a02c"),
  sex            = c("male" = "#4E589F", "female" = "#F5A4A7")
)

df_outcome_only <- dplyr::filter(df, outcome == "Outcome")
plots <- lapply(names(group_colors), function(gv) {
  plot_binned_group_trend(df_outcome_only, "cog_score", "time_to_outcome", gv,
                          colors = group_colors[[gv]], title = paste("By", gv))
})

combined <- patchwork::wrap_plots(plots, ncol = length(plots))
ggplot2::ggsave(file.path(config$output_dir, "trajectory_plots.svg"),
               combined, width = 4 * length(plots), height = 4)
