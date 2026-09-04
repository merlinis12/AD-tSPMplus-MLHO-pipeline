## =============================================================================
## utils_stats.R
## Generic, group-agnostic statistical helpers used by the longitudinal
## cognitive modeling and biomarker validation scripts. A single
## implementation replaces the many near-duplicate functions in the original
## analysis scripts (one per grouping variable).
## =============================================================================

#' Fit a linear mixed-effects model of a longitudinal outcome as a function
#' of time interacted with a grouping variable, with a per-patient random
#' intercept and slope, and return group-specific slope estimates.
#'
#' @param df long-format data with one row per repeated measurement
#' @param outcome_var name of the repeated-measures outcome (e.g., cognitive score)
#' @param group_var name of the (categorical) grouping variable of interest
#' @param time_var name of the time variable (e.g., years to/from index date)
#' @param id_var name of the subject identifier (random-effects grouping factor)
#' @param covariates optional character vector of additional fixed-effect covariates
fit_lme_group_trends <- function(df, outcome_var, group_var, time_var, id_var,
                                 covariates = character()) {

  df <- df %>%
    dplyr::filter(!is.na(.data[[group_var]]), !is.na(.data[[outcome_var]])) %>%
    dplyr::mutate(time_c = .data[[time_var]])

  extra_terms <- if (length(covariates) > 0) paste("+", paste(covariates, collapse = " + ")) else ""
  form <- stats::as.formula(
    paste0(outcome_var, " ~ time_c * ", group_var, extra_terms,
          " + (1 + time_c | ", id_var, ")")
  )

  model <- lme4::lmer(form, data = df, REML = FALSE,
                      control = lme4::lmerControl(optimizer = "bobyqa"))

  em_trends <- emmeans::emtrends(
    model, specs = stats::as.formula(paste("~", group_var)),
    var = "time_c", mode = "satterthwaite"
  )

  slopes <- as.data.frame(em_trends) %>% dplyr::rename(slope = time_c.trend)
  slope_pairs <- as.data.frame(emmeans::contrast(em_trends, method = "pairwise", adjust = "fdr"))
  interaction_test <- stats::anova(model)

  list(
    model = model,
    slopes = slopes,
    slope_pairs = slope_pairs,
    interaction_test = interaction_test,
    tidy_fixed = broom.mixed::tidy(model, effects = "fixed")
  )
}

#' Format group-level slope estimates into a compact publication-style table.
format_slope_table <- function(slopes_df, group_label) {
  slopes_df %>%
    dplyr::mutate(CI = sprintf("%.3f, %.3f", slope - 1.96 * SE, slope + 1.96 * SE)) %>%
    dplyr::transmute(
      !!group_label := .data[[names(slopes_df)[1]]],
      `Slope` = round(slope, 3),
      `SE` = round(SE, 3),
      `95% CI` = CI
    )
}

#' Format pairwise slope-contrast results into a compact publication-style table.
format_contrast_table <- function(slope_pairs_df) {
  slope_pairs_df %>%
    dplyr::transmute(
      Comparison = contrast,
      `Delta slope` = round(estimate, 3),
      SE = round(SE, 3),
      df = round(df, 1),
      `p (FDR)` = signif(p.value, 3)
    )
}

## ---------------------------------------------------------------------------
## Cross-sectional group comparisons (e.g., biomarker validation)
## ---------------------------------------------------------------------------

#' Winsorize-by-removal: replace values outside 1.5*IQR with NA.
remove_outliers_iqr <- function(x) {
  q1 <- stats::quantile(x, 0.25, na.rm = TRUE)
  q3 <- stats::quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  x[x < q1 - 1.5 * iqr | x > q3 + 1.5 * iqr] <- NA
  x
}

#' Descriptive summary of a continuous variable by group.
summarize_by_group <- function(df, group_var, value_var) {
  df %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean = mean(.data[[value_var]], na.rm = TRUE),
      sd = stats::sd(.data[[value_var]], na.rm = TRUE),
      median = stats::median(.data[[value_var]], na.rm = TRUE),
      IQR = stats::IQR(.data[[value_var]], na.rm = TRUE),
      .groups = "drop"
    )
}

#' t-test (2 groups) or one-way ANOVA (>2 groups) for a continuous variable
#' across a categorical grouping variable.
test_by_group <- function(df, group_var, value_var) {
  form <- stats::as.formula(paste(value_var, "~", group_var))
  if (dplyr::n_distinct(df[[group_var]]) == 2) {
    broom::tidy(stats::t.test(form, data = df))
  } else {
    broom::tidy(stats::aov(form, data = df))
  }
}
