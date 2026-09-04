## =============================================================================
## utils_plotting.R
## Generic plotting helpers. Each function accepts an optional named vector
## of colors so a single implementation can be reused across every grouping
## variable (sex, genetic risk, trajectory group, comorbidity-burden group,
## ...) instead of duplicating the plotting code per group as in the
## original analysis scripts.
## =============================================================================

#' Horizontal lollipop plot of feature/pathway importance scores.
plot_feature_importance <- function(df, label_col, score_col, top_n = NULL, title = NULL) {
  if (!is.null(top_n)) df <- df[order(-df[[score_col]]), ][seq_len(min(top_n, nrow(df))), ]

  ggplot2::ggplot(df, ggplot2::aes(
    x = reorder(.data[[label_col]], .data[[score_col]]),
    y = .data[[score_col]]
  )) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data[[label_col]], y = 0, yend = .data[[score_col]]),
                          alpha = 0.5) +
    ggplot2::geom_point(color = "red", alpha = 0.6, size = 2) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Importance", title = title) +
    ggplot2::theme_minimal()
}

#' Boxplot of a continuous value across a grouping variable, with an optional
#' pairwise-comparison test annotation (via ggpubr).
plot_group_boxplot <- function(df, group_var, value_var, colors = NULL, title = NULL) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[group_var]], y = .data[[value_var]],
                                        fill = .data[[group_var]])) +
    ggplot2::geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 21, size = 3,
                          fill = "white", stroke = 1) +
    ggplot2::labs(title = title, x = group_var, y = value_var) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none",
                  axis.text.x = ggplot2::element_text(angle = 15, hjust = 1))

  if (!is.null(colors)) p <- p + ggplot2::scale_fill_manual(values = colors)
  p
}

#' Run a group comparison test and only render/save the boxplot if the result
#' is significant at `alpha`. Outliers are removed (IQR rule) before both the
#' test and the plot.
plot_if_significant <- function(df, group_var, value_var, title, save_path,
                                colors = NULL, alpha = 0.05) {
  df[[value_var]] <- remove_outliers_iqr(df[[value_var]])
  test_result <- test_by_group(df, group_var, value_var)
  p_value <- test_result$p.value[1]

  if (is.na(p_value) || p_value >= alpha) {
    message(sprintf("Not significant (p = %.4f) - skipping plot for %s", p_value, value_var))
    return(invisible(NULL))
  }

  message(sprintf("Significant (p = %.4f) - saving plot for %s", p_value, value_var))
  p <- plot_group_boxplot(df, group_var, value_var, colors = colors, title = title)
  ggplot2::ggsave(save_path, plot = p, width = 5, height = 5)
  p
}

#' Mean +/- SE trend of a continuous outcome over binned time, by group.
plot_binned_group_trend <- function(df, outcome_var, time_var, group_var,
                                    colors = NULL, bin_width = 5, title = NULL) {
  df_plot <- df %>%
    dplyr::filter(!is.na(.data[[group_var]])) %>%
    dplyr::mutate(time_bin = floor(.data[[time_var]] / bin_width) * bin_width) %>%
    dplyr::group_by(time_bin, .data[[group_var]]) %>%
    dplyr::summarise(
      mean_val = mean(.data[[outcome_var]], na.rm = TRUE),
      sd_val = stats::sd(.data[[outcome_var]], na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(se_val = sd_val / sqrt(n))

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = time_bin, y = mean_val, color = .data[[group_var]])) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = mean_val - se_val, ymax = mean_val + se_val,
                                      fill = .data[[group_var]]), alpha = 0.2, color = NA) +
    ggplot2::labs(y = outcome_var, title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  if (!is.null(colors)) {
    p <- p + ggplot2::scale_color_manual(values = colors) + ggplot2::scale_fill_manual(values = colors)
  }
  p
}

#' Bar chart of estimated per-year slopes with 95% CI error bars (for
#' visualizing fit_lme_group_trends() output).
plot_slope_estimates <- function(slopes_df, group_col, title = NULL) {
  df <- slopes_df %>%
    dplyr::rename(group = dplyr::all_of(group_col), lower = lower.CL, upper = upper.CL)

  ggplot2::ggplot(df, ggplot2::aes(x = reorder(group, slope), y = slope, fill = group)) +
    ggplot2::geom_col(width = 0.5, alpha = 0.8) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(x = NULL, y = "Estimated change per year", title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}

#' Generic tile heatmap, e.g. prevalence (%) of a condition/sequence across a
#' binned time axis, or a signed contrast between two subgroups.
plot_heatmap <- function(df, x_var, y_var, fill_var, title = NULL,
                         diverging = FALSE, fill_label = fill_var) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[fill_var]])) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::labs(x = x_var, y = y_var, fill = fill_label, title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (diverging) {
    p + ggplot2::scale_fill_gradient2(low = "#1a9850", mid = "white", high = "#d73027", midpoint = 0)
  } else {
    p + ggplot2::scale_fill_gradient(low = "#deebf7", high = "#08519c")
  }
}
