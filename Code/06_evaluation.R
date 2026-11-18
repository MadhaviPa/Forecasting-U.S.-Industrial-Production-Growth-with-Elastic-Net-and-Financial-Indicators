# ============================================================================
# 06_evaluation.R
# Derived from the uploaded Madhavi_panchangum_fred.R; original methods preserved.
# Run from the repository root via Rscript run_all.R.
# ============================================================================


# 16. POINT FORECAST METRICS
# ==============================================================================

point_metrics <- point_forecasts |>
  dplyr::group_by(window, horizon, model) |>
  dplyr::summarise(
    n_forecasts = dplyr::n(),
    MSFE = mean(error^2, na.rm = TRUE),
    RMSFE = sqrt(MSFE),
    MAE = mean(abs(error), na.rm = TRUE),
    .groups = "drop"
  )

point_ar_metrics <- point_metrics |>
  dplyr::filter(model == "AR(1)") |>
  dplyr::select(window, horizon, AR1_MSFE = MSFE, AR1_RMSFE = RMSFE)

point_metrics <- point_metrics |>
  dplyr::left_join(point_ar_metrics, by = c("window", "horizon")) |>
  dplyr::mutate(
    relative_MSFE_to_AR1 = MSFE / AR1_MSFE,
    Theils_U_to_AR1 = RMSFE / AR1_RMSFE
  ) |>
  dplyr::arrange(window, horizon, relative_MSFE_to_AR1)


# ==============================================================================
# 17. POINT VARIABLE-SELECTION SUMMARY
# ==============================================================================

point_fit_counts <- point_fit_log |>
  dplyr::count(window, horizon, model, name = "number_of_successful_fits")

if (nrow(point_selection_log) == 0L) {

  point_selection_summary <- tibble::tibble(
    window = integer(),
    horizon = integer(),
    model = character(),
    variable = character(),
    times_selected = integer(),
    number_of_successful_fits = integer(),
    selection_percentage = numeric()
  )

} else {

  point_selection_summary <- point_selection_log |>
    dplyr::count(window, horizon, model, variable, name = "times_selected") |>
    dplyr::left_join(
      point_fit_counts,
      by = c("window", "horizon", "model")
    ) |>
    dplyr::mutate(
      selection_percentage =
        100 * times_selected / number_of_successful_fits,
      variable = strip_lag_suffix(variable)
    ) |>
    dplyr::arrange(window, horizon, model, dplyr::desc(selection_percentage))
}


# ==============================================================================
# 18. POINT DM TESTS: EACH EN MODEL AGAINST AR(1) ONLY
#     ALL ROLLING WINDOWS x ALL HORIZONS
# ==============================================================================

point_dm_rows <- list()
point_dm_id <- 0L

for (w in point_rolling_windows) {
  for (h in forecast_horizons) {
    for (m in names(point_specs)) {

      pair <- point_forecasts |>
        dplyr::filter(
          window == w,
          horizon == h,
          model %in% c("AR(1)", m)
        ) |>
        dplyr::select(target_date, actual, model, forecast) |>
        tidyr::pivot_wider(
          names_from = model,
          values_from = forecast
        )

      if (!all(c("AR(1)", m) %in% names(pair))) next

      pair <- pair |>
        dplyr::filter(
          is.finite(actual),
          is.finite(.data[["AR(1)"]]),
          is.finite(.data[[m]])
        )

      if (nrow(pair) == 0L) next

      test <- dm_test_benchmark_vs_model(
        benchmark_loss = (pair$actual - pair[["AR(1)"]])^2,
        model_loss = (pair$actual - pair[[m]])^2,
        horizon = h
      )

      point_dm_id <- point_dm_id + 1L
      point_dm_rows[[point_dm_id]] <- test |>
        dplyr::mutate(
          window = w,
          horizon = h,
          benchmark_model = "AR(1)",
          model = m,
          model_better_5pct =
            is.finite(p_value_model_better) & p_value_model_better < 0.05,
          model_better_10pct =
            is.finite(p_value_model_better) & p_value_model_better < 0.10,
          .before = 1
        )
    }
  }
}

point_dm_vs_ar1 <- dplyr::bind_rows(point_dm_rows) |>
  dplyr::group_by(window, horizon) |>
  dplyr::mutate(
    p_value_model_better_holm =
      stats::p.adjust(p_value_model_better, method = "holm"),
    model_better_holm_5pct =
      is.finite(p_value_model_better_holm) &
      p_value_model_better_holm < 0.05
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(window, horizon, p_value_model_better_holm)


# ==============================================================================
# 19. CLARK-WEST (2007) MSPE-ADJUSTED TEST
#     GENUINELY NESTED EN MODELS, WINDOW = 180, ALL HORIZONS
# ==============================================================================

point_nested_pairs <- build_nested_pairs(point_specs)

cat("\nNested EN pairs used by Clark-West MSPE-adjusted test:\n")
print(point_nested_pairs)

point_cw_rows <- list()
point_cw_id <- 0L

for (h in forecast_horizons) {

  for (j in seq_len(nrow(point_nested_pairs))) {

    small <- point_nested_pairs$smaller_model[j]
    large <- point_nested_pairs$larger_model[j]

    pair <- point_forecasts |>
      dplyr::filter(
        window == clark_west_window,
        horizon == h,
        model %in% c(small, large)
      ) |>
      dplyr::select(target_date, actual, model, forecast) |>
      tidyr::pivot_wider(
        names_from = model,
        values_from = forecast
      )

    if (!all(c(small, large) %in% names(pair))) next

    pair <- pair |>
      dplyr::filter(
        is.finite(actual),
        is.finite(.data[[small]]),
        is.finite(.data[[large]])
      )

    if (nrow(pair) == 0L) next

    test <- clark_west_mspe_adjusted(
      actual = pair$actual,
      forecast_small = pair[[small]],
      forecast_large = pair[[large]],
      horizon = h
    )

    point_cw_id <- point_cw_id + 1L
    point_cw_rows[[point_cw_id]] <- test |>
      dplyr::mutate(
        window = clark_west_window,
        horizon = h,
        smaller_model = small,
        larger_model = large,
        extra_predictors =
          point_nested_pairs$extra_predictors[j],
        larger_model_better_5pct = reject_5pct,
        larger_model_better_10pct = reject_10pct,
        .before = 1
      )
  }
}

point_clark_west <- dplyr::bind_rows(point_cw_rows) |>
  dplyr::arrange(
    horizon,
    smaller_model,
    larger_model
  )






# ==============================================================================
# 23. QUANTILE FORECAST METRICS, WINDOW 180
# ==============================================================================

quantile_metrics <- quantile_forecasts |>
  dplyr::group_by(window, horizon, tau, model) |>
  dplyr::summarise(
    n_forecasts = dplyr::n(),
    mean_pinball_loss = mean(pinball_loss, na.rm = TRUE),
    MAE = mean(abs(error), na.rm = TRUE),
    empirical_coverage = mean(actual <= forecast, na.rm = TRUE),
    calibration_error = empirical_coverage - tau,
    .groups = "drop"
  )

qar_metrics <- quantile_metrics |>
  dplyr::filter(model == "QAR(1)") |>
  dplyr::select(
    window,
    horizon,
    tau,
    QAR1_pinball = mean_pinball_loss
  )

quantile_metrics <- quantile_metrics |>
  dplyr::left_join(
    qar_metrics,
    by = c("window", "horizon", "tau")
  ) |>
  dplyr::mutate(
    relative_pinball_loss_to_QAR1 =
      mean_pinball_loss / QAR1_pinball
  ) |>
  dplyr::arrange(horizon, tau, relative_pinball_loss_to_QAR1)


# ==============================================================================
# 24. QUANTILE VARIABLE-SELECTION SUMMARY, WINDOW 180
# ==============================================================================

quantile_fit_counts <- quantile_fit_log |>
  dplyr::count(
    window, horizon, tau, model,
    name = "number_of_successful_fits"
  )

if (nrow(quantile_selection_log) == 0L) {

  quantile_selection_summary <- tibble::tibble(
    window = integer(),
    horizon = integer(),
    tau = numeric(),
    model = character(),
    variable = character(),
    times_selected = integer(),
    number_of_successful_fits = integer(),
    selection_percentage = numeric()
  )

} else {

  quantile_selection_summary <- quantile_selection_log |>
    dplyr::count(
      window, horizon, tau, model, variable,
      name = "times_selected"
    ) |>
    dplyr::left_join(
      quantile_fit_counts,
      by = c("window", "horizon", "tau", "model")
    ) |>
    dplyr::mutate(
      selection_percentage =
        100 * times_selected / number_of_successful_fits,
      variable = strip_lag_suffix(variable)
    ) |>
    dplyr::arrange(
      horizon, tau, model,
      dplyr::desc(selection_percentage)
    )
}


# ==============================================================================
# 25. QUANTILE DM TESTS:
#     EACH Q-EN MODEL AGAINST QAR(1) ONLY
#     WINDOW 180 x ALL HORIZONS x ALL TAUS
# ==============================================================================

quantile_dm_rows <- list()
quantile_dm_id <- 0L

for (h in forecast_horizons) {
  for (tau_value in quantile_levels) {
    for (m in names(quantile_specs)) {

      pair <- quantile_forecasts |>
        dplyr::filter(
          window == 180L,
          horizon == h,
          tau == tau_value,
          model %in% c("QAR(1)", m)
        ) |>
        dplyr::select(
          target_date, actual, model, forecast, pinball_loss
        ) |>
        dplyr::select(-forecast, -actual) |>
        tidyr::pivot_wider(
          names_from = model,
          values_from = pinball_loss
        )

      if (!all(c("QAR(1)", m) %in% names(pair))) next

      pair <- pair |>
        dplyr::filter(
          is.finite(.data[["QAR(1)"]]),
          is.finite(.data[[m]])
        )

      if (nrow(pair) == 0L) next

      test <- dm_test_benchmark_vs_model(
        benchmark_loss = pair[["QAR(1)"]],
        model_loss = pair[[m]],
        horizon = h
      )

      quantile_dm_id <- quantile_dm_id + 1L
      quantile_dm_rows[[quantile_dm_id]] <- test |>
        dplyr::mutate(
          window = 180L,
          horizon = h,
          tau = tau_value,
          benchmark_model = "QAR(1)",
          model = m,
          model_better_5pct =
            is.finite(p_value_model_better) &
            p_value_model_better < 0.05,
          model_better_10pct =
            is.finite(p_value_model_better) &
            p_value_model_better < 0.10,
          .before = 1
        )
    }
  }
}

quantile_dm_vs_qar1 <- dplyr::bind_rows(quantile_dm_rows) |>
  dplyr::group_by(horizon, tau) |>
  dplyr::mutate(
    p_value_model_better_holm =
      stats::p.adjust(p_value_model_better, method = "holm"),
    model_better_holm_5pct =
      is.finite(p_value_model_better_holm) &
      p_value_model_better_holm < 0.05
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(horizon, tau, p_value_model_better_holm)


# ==============================================================================
# 26. NESTED Q-EN COMPARISONS, WINDOW 180
# ==============================================================================



quantile_nested_pairs <- build_nested_pairs(quantile_specs)

quantile_nested_dm_rows <- list()
quantile_nested_dm_id <- 0L

for (h in forecast_horizons) {
  for (tau_value in quantile_levels) {
    for (j in seq_len(nrow(quantile_nested_pairs))) {

      small <- quantile_nested_pairs$smaller_model[j]
      large <- quantile_nested_pairs$larger_model[j]

      pair <- quantile_forecasts |>
        dplyr::filter(
          window == 180L,
          horizon == h,
          tau == tau_value,
          model %in% c(small, large)
        ) |>
        dplyr::select(target_date, model, pinball_loss) |>
        tidyr::pivot_wider(
          names_from = model,
          values_from = pinball_loss
        )

      if (!all(c(small, large) %in% names(pair))) next

      pair <- pair |>
        dplyr::filter(
          is.finite(.data[[small]]),
          is.finite(.data[[large]])
        )

      if (nrow(pair) == 0L) next

      # Smaller model is the benchmark; H1 says the larger nested model
      # has a lower expected pinball loss.
      test <- dm_test_benchmark_vs_model(
        benchmark_loss = pair[[small]],
        model_loss = pair[[large]],
        horizon = h
      )

      quantile_nested_dm_id <- quantile_nested_dm_id + 1L
      quantile_nested_dm_rows[[quantile_nested_dm_id]] <- test |>
        dplyr::mutate(
          window = 180L,
          horizon = h,
          tau = tau_value,
          smaller_model = small,
          larger_model = large,
          extra_predictors =
            quantile_nested_pairs$extra_predictors[j],
          larger_model_better_5pct =
            is.finite(p_value_model_better) &
            p_value_model_better < 0.05,
          .before = 1
        )
    }
  }
}

quantile_nested_model_tests <- dplyr::bind_rows(
  quantile_nested_dm_rows
) |>
  dplyr::group_by(horizon, tau) |>
  dplyr::mutate(
    p_value_holm =
      stats::p.adjust(p_value_model_better, method = "holm"),
    larger_model_better_holm_5pct =
      is.finite(p_value_holm) & p_value_holm < 0.05
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(horizon, tau, p_value_holm)


# ==============================================================================
# 27. QUANTILE-CROSSING DIAGNOSTICS, WINDOW 180
# ==============================================================================

quantile_crossing_diagnostics <- quantile_forecasts |>
  dplyr::select(
    window, horizon, target_date, model, tau, forecast
  ) |>
  tidyr::pivot_wider(
    names_from = tau,
    values_from = forecast,
    names_prefix = "tau_"
  )

tau_columns <- paste0("tau_", quantile_levels)

quantile_crossing_diagnostics <- quantile_crossing_diagnostics |>
  dplyr::rowwise() |>
  dplyr::mutate(
    quantile_crossing = {
      vals <- c_across(dplyr::all_of(tau_columns))
      vals <- vals[is.finite(vals)]
      length(vals) >= 2L && any(diff(vals) < 0)
    }
  ) |>
  dplyr::ungroup()

quantile_crossing_summary <- quantile_crossing_diagnostics |>
  dplyr::group_by(window, horizon, model) |>
  dplyr::summarise(
    n_complete_or_partial_quantile_curves = dplyr::n(),
    n_crossing = sum(quantile_crossing, na.rm = TRUE),
    pct_crossing = 100 * mean(quantile_crossing, na.rm = TRUE),
    .groups = "drop"
  )





print(point_metrics)
print(point_dm_vs_ar1)
print(point_clark_west)
print(quantile_metrics)
print(quantile_dm_vs_qar1)


# Export reproducible numerical outputs from the CURRENT run.
# Original detailed objects are kept within the R session.
dir.create("results/generated", recursive = TRUE, showWarnings = FALSE)
readr::write_csv(point_metrics, "results/generated/msfe_results.csv")
readr::write_csv(point_dm_vs_ar1, "results/generated/dm_results.csv")
readr::write_csv(quantile_metrics, "results/generated/quantile_results.csv")
readr::write_csv(point_selection_summary, "results/generated/variable_selection_frequency.csv")
