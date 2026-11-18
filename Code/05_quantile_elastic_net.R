# ============================================================================
# 05_quantile_elastic_net.R
# Derived from the uploaded Madhavi_panchangum_fred.R; original methods preserved.
# Run from the repository root via Rscript run_all.R.
# ============================================================================


# PART II. QUANTILE FORECASTING -- MAIN WINDOW ONLY (W = 180)
# ==============================================================================



quantile_forecast_rows <- list()
quantile_selection_rows <- list()
quantile_fit_rows <- list()

quantile_forecast_id <- 0L
quantile_selection_id <- 0L
quantile_fit_id <- 0L

quantile_tuning_cache <- list()


# ==============================================================================
# 21. QAR(1) BENCHMARK, WINDOW 180 ONLY
# ==============================================================================

for (h in forecast_horizons) {

  horizon_data <- master_panel |>
    dplyr::mutate(
      target_future = dplyr::lead(INDPRO_growth, h),
      target_date = dplyr::lead(date, h)
    )

  for (window_size in quantile_rolling_windows) {

    first_origin <- if (use_common_oos_start) {
      # Keep the same common OOS start used by the point exercise so that the
      # main 180-month quantile results are directly aligned with point W=180.
      max(point_rolling_windows) + max(forecast_horizons)
    } else {
      window_size + h
    }

    last_origin <- nrow(horizon_data) - h
    if (first_origin > last_origin) next

    for (forecast_origin in seq.int(first_origin, last_origin)) {

      sample_i <- get_rolling_sample(
        horizon_data, h, window_size, forecast_origin
      )
      if (is.null(sample_i)) next

      qar_train <- sample_i$training_df |>
        dplyr::filter(is.finite(AR1))

      if (nrow(qar_train) < 30L || !is.finite(sample_i$current_df$AR1)) next

      for (tau_value in quantile_levels) {

        qar_model <- tryCatch(
          quantreg::rq(
            target_future ~ AR1,
            tau = tau_value,
            data = qar_train,
            method = "fn"
          ),
          error = function(e) NULL
        )
        if (is.null(qar_model)) next

        qar_forecast <- tryCatch(
          as.numeric(
            stats::predict(
              qar_model,
              newdata = sample_i$current_df
            )
          ),
          error = function(e) NA_real_
        )
        if (!is.finite(qar_forecast)) next

        quantile_forecast_id <- quantile_forecast_id + 1L
        quantile_forecast_rows[[quantile_forecast_id]] <- tibble::tibble(
          window = window_size,
          horizon = h,
          forecast_origin = sample_i$origin_date,
          target_date = sample_i$target_date,
          model = "QAR(1)",
          tau = tau_value,
          actual = sample_i$actual,
          forecast = qar_forecast,
          error = sample_i$actual - qar_forecast,
          pinball_loss =
            pinball_loss(sample_i$actual, qar_forecast, tau_value),
          alpha = NA_real_,
          lambda = NA_real_
        )
      }
    }
  }
}


# ==============================================================================
# 22. QUANTILE ELASTIC NET Q-EN-1 ... Q-EN-5, WINDOW 180 ONLY
# ==============================================================================

run_quantile_en_model <- function(model_name) {

  predictors <- quantile_specs[[model_name]]

  cat("\nRunning quantile model:", model_name, "\n")

  successes <- 0L
  failures <- 0L

  for (h in forecast_horizons) {

    horizon_data <- master_panel |>
      dplyr::mutate(
        target_future = dplyr::lead(INDPRO_growth, h),
        target_date = dplyr::lead(date, h)
      )

    for (window_size in quantile_rolling_windows) {

      first_origin <- if (use_common_oos_start) {
        max(point_rolling_windows) + max(forecast_horizons)
      } else {
        window_size + h
      }

      last_origin <- nrow(horizon_data) - h
      if (first_origin > last_origin) next

      for (forecast_origin in seq.int(first_origin, last_origin)) {

        sample_i <- get_rolling_sample(
          horizon_data, h, window_size, forecast_origin
        )
        if (is.null(sample_i)) next

        training_df <- sample_i$training_df
        current_df <- sample_i$current_df

        for (tau_value in quantile_levels) {

          cache_key <- paste(
            window_size, h, model_name, tau_value, sep = "|"
          )
          cached <- quantile_tuning_cache[[cache_key]]

          need_retune <- is.null(cached) ||
            (forecast_origin - cached$forecast_origin >= retune_every)

          if (need_retune) {

            tuning <- tune_quantile_en(
              training_df = training_df,
              predictors = predictors,
              tau_value = tau_value
            )

            if (!is.null(tuning)) {
              quantile_tuning_cache[[cache_key]] <<- list(
                forecast_origin = forecast_origin,
                tuning = tuning
              )
            }

          } else {
            tuning <- cached$tuning
          }

          if (is.null(tuning)) {
            failures <- failures + 1L
            next
          }

          best_alpha <- tuning$alpha
          best_lambda <- tuning$lambda

          lambda_path <- as.numeric(tuning$lambda_path)
          lambda_path <- unique(lambda_path[is.finite(lambda_path)])

          if (length(lambda_path) < 2L) {
            failures <- failures + 1L
            next
          }

          matrices <- prepare_predictor_matrices(
            fit_df = training_df,
            apply_df = current_df,
            predictors = predictors
          )

          if (is.null(matrices)) {
            failures <- failures + 1L
            next
          }

          y_train <- training_df$target_future
          valid_y <- is.finite(y_train)

          x_train <- matrices$x_fit[valid_y, , drop = FALSE]
          y_train <- y_train[valid_y]
          x_current <- matrices$x_apply

          if (
            length(y_train) < 30L ||
            ncol(x_train) == 0L ||
            nrow(x_current) != 1L
          ) {
            failures <- failures + 1L
            next
          }

          final_q_model <- tryCatch(
            rqPen::rq.pen(
              x = x_train,
              y = y_train,
              tau = tau_value,
              penalty = "ENet",
              a = best_alpha,
              lambda = lambda_path,
              scalex = TRUE,
              lambda.discard = FALSE,
              alg = "huber"
            ),
            error = function(e) NULL
          )

          if (is.null(final_q_model)) {
            failures <- failures + 1L
            next
          }

          final_lambda_path <- as.numeric(final_q_model$lambda)
          final_lambda_path <- final_lambda_path[is.finite(final_lambda_path)]

          if (length(final_lambda_path) < 1L) {
            failures <- failures + 1L
            next
          }

          selected_index <- which.min(
            abs(final_lambda_path - best_lambda)
          )

          qen_pred <- tryCatch(
            stats::predict(
              final_q_model,
              newx = x_current,
              tau = tau_value,
              a = best_alpha,
              lambdaIndex = selected_index,
              sort = FALSE
            ),
            error = function(e) NULL
          )

          if (is.null(qen_pred)) {
            failures <- failures + 1L
            next
          }

          values <- as.numeric(qen_pred)
          values <- values[is.finite(values)]

          if (length(values) < 1L) {
            failures <- failures + 1L
            next
          }

          qen_forecast <- values[[1]]

          quantile_forecast_id <<- quantile_forecast_id + 1L
          quantile_forecast_rows[[quantile_forecast_id]] <<-
            tibble::tibble(
              window = window_size,
              horizon = h,
              forecast_origin = sample_i$origin_date,
              target_date = sample_i$target_date,
              model = model_name,
              tau = tau_value,
              actual = sample_i$actual,
              forecast = qen_forecast,
              error = sample_i$actual - qen_forecast,
              pinball_loss = pinball_loss(
                sample_i$actual,
                qen_forecast,
                tau_value
              ),
              alpha = best_alpha,
              lambda = final_lambda_path[selected_index]
            )

          quantile_fit_id <<- quantile_fit_id + 1L
          quantile_fit_rows[[quantile_fit_id]] <<-
            tibble::tibble(
              window = window_size,
              horizon = h,
              forecast_origin = sample_i$origin_date,
              model = model_name,
              tau = tau_value
            )

          successes <- successes + 1L

          qcoef <- extract_qen_coefficients_by_index(
            model = final_q_model,
            tau_value = tau_value,
            alpha_value = best_alpha,
            lambda_index = selected_index,
            predictor_names = colnames(x_train)
          )

          if (!is.null(qcoef)) {

            selected_q <- names(qcoef)[
              !is.na(names(qcoef)) &
                !(names(qcoef) %in% c("(Intercept)", "Intercept")) &
                is.finite(qcoef) &
                abs(qcoef) > selection_tolerance
            ]

            if (length(selected_q) > 0L) {
              quantile_selection_id <<- quantile_selection_id + 1L
              quantile_selection_rows[[quantile_selection_id]] <<-
                tibble::tibble(
                  window = window_size,
                  horizon = h,
                  forecast_origin = sample_i$origin_date,
                  model = model_name,
                  tau = tau_value,
                  variable = selected_q
                )
            }
          }
        }
      }
    }
  }

  cat(
    "Finished", model_name,
    "| successful forecasts:", successes,
    "| skipped/failed:", failures, "\n"
  )

  invisible(NULL)
}


purrr::walk(names(quantile_specs), run_quantile_en_model)

quantile_forecasts <- dplyr::bind_rows(quantile_forecast_rows) |>
  dplyr::arrange(window, horizon, tau, target_date, model)

quantile_fit_log <- dplyr::bind_rows(quantile_fit_rows)
quantile_selection_log <- dplyr::bind_rows(quantile_selection_rows)

if (nrow(quantile_forecasts) == 0L) {
  stop("No quantile forecasts were produced.")
}

stopifnot(all(quantile_forecasts$window == 180L))


# ==============================================================================
