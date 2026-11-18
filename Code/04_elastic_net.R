# ============================================================================
# 04_elastic_net.R
# Derived from the uploaded Madhavi_panchangum_fred.R; original methods preserved.
# Run from the repository root via Rscript run_all.R.
# ============================================================================


# 11. HELPER FUNCTIONS
# ==============================================================================

pinball_loss <- function(actual, forecast, tau) {
  error <- actual - forecast
  ifelse(error >= 0, tau * error, (tau - 1) * error)
}

strip_lag_suffix <- function(x) {
  sub("__L1$", "", x)
}

# Chronological validation rows from an outer training sample.
chronological_split_index <- function(n) {
  validation_n <- max(
    minimum_inner_validation,
    floor(inner_validation_share * n)
  )

  # Keep at least 36 observations in the inner fitting sample whenever possible.
  validation_n <- min(validation_n, max(1L, n - 36L))

  fit_end <- n - validation_n

  list(
    fit = seq_len(fit_end),
    validation = seq.int(fit_end + 1L, n)
  )
}

# ------------------------------------------------------------------------------
# Leakage-free Stock-Watson PCA-EM for FRED-MD predictors
# ------------------------------------------------------------------------------

# Fit PCA-EM using ONLY the supplied training matrix. Observed cells are never
# changed. Missing cells are initialized to the training-series mean, the panel
# is standardized, a rank-r PCA common component is estimated, and only the
# originally missing cells are updated from that common component. Means and
# standard deviations are recomputed after each update, as described by
# McCracken-Ng / Stock-Watson.
fit_stock_watson_pca_em <- function(x, r = pca_em_factors) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x[!is.finite(x)] <- NA_real_

  if (nrow(x) < 3L || ncol(x) < 2L) return(NULL)

  observed_n <- colSums(is.finite(x))
  observed_sd <- apply(
    x, 2,
    function(z) {
      z <- z[is.finite(z)]
      if (length(z) < 2L) NA_real_ else stats::sd(z)
    }
  )

  usable <- observed_n >= pca_em_min_observed &
    is.finite(observed_sd) & observed_sd > 0

  x <- x[, usable, drop = FALSE]
  if (ncol(x) < 2L) return(NULL)

  missing_mask <- !is.finite(x)
  initial_means <- colMeans(x, na.rm = TRUE)
  completed <- x
  for (j in seq_len(ncol(completed))) {
    completed[missing_mask[, j], j] <- initial_means[j]
  }

  previous_missing <- completed[missing_mask]
  final_v <- NULL
  final_means <- NULL
  final_sds <- NULL
  converged <- FALSE

  for (iter in seq_len(pca_em_max_iter)) {
    means <- colMeans(completed)
    sds <- apply(completed, 2, stats::sd)
    good <- is.finite(sds) & sds > 0
    if (!all(good)) {
      # This is uncommon after the initial usability screen; dropping here is
      # safer than manufacturing variation.
      completed <- completed[, good, drop = FALSE]
      x <- x[, good, drop = FALSE]
      missing_mask <- missing_mask[, good, drop = FALSE]
      if (ncol(completed) < 2L) return(NULL)
      means <- colMeans(completed)
      sds <- apply(completed, 2, stats::sd)
    }

    z <- sweep(completed, 2, means, "-")
    z <- sweep(z, 2, sds, "/")

    r_eff <- min(
      as.integer(r),
      nrow(z) - 1L,
      ncol(z) - 1L
    )
    if (!is.finite(r_eff) || r_eff < 1L) return(NULL)

    sv <- tryCatch(
      base::svd(z, nu = r_eff, nv = r_eff),
      error = function(e) NULL
    )
    if (is.null(sv)) return(NULL)

    v_r <- sv$v[, seq_len(r_eff), drop = FALSE]
    scores <- z %*% v_r
    z_common <- scores %*% t(v_r)
    reconstructed <- sweep(z_common, 2, sds, "*")
    reconstructed <- sweep(reconstructed, 2, means, "+")

    # EM update: only missing training cells are replaced. Observed values stay
    # exactly equal to the data available at this forecast origin.
    completed[missing_mask] <- reconstructed[missing_mask]

    current_missing <- completed[missing_mask]
    change <- if (length(current_missing) == 0L) {
      0
    } else {
      max(abs(current_missing - previous_missing), na.rm = TRUE)
    }

    previous_missing <- current_missing
    final_v <- v_r
    final_means <- colMeans(completed)
    final_sds <- apply(completed, 2, stats::sd)

    if (is.finite(change) && change < pca_em_tolerance) {
      converged <- TRUE
      break
    }
  }

  # Re-estimate final loadings from the final completed and standardized
  # TRAINING panel only. These are frozen when validation/current rows are
  # imputed below.
  z_final <- sweep(completed, 2, final_means, "-")
  z_final <- sweep(z_final, 2, final_sds, "/")
  r_eff <- min(as.integer(r), nrow(z_final) - 1L, ncol(z_final) - 1L)
  sv_final <- tryCatch(
    base::svd(z_final, nu = r_eff, nv = r_eff),
    error = function(e) NULL
  )
  if (is.null(sv_final)) return(NULL)

  list(
    completed = completed,
    variables = colnames(completed),
    means = final_means,
    sds = final_sds,
    loadings = sv_final$v[, seq_len(r_eff), drop = FALSE],
    r = r_eff,
    iterations = iter,
    converged = converged,
    n_imputed_training = sum(missing_mask)
  )
}

# Impute validation/current rows using a PCA-EM model estimated on the training
# sample. Crucially, apply rows NEVER update the training means, SDs or loadings.
# A factor score for each apply row is estimated from that row's observed FRED-MD
# entries and the fixed training loadings; missing entries are then reconstructed.
apply_stock_watson_pca_em <- function(model, new_x) {
  new_x <- as.matrix(new_x[, model$variables, drop = FALSE])
  storage.mode(new_x) <- "double"
  new_x[!is.finite(new_x)] <- NA_real_

  out <- new_x
  L <- model$loadings

  for (i in seq_len(nrow(new_x))) {
    row <- new_x[i, ]
    observed <- is.finite(row)
    missing <- !observed

    if (!any(missing)) next

    # Standardize only with training-sample parameters.
    z_obs <- (row[observed] - model$means[observed]) / model$sds[observed]
    L_obs <- L[observed, , drop = FALSE]

    if (sum(observed) >= model$r && nrow(L_obs) > 0L) {
      gram <- crossprod(L_obs)
      ridge <- diag(1e-8, nrow = ncol(gram))
      score <- tryCatch(
        solve(gram + ridge, crossprod(L_obs, z_obs)),
        error = function(e) NULL
      )
    } else {
      score <- NULL
    }

    if (is.null(score) || any(!is.finite(score))) {
      # If too little information exists in this row, the leakage-free fallback
      # is the TRAINING mean, not a value computed from validation/future rows.
      out[i, missing] <- model$means[missing]
    } else {
      z_missing <- as.numeric(L[missing, , drop = FALSE] %*% score)
      out[i, missing] <- model$means[missing] + model$sds[missing] * z_missing
    }
  }

  out
}

# Prepare model matrices without data leakage.
# FRED-MD predictors are completed by training-only Stock-Watson PCA-EM.
# Non-FRED predictors (the separately constructed S&P/Treasury variables) use
# training-sample medians only if a residual missing value remains.
prepare_predictor_matrices <- function(fit_df, apply_df, predictors) {
  predictors <- intersect(predictors, names(fit_df))
  if (length(predictors) == 0L) return(NULL)

  x_fit_all <- fit_df |>
    dplyr::select(dplyr::all_of(predictors)) |>
    as.data.frame()
  x_apply_all <- apply_df |>
    dplyr::select(dplyr::all_of(predictors)) |>
    as.data.frame()

  # fred_lagged is defined from the transformed FRED-MD predictor set before
  # forecasting starts. Restrict PCA-EM to that block only.
  fred_cols <- intersect(predictors, fred_lagged)
  other_cols <- setdiff(predictors, fred_cols)

  completed_fit <- list()
  completed_apply <- list()

  if (length(fred_cols) >= 2L) {
    fred_fit <- x_fit_all[, fred_cols, drop = FALSE]
    em_model <- fit_stock_watson_pca_em(fred_fit)

    if (!is.null(em_model)) {
      fred_apply <- x_apply_all[, em_model$variables, drop = FALSE]
      completed_fit[["fred"]] <- as.data.frame(em_model$completed)
      completed_apply[["fred"]] <- as.data.frame(
        apply_stock_watson_pca_em(em_model, fred_apply)
      )
    }
  }

  # Any FRED columns not usable by PCA-EM, and all non-FRED predictors, are
  # handled from the fitting sample alone. This preserves predictors where
  # possible without importing information from apply/future observations.
  used_fred <- if ("fred" %in% names(completed_fit)) {
    names(completed_fit[["fred"]])
  } else character(0)

  fallback_cols <- unique(c(setdiff(fred_cols, used_fred), other_cols))

  if (length(fallback_cols) > 0L) {
    fb_fit <- x_fit_all[, fallback_cols, drop = FALSE]
    fb_apply <- x_apply_all[, fallback_cols, drop = FALSE]

    medians <- vapply(
      fb_fit,
      function(x) {
        z <- as.numeric(x)
        z <- z[is.finite(z)]
        if (length(z) == 0L) NA_real_ else stats::median(z)
      },
      numeric(1)
    )
    usable_fb <- is.finite(medians)
    fb_fit <- fb_fit[, usable_fb, drop = FALSE]
    fb_apply <- fb_apply[, names(fb_fit), drop = FALSE]
    medians <- medians[names(fb_fit)]

    for (v in names(fb_fit)) {
      a <- as.numeric(fb_fit[[v]])
      b <- as.numeric(fb_apply[[v]])
      a[!is.finite(a)] <- medians[[v]]
      b[!is.finite(b)] <- medians[[v]]
      fb_fit[[v]] <- a
      fb_apply[[v]] <- b
    }

    completed_fit[["fallback"]] <- fb_fit
    completed_apply[["fallback"]] <- fb_apply
  }

  if (length(completed_fit) == 0L) return(NULL)

  x_fit <- dplyr::bind_cols(!!!unname(completed_fit)) |> as.data.frame()
  x_apply <- dplyr::bind_cols(!!!unname(completed_apply)) |> as.data.frame()

  # Restore the requested predictor order where possible.
  keep_order <- predictors[predictors %in% names(x_fit)]
  x_fit <- x_fit[, keep_order, drop = FALSE]
  x_apply <- x_apply[, keep_order, drop = FALSE]

  valid_variance <- vapply(
    x_fit,
    function(x) {
      x <- as.numeric(x)
      all(is.finite(x)) && is.finite(stats::sd(x)) && stats::sd(x) > 0
    },
    logical(1)
  )

  x_fit <- x_fit[, valid_variance, drop = FALSE]
  x_apply <- x_apply[, names(x_fit), drop = FALSE]
  if (ncol(x_fit) == 0L) return(NULL)

  list(
    x_fit = as.matrix(x_fit),
    x_apply = as.matrix(x_apply),
    variables = colnames(x_fit)
  )
}

# Extract a coefficient vector robustly from rqPen output.
extract_rqpen_coefficients <- function(model, tau_value) {
  co <- tryCatch(
    stats::coef(model, tau = tau_value, cvmin = TRUE),
    error = function(e) NULL
  )

  if (is.null(co)) {
    return(numeric(0))
  }

  if (is.vector(co)) {
    out <- as.numeric(co)
    names(out) <- names(co)
    return(out)
  }

  mat <- as.matrix(co)
  out <- as.numeric(mat[, 1])
  names(out) <- rownames(mat)
  out
}

# ------------------------------------------------------------------------------
# Point Elastic Net tuning: chronological inner holdout.
# ------------------------------------------------------------------------------

tune_point_en <- function(training_df, predictors) {
  n <- nrow(training_df)
  split <- chronological_split_index(n)

  fit_df <- training_df[split$fit, , drop = FALSE]
  val_df <- training_df[split$validation, , drop = FALSE]

  matrices <- prepare_predictor_matrices(
    fit_df = fit_df,
    apply_df = val_df,
    predictors = predictors
  )

  if (is.null(matrices)) {
    return(NULL)
  }

  y_fit <- fit_df$target_future
  y_val <- val_df$target_future

  valid_y_fit <- is.finite(y_fit)
  valid_y_val <- is.finite(y_val)

  x_fit <- matrices$x_fit[valid_y_fit, , drop = FALSE]
  y_fit <- y_fit[valid_y_fit]

  x_val <- matrices$x_apply[valid_y_val, , drop = FALSE]
  y_val <- y_val[valid_y_val]

  if (length(y_fit) < 30L || length(y_val) < 5L || ncol(x_fit) == 0L) {
    return(NULL)
  }

  tuning_rows <- list()
  id <- 0L

  for (a in alpha_grid) {
    fit <- tryCatch(
      glmnet::glmnet(
        x = x_fit,
        y = y_fit,
        family = "gaussian",
        alpha = a,
        standardize = TRUE,
        intercept = TRUE,
        nlambda = point_nlambda
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) next

    preds <- tryCatch(
      predict(fit, newx = x_val),
      error = function(e) NULL
    )

    if (is.null(preds)) next

    lambda_values <- fit$lambda
    k <- min(ncol(preds), length(lambda_values))

    if (k == 0L) next

    mse <- vapply(
      seq_len(k),
      function(j) mean((y_val - as.numeric(preds[, j]))^2, na.rm = TRUE),
      numeric(1)
    )

    id <- id + 1L
    tuning_rows[[id]] <- tibble::tibble(
      alpha = a,
      lambda = lambda_values[seq_len(k)],
      validation_loss = mse
    )
  }

  tuning_table <- dplyr::bind_rows(tuning_rows)

  if (nrow(tuning_table) == 0L || !any(is.finite(tuning_table$validation_loss))) {
    return(NULL)
  }

  tuning_table |>
    dplyr::filter(is.finite(validation_loss)) |>
    dplyr::slice_min(validation_loss, n = 1, with_ties = FALSE)
}

# ------------------------------------------------------------------------------
# Quantile Elastic Net tuning: chronological inner holdout.
# ------------------------------------------------------------------------------

tune_quantile_en <- function(training_df, predictors, tau_value) {
  n <- nrow(training_df)
  split <- chronological_split_index(n)

  fit_df <- training_df[split$fit, , drop = FALSE]
  val_df <- training_df[split$validation, , drop = FALSE]

  matrices <- prepare_predictor_matrices(
    fit_df = fit_df,
    apply_df = val_df,
    predictors = predictors
  )
  if (is.null(matrices)) return(NULL)

  y_fit <- fit_df$target_future
  y_val <- val_df$target_future

  valid_y_fit <- is.finite(y_fit)
  valid_y_val <- is.finite(y_val)

  x_fit <- matrices$x_fit[valid_y_fit, , drop = FALSE]
  y_fit <- y_fit[valid_y_fit]

  x_val <- matrices$x_apply[valid_y_val, , drop = FALSE]
  y_val <- y_val[valid_y_val]

  if (length(y_fit) < 30L || length(y_val) < 5L || ncol(x_fit) == 0L) {
    return(NULL)
  }

  best_result <- NULL
  best_loss <- Inf

  for (a_value in alpha_grid) {

    qfit <- tryCatch(
      rqPen::rq.pen(
        x = x_fit,
        y = y_fit,
        tau = tau_value,
        penalty = "ENet",
        a = a_value,
        nlambda = quantile_nlambda,
        scalex = TRUE,
        lambda.discard = FALSE,
        alg = "huber"
      ),
      error = function(e) {
        message(
          "Q-EN tuning fit failed: tau=", tau_value,
          ", alpha=", a_value,
          " | ", conditionMessage(e)
        )
        NULL
      }
    )

    if (is.null(qfit)) next

    lambda_path <- as.numeric(qfit$lambda)
    lambda_path <- lambda_path[is.finite(lambda_path)]

    if (length(lambda_path) < 2L) {
      message(
        "Q-EN tuning skipped because rqPen returned fewer than 2 lambdas: ",
        "tau=", tau_value, ", alpha=", a_value
      )
      next
    }

    # Use rqPen's documented lambdaIndex interface instead of matching
    # floating-point lambda values. One prediction column is requested for
    # each fitted lambda on the path.
    pred_path <- tryCatch(
      stats::predict(
        qfit,
        newx = x_val,
        tau = tau_value,
        a = a_value,
        lambdaIndex = seq_along(lambda_path),
        sort = FALSE
      ),
      error = function(e) {
        message(
          "Q-EN tuning prediction failed: tau=", tau_value,
          ", alpha=", a_value,
          " | ", conditionMessage(e)
        )
        NULL
      }
    )

    if (is.null(pred_path)) next

    pred_path <- as.matrix(pred_path)

    # rqPen should return validation observations by lambda. If a package
    # version returns the transpose, correct that deterministically.
    if (nrow(pred_path) != length(y_val) &&
        ncol(pred_path) == length(y_val)) {
      pred_path <- t(pred_path)
    }

    if (nrow(pred_path) != length(y_val)) {
      message(
        "Q-EN tuning prediction dimension mismatch: tau=", tau_value,
        ", alpha=", a_value,
        ", expected rows=", length(y_val),
        ", received=", nrow(pred_path)
      )
      next
    }

    k <- min(ncol(pred_path), length(lambda_path))
    if (k < 1L) next

    pred_path <- pred_path[, seq_len(k), drop = FALSE]
    lambda_used <- lambda_path[seq_len(k)]

    losses <- vapply(
      seq_len(k),
      function(j) {
        pred_j <- as.numeric(pred_path[, j])
        if (!all(is.finite(pred_j))) return(Inf)

        loss_j <- mean(
          pinball_loss(
            actual = y_val,
            forecast = pred_j,
            tau = tau_value
          ),
          na.rm = TRUE
        )

        if (is.finite(loss_j)) loss_j else Inf
      },
      numeric(1)
    )

    if (!any(is.finite(losses))) next

    idx <- which.min(losses)

    if (is.finite(losses[idx]) && losses[idx] < best_loss) {
      best_loss <- losses[idx]

      best_result <- list(
        alpha = as.numeric(a_value),
        lambda = as.numeric(lambda_used[idx]),
        lambda_index = as.integer(idx),
        lambda_path = as.numeric(lambda_used),
        validation_loss = as.numeric(losses[idx])
      )
    }
  }

  best_result
}

# Robust coefficient extraction using the documented rqPen lambdaIndex API.
extract_qen_coefficients_by_index <- function(model, tau_value, alpha_value,
                                              lambda_index, predictor_names) {
  co <- tryCatch(
    stats::coef(
      model,
      tau = tau_value,
      a = alpha_value,
      lambdaIndex = lambda_index
    ),
    error = function(e) {
      message(
        "Q-EN coefficient extraction failed: tau=", tau_value,
        ", alpha=", alpha_value,
        ", lambdaIndex=", lambda_index,
        " | ", conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(co)) return(NULL)

  # coef.rq.pen.seq can return a list containing the requested matrix.
  while (is.list(co) && length(co) == 1L) {
    co <- co[[1]]
  }

  if (is.list(co)) {
    # If more than one component remains, use the first non-empty numeric/matrix
    # component because this call requests exactly one tau, a, and lambdaIndex.
    candidate <- NULL
    for (item in co) {
      if (is.numeric(item) || is.matrix(item) || inherits(item, "Matrix")) {
        candidate <- item
        break
      }
    }
    co <- candidate
  }

  if (is.null(co)) return(NULL)

  if (is.vector(co)) {
    values <- as.numeric(co)
    nm <- names(co)
  } else {
    mat <- as.matrix(co)
    if (ncol(mat) < 1L) return(NULL)
    values <- as.numeric(mat[, 1])
    nm <- rownames(mat)
  }

  expected_n <- length(predictor_names) + 1L
  if (length(values) != expected_n) {
    message(
      "Q-EN coefficient length mismatch: expected ", expected_n,
      ", received ", length(values)
    )
    return(NULL)
  }

  if (is.null(nm) || length(nm) != expected_n) {
    nm <- c("(Intercept)", predictor_names)
  }

  names(values) <- nm
  values
}


# ==============================================================================
# 12. FORECAST / TEST SETTINGS
# ==============================================================================

# POINT robustness: all requested rolling windows.
point_rolling_windows <- rolling_windows

# QUANTILES: main specification only, exactly as requested.
quantile_rolling_windows <- 180L

# Clark-West is only used for nested POINT models and only for W = 180.
clark_west_window <- 180L

# Minimum number of matched OOS observations required by the forecast tests.
minimum_test_observations <- 20L


# ==============================================================================
# 13. COMMON FORECAST AND TEST HELPERS
# ==============================================================================

get_rolling_sample <- function(horizon_data, h, window_size, forecast_origin) {

  actual <- horizon_data$target_future[forecast_origin]
  target_date <- horizon_data$target_date[forecast_origin]
  origin_date <- horizon_data$date[forecast_origin]

  if (!is.finite(actual) || is.na(target_date)) return(NULL)

  # Direct h-step forecasting:
  # the last training origin s must satisfy s + h <= current forecast origin t.
  last_training_origin <- forecast_origin - h
  first_training_origin <- last_training_origin - window_size + 1L

  if (first_training_origin < 1L) return(NULL)

  training_rows <- seq.int(first_training_origin, last_training_origin)

  training_df <- horizon_data[
    training_rows,
    ,
    drop = FALSE
  ] |>
    dplyr::filter(is.finite(target_future))

  if (nrow(training_df) < floor(0.80 * window_size)) return(NULL)

  current_df <- horizon_data[
    forecast_origin,
    ,
    drop = FALSE
  ]

  list(
    actual = actual,
    target_date = target_date,
    origin_date = origin_date,
    training_df = training_df,
    current_df = current_df
  )
}


# Newey-West/Bartlett long-run variance of the sample mean.
hac_mean_test <- function(d, horizon, alternative = c("two.sided", "greater", "less")) {

  alternative <- match.arg(alternative)

  d <- as.numeric(d)
  d <- d[is.finite(d)]
  n <- length(d)

  if (n < minimum_test_observations) {
    return(tibble::tibble(
      n = n,
      statistic = NA_real_,
      p_value = NA_real_,
      mean_difference = if (n > 0L) mean(d) else NA_real_,
      hac_lag = NA_integer_,
      test_status = "Too few matched OOS observations"
    ))
  }

  # h-step direct forecasts can induce serial correlation through h-1.
  hac_lag <- min(max(0L, as.integer(horizon) - 1L), n - 1L)

  d_bar <- mean(d)
  z <- d - d_bar

  lrv <- sum(z^2) / n

  if (hac_lag > 0L) {
    for (j in seq_len(hac_lag)) {
      gamma_j <- sum(z[(j + 1L):n] * z[1L:(n - j)]) / n
      weight_j <- 1 - j / (hac_lag + 1L)
      lrv <- lrv + 2 * weight_j * gamma_j
    }
  }

  if (!is.finite(lrv) || lrv <= 0) {
    return(tibble::tibble(
      n = n,
      statistic = NA_real_,
      p_value = NA_real_,
      mean_difference = d_bar,
      hac_lag = hac_lag,
      test_status = "Non-positive HAC variance"
    ))
  }

  statistic <- d_bar / sqrt(lrv / n)

  p_value <- switch(
    alternative,
    "two.sided" = 2 * stats::pnorm(-abs(statistic)),
    "greater" = stats::pnorm(statistic, lower.tail = FALSE),
    "less" = stats::pnorm(statistic)
  )

  tibble::tibble(
    n = n,
    statistic = statistic,
    p_value = p_value,
    mean_difference = d_bar,
    hac_lag = hac_lag,
    test_status = "OK"
  )
}


# ------------------------------------------------------------------------------
# Diebold-Mariano test used in this file
# ------------------------------------------------------------------------------

# Convention:
#   d_t = loss_benchmark - loss_competitor
# so a POSITIVE mean difference means the competitor has smaller loss.
#
# We use the Harvey-Leybourne-Newbold small-sample correction and an h-1
# Bartlett HAC bandwidth.
dm_test_benchmark_vs_model <- function(
    benchmark_loss,
    model_loss,
    horizon
) {

  ok <- is.finite(benchmark_loss) & is.finite(model_loss)
  l0 <- as.numeric(benchmark_loss[ok])
  l1 <- as.numeric(model_loss[ok])
  n <- length(l0)

  if (n < minimum_test_observations) {
    return(tibble::tibble(
      n = n,
      mean_loss_benchmark = if (n > 0L) mean(l0) else NA_real_,
      mean_loss_model = if (n > 0L) mean(l1) else NA_real_,
      mean_loss_difference = if (n > 0L) mean(l0 - l1) else NA_real_,
      dm_statistic = NA_real_,
      hln_statistic = NA_real_,
      p_value_two_sided = NA_real_,
      p_value_model_better = NA_real_,
      hac_lag = NA_integer_,
      test_status = "Too few matched OOS observations"
    ))
  }

  d <- l0 - l1
  d_bar <- mean(d)
  z <- d - d_bar
  hac_lag <- min(max(0L, as.integer(horizon) - 1L), n - 1L)

  lrv <- sum(z^2) / n

  if (hac_lag > 0L) {
    for (j in seq_len(hac_lag)) {
      gamma_j <- sum(z[(j + 1L):n] * z[1L:(n - j)]) / n
      weight_j <- 1 - j / (hac_lag + 1L)
      lrv <- lrv + 2 * weight_j * gamma_j
    }
  }

  if (!is.finite(lrv) || lrv <= 0) {
    return(tibble::tibble(
      n = n,
      mean_loss_benchmark = mean(l0),
      mean_loss_model = mean(l1),
      mean_loss_difference = d_bar,
      dm_statistic = NA_real_,
      hln_statistic = NA_real_,
      p_value_two_sided = NA_real_,
      p_value_model_better = NA_real_,
      hac_lag = hac_lag,
      test_status = "Non-positive HAC variance"
    ))
  }

  dm_stat <- d_bar / sqrt(lrv / n)
  h <- as.integer(horizon)

  hln_inside <- (n + 1 - 2 * h + h * (h - 1) / n) / n

  if (!is.finite(hln_inside) || hln_inside <= 0) {
    hln_stat <- NA_real_
    p_two <- NA_real_
    p_better <- NA_real_
    status <- "HLN correction unavailable"
  } else {
    hln_stat <- dm_stat * sqrt(hln_inside)
    p_two <- 2 * stats::pt(-abs(hln_stat), df = n - 1L)
    # H1: benchmark loss - model loss > 0, i.e. model is better.
    p_better <- stats::pt(hln_stat, df = n - 1L, lower.tail = FALSE)
    status <- "OK"
  }

  tibble::tibble(
    n = n,
    mean_loss_benchmark = mean(l0),
    mean_loss_model = mean(l1),
    mean_loss_difference = d_bar,
    dm_statistic = dm_stat,
    hln_statistic = hln_stat,
    p_value_two_sided = p_two,
    p_value_model_better = p_better,
    hac_lag = hac_lag,
    test_status = status
  )
}


# ------------------------------------------------------------------------------
# CLARK-WEST (2007) MSPE-ADJUSTED TEST -- PROFESSOR'S IMPLEMENTATION
# ------------------------------------------------------------------------------



professor_newey_west_variance <- function(y, qn) {

  y <- as.numeric(y)
  y <- y[is.finite(y)]

  Tn <- length(y)

  if (Tn < 2L) {
    return(NA_real_)
  }

  qn <- max(1L, as.integer(qn))
  qn <- min(qn, Tn)

  ybar <- mean(y)
  dy <- y - ybar

  # MATLAB: G0 = dy' * dy / T
  G0 <- sum(dy^2) / Tn

  # MATLAB:
  # for j = 1:qn-1
  #   gamma = dy(j+1:T)' * dy(1:T-j) / (T-1)
  #   G0 = G0 + (gamma + gamma') * (1 - abs(j/qn))
  # end
  #
  # For a scalar series gamma == gamma', hence add 2*gamma*weight.
  if (qn > 1L) {

    for (j in seq_len(qn - 1L)) {

      gamma_j <-
        sum(
          dy[(j + 1L):Tn] *
            dy[1L:(Tn - j)]
        ) / (Tn - 1L)

      bartlett_weight <- 1 - abs(j / qn)

      G0 <- G0 +
        2 * gamma_j * bartlett_weight
    }
  }

  as.numeric(G0)
}


clark_west_mspe_adjusted <- function(
    actual,
    forecast_small,
    forecast_large,
    horizon
) {

  ok <-
    is.finite(actual) &
    is.finite(forecast_small) &
    is.finite(forecast_large)

  y <- as.numeric(actual[ok])
  yf1 <- as.numeric(forecast_small[ok])
  yf2 <- as.numeric(forecast_large[ok])

  P <- length(y)

  if (P < minimum_test_observations) {

    return(
      tibble::tibble(
        n_predictions = P,
        MSPE_small = NA_real_,
        MSPE_large = NA_real_,
        CW_adjustment = NA_real_,
        MSPE_large_adjusted = NA_real_,
        MSPE_adjusted_difference = NA_real_,
        NW_qn = as.integer(horizon),
        NW_max_autocovariance_lag =
          max(0L, as.integer(horizon) - 1L),
        CW_MSPE_adjusted_statistic = NA_real_,
        p_value_one_sided = NA_real_,
        reject_10pct = FALSE,
        reject_5pct = FALSE,
        test_status = "Too few matched OOS observations"
      )
    )
  }

  e1 <- y - yf1
  e2 <- y - yf2

  # Eq. (2.1) / supplied CW_test.m.
  froll_adj <-
    e1^2 -
    (
      e2^2 -
        (yf1 - yf2)^2
    )

  MSPE_1 <- mean(e1^2)
  MSPE_2 <- mean(e2^2)

  adjustment <-
    mean((yf1 - yf2)^2)

  MSPE_2_adjusted <-
    MSPE_2 - adjustment

  MSPE_adjusted_difference <-
    MSPE_1 - MSPE_2_adjusted

  # For h-step forecasts use qn=h in the professor's NW implementation,
  # which includes autocovariances j=1,...,h-1.
  nw_qn <- max(1L, as.integer(horizon))

  var_froll_adj <-
    professor_newey_west_variance(
      y = froll_adj,
      qn = nw_qn
    )

  if (
    !is.finite(var_froll_adj) ||
    var_froll_adj <= 0
  ) {

    return(
      tibble::tibble(
        n_predictions = P,
        MSPE_small = MSPE_1,
        MSPE_large = MSPE_2,
        CW_adjustment = adjustment,
        MSPE_large_adjusted = MSPE_2_adjusted,
        MSPE_adjusted_difference =
          MSPE_adjusted_difference,
        NW_qn = nw_qn,
        NW_max_autocovariance_lag = nw_qn - 1L,
        CW_MSPE_adjusted_statistic = NA_real_,
        p_value_one_sided = NA_real_,
        reject_10pct = FALSE,
        reject_5pct = FALSE,
        test_status = "Non-positive Newey-West variance"
      )
    )
  }

  # Exactly mirrors:
  # CWstat = sqrt(P) * mean(froll_adj) / sqrt(varfroll_adj)
  CWstat <-
    sqrt(P) *
    mean(froll_adj) /
    sqrt(var_froll_adj)

  # Clark-West uses a one-sided test; positive values favor the larger model.
  p_one_sided <-
    stats::pnorm(
      CWstat,
      lower.tail = FALSE
    )

  tibble::tibble(
    n_predictions = P,
    MSPE_small = MSPE_1,
    MSPE_large = MSPE_2,
    CW_adjustment = adjustment,
    MSPE_large_adjusted = MSPE_2_adjusted,
    MSPE_adjusted_difference =
      MSPE_adjusted_difference,
    NW_qn = nw_qn,
    NW_max_autocovariance_lag = nw_qn - 1L,
    CW_MSPE_adjusted_statistic = CWstat,
    p_value_one_sided = p_one_sided,
    reject_10pct =
      is.finite(CWstat) && CWstat > 1.282,
    reject_5pct =
      is.finite(CWstat) && CWstat > 1.645,
    test_status = "OK"
  )
}


# Determine the genuinely nested EN model pairs from the model specifications.
# This deliberately does NOT assume EN-1 -> EN-2 -> EN-3 -> EN-4 -> EN-5 is
# sequentially nested, because EN-2 and EN-3 contain different added variables.
build_nested_pairs <- function(specs) {

  model_names <- names(specs)
  rows <- list()
  id <- 0L

  for (small in model_names) {
    for (large in model_names) {

      if (identical(small, large)) next

      small_set <- unique(specs[[small]])
      large_set <- unique(specs[[large]])

      is_nested <- length(small_set) < length(large_set) &&
        all(small_set %in% large_set)

      if (!is_nested) next

      id <- id + 1L
      rows[[id]] <- tibble::tibble(
        smaller_model = small,
        larger_model = large,
        n_predictors_small = length(small_set),
        n_predictors_large = length(large_set),
        extra_predictors = paste(
          strip_lag_suffix(setdiff(large_set, small_set)),
          collapse = " + "
        )
      )
    }
  }

  dplyr::bind_rows(rows) |>
    dplyr::distinct(smaller_model, larger_model, .keep_all = TRUE) |>
    dplyr::arrange(n_predictors_small, n_predictors_large,
                   smaller_model, larger_model)
}


# ==============================================================================
# PART I. POINT FORECASTING

point_forecast_rows <- list()
point_selection_rows <- list()
point_fit_rows <- list()

point_forecast_id <- 0L
point_selection_id <- 0L
point_fit_id <- 0L

point_tuning_cache <- list()


# ==============================================================================
# 14. POINT BENCHMARK: AR(1), ALL WINDOWS x ALL HORIZONS
# ==============================================================================

cat("\nRunning point benchmark AR(1)\n")

for (h in forecast_horizons) {

  horizon_data <- master_panel |>
    dplyr::mutate(
      target_future = dplyr::lead(INDPRO_growth, h),
      target_date = dplyr::lead(date, h)
    )

  for (window_size in point_rolling_windows) {

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

      ar_train <- sample_i$training_df |>
        dplyr::filter(is.finite(AR1))

      if (nrow(ar_train) < 30L || !is.finite(sample_i$current_df$AR1)) next

      ar_model <- tryCatch(
        stats::lm(target_future ~ AR1, data = ar_train),
        error = function(e) NULL
      )
      if (is.null(ar_model)) next

      ar_forecast <- tryCatch(
        as.numeric(stats::predict(ar_model, newdata = sample_i$current_df)),
        error = function(e) NA_real_
      )
      if (!is.finite(ar_forecast)) next

      point_forecast_id <- point_forecast_id + 1L
      point_forecast_rows[[point_forecast_id]] <- tibble::tibble(
        window = window_size,
        horizon = h,
        forecast_origin = sample_i$origin_date,
        target_date = sample_i$target_date,
        model = "AR(1)",
        actual = sample_i$actual,
        forecast = ar_forecast,
        error = sample_i$actual - ar_forecast,
        alpha = NA_real_,
        lambda = NA_real_
      )
    }
  }
}


# ==============================================================================
# 15. POINT ELASTIC-NET MODELS EN-1 ... EN-5, ALL WINDOWS x ALL HORIZONS
# ==============================================================================

run_point_en_model <- function(model_name) {

  predictors <- point_specs[[model_name]]

  cat("\nRunning point model:", model_name, "\n")

  for (h in forecast_horizons) {

    horizon_data <- master_panel |>
      dplyr::mutate(
        target_future = dplyr::lead(INDPRO_growth, h),
        target_date = dplyr::lead(date, h)
      )

    for (window_size in point_rolling_windows) {

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

        cache_key <- paste(window_size, h, model_name, sep = "|")
        cached <- point_tuning_cache[[cache_key]]

        need_retune <- is.null(cached) ||
          (forecast_origin - cached$forecast_origin >= retune_every)

        if (need_retune) {
          tuning <- tune_point_en(training_df, predictors)

          if (!is.null(tuning)) {
            point_tuning_cache[[cache_key]] <<- list(
              forecast_origin = forecast_origin,
              tuning = tuning
            )
          }
        } else {
          tuning <- cached$tuning
        }

        if (is.null(tuning)) next

        best_alpha <- tuning$alpha[[1]]
        best_lambda <- tuning$lambda[[1]]

        matrices <- prepare_predictor_matrices(
          fit_df = training_df,
          apply_df = current_df,
          predictors = predictors
        )
        if (is.null(matrices)) next

        y_train <- training_df$target_future
        valid_y <- is.finite(y_train)

        x_train <- matrices$x_fit[valid_y, , drop = FALSE]
        y_train <- y_train[valid_y]
        x_current <- matrices$x_apply

        if (length(y_train) < 30L || ncol(x_train) == 0L) next

        final_model <- tryCatch(
          glmnet::glmnet(
            x = x_train,
            y = y_train,
            family = "gaussian",
            alpha = best_alpha,
            lambda = best_lambda,
            standardize = TRUE,
            intercept = TRUE
          ),
          error = function(e) NULL
        )
        if (is.null(final_model)) next

        en_forecast <- tryCatch(
          as.numeric(
            stats::predict(
              final_model,
              newx = x_current,
              s = best_lambda
            )
          ),
          error = function(e) NA_real_
        )
        if (!is.finite(en_forecast)) next

        point_forecast_id <<- point_forecast_id + 1L
        point_forecast_rows[[point_forecast_id]] <<- tibble::tibble(
          window = window_size,
          horizon = h,
          forecast_origin = sample_i$origin_date,
          target_date = sample_i$target_date,
          model = model_name,
          actual = sample_i$actual,
          forecast = en_forecast,
          error = sample_i$actual - en_forecast,
          alpha = best_alpha,
          lambda = best_lambda
        )

        point_fit_id <<- point_fit_id + 1L
        point_fit_rows[[point_fit_id]] <<- tibble::tibble(
          window = window_size,
          horizon = h,
          forecast_origin = sample_i$origin_date,
          model = model_name
        )

        co <- as.matrix(stats::coef(final_model, s = best_lambda))
        co_values <- as.numeric(co[, 1])
        co_names <- rownames(co)

        selected <- co_names[
          co_names != "(Intercept)" &
            is.finite(co_values) &
            abs(co_values) > selection_tolerance
        ]

        if (length(selected) > 0L) {
          point_selection_id <<- point_selection_id + 1L
          point_selection_rows[[point_selection_id]] <<- tibble::tibble(
            window = window_size,
            horizon = h,
            forecast_origin = sample_i$origin_date,
            model = model_name,
            variable = selected
          )
        }
      }
    }
  }

  invisible(NULL)
}


purrr::walk(names(point_specs), run_point_en_model)

point_forecasts <- dplyr::bind_rows(point_forecast_rows) |>
  dplyr::arrange(window, horizon, target_date, model)

point_fit_log <- dplyr::bind_rows(point_fit_rows)
point_selection_log <- dplyr::bind_rows(point_selection_rows)

if (nrow(point_forecasts) == 0L) {
  stop("No point forecasts were produced.")
}


# ==============================================================================
