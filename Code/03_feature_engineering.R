# ============================================================================
# 03_feature_engineering.R
# Derived from the uploaded Madhavi_panchangum_fred.R; original methods preserved.
# Run from the repository root via Rscript run_all.R.
# ============================================================================


# 7. KEEP ACTUAL TRADING-DAY FINANCIAL VALUES AS ELASTIC-NET PREDICTORS
# ==============================================================================



financial_daily_model_data <- financial_daily |>
  dplyr::filter(
    is.finite(SP500_DAILY_LOG_RETURN),
    is.finite(DGS10_DAILY_LOG_CHANGE),
    is.finite(SP500_X_DGS10_MODEL)
  ) |>
  dplyr::mutate(
    month = lubridate::floor_date(date, unit = "month")
  ) |>
  dplyr::arrange(month, date) |>
  dplyr::group_by(month) |>
  dplyr::mutate(
    trading_day_in_month = dplyr::row_number()
  ) |>
  dplyr::ungroup()

trading_days_per_month <- financial_daily_model_data |>
  dplyr::count(month, name = "actual_trading_days") |>
  dplyr::arrange(month)

max_trading_days <- max(
  trading_days_per_month$actual_trading_days,
  na.rm = TRUE
)

min_trading_days <- min(
  trading_days_per_month$actual_trading_days,
  na.rm = TRUE
)


make_trading_day_block <- function(
    data,
    value_col,
    prefix,
    max_slots
) {

  value_col <- rlang::ensym(value_col)

  long <- data |>
    dplyr::select(
      month,
      trading_day_in_month,
      value = !!value_col
    ) |>
    dplyr::filter(
      is.finite(value),
      trading_day_in_month <= max_slots
    ) |>
    dplyr::mutate(
      value_name = sprintf(
        "%s_D%02d",
        prefix,
        trading_day_in_month
      ),
      availability_name = sprintf(
        "%s_AVAILABLE_D%02d",
        prefix,
        trading_day_in_month
      )
    )

  value_wide <- long |>
    dplyr::select(month, value_name, value) |>
    tidyr::pivot_wider(
      names_from = value_name,
      values_from = value
    )

  availability_wide <- long |>
    dplyr::transmute(
      month,
      availability_name,
      available = 1
    ) |>
    tidyr::pivot_wider(
      names_from = availability_name,
      values_from = available
    )

  all_months <- tibble::tibble(
    month = sort(unique(data$month))
  )

  out <- all_months |>
    dplyr::left_join(value_wide, by = "month") |>
    dplyr::left_join(availability_wide, by = "month") |>
    dplyr::rename(date = month)

  value_names <- sprintf(
    "%s_D%02d",
    prefix,
    seq_len(max_slots)
  )

  availability_names <- sprintf(
    "%s_AVAILABLE_D%02d",
    prefix,
    seq_len(max_slots)
  )

  for (nm in setdiff(value_names, names(out))) {
    out[[nm]] <- NA_real_
  }

  for (nm in setdiff(availability_names, names(out))) {
    out[[nm]] <- NA_real_
  }

  # Structural absent slots are encoded, never LOCF-filled.
  out[value_names] <- lapply(
    out[value_names],
    function(x) {
      x[!is.finite(x)] <- 0
      as.numeric(x)
    }
  )

  out[availability_names] <- lapply(
    out[availability_names],
    function(x) {
      x[!is.finite(x)] <- 0
      as.numeric(x)
    }
  )

  out <- out |>
    dplyr::select(
      date,
      dplyr::all_of(value_names),
      dplyr::all_of(availability_names)
    ) |>
    dplyr::arrange(date)

  list(
    data = out,
    value_names = value_names,
    availability_names = availability_names
  )
}


sp500_daily_block <- make_trading_day_block(
  data = financial_daily_model_data,
  value_col = SP500_DAILY_LOG_RETURN,
  prefix = "SP500_RET",
  max_slots = max_trading_days
)

dgs10_daily_block <- make_trading_day_block(
  data = financial_daily_model_data,
  value_col = DGS10_DAILY_LOG_CHANGE,
  prefix = "DGS10_LOG",
  max_slots = max_trading_days
)

interaction_daily_block <- make_trading_day_block(
  data = financial_daily_model_data,
  value_col = SP500_X_DGS10_MODEL,
  prefix = "INTERACTION",
  max_slots = max_trading_days
)


financial_daily_predictors <- sp500_daily_block$data |>
  dplyr::inner_join(
    dgs10_daily_block$data,
    by = "date"
  ) |>
  dplyr::inner_join(
    interaction_daily_block$data,
    by = "date"
  ) |>
  dplyr::arrange(date)

sp500_daily_names <- sp500_daily_block$value_names
dgs10_daily_names <- dgs10_daily_block$value_names
interaction_daily_names <- interaction_daily_block$value_names

sp500_availability_names <- sp500_daily_block$availability_names
dgs10_availability_names <- dgs10_daily_block$availability_names
interaction_availability_names <- interaction_daily_block$availability_names


financial_daily_block_diagnostics <- trading_days_per_month |>
  dplyr::mutate(
    max_slots_in_rectangular_model_matrix = max_trading_days,
    structural_slots_not_observed =
      max_slots_in_rectangular_model_matrix - actual_trading_days,
    locf_positions_added_for_nontrading_days = 0L
  )


stopifnot(
  length(sp500_daily_names) == max_trading_days,
  length(dgs10_daily_names) == max_trading_days,
  length(interaction_daily_names) == max_trading_days,
  all(
    financial_daily_block_diagnostics$
      locf_positions_added_for_nontrading_days == 0L
  )
)


# ------------------------------------------------------------------------------
# 7A. CORRELATION DIAGNOSTIC USING RAW LEVELS AFTER LOCF
# ------------------------------------------------------------------------------

# This plot remains a descriptive monthly correlation diagnostic only; it is not
# used to construct the elastic-net regressors. The elastic nets receive the
# actual-trading-day columns created above.

financial_raw_month_end <- financial_daily_locf |>
  dplyr::mutate(date = lubridate::floor_date(date, unit = "month")) |>
  dplyr::group_by(date) |>
  dplyr::summarise(
    SP500_RAW = dplyr::last(SP500),
    DGS10_RAW = dplyr::last(DGS10),
    .groups = "drop"
  )

indpro_raw_monthly <- indpro_data |>
  dplyr::transmute(
    date = lubridate::floor_date(date, unit = "month"),
    INDPRO_RAW = INDPRO
  ) |>
  dplyr::filter(is.finite(INDPRO_RAW))

correlation_panel_raw <- indpro_raw_monthly |>
  dplyr::inner_join(financial_raw_month_end, by = "date") |>
  dplyr::arrange(date)

correlation_matrix_raw <- stats::cor(
  correlation_panel_raw |>
    dplyr::select(INDPRO_RAW, SP500_RAW, DGS10_RAW),
  use = "pairwise.complete.obs",
  method = "pearson"
)

print(round(correlation_matrix_raw, 3))


# ==============================================================================
# 8. DEFINE THE FRED-MD PREDICTOR SET
# ==============================================================================

# Keep the original exclusion rule so that the separately downloaded daily S&P
# 500 and DGS10 blocks are not duplicated by stock/interest-rate FRED-MD series.

fredmd_interest_rate_variables <- c(
  "FEDFUNDS",
  "CP3Mx",
  "TB3MS",
  "TB6MS",
  "GS1",
  "GS5",
  "GS10",
  "AAA",
  "BAA",
  "COMPAPFFx",
  "TB3SMFFM",
  "TB6SMFFM",
  "T1YFFM",
  "T5YFFM",
  "T10YFFM",
  "AAAFFM",
  "BAAFFM"
)

fredmd_stock_market_variables <- c(
  "S&P 500",
  "S&P: indust",
  "S&P div yield",
  "S&P PE ratio",
  "SP500",
  "SP500x",
  "S&P500"
)

all_macro_variables <- setdiff(
  names(macro_transformed),
  "date"
)

stock_regex_matches <- all_macro_variables[
  stringr::str_detect(
    all_macro_variables,
    regex("(^S&P)|(^SP500)|stock|equity", ignore_case = TRUE)
  )
]

fredmd_excluded_variables <- unique(c(
  "INDPRO",
  intersect(fredmd_interest_rate_variables, all_macro_variables),
  intersect(fredmd_stock_market_variables, all_macro_variables),
  stock_regex_matches
))

fredmd_predictors <- setdiff(
  all_macro_variables,
  fredmd_excluded_variables
)

cat("\nFRED-MD predictor definition\n")
cat("All transformed macro variables:", length(all_macro_variables), "\n")
cat("Excluded target/interest/stock variables:", length(fredmd_excluded_variables), "\n")
cat("FRED-MD predictors retained:", length(fredmd_predictors), "\n\n")


# ==============================================================================
# 9. BUILD MASTER MONTHLY PANEL WITH DAILY FINANCIAL PREDICTOR COLUMNS
# ==============================================================================

macro_panel <- macro_transformed |>
  dplyr::mutate(date = lubridate::floor_date(date, unit = "month")) |>
  dplyr::select(date, dplyr::all_of(fredmd_predictors))

master_panel <- indpro_growth_data |>
  dplyr::mutate(date = lubridate::floor_date(date, unit = "month")) |>
  dplyr::left_join(macro_panel, by = "date") |>
  dplyr::left_join(financial_daily_predictors, by = "date") |>
  dplyr::filter(
    date >= lubridate::floor_date(financial_start_date, unit = "month")
  ) |>
  dplyr::arrange(date)



all_daily_financial_predictors <- c(
  sp500_daily_names,
  sp500_availability_names,
  dgs10_daily_names,
  dgs10_availability_names,
  interaction_daily_names,
  interaction_availability_names
)

all_exogenous_predictors <- c(
  fredmd_predictors,
  all_daily_financial_predictors
)

master_panel <- master_panel |>
  dplyr::mutate(AR1 = INDPRO_growth)

cat("\nMASTER PANEL DATE RANGE\n")
cat("Start:", format(min(master_panel$date, na.rm = TRUE)), "\n")
cat("End:", format(max(master_panel$date, na.rm = TRUE)), "\n")
cat("Number of rows:", nrow(master_panel), "\n")
cat("Maximum S&P trading-day slots:", length(sp500_daily_names), "\n")
cat("Maximum DGS10 trading-day slots:", length(dgs10_daily_names), "\n")
cat("Maximum interaction trading-day slots:", length(interaction_daily_names), "\n")
cat("Actual number of trading days for each month is saved separately.\n\n")


# ==============================================================================
# 10. MODEL SPECIFICATIONS USING DAILY FINANCIAL VALUES
# ==============================================================================

# Current-month predictor names. The historical object names are retained so that
# the helper functions below do not need to be rewritten, but these variables are
# NOT lagged: they refer to information available at forecast origin month t.
fred_lagged <- fredmd_predictors

sp_lagged <- c(
  sp500_daily_names,
  sp500_availability_names
)

dgs_lagged <- c(
  dgs10_daily_names,
  dgs10_availability_names
)

interaction_lagged <- c(
  interaction_daily_names,
  interaction_availability_names
)

# The economic EN definitions are preserved, but every financial component is
# now a BLOCK OF ACTUAL TRADING-DAY VALUES rather than one monthly aggregate.
#
# EN-1 : FRED-MD only
# EN-2 : FRED-MD + actual trading-day S&P return values
# EN-3 : FRED-MD + actual trading-day DGS10 values
# EN-4 : FRED-MD + both actual-trading-day financial blocks
# EN-5 : EN-4 + actual trading-day interaction values

point_specs <- list(
  "EN-1" = fred_lagged,
  "EN-2" = c(fred_lagged, sp_lagged),
  "EN-3" = c(fred_lagged, dgs_lagged),
  "EN-4" = c(fred_lagged, sp_lagged, dgs_lagged),
  "EN-5" = c(fred_lagged, sp_lagged, dgs_lagged, interaction_lagged)
)

quantile_specs <- list(
  "Q-EN-1" = fred_lagged,
  "Q-EN-2" = c(fred_lagged, sp_lagged),
  "Q-EN-3" = c(fred_lagged, dgs_lagged),
  "Q-EN-4" = c(fred_lagged, sp_lagged, dgs_lagged),
  "Q-EN-5" = c(fred_lagged, sp_lagged, dgs_lagged, interaction_lagged)
)

point_specs <- purrr::map(
  point_specs,
  ~ intersect(.x, names(master_panel))
)

quantile_specs <- purrr::map(
  quantile_specs,
  ~ intersect(.x, names(master_panel))
)

point_predictor_map <- purrr::imap_dfr(
  point_specs,
  ~ tibble::tibble(model = .y, variable = .x)
)

quantile_predictor_map <- purrr::imap_dfr(
  quantile_specs,
  ~ tibble::tibble(model = .y, variable = .x)
)



# ==============================================================================
