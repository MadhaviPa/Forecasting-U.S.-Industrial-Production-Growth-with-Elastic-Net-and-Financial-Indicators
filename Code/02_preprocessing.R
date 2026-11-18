# ============================================================================
# 02_preprocessing.R
# Derived from the uploaded Madhavi_panchangum_fred.R; original methods preserved.
# Run from the repository root via Rscript run_all.R.
# ============================================================================


# 3A. CHECK MISSING VALUES IN THE RAW FRED-MD DATA
# ==============================================================================

fredmd_numeric_for_diagnostics <- raw_data |>
  dplyr::select(-sasdate) |>
  dplyr::mutate(
    dplyr::across(
      dplyr::everything(),
      as.numeric
    )
  )

missing_before_by_series <- vapply(
  fredmd_numeric_for_diagnostics,
  function(x) sum(!is.finite(x)),
  integer(1)
)

fredmd_missing_diagnostics <- tibble::tibble(
  series = names(missing_before_by_series),
  missing_in_raw_fredmd = as.integer(missing_before_by_series),
  pct_missing_in_raw_fredmd =
    100 * missing_in_raw_fredmd /
    nrow(fredmd_numeric_for_diagnostics)
) |>
  dplyr::arrange(
    dplyr::desc(missing_in_raw_fredmd)
  )

cat("FRED-MD missing-value diagnostics\n")
cat(
  "Total missing values:",
  sum(missing_before_by_series),
  "\n"
)
cat(
  "Series containing missing values:",
  sum(missing_before_by_series > 0L),
  "\n\n"
)


# ==============================================================================
# 4. FRED-MD TRANSFORMATION FUNCTIONS
# ==============================================================================

# FRED-MD transformation codes:
#
# 1 = x_t
# 2 = Delta x_t
# 3 = Delta^2 x_t
# 4 = log(x_t)
# 5 = Delta log(x_t)
# 6 = Delta^2 log(x_t)
# 7 = Delta(x_t / x_{t-1} - 1)
#
# These transformations follow the transformation codes contained in
# the first row of the FRED-MD dataset.

transformation_name <- function(tcode) {
  
  dplyr::case_when(
    tcode == 1 ~ "Level",
    tcode == 2 ~ "First difference",
    tcode == 3 ~ "Second difference",
    tcode == 4 ~ "Log",
    tcode == 5 ~ "First difference of log",
    tcode == 6 ~ "Second difference of log",
    tcode == 7 ~ "Difference of percentage change",
    TRUE ~ "Unknown"
  )
}


apply_fredmd_transform <- function(x, tcode) {
  
  x <- as.numeric(x)
  
  # If the transformation code is missing, return NA.
  if (!is.finite(tcode)) {
    return(rep(NA_real_, length(x)))
  }
  
  
  # ---------------------------------------------------------------------------
  # Code 1: No transformation
  # x_t
  # ---------------------------------------------------------------------------
  
  if (tcode == 1) {
    
    out <- x
    
    
    # ---------------------------------------------------------------------------
    # Code 2: First difference
    # Delta x_t = x_t - x_{t-1}
    # ---------------------------------------------------------------------------
    
  } else if (tcode == 2) {
    
    out <- c(
      NA_real_,
      diff(x)
    )
    
    
    # ---------------------------------------------------------------------------
    # Code 3: Second difference
    # Delta^2 x_t
    # ---------------------------------------------------------------------------
    
  } else if (tcode == 3) {
    
    out <- c(
      NA_real_,
      NA_real_,
      diff(x, differences = 2)
    )
    
    
    # ---------------------------------------------------------------------------
    # Code 4: Log level
    # log(x_t)
    # ---------------------------------------------------------------------------
    
  } else if (tcode == 4) {
    
    log_x <- rep(NA_real_, length(x))
    
    positive <- is.finite(x) & x > 0
    
    log_x[positive] <- log(x[positive])
    
    out <- log_x
    
    
    # ---------------------------------------------------------------------------
    # Code 5: First difference of log
    # Delta log(x_t)
    # ---------------------------------------------------------------------------
    
  } else if (tcode == 5) {
    
    log_x <- rep(NA_real_, length(x))
    
    positive <- is.finite(x) & x > 0
    
    log_x[positive] <- log(x[positive])
    
    out <- c(
      NA_real_,
      diff(log_x)
    )
    
    
    # ---------------------------------------------------------------------------
    # Code 6: Second difference of log
    # Delta^2 log(x_t)
    # ---------------------------------------------------------------------------
    
  } else if (tcode == 6) {
    
    log_x <- rep(NA_real_, length(x))
    
    positive <- is.finite(x) & x > 0
    
    log_x[positive] <- log(x[positive])
    
    out <- c(
      NA_real_,
      NA_real_,
      diff(log_x, differences = 2)
    )
    
    
    # ---------------------------------------------------------------------------
    # Code 7:
    # Delta(x_t / x_{t-1} - 1)
    # ---------------------------------------------------------------------------
    
  } else if (tcode == 7) {
    
    growth <- x / dplyr::lag(x) - 1
    
    out <- c(
      NA_real_,
      diff(growth)
    )
    
    
  } else {
    
    stop(
      paste(
        "Unknown FRED-MD transformation code:",
        tcode
      )
    )
  }
  
  as.numeric(out)
}


# ==============================================================================
# 4A. AUGMENTED DICKEY-FULLER TEST
# ==============================================================================

# H0: the series contains a unit root.
# p < 0.05 is treated as evidence of stationarity.

adf_stationarity <- function(x, alpha = 0.05) {
  
  x <- as.numeric(x)
  
  x <- x[is.finite(x)]
  
  if (
    length(x) < 30L ||
    !is.finite(stats::sd(x)) ||
    stats::sd(x) == 0
  ) {
    
    return(
      list(
        stationary = FALSE,
        p_value = NA_real_
      )
    )
  }
  
  test <- tryCatch(
    suppressWarnings(
      tseries::adf.test(
        x,
        alternative = "stationary"
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(test)) {
    
    return(
      list(
        stationary = FALSE,
        p_value = NA_real_
      )
    )
  }
  
  list(
    stationary =
      is.finite(test$p.value) &&
      test$p.value < alpha,
    
    p_value =
      as.numeric(test$p.value)
  )
}


# ==============================================================================
# 5. INDUSTRIAL PRODUCTION TARGET
# ==============================================================================

if (!"INDPRO" %in% names(raw_data)) {
  
  stop(
    "INDPRO is not present in the FRED-MD dataset."
  )
}


indpro_data <- raw_data |>
  dplyr::transmute(
    date = sasdate,
    INDPRO = as.numeric(INDPRO)
  ) |>
  dplyr::arrange(date)


# ------------------------------------------------------------------------------
# INDPRO has FRED-MD transformation code 5:
#
# Delta log(IP_t)
#
# For interpretation as monthly percentage log growth, multiply the
# transformation by 100.
# ------------------------------------------------------------------------------

indpro_tcode <- transformation_codes[["INDPRO"]]

if (
  is.null(indpro_tcode) ||
  !is.finite(indpro_tcode)
) {
  
  stop(
    "No valid FRED-MD transformation code was found for INDPRO."
  )
}


indpro_transformed <- apply_fredmd_transform(
  x = indpro_data$INDPRO,
  tcode = indpro_tcode
)


indpro_growth_data <- indpro_data |>
  dplyr::mutate(
    INDPRO_growth =
      100 * indpro_transformed
  ) |>
  dplyr::select(
    date,
    INDPRO,
    INDPRO_growth
  )


# Check stationarity of the transformed target.

indpro_test <- adf_stationarity(
  indpro_growth_data$INDPRO_growth
)


indpro_stationarity_diagnostics <- tibble::tibble(
  series = "INDPRO",
  frequency = "monthly",
  transformation_code = indpro_tcode,
  transformation =
    transformation_name(indpro_tcode),
  scale = "100 x transformed value",
  adf_p_value = indpro_test$p_value,
  stationary_after_transform =
    indpro_test$stationary
)


cat("\nINDPRO transformation\n")
cat(
  "Transformation code:",
  indpro_tcode,
  "\n"
)
cat(
  "Transformation:",
  transformation_name(indpro_tcode),
  "\n"
)
cat(
  "ADF p-value:",
  indpro_test$p_value,
  "\n"
)
cat(
  "Stationary:",
  indpro_test$stationary,
  "\n\n"
)


# ==============================================================================
# 5A. APPLY FRED-MD TRANSFORMATION CODES TO ALL MACRO VARIABLES
# ==============================================================================

# Keep the original raw data unchanged.
# macro_transformed will contain the transformed variables used in modelling.

macro_transformed <- raw_data


transformation_diagnostics <- list()

diagnostic_id <- 0L


for (variable in names(transformation_codes)) {
  
  # Skip variables not found in the current vintage.
  if (!variable %in% names(macro_transformed)) {
    next
  }
  
  # Convert the series safely to numeric.
  x <- suppressWarnings(
    as.numeric(
      macro_transformed[[variable]]
    )
  )
  
  tcode <- transformation_codes[[variable]]
  
  # Apply the transformation specified by FRED-MD.
  transformed_values <- apply_fredmd_transform(
    x = x,
    tcode = tcode
  )
  
  macro_transformed[[variable]] <-
    transformed_values
  
  
  # Test the transformed series for stationarity.
  stationarity_test <- adf_stationarity(
    transformed_values
  )
  
  
  diagnostic_id <- diagnostic_id + 1L
  
  
  transformation_diagnostics[[diagnostic_id]] <-
    tibble::tibble(
      
      series = variable,
      
      transformation_code =
        as.integer(tcode),
      
      transformation =
        transformation_name(tcode),
      
      adf_p_value =
        stationarity_test$p_value,
      
      stationary_after_transform =
        stationarity_test$stationary
    )
}


transformation_diagnostics <-
  dplyr::bind_rows(
    transformation_diagnostics
  )


# Rename the date variable after all transformations are applied.

macro_transformed <- macro_transformed |>
  dplyr::rename(
    date = sasdate
  ) |>
  dplyr::arrange(date)


# ==============================================================================
# 5B. TRANSFORMATION SUMMARY
# ==============================================================================

cat("\nFRED-MD transformation summary\n")
cat(
  "Number of transformed series:",
  nrow(transformation_diagnostics),
  "\n\n"
)


print(
  transformation_diagnostics |>
    dplyr::count(
      transformation_code,
      transformation,
      name = "number_of_series"
    ) |>
    dplyr::arrange(
      transformation_code
    )
)


# Optional: save transformation diagnostics.

readr::write_csv(
  transformation_diagnostics,
  file.path(
    output_dir,
    "fredmd_transformation_diagnostics.csv"
  )
)
# ==============================================================================
# 6B. CHECK EVERY CALENDAR DATE, MONTH BY MONTH:
#     NON-TRADING DAY OR TRUE MISSING VALUE ON A TRADING DAY?
# ------------------------------------------------------------------------------

# The purpose of this section is to distinguish:
#
#   1. NON-TRADING DATE
#      Weekend / exchange holiday / other date on which the S&P 500 did not
#      trade. These dates are NOT missing observations, are NOT counted as
#      missing, and are NEVER filled.
#
#   2. ACTUAL TRADING DATE
#      A date appearing in Yahoo's S&P 500 daily date index. On these dates we
#      inspect SP500 and DGS10 separately. If either required raw value is
#      missing, it is a TRUE missing observation on a trading date and is
#      eligible for LOCF.
#
# We explicitly build a complete CALENDAR-DATE audit panel only for diagnostics.
# The modelling panel later drops all non-trading dates before LOCF.

financial_calendar_start <- max(
  min(sp500_raw_daily$date, na.rm = TRUE),
  min(treasury10_raw_daily$date, na.rm = TRUE)
)

financial_calendar_end <- min(
  max(sp500_raw_daily$date, na.rm = TRUE),
  max(treasury10_raw_daily$date, na.rm = TRUE)
)

all_calendar_dates <- tibble::tibble(
  date = seq.Date(
    financial_calendar_start,
    financial_calendar_end,
    by = "day"
  )
) |>
  dplyr::mutate(
    month = lubridate::floor_date(date, unit = "month")
  )


# A date is classified as a TRADING DATE from the actual Yahoo S&P 500 date
# index, not from day-of-week rules. Therefore exchange holidays are correctly
# treated as non-trading dates even when they fall Monday-Friday.
sp500_trading_date_index <- sp500_raw_daily |>
  dplyr::distinct(date) |>
  dplyr::mutate(is_trading_day = TRUE)


financial_date_audit <- all_calendar_dates |>
  dplyr::left_join(
    sp500_trading_date_index,
    by = "date"
  ) |>
  dplyr::mutate(
    is_trading_day = dplyr::coalesce(is_trading_day, FALSE)
  ) |>
  dplyr::left_join(
    sp500_raw_daily |>
      dplyr::select(date, SP500),
    by = "date"
  ) |>
  dplyr::left_join(
    treasury10_raw_daily |>
      dplyr::select(date, DGS10),
    by = "date"
  ) |>
  dplyr::mutate(

    SP500_missing = !is.finite(SP500),
    DGS10_missing = !is.finite(DGS10),

    date_status = dplyr::case_when(

      !is_trading_day ~
        "NON_TRADING_DAY_NOT_MISSING_NOT_FILLED",

      is_trading_day &
        !SP500_missing &
        !DGS10_missing ~
        "TRADING_DAY_COMPLETE",

      is_trading_day &
        SP500_missing &
        !DGS10_missing ~
        "TRADING_DAY_TRUE_SP500_MISSING",

      is_trading_day &
        !SP500_missing &
        DGS10_missing ~
        "TRADING_DAY_TRUE_DGS10_MISSING",

      is_trading_day &
        SP500_missing &
        DGS10_missing ~
        "TRADING_DAY_BOTH_TRUE_MISSING",

      TRUE ~
        "UNCLASSIFIED"
    ),

    count_as_missing_SP500 =
      is_trading_day & SP500_missing,

    count_as_missing_DGS10 =
      is_trading_day & DGS10_missing,

    eligible_for_locf_SP500 =
      is_trading_day & SP500_missing,

    eligible_for_locf_DGS10 =
      is_trading_day & DGS10_missing
  ) |>
  dplyr::arrange(date)


# Save the complete date-by-date audit.
readr::write_csv(
  financial_date_audit,
  file.path(
    output_dir,
    "daily_financial_date_audit_trading_vs_nontrading.csv"
  )
)


# Monthly diagnostics: this is the explicit month-by-month check requested.
financial_monthly_missing_audit <- financial_date_audit |>
  dplyr::group_by(month) |>
  dplyr::summarise(

    calendar_dates_in_month = dplyr::n(),

    actual_trading_dates =
      sum(is_trading_day),

    non_trading_dates =
      sum(!is_trading_day),

    SP500_true_missing_on_trading_days =
      sum(count_as_missing_SP500),

    DGS10_true_missing_on_trading_days =
      sum(count_as_missing_DGS10),

    trading_dates_with_any_true_missing =
      sum(
        is_trading_day &
          (SP500_missing | DGS10_missing)
      ),

    trading_dates_complete_before_locf =
      sum(
        is_trading_day &
          !SP500_missing &
          !DGS10_missing
      ),

    .groups = "drop"
  ) |>
  dplyr::arrange(month)


cat("\n============================================================\n")
cat("MONTH-BY-MONTH FINANCIAL DATE AUDIT\n")
cat("============================================================\n")
cat(
  "Non-trading dates are classified separately and are NOT counted as missing.\n"
)
cat(
  "Only missing values on actual Yahoo S&P trading dates are true missing values.\n\n"
)

print(financial_monthly_missing_audit)


readr::write_csv(
  financial_monthly_missing_audit,
  file.path(
    output_dir,
    "monthly_financial_missing_audit_before_LOCF.csv"
  )
)


# Overall missing summary based ONLY on actual trading dates.
financial_missing_before_locf <- tibble::tibble(

  series = c("SP500", "DGS10"),

  n_actual_trading_dates =
    sum(financial_date_audit$is_trading_day),

  missing_before_locf = c(
    sum(financial_date_audit$count_as_missing_SP500),
    sum(financial_date_audit$count_as_missing_DGS10)
  )
) |>
  dplyr::mutate(
    pct_missing_before_locf =
      100 * missing_before_locf /
      n_actual_trading_dates
  )


cat("\n============================================================\n")
cat("TRUE MISSING VALUES ON ACTUAL TRADING DATES BEFORE LOCF\n")
cat("============================================================\n")
print(financial_missing_before_locf)


readr::write_csv(
  financial_missing_before_locf,
  file.path(
    output_dir,
    "daily_financial_missing_before_LOCF.csv"
  )
)


# ------------------------------------------------------------------------------
# 6C. REMOVE NON-TRADING DATES, THEN LOCF TRUE TRADING-DAY MISSING VALUES ONLY
# ------------------------------------------------------------------------------

# Keep ONLY actual trading dates before applying LOCF.
#
# Thus:
#   * weekends never enter LOCF;
#   * exchange holidays never enter LOCF;
#   * non-trading dates are never manufactured as observations;
#   * only true gaps on actual trading dates can be filled.
#
# LOCF runs through the ordered trading-day sequence, not separately within each
# month. This remains leakage-free and permits, for example, a genuinely missing
# value on the first trading day of a month to use the last observed value from
# the final trading day of the previous month. Future observations are never used.

financial_trading_days_before_locf <- financial_date_audit |>
  dplyr::filter(is_trading_day) |>
  dplyr::select(
    date,
    month,
    SP500,
    DGS10,
    SP500_missing_before_locf = SP500_missing,
    DGS10_missing_before_locf = DGS10_missing
  ) |>
  dplyr::arrange(date)


financial_daily_locf <- financial_trading_days_before_locf |>
  dplyr::mutate(

    SP500_before_locf = SP500,
    DGS10_before_locf = DGS10,

    SP500 = zoo::na.locf(
      SP500,
      na.rm = FALSE
    ),

    DGS10 = zoo::na.locf(
      DGS10,
      na.rm = FALSE
    ),

    SP500_filled_by_locf =
      SP500_missing_before_locf &
      is.finite(SP500),

    DGS10_filled_by_locf =
      DGS10_missing_before_locf &
      is.finite(DGS10)
  ) |>
  dplyr::arrange(date)


# Month-by-month audit AFTER LOCF.
financial_monthly_locf_audit <- financial_daily_locf |>
  dplyr::group_by(month) |>
  dplyr::summarise(

    actual_trading_dates =
      dplyr::n(),

    SP500_true_missing_before_locf =
      sum(SP500_missing_before_locf),

    SP500_values_filled_by_locf =
      sum(SP500_filled_by_locf),

    SP500_missing_after_locf =
      sum(!is.finite(SP500)),

    DGS10_true_missing_before_locf =
      sum(DGS10_missing_before_locf),

    DGS10_values_filled_by_locf =
      sum(DGS10_filled_by_locf),

    DGS10_missing_after_locf =
      sum(!is.finite(DGS10)),

    .groups = "drop"
  ) |>
  dplyr::arrange(month)


cat("\n============================================================\n")
cat("MONTH-BY-MONTH TRUE TRADING-DAY MISSING VALUES AFTER LOCF\n")
cat("============================================================\n")
print(financial_monthly_locf_audit)


readr::write_csv(
  financial_monthly_locf_audit,
  file.path(
    output_dir,
    "monthly_financial_LOCF_audit.csv"
  )
)


# Save every actual trading date with before/after values and LOCF flags.
readr::write_csv(
  financial_daily_locf,
  file.path(
    output_dir,
    "daily_financial_trading_days_before_after_LOCF.csv"
  )
)


financial_missing_after_locf <- tibble::tibble(

  series = c("SP500", "DGS10"),

  n_actual_trading_dates =
    nrow(financial_daily_locf),

  true_missing_before_locf = c(
    sum(financial_daily_locf$SP500_missing_before_locf),
    sum(financial_daily_locf$DGS10_missing_before_locf)
  ),

  filled_by_locf = c(
    sum(financial_daily_locf$SP500_filled_by_locf),
    sum(financial_daily_locf$DGS10_filled_by_locf)
  ),

  missing_after_locf = c(
    sum(!is.finite(financial_daily_locf$SP500)),
    sum(!is.finite(financial_daily_locf$DGS10))
  )
)


cat("\n============================================================\n")
cat("OVERALL TRADING-DAY MISSING VALUES AFTER LOCF\n")
cat("============================================================\n")
print(financial_missing_after_locf)


readr::write_csv(
  financial_missing_after_locf,
  file.path(
    output_dir,
    "daily_financial_missing_after_LOCF.csv"
  )
)


# Any unresolved leading missing value has no earlier observation available for
# leakage-free LOCF. Remove only such unresolved leading trading-date rows before
# transformations and modelling.
financial_daily_locf <- financial_daily_locf |>
  dplyr::filter(
    is.finite(SP500),
    is.finite(DGS10)
  ) |>
  dplyr::select(
    date,
    SP500,
    DGS10
  ) |>
  dplyr::arrange(date)


stopifnot(
  all(is.finite(financial_daily_locf$SP500)),
  all(is.finite(financial_daily_locf$DGS10))
)


# ------------------------------------------------------------------------------
# 6D. COMPUTE DAILY STATIONARY FINANCIAL VARIABLES AFTER LOCF
# ------------------------------------------------------------------------------

# S&P 500 standalone predictor:
#   100 * Delta log(SP500_t)
#
# Treasury standalone predictor:
#   100 * Delta log(DGS10_t)
#
# Interaction component for Treasury:
#   Delta DGS10_t

financial_daily <- financial_daily_locf |>
  dplyr::mutate(
    SP500_DAILY_LOG_RETURN = dplyr::if_else(
      SP500 > 0 & dplyr::lag(SP500) > 0,
      100 * (log(SP500) - log(dplyr::lag(SP500))),
      NA_real_
    ),
    DGS10_DAILY_LOG_CHANGE = dplyr::if_else(
      DGS10 > 0 & dplyr::lag(DGS10) > 0,
      100 * (log(DGS10) - log(dplyr::lag(DGS10))),
      NA_real_
    ),
    DGS10_DAILY_CHANGE = DGS10 - dplyr::lag(DGS10)
  )


# ------------------------------------------------------------------------------
# 6E. DAILY S&P 500 x 10-YEAR TREASURY INTERACTION
# ------------------------------------------------------------------------------

# Construct the interaction only AFTER the raw levels have been completed by
# LOCF and the daily transformations have been formed:
#
#   Interaction_t = S&P500 daily log return_t * Delta DGS10_t

# 1. Construct interaction

financial_daily <- financial_daily |>
  dplyr::mutate(
    SP500_X_DGS10_DAILY = dplyr::if_else(
      is.finite(SP500_DAILY_LOG_RETURN) &
        is.finite(DGS10_DAILY_CHANGE),
      SP500_DAILY_LOG_RETURN * DGS10_DAILY_CHANGE,
      NA_real_
    )
  ) |>
  dplyr::arrange(date)


# 2. Test interaction stationarity
interaction_level_test <- adf_stationarity(
  financial_daily$SP500_X_DGS10_DAILY
)


# 3. Decide level versus first difference
if (isTRUE(interaction_level_test$stationary)) {
  
  interaction_transformation_used <- "Level interaction"
  interaction_was_differenced <- FALSE
  
  financial_daily <- financial_daily |>
    dplyr::mutate(
      SP500_X_DGS10_FIRST_DIFF = NA_real_,
      SP500_X_DGS10_MODEL = SP500_X_DGS10_DAILY
    )
  
  interaction_final_test <- interaction_level_test
  
} else {
  
  interaction_transformation_used <- "First difference of interaction"
  interaction_was_differenced <- TRUE
  
  financial_daily <- financial_daily |>
    dplyr::mutate(
      SP500_X_DGS10_FIRST_DIFF =
        SP500_X_DGS10_DAILY -
        dplyr::lag(SP500_X_DGS10_DAILY),
      
      SP500_X_DGS10_MODEL =
        SP500_X_DGS10_FIRST_DIFF
    )
  
  # Re-test after first differencing
  interaction_final_test <- adf_stationarity(
    financial_daily$SP500_X_DGS10_MODEL
  )
}


# 4. NOW create the final interaction object
daily_financial_interaction <- financial_daily |>
  dplyr::transmute(
    date,
    SP500_DAILY_LOG_RETURN,
    DGS10_DAILY_CHANGE,
    SP500_X_DGS10_LEVEL = SP500_X_DGS10_DAILY,
    SP500_X_DGS10 = SP500_X_DGS10_MODEL
  )

# Compatibility objects retained for any downstream diagnostic code.
sp500_daily <- financial_daily |>
  dplyr::transmute(
    date,
    SP500,
    SP500_DAILY_LOG_RETURN,
    SP500_STATIONARY = SP500_DAILY_LOG_RETURN
  )

treasury10_daily <- financial_daily |>
  dplyr::transmute(
    date,
    DGS10,
    DGS10_STATIONARY = DGS10_DAILY_LOG_CHANGE,
    DGS10_DAILY_CHANGE
  )

daily_financial_interaction <- financial_daily |>
  dplyr::transmute(
    date,
    SP500_DAILY_LOG_RETURN,
    DGS10_DAILY_CHANGE,
    SP500_X_DGS10_LEVEL = SP500_X_DGS10_DAILY,
    SP500_X_DGS10 = SP500_X_DGS10_MODEL
  )


# ------------------------------------------------------------------------------
# 6F. STATIONARITY CHECKS AFTER LOCF AND AFTER INTERACTION CONSTRUCTION
# ------------------------------------------------------------------------------

sp500_test <- adf_stationarity(
  financial_daily$SP500_DAILY_LOG_RETURN
)

dgs10_test <- adf_stationarity(
  financial_daily$DGS10_DAILY_LOG_CHANGE
)


# ------------------------------------------------------------------------------
# INTERACTION STATIONARITY RULE USED BY THE FORECAST MODELS
# ------------------------------------------------------------------------------


interaction_level_test <- adf_stationarity(
  financial_daily$SP500_X_DGS10_DAILY
)

if (isTRUE(interaction_level_test$stationary)) {

  interaction_transformation_used <- "Level interaction"
  interaction_was_differenced <- FALSE

  financial_daily <- financial_daily |>
    dplyr::mutate(
      SP500_X_DGS10_FIRST_DIFF = NA_real_,
      SP500_X_DGS10_MODEL = SP500_X_DGS10_DAILY
    )

  interaction_final_test <- interaction_level_test

} else {

  interaction_transformation_used <- "First difference of interaction"
  interaction_was_differenced <- TRUE

  financial_daily <- financial_daily |>
    dplyr::mutate(
      SP500_X_DGS10_FIRST_DIFF =
        SP500_X_DGS10_DAILY -
        dplyr::lag(SP500_X_DGS10_DAILY),

      SP500_X_DGS10_MODEL =
        SP500_X_DGS10_FIRST_DIFF
    )

  interaction_final_test <- adf_stationarity(
    financial_daily$SP500_X_DGS10_MODEL
  )
}





financial_stationarity_diagnostics <- dplyr::bind_rows(

  tibble::tibble(
    series = "SP500_DAILY_LOG_RETURN",
    frequency = "daily",
    transformation =
      "100 x first difference of log after raw-level LOCF",
    adf_p_value = sp500_test$p_value,
    stationary_after_transform = sp500_test$stationary,
    used_in_model = TRUE
  ),

  tibble::tibble(
    series = "DGS10_DAILY_LOG_CHANGE",
    frequency = "daily",
    transformation =
      "100 x first difference of log after raw-level LOCF",
    adf_p_value = dgs10_test$p_value,
    stationary_after_transform = dgs10_test$stationary,
    used_in_model = TRUE
  ),

  tibble::tibble(
    series = "SP500_X_DGS10_DAILY_LEVEL",
    frequency = "daily",
    transformation =
      "SP500 daily log return x daily level change in DGS10",
    adf_p_value = interaction_level_test$p_value,
    stationary_after_transform =
      interaction_level_test$stationary,
    used_in_model =
      !interaction_was_differenced
  ),

  tibble::tibble(
    series = "SP500_X_DGS10_MODEL",
    frequency = "daily",
    transformation =
      interaction_transformation_used,
    adf_p_value = interaction_final_test$p_value,
    stationary_after_transform =
      interaction_final_test$stationary,
    used_in_model = TRUE
  )
)


interaction_stationarity_decision <- tibble::tibble(
  interaction_level_adf_p_value =
    interaction_level_test$p_value,

  interaction_level_stationary =
    interaction_level_test$stationary,

  first_difference_applied =
    interaction_was_differenced,

  transformation_used_in_EN5_and_QEN5 =
    interaction_transformation_used,

  final_interaction_adf_p_value =
    interaction_final_test$p_value,

  final_interaction_stationary =
    interaction_final_test$stationary
)


print(financial_stationarity_diagnostics)

readr::write_csv(
  financial_stationarity_diagnostics,
  file.path(
    output_dir,
    "daily_financial_stationarity_after_LOCF.csv"
  )
)

readr::write_csv(
  interaction_stationarity_decision,
  file.path(
    output_dir,
    "interaction_stationarity_decision.csv"
  )
)


# Save the completed daily data used to construct the actual model predictors.
readr::write_csv(
  financial_daily,
  file.path(output_dir, "daily_financial_data_after_LOCF_and_transformations.csv")
)


# ==============================================================================
