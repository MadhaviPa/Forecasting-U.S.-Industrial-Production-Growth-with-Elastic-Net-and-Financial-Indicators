# Original uploaded analysis, with personal path and API key redacted.



# ==============================================================================
# 1. PACKAGES
# ==============================================================================

packages <- c(
  "tidyverse",
  "lubridate",
  "fredr",
  "quantmod",
  "zoo",
  "glmnet",
  "quantreg",
  "rqPen",
  "scales",
  "tseries"
)

new_packages <- setdiff(packages, rownames(installed.packages()))

if (length(new_packages) > 0L) {
  install.packages(new_packages, dependencies = TRUE)
}

invisible(lapply(packages, library, character.only = TRUE))


# ==============================================================================
# 2. USER SETTINGS
# ==============================================================================

# Change this path if Current.csv is stored elsewhere.
file_path <- Sys.getenv("FRED_MD_FILE", "data/raw/2026-07-MD.csv")


# Financial series begin in 1992 in the original setup.
financial_start_date <- as.Date("1990-01-01")
end_date <- Sys.Date()



# Forecast design from the report.
forecast_horizons <- c(1L, 3L, 6L, 12L, 24L)
rolling_windows <- c(120L, 180L, 240L)
quantile_levels <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
alpha_grid <- seq(0.1, 0.9, by = 0.1)

# Chronological inner tuning split.

inner_validation_share <- 0.20
minimum_inner_validation <- 18L

# Number of lambda values used when glmnet / rqPen automatically builds a path.
point_nlambda <- 50L
quantile_nlambda <- 30L

pca_em_factors <- 8L
pca_em_max_iter <- 200L
pca_em_tolerance <- 1e-6
pca_em_min_observed <- 12L


retune_every <- 12L

use_common_oos_start <- TRUE

selection_tolerance <- 1e-8

set.seed(123)

# Output folders.
output_dir <- "clarkwest"
plot_dir <- file.path(output_dir, "plots")
error_plot_dir <- file.path(plot_dir, "quantile_error_distributions")

for (d in c(output_dir, plot_dir, error_plot_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

fred_key <- Sys.getenv("FRED_API_KEY")

fredr::fredr_set_key(fred_key)


# ==============================================================================
# 3. LOAD FRED-MD DATA AND TRANSFORMATION CODES
# ==============================================================================

fred_md_raw <- readr::read_csv(
  file_path,
  na = c("", "NA", "N/A", ".", "null", "NULL"),
  show_col_types = FALSE
)

# ------------------------------------------------------------------------------
# The first row of the FRED-MD CSV contains the recommended transformation code
# for each economic series.
# ------------------------------------------------------------------------------

transformation_codes <- fred_md_raw |>
  dplyr::slice(1) |>
  dplyr::select(-sasdate) |>
  unlist(use.names = TRUE) |>
  as.numeric()

names(transformation_codes) <- names(
  fred_md_raw |> dplyr::select(-sasdate)
)

# ------------------------------------------------------------------------------
# Remove the transformation-code row and prepare the actual monthly observations.
# ------------------------------------------------------------------------------

raw_data <- fred_md_raw |>
  dplyr::filter(sasdate != "Transform:") |>
  dplyr::mutate(
    sasdate = lubridate::mdy(sasdate)
  ) |>
  dplyr::filter(!is.na(sasdate)) |>
  dplyr::arrange(sasdate)

cat("\nFRED-MD data\n")
cat("Rows:", nrow(raw_data), "\n")
cat("Variables:", ncol(raw_data) - 1L, "\n")
cat(
  "Period:",
  format(min(raw_data$sasdate)),
  "to",
  format(max(raw_data$sasdate)),
  "\n\n"
)


# ==============================================================================
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
# 6. EXTERNAL DAILY FINANCIAL DATA: MISSING VALUES -> LOCF -> TRANSFORMATIONS
# ==============================================================================

# The financial predictors remain DAILY. They are NOT summed or averaged to one
# monthly number. Non-trading days are explicitly identified in a date-audit table. The S&P 500 Yahoo
# date index defines actual trading dates; weekends and exchange holidays are
# classified as non-trading, not counted as missing, and never filled.
# Section 7 reshapes each month's ACTUAL trading-day observations into D01...Dmax
# predictor columns, where Dmax is determined from the observed data.

# ------------------------------------------------------------------------------
# 6A. DOWNLOAD RAW S&P 500 AND 10-YEAR TREASURY YIELD
# ------------------------------------------------------------------------------

sp500_xts <- quantmod::getSymbols(
  "^GSPC",
  src = "yahoo",
  from = financial_start_date,
  to = end_date,
  auto.assign = FALSE
)

sp500_raw_daily <- tibble::tibble(
  date = as.Date(zoo::index(sp500_xts)),
  SP500 = as.numeric(quantmod::Cl(sp500_xts))
) |>
  dplyr::arrange(date)

treasury10_raw_daily <- fredr::fredr(
  series_id = "DGS10",
  observation_start = financial_start_date,
  observation_end = end_date
) |>
  dplyr::transmute(
    date = as.Date(date),
    DGS10 = as.numeric(value)
  ) |>
  dplyr::arrange(date)


# ------------------------------------------------------------------------------
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
