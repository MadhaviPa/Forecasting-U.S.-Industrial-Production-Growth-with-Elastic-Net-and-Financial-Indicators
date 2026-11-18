# ============================================================================
# 01_download_data.R
# Derived from the uploaded Madhavi_panchangum_fred.R; original methods preserved.
# Run from the repository root via Rscript run_all.R.
# ============================================================================





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
if (!file.exists(file_path)) stop("FRED-MD CSV not found: ", file_path, "
Download the FRED-MD vintage and save it under data/raw/ or set FRED_MD_FILE.")


# Financial series begin in 1992 in the original setup.
financial_start_date <- as.Date("1990-01-01")
end_date <- as.Date(Sys.getenv("FORECAST_END_DATE", "2026-06-30"))



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
output_dir <- "results/details"
plot_dir <- file.path(output_dir, "plots")
error_plot_dir <- file.path(plot_dir, "quantile_error_distributions")

for (d in c(output_dir, plot_dir, error_plot_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

fred_key <- Sys.getenv("FRED_API_KEY", "")
if (!nzchar(fred_key)) stop("Set FRED_API_KEY in your environment before running the analysis.")
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
