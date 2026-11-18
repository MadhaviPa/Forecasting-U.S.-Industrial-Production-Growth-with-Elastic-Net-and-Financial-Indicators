# Run all stages in one R session, preserving intermediate objects.
# Start from the repository root: Rscript run_all.R
steps <- c(
  "01_download_data.R",
  "02_preprocessing.R",
  "03_feature_engineering.R",
  "04_elastic_net.R",
  "05_quantile_elastic_net.R",
  "06_evaluation.R"
)
for (step in steps) {
  message("\n=== Running ", step, " ===")
  source(file.path("code", step), local = FALSE, encoding = "UTF-8")
}
message("All analysis stages finished.")
