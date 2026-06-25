# =============================================================================
# 06_normalization_filters.R
# Compatibility loader for normalization and filter modules
# Keep this file as the public source() entry point for app, scripts, and run_pipeline.R.
# =============================================================================

normalization_module_files <- c(
  "pipeline/R/06a_normalization_core.R",
  "pipeline/R/06b_normalization_plots.R",
  "pipeline/R/06c_filter_helpers.R"
)

for (normalization_module_file in normalization_module_files) {
  source(normalization_module_file, local = .GlobalEnv)
}

rm(normalization_module_file, normalization_module_files)
