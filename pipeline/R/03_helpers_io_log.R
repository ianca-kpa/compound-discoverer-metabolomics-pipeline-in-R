# =============================================================================
# 03_helpers_io_log.R
# Compatibility loader for helper modules
# Keep this file as the public source() entry point for app, scripts, and run_pipeline.R.
# =============================================================================

helper_module_files <- c(
  "pipeline/R/03a_logging_helpers.R",
  "pipeline/R/03b_io_reference_helpers.R",
  "pipeline/R/03c_text_config_helpers.R",
  "pipeline/R/03d_setting_value_helpers.R",
  "pipeline/R/03e_summary_path_helpers.R",
  "pipeline/R/03f_runtime_console_helpers.R"
)

for (helper_module_file in helper_module_files) {
  source(helper_module_file, local = .GlobalEnv)
}

rm(helper_module_file, helper_module_files)
