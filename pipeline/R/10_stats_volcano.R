# =============================================================================
# 10_stats_volcano.R
# Compatibility loader for statistics, volcano, and export modules
# Keep this file as the public source() entry point for app, scripts, and run_pipeline.R.
# MULTIGROUP_GLOBAL is non-directional by contract: its FC/log2FC values remain
# NA and the volcano module must only receive pairwise comparison results.
# =============================================================================

stats_module_files <- c(
  "pipeline/R/10a_stats_core.R",
  "pipeline/R/10b_volcano_plots.R",
  "pipeline/R/10c_stats_driver.R",
  "pipeline/R/10d_stats_exports.R"
)

for (stats_module_file in stats_module_files) {
  source(stats_module_file, local = .GlobalEnv)
}

rm(stats_module_file, stats_module_files)
