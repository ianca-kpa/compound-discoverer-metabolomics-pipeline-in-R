# =============================================================================
# run_pipeline.R
# Entry point for the untargeted metabolomics pipeline
# =============================================================================

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(scipen = 999)

cat("
==================================================
")
cat("Untargeted metabolomics pipeline
")
cat("Runner started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "
")
cat("==================================================
")

tryCatch(
  {
    this_file <- if (!is.null(sys.frames()[[1]]$ofile)) sys.frames()[[1]]$ofile else "pipeline/run_pipeline.R"
    script_dir <- normalizePath(dirname(this_file), winslash = "/", mustWork = FALSE)
    project_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
    if (!dir.exists(project_dir)) project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    setwd(project_dir)

    cat("
Current working directory:
")
    cat(getwd(), "
")

        if (!file.exists("pipeline/config/settings.R")) {
            stop("pipeline/config/settings.R not found. Copy pipeline/config/settings.example.R to pipeline/config/settings.R and edit it before running the pipeline.")
    }

    cat("
[1/4] Loading settings.R ...
")
    source("pipeline/config/settings.R", local = .GlobalEnv)
    cat("settings.R loaded successfully.
")

    cat("
[2/4] Loading pipeline modules ...
")

    module_files <- c(
            "pipeline/R/00_packages.R",
            "pipeline/R/01_validation.R",
            "pipeline/R/02_comparisons.R",
            "pipeline/R/03_helpers_io_log.R",
            "pipeline/R/04_metadata.R",
            "pipeline/R/05_features_assay.R",
            "pipeline/R/06_normalization_filters.R",
            "pipeline/R/07_duplicates.R",
            "pipeline/R/08_exports.R",
            "pipeline/R/09_pca.R",
            "pipeline/R/10_stats_volcano.R",
            "pipeline/R/11_heatmaps.R",
            "pipeline/R/12_main_pipeline.R"
    )

    for (f in module_files) {
      cat("Loading:", f, "
")
      source(f, local = .GlobalEnv)
    }

    cat("
All modules loaded successfully.
")

    cat("
[3/4] Validating configuration ...
")
    validate_settings()
    cat("Configuration OK.
")

    cat("[4/4] Running pipeline ...
")
    pipeline_result <- run_untargeted_pipeline()

#     cat("
# ==================================================
# ")
    cat("Pipeline finished successfully.
")
    cat("Finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "")
    cat("
==================================================
")

    if (is.list(pipeline_result)) {
      cat("
Returned objects:
")
      cat(paste(names(pipeline_result), collapse = "
"), "
")
    }

    if (exists("output_dir")) {
      cat("
Output directory:
")
      cat(output_dir, "
")
    }

    if (exists("pipeline_result") && is.list(pipeline_result) && "log_path" %in% names(pipeline_result)) {
      cat("
Log file:
")
      cat(pipeline_result$log_path, "
")
    }
  },
  error = function(e) {
    cat("
==================================================
")
    cat("PIPELINE FAILED
")
    cat("Message:", conditionMessage(e), "
")
    cat("Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "
")
    cat("
==================================================
")
    stop(e)
  }
)
