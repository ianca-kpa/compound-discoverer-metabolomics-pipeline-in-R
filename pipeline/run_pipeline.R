# =============================================================================
# run_pipeline.R
# Entry point for the untargeted metabolomics pipeline
# =============================================================================

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(scipen = 999)

cat("\n==================================================\n")
cat("Untargeted metabolomics pipeline\n")
cat("Runner started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==================================================\n")

tryCatch(
  {
    locate_project_root <- function(start_paths) {
      for (start_path in start_paths) {
        if (is.null(start_path) || !nzchar(start_path)) {
          next
        }

        current_path <- normalizePath(start_path, winslash = "/", mustWork = FALSE)

        repeat {
          settings_path <- file.path(current_path, "pipeline", "config", "settings.R")
          if (file.exists(settings_path)) {
            return(current_path)
          }

          parent_path <- normalizePath(file.path(current_path, ".."), winslash = "/", mustWork = FALSE)
          if (identical(parent_path, current_path)) {
            break
          }

          current_path <- parent_path
        }
      }

      stop(
        "Could not locate the project root. Open the repository root in RStudio, or make sure pipeline/config/settings.R exists."
      )
    }

    this_file <- if (!is.null(sys.frames()[[1]]$ofile)) sys.frames()[[1]]$ofile else NA_character_
    script_dir <- if (!is.na(this_file) && nzchar(this_file)) dirname(this_file) else NA_character_
    project_dir <- locate_project_root(c(script_dir, getwd()))
    setwd(project_dir)

    cat("\nCurrent working directory:\n")
    cat(getwd(), "\n")

    if (!file.exists("pipeline/config/settings.R")) {
      stop("pipeline/config/settings.R not found. Copy pipeline/config/settings.example.R to pipeline/config/settings.R and edit it before running the pipeline.")
    }

    cat("\n[1/4] Loading settings.R ...\n")
    source("pipeline/config/settings.R", local = .GlobalEnv)
    cat("settings.R loaded successfully.\n")

    cat("\n[2/4] Loading pipeline modules ...\n")

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
      cat("Loading:", f, "\n")
      source(f, local = .GlobalEnv)
    }

    cat("\nAll modules loaded successfully.\n")

    cat("\n[3/4] Validating configuration ...\n")
    validate_settings()
    cat("Configuration OK.\n")

    cat("[4/4] Running pipeline ...\n")
    pipeline_result <- run_untargeted_pipeline()

    cat("Pipeline finished successfully.\n")
    cat("Finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "")
    cat("\n==================================================\n")

    if (is.list(pipeline_result)) {
      cat("\nReturned objects:\n")
      cat(paste(names(pipeline_result), collapse = "\n"), "\n")
    }

    if (exists("output_dir")) {
      cat("\nOutput directory:\n")
      cat(output_dir, "\n")
    }

    if (exists("pipeline_result") && is.list(pipeline_result) && "log_path" %in% names(pipeline_result)) {
      cat("\nLog file:\n")
      cat(pipeline_result$log_path, "\n")
    }
  },
  error = function(e) {
    cat("\n==================================================\n")
    cat("PIPELINE FAILED\n")
    cat("Message:", conditionMessage(e), "\n")
    cat("Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    cat("\n==================================================\n")
    stop(e)
  }
)
