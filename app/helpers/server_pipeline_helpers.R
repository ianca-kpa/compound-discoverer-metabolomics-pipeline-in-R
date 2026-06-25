# Helpers for launching and stopping the pipeline process from the Shiny app.

pipeline_rscript_path <- function() {
  rscript_cmd <- file.path(R.home("bin"), "Rscript")
  if (.Platform$OS.type == "windows") {
    rscript_cmd <- paste0(rscript_cmd, ".exe")
  }

  rscript_cmd
}

launch_pipeline_process <- function(project_root, output_log_file, status_message, pipeline_log_text, process_state) {
  rscript_cmd <- pipeline_rscript_path()

  if (!file.exists(rscript_cmd)) {
    status_message("R script executable not found. Cannot run pipeline from UI.")
    return(invisible(FALSE))
  }

  log_file <- tempfile(pattern = "pipeline_run_", fileext = ".log")
  pipeline_log_text(paste0(
    "Starting pipeline at ",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "...\nWaiting for Rscript output."
  ))

  if (requireNamespace("processx", quietly = TRUE)) {
    proc <- tryCatch(
      processx::process$new(
        command = rscript_cmd,
        args = c("pipeline/run_pipeline.R"),
        wd = project_root,
        stdout = log_file,
        stderr = log_file,
        cleanup = TRUE
      ),
      error = function(e) e
    )

    if (inherits(proc, "error")) {
      status_message(paste("Failed to start pipeline:", conditionMessage(proc)))
      return(invisible(FALSE))
    }

    process_state$proc <- proc
    process_state$log_file <- log_file
    process_state$pipeline_log_file <- output_log_file
    process_state$running <- TRUE
    process_state$started_at <- Sys.time()
    process_state$pid <- tryCatch(proc$get_pid(), error = function(e) NULL)
    status_message("Pipeline is running in background. Open 'Pipeline Log' tab for live output.")
    return(invisible(TRUE))
  }

  status_message("Package 'processx' not installed. Running pipeline synchronously without live streaming.")

  exit_code <- tryCatch(
    system2(
      rscript_cmd,
      args = c("pipeline/run_pipeline.R"),
      stdout = log_file,
      stderr = log_file,
      wait = TRUE
    ),
    error = function(e) {
      writeLines(paste("Pipeline execution failed:", conditionMessage(e)), log_file)
      1L
    }
  )

  if (file.exists(output_log_file)) {
    pipeline_log_text(read_log_preview(output_log_file))
  } else if (file.exists(log_file)) {
    pipeline_log_text(read_log_preview(log_file))
  } else {
    pipeline_log_text("No log file generated.")
  }

  if (isTRUE(as.integer(exit_code) == 0L)) {
    status_message("Pipeline run completed successfully.")
  } else {
    status_message(
      paste(
        "Pipeline run failed with exit code",
        as.integer(exit_code),
        ". Check Pipeline Log tab."
      )
    )
  }

  invisible(isTRUE(as.integer(exit_code) == 0L))
}

stop_pipeline_process <- function(process_state, status_message) {
  if (!isTRUE(process_state$running) || is.null(process_state$proc)) {
    status_message("No running pipeline process to stop.")
    return(invisible(FALSE))
  }

  if (process_state$proc$is_alive()) {
    process_state$proc$kill()
    status_message("Pipeline process stopped by user.")
  } else {
    status_message("Pipeline process is not running.")
  }

  process_state$running <- FALSE
  invisible(TRUE)
}
