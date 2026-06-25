# Helpers for uploaded files, config paths, and output directory resolution.

resolve_config_input_path <- function(uploaded, external_path, kind) {
  if (!is.null(uploaded)) {
    return(list(path = file.path("data", basename(uploaded$name)), missing = FALSE))
  }

  ext <- safe_trimws(external_path)
  if (nzchar(ext)) {
    return(list(path = ext, missing = FALSE))
  }

  list(
    path = NULL,
    missing = TRUE,
    message = paste(
      "No",
      kind,
      "path was provided. Keep existing config value or provide one."
    )
  )
}

copy_uploaded_file_to_data <- function(uploaded, project_root) {
  if (is.null(uploaded)) {
    return(invisible(FALSE))
  }

  dir.create(file.path(project_root, "data"), recursive = TRUE, showWarnings = FALSE)
  file.copy(
    uploaded$datapath,
    file.path(project_root, "data", basename(uploaded$name)),
    overwrite = TRUE
  )

  invisible(TRUE)
}

persist_injection_order_upload_to_data <- function(uploaded, project_root) {
  if (is.null(uploaded)) {
    return(invisible(FALSE))
  }

  dir.create(file.path(project_root, "data"), recursive = TRUE, showWarnings = FALSE)
  injection_order_canonical_path <- file.path(project_root, "data", "Input Files.xlsx")
  file.copy(uploaded$datapath, injection_order_canonical_path, overwrite = TRUE)

  injection_order_original_path <- file.path(project_root, "data", basename(uploaded$name))
  canonical_norm <- tolower(normalizePath(injection_order_canonical_path, winslash = "/", mustWork = FALSE))
  original_norm <- tolower(normalizePath(injection_order_original_path, winslash = "/", mustWork = FALSE))
  if (!identical(original_norm, canonical_norm)) {
    file.copy(uploaded$datapath, injection_order_original_path, overwrite = TRUE)
  }

  invisible(TRUE)
}

resolve_output_dir_abs_value <- function(output_dir, project_root) {
  out <- strip_outer_quotes(output_dir)
  if (!nzchar(out)) {
    out <- "output"
  }

  if (is_absolute_path(out)) {
    return(normalizePath(out, winslash = "/", mustWork = FALSE))
  }

  normalizePath(file.path(project_root, out), winslash = "/", mustWork = FALSE)
}

is_current_file_path <- function(file_path, clear_ts = NULL, session_started_at = NULL) {
  if (is.null(file_path) || !file.exists(file_path)) {
    return(FALSE)
  }

  file_mtime <- file.mtime(file_path)

  if (!is.null(clear_ts)) {
    return(file_mtime >= (clear_ts - 1))
  }

  if (!is.null(session_started_at)) {
    return(file_mtime >= (session_started_at - 2))
  }

  TRUE
}
