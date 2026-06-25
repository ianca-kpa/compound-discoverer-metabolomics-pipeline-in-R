# Helpers for Results Gallery filtering and Pipeline Log summaries.

normalize_result_gallery_filter <- function(filter) {
  filter <- safe_trimws(filter)
  valid_filters <- c("all", "pca", "volcano", "heatmap", "qc_norm", "other")

  if (!nzchar(filter) || !filter %in% valid_filters) {
    return("all")
  }

  filter
}

classify_result_image_path <- function(path, out_dir) {
  rel <- tolower(rel_path_from_output(path, out_dir))

  if (grepl("pca", rel, fixed = TRUE)) {
    return("pca")
  }
  if (grepl("volcano", rel, fixed = TRUE)) {
    return("volcano")
  }
  if (grepl("heatmap", rel, fixed = TRUE)) {
    return("heatmap")
  }
  if (grepl("qc|qcrsc|loess|normalization|drift|rsd", rel, perl = TRUE)) {
    return("qc_norm")
  }

  "other"
}

filter_result_image_files <- function(files, mtime, out_dir, filter = "all") {
  filter <- normalize_result_gallery_filter(filter)

  if (identical(filter, "all") || length(files) == 0) {
    return(list(files = files, mtime = mtime))
  }

  categories <- vapply(files, classify_result_image_path, character(1), out_dir = out_dir)
  keep <- categories == filter

  list(files = files[keep], mtime = mtime[keep])
}

build_pipeline_log_summary_ui <- function(log_text) {
  if (is.null(log_text) || !nzchar(log_text) || identical(log_text, "No run executed yet.")) {
    return(tags$div(
      class = "pipeline-log-summary pipeline-log-summary-idle",
      tags$strong("Log status"),
      tags$span("No run executed yet.")
    ))
  }

  lines <- unlist(strsplit(log_text, "\\r?\\n", perl = TRUE), use.names = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]

  error_lines <- grep("\\b(error|failed|failure|fatal|cannot|did not start)\\b", lines, ignore.case = TRUE, value = TRUE)
  warning_lines <- grep("\\b(warning|warn|missing|skipped)\\b", lines, ignore.case = TRUE, value = TRUE)
  warning_lines <- setdiff(warning_lines, error_lines)

  if (length(error_lines) > 0) {
    examples <- tail(error_lines, 3)
    return(tags$div(
      class = "pipeline-log-summary pipeline-log-summary-error",
      tags$strong(paste(length(error_lines), "error-related log line(s) found")),
      tags$ul(lapply(examples, tags$li))
    ))
  }

  if (length(warning_lines) > 0) {
    examples <- tail(warning_lines, 3)
    return(tags$div(
      class = "pipeline-log-summary pipeline-log-summary-warning",
      tags$strong(paste(length(warning_lines), "warning/notice log line(s) found")),
      tags$ul(lapply(examples, tags$li))
    ))
  }

  tags$div(
    class = "pipeline-log-summary pipeline-log-summary-ok",
    tags$strong("Log status"),
    tags$span("No obvious errors or warnings detected in the current log preview.")
  )
}
