# =============================================================================
# 03a_logging_helpers.R
# Logging helpers
# =============================================================================

# =============================================================================
# 03_helpers_io_log.R
# General helpers: logging / IO / text / directories
# =============================================================================

# -----------------------------------------------------------------------------
# Logging helpers
# -----------------------------------------------------------------------------
append_log_line <- function(log_path, line) {
  write(line, file = log_path, append = TRUE)
}

# Format dimensions of a data frame, matrix, or vector for logging
fmt_dims <- function(x) {
  if (is.null(x)) return("rows=NA cols=NA")
  if (is.data.frame(x)) return(paste0("rows=", nrow(x), " cols=", ncol(x)))
  if (is.matrix(x)) return(paste0("rows=", nrow(x), " cols=", ncol(x)))
  if (is.vector(x)) return(paste0("len=", length(x)))
  
  tryCatch(
    paste0("rows=", nrow(x), " cols=", ncol(x)),
    error = function(e) "rows=NA cols=NA"
  )
}

# Log a message about a written object, including its dimensions and file path
log_written_object <- function(log_path, file_path, object, note = NULL) {
  msg <- paste0("- ", basename(file_path), " -> ", fmt_dims(object), " | path: ", file_path)
  
  if (!is.null(note) && nzchar(note)) {
    msg <- paste0(msg, " | note: ", note)
  }
  
  append_log_line(log_path, msg)
}

# Run an expression while capturing all console output (messages, warnings, errors) to a log file
# Capture console output into the pipeline log while still mirroring it to the console.
with_console_capture_to_file <- function(log_path, expr) {
  con <- file(log_path, open = "a", encoding = "UTF-8")
  
  on.exit({
    while (sink.number() > 0) {
      try(sink(), silent = TRUE)
    }
    try(close(con), silent = TRUE)
  }, add = TRUE)
  
  sink(con, append = TRUE, split = TRUE)
  
  tryCatch(
    withCallingHandlers(
      expr,
      message = function(m) {
        # Emit messages through stdout so sink(split=TRUE) mirrors to console + log once.
        cat(conditionMessage(m), sep = "")
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        cat("[WARNING] ", conditionMessage(w), sep = "")
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      cat("[ERROR] ", conditionMessage(e), sep = "")
      stop(e)
    }
  )
}
