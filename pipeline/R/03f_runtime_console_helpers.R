# =============================================================================
# 03f_runtime_console_helpers.R
# Runtime and console helpers
# =============================================================================

# -----------------------------------------------------------------------------
# Runtime helper
# -----------------------------------------------------------------------------
run_step <- function(step_name, expr) {
  message("\n==================================================")
  message("START: ", step_name)
  message("TIME : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  
  t0 <- Sys.time()
  result <- force(expr)
  t1 <- Sys.time()
  
  elapsed <- round(as.numeric(difftime(t1, t0, units = "secs")), 2)
  
  message("DONE : ", step_name)
  message("ELAPSED: ", elapsed, " sec")
  message("==================================================")
  
  invisible(result)
}

# -----------------------------------------------------------------------------
# Console helpers
# -----------------------------------------------------------------------------
fmt_time_sec <- function(t0, t1 = Sys.time()) {
  round(as.numeric(difftime(t1, t0, units = "secs")), 2)
}

runtime_profile_reset <- function(output_dir = NULL) {
  assign(
    ".runtime_profile",
    list(
      rows = list(),
      output_dir = output_dir,
      started_at = Sys.time()
    ),
    envir = .GlobalEnv
  )
  invisible(TRUE)
}

runtime_profile_exists <- function() {
  exists(".runtime_profile", envir = .GlobalEnv, inherits = FALSE)
}

runtime_profile_set_output_dir <- function(output_dir) {
  if (!runtime_profile_exists()) {
    runtime_profile_reset(output_dir)
  } else {
    profile <- get(".runtime_profile", envir = .GlobalEnv, inherits = FALSE)
    profile$output_dir <- output_dir
    assign(".runtime_profile", profile, envir = .GlobalEnv)
  }
  invisible(TRUE)
}

runtime_profile_record <- function(label,
                                   t0,
                                   t1 = Sys.time(),
                                   category = "step",
                                   status = "ok",
                                   detail = "") {
  if (is.null(t0) || !inherits(t0, "POSIXt")) {
    return(invisible(FALSE))
  }

  if (!runtime_profile_exists()) {
    runtime_profile_reset()
  }

  profile <- get(".runtime_profile", envir = .GlobalEnv, inherits = FALSE)
  elapsed_sec <- fmt_time_sec(t0, t1)

  row <- data.frame(
    category = as.character(category)[1],
    label = as.character(label)[1],
    status = as.character(status)[1],
    start_time = format(t0, "%Y-%m-%d %H:%M:%S"),
    end_time = format(t1, "%Y-%m-%d %H:%M:%S"),
    elapsed_sec = elapsed_sec,
    detail = as.character(detail)[1],
    stringsAsFactors = FALSE
  )

  profile$rows[[length(profile$rows) + 1]] <- row
  assign(".runtime_profile", profile, envir = .GlobalEnv)
  invisible(TRUE)
}

runtime_profile_data <- function() {
  if (!runtime_profile_exists()) {
    return(data.frame())
  }

  profile <- get(".runtime_profile", envir = .GlobalEnv, inherits = FALSE)
  if (length(profile$rows) == 0) {
    return(data.frame())
  }

  do.call(rbind, profile$rows)
}

runtime_profile_write <- function(path = NULL, top_n = 10) {
  df <- runtime_profile_data()
  if (nrow(df) == 0) {
    step_info("Runtime profile: no timing rows collected.")
    return(invisible(df))
  }

  if (is.null(path)) {
    profile <- get(".runtime_profile", envir = .GlobalEnv, inherits = FALSE)
    if (!is.null(profile$output_dir) && nzchar(profile$output_dir)) {
      path <- file.path(profile$output_dir, "global", "audits_global", "runtime_profile.csv")
    }
  }

  if (!is.null(path) && nzchar(path)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    if (exists("write_csv_safe", mode = "function")) {
      write_csv_safe(df, path)
    } else {
      write.csv(df, path, row.names = FALSE)
    }
    step_info("Runtime profile exported: ", path)
  }

  ord <- order(df$elapsed_sec, decreasing = TRUE, na.last = NA)
  top <- df[head(ord, min(top_n, length(ord))), , drop = FALSE]
  step_info("Slowest runtime blocks:")
  for (i in seq_len(nrow(top))) {
    step_info(
      sprintf(
        "  %02d. [%s] %s: %.2f sec",
        i,
        top$category[i],
        top$label[i],
        top$elapsed_sec[i]
      )
    )
  }

  invisible(df)
}

profile_expr <- function(label, expr, category = "block", detail = "") {
  t0 <- Sys.time()
  status <- "ok"
  on.exit({
    runtime_profile_record(
      label = label,
      t0 = t0,
      category = category,
      status = status,
      detail = detail
    )
  }, add = TRUE)

  tryCatch(
    force(expr),
    error = function(e) {
      status <<- "error"
      stop(e)
    }
  )
}

console_rule <- function(char = "=", n = 60) {
  paste(rep(char, n), collapse = "")
}

step_start <- function(step_no, step_total, title) {
  message(console_rule())
  message(sprintf("[STEP %02d/%02d] %s", step_no, step_total, title))
  message(console_rule())
  invisible(Sys.time())
}

step_info <- function(...) {
  message("[INFO] ", paste0(..., collapse = ""))
}

step_ok <- function(title, t0 = NULL) {
  if (is.null(t0)) {
    message("[OK] ", title)
  } else {
    elapsed <- fmt_time_sec(t0)
    runtime_profile_record(
      label = title,
      t0 = t0,
      category = "step",
      status = "ok"
    )
    message("[OK] ", title, " finished in ", elapsed, " sec")
  }
}

step_warn <- function(...) {
  message("[WARN] ", paste0(..., collapse = ""))
}

step_fail <- function(...) {
  message("[ERROR] ", paste0(..., collapse = ""))
}

# Replace or append a key assignment in a settings/config text blob.
replace_or_append <- function(text, key, value_expr) {
  pattern <- paste0("^\\s*", key, "\\s*<-")
  replacement <- paste0(key, " <- ", value_expr)
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]

  idx <- grep(pattern, lines)
  if (length(idx) > 0) {
    lines[idx[1]] <- replacement
  } else {
    lines <- c(lines, replacement)
  }

  paste(lines, collapse = "\n")
}
