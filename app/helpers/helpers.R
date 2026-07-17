# App-only helpers for reading config text and building Settings Builder input IDs.

safe_read_file <- function(path) {
  if (!file.exists(path)) {
    return(paste("File not found:", path))
  }
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

extract_config_value <- function(config_text, key) {
  pattern <- paste0("^\\s*", key, "\\s*<-\\s*(.+)$")
  lines <- strsplit(config_text, "\\n", fixed = FALSE)[[1]]
  hit <- grep(pattern, lines, perl = TRUE, value = TRUE)

  if (length(hit) == 0) {
    return(NULL)
  }

  raw <- sub(pattern, "\\1", hit[length(hit)], perl = TRUE)
  raw <- sub("\\s*#.*$", "", trimws(raw))

  while (grepl("^\".*\"$|^'.*'$|^“.*”$|^‘.*’$", raw)) {
    raw <- gsub("^\"|\"$|^'|'$|^“|”$|^‘|’$", "", raw)
  }

  raw
}

is_absolute_path <- function(path) {
  grepl("^[A-Za-z]:[/\\]|^/|^~", path)
}

setting_input_id <- function(key) {
  paste0("settings_", key)
}

setting_has_value <- function(value) {
  !is.null(value) && length(value) > 0 && !all(is.na(value)) && nzchar(trimws(as.character(value)[1]))
}

setting_default_assignment <- function(key, value, type = "text") {
  encoded <- switch(type,
    checkbox = setting_value_logical(value),
    logical_select = setting_value_logical(value),
    numeric = setting_value_numeric(value, default = 0),
    numeric_or_inf = setting_value_numeric(value, default = Inf),
    integer = setting_value_integer(value, default = 0L),
    multiselect = setting_value_vector_text(value),
    vector_numeric = setting_value_vector_numeric(value),
    vector_text = setting_value_vector_text(value),
    detected_multiselect = setting_value_vector_text(value),
    nullable_vector_text = setting_value_vector_text(value, allow_null = TRUE),
    sheet = setting_value_sheet(value),
    setting_value_text(value)
  )

  paste0(key, " <- ", encoded)
}

build_default_config_from_global <- function() {
  if (!exists("settings_form_sections", inherits = TRUE)) {
    return("# App defaults are not available yet.")
  }

  lines <- c(
    "# Generated from app/global.R defaults.",
    "# This is the app startup state; settings.R is only updated when you save or run.",
    ""
  )

  for (section in settings_form_sections) {
    lines <- c(lines, paste0("# ", section$title))
    for (spec in section$fields) {
      lines <- c(lines, setting_default_assignment(spec$key, spec$default, spec$type))
    }
    lines <- c(lines, "")
  }

  paste(lines, collapse = "\n")
}

read_initial_config <- function() {
  build_default_config_from_global()
}
