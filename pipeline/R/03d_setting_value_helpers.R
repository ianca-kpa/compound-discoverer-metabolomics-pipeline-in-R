# =============================================================================
# 03d_setting_value_helpers.R
# Setting parsing and encoding helpers
# =============================================================================

# -----------------------------------------------------------------------------
# Helper: detect missing-like values
# -----------------------------------------------------------------------------
is_missing_like <- function(x) {
  xn <- toupper(trimws(as.character(x)))
  xn %in% c("", "NA", "N/A", "NULL")
}

# -----------------------------------------------------------------------------
# General helpers moved from app/global.R for sharing: trimming, parsing, encoding
# -----------------------------------------------------------------------------
safe_trimws <- function(value) {
  if (is.null(value) || length(value) == 0 || all(is.na(value))) {
    return("")
  }
  trimws(as.character(value)[1])
}

# Normalize the public output-level setting. Kept in the shared helper loader so
# it is available both to the pipeline and during Shiny server initialization.
normalize_output_level <- function(value = NULL, legacy_minimal = FALSE) {
  if (is.null(value) || length(value) == 0 || is.na(value[1]) || !nzchar(trimws(as.character(value)[1]))) {
    return(if (isTRUE(legacy_minimal)) "minimal" else "standard")
  }
  level <- tolower(trimws(as.character(value)[1]))
  level <- gsub("[- /]+", "_", level)
  if (level %in% c("full", "debug", "full_debug")) return("full_debug")
  if (level %in% c("minimal", "standard")) return(level)
  "standard"
}

derive_output_flags <- function(output_level) {
  output_level <- normalize_output_level(output_level)
  list(
    output_level = output_level,
    export_debug_outputs = identical(output_level, "full_debug"),
    export_intermediate_tables = identical(output_level, "full_debug"),
    export_all_plots = identical(output_level, "full_debug"),
    export_normalization_audit = output_level %in% c("standard", "full_debug"),
    export_qc_summary = output_level %in% c("standard", "full_debug"),
    export_multigroup_outputs = output_level %in% c("standard", "full_debug"),
    export_all_pairwise_multigroup = identical(output_level, "full_debug")
  )
}

get_optional_setting_value <- function(var_name, env = parent.frame()) {
  if (!exists(var_name, envir = env, inherits = TRUE)) {
    return(NULL)
  }

  value <- safe_trimws(get(var_name, envir = env, inherits = TRUE))
  if (!nzchar(value)) {
    return(NULL)
  }

  value
}

replace_missing_matrix_values <- function(mat, replacement = 0) {
  out <- as.matrix(mat)
  out[is.na(out)] <- replacement
  out
}

scale_matrix_columns_zscore <- function(mat) {
  out <- scale(mat, center = TRUE, scale = TRUE)
  replace_missing_matrix_values(out)
}

scale_matrix_columns_pareto <- function(mat) {
  mu <- colMeans(mat, na.rm = TRUE)
  sdv <- apply(mat, 2, stats::sd, na.rm = TRUE)
  denom <- sqrt(sdv)
  denom[is.na(denom) | denom == 0] <- 1

  out <- sweep(sweep(mat, 2, mu, "-"), 2, denom, "/")
  replace_missing_matrix_values(out)
}

safe_read_table <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("csv")) {
    return(utils::read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  if (ext %in% c("tsv", "txt")) {
    return(utils::read.delim(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Package 'readxl' is required to read Excel metadata files.")
    }
    return(as.data.frame(readxl::read_excel(path), stringsAsFactors = FALSE))
  }
  stop("Unsupported metadata extension: ", ext)
}

# Parsing helpers
parse_as_text <- function(value) {
  value <- trimws(as.character(value)[1])
  value <- sub("^c\\((.*)\\)$", "\\1", value, perl = TRUE)
  gsub("^\"|\"$|^'|'$", "", value)
}

parse_as_logical <- function(value) {
  tolower(trimws(as.character(value)[1])) %in% c("true", "t", "1", "yes")
}

parse_as_numeric <- function(value) {
  value <- trimws(as.character(value)[1])
  value <- sub("^c\\((.*)\\)$", "\\1", value, perl = TRUE)
  value <- strsplit(value, ",", fixed = TRUE)[[1]][1]
  suppressWarnings(as.numeric(trimws(value)))
}

parse_as_vector <- function(value) {
  if (length(value) == 1 && toupper(value) == "NULL") {
    return(character(0))
  }
  value <- sub("^c\\((.*)\\)$", "\\1", value, perl = TRUE)
  parts <- strsplit(value, ",", fixed = TRUE)[[1]]
  parts <- gsub("^\"|\"$|^'|'$", "", trimws(parts))
  parts[nzchar(parts)]
}

extract_and_parse_setting <- function(config_text, key, parser_fn, default = NULL) {
  value <- extract_config_value(config_text, key)
  if (!setting_has_value(value)) {
    return(default)
  }
  parser_fn(value)
}

# Setting display helpers
setting_display_value <- function(config_text, key, default = "") {
  extract_and_parse_setting(config_text, key, parse_as_text, default)
}

setting_display_logical <- function(config_text, key, default = FALSE) {
  extract_and_parse_setting(config_text, key, parse_as_logical, default)
}

setting_display_numeric <- function(config_text, key, default = NA_real_) {
  numeric_value <- extract_and_parse_setting(config_text, key, parse_as_numeric, NA_real_)
  if (length(numeric_value) != 1 || is.na(numeric_value)) default else numeric_value
}

setting_display_integer <- function(config_text, key, default = NA_integer_) {
  value <- setting_display_numeric(config_text, key, default = NA_real_)
  if (length(value) != 1 || is.na(value)) default else as.integer(round(value))
}

setting_display_vector <- function(config_text, key) {
  extract_and_parse_setting(config_text, key, parse_as_vector, character(0))
}

setting_display_csv <- function(config_text, key, default = "") {
  vec <- setting_display_vector(config_text, key)
  if (length(vec) == 0) default else paste(vec, collapse = ", ")
}

setting_display_sheet <- function(config_text, key, default = "") {
  value <- setting_display_value(config_text, key, default = default)
  if (!setting_has_value(value)) {
    return(default)
  }
  suppressWarnings({
    numeric_value <- as.numeric(value)
  })
  if (!is.na(numeric_value) && abs(numeric_value - round(numeric_value)) < .Machine$double.eps^0.5) {
    as.character(as.integer(round(numeric_value)))
  } else {
    value
  }
}

# Normalize vector defaults used by Shiny controls. These accept either an
# already parsed vector or an R-style value read from generated config text.
setting_default_vector <- function(value) {
  if (is.null(value) || length(value) == 0 || all(is.na(value))) {
    return(character(0))
  }

  if (length(value) > 1) {
    out <- trimws(as.character(value))
    return(unique(out[!is.na(out) & nzchar(out)]))
  }

  raw <- trimws(as.character(value)[1])
  if (!nzchar(raw) || toupper(raw) == "NULL" || raw %in% c("c()", "character(0)")) {
    return(character(0))
  }

  unique(parse_as_vector(raw))
}

setting_default_numeric_vector <- function(value) {
  parsed <- setting_default_vector(value)
  if (length(parsed) == 0) {
    return(numeric(0))
  }

  out <- suppressWarnings(as.numeric(parsed))
  unique(out[is.finite(out)])
}

# Encoding helpers
encode_r_string <- function(value) {
  value <- gsub("\\\\", "\\\\\\\\", as.character(value))
  value <- gsub("\"", "\\\\\"", value, fixed = TRUE)
  paste0("\"", value, "\"")
}

encode_text <- function(value, nullable = FALSE) {
  value <- safe_trimws(value)
  if (nullable && (!nzchar(value) || toupper(value) == "NULL")) {
    return("NULL")
  }
  if (!nzchar(value)) {
    return('""')
  }
  encode_r_string(value)
}

encode_logical <- function(value) {
  val <- toupper(safe_trimws(value))
  if (val %in% c("TRUE", "T", "1", "YES")) "TRUE" else "FALSE"
}

encode_numeric <- function(value, default = 0) {
  numeric_value <- suppressWarnings(as.numeric(value))
  if (length(numeric_value) == 0 || is.na(numeric_value)) {
    numeric_value <- default
  }
  format(numeric_value, scientific = FALSE, trim = TRUE)
}

encode_integer <- function(value, default = 0L) {
  integer_value <- suppressWarnings(as.integer(round(as.numeric(value))))
  if (length(integer_value) == 0 || is.na(integer_value)) {
    integer_value <- default
  }
  as.character(integer_value)
}

encode_sheet <- function(value) {
  value <- safe_trimws(value)
  if (!nzchar(value)) {
    return('""')
  }
  if (grepl("^-?[0-9]+(\\.[0]+)?$", value)) {
    as.character(as.integer(round(as.numeric(value))))
  } else {
    encode_r_string(value)
  }
}

encode_vector_text <- function(value, allow_null = FALSE) {
  items <- trimws(unlist(strsplit(as.character(value), ",", fixed = TRUE)))
  items <- items[nzchar(items)]
  if (length(items) == 0) {
    return(if (allow_null) "NULL" else "c()")
  }
  paste0("c(", paste(vapply(items, encode_r_string, character(1)), collapse = ", "), ")")
}

encode_vector_numeric <- function(value) {
  items <- trimws(unlist(strsplit(as.character(value), ",", fixed = TRUE)))
  items <- items[nzchar(items)]
  if (length(items) == 0) {
    return("c()")
  }
  numeric_items <- suppressWarnings(as.numeric(items))
  numeric_items <- numeric_items[!is.na(numeric_items)]
  if (length(numeric_items) == 0) {
    return("c()")
  }
  paste0("c(", paste(format(numeric_items, scientific = FALSE, trim = TRUE), collapse = ", "), ")")
}

# Setting value encoders
setting_value_text <- function(value) { encode_text(value, nullable = FALSE) }
setting_value_nullable_text <- function(value) { encode_text(value, nullable = TRUE) }
setting_value_logical <- function(value) { encode_logical(value) }
setting_value_numeric <- function(value, default = 0) { encode_numeric(value, default) }
setting_value_integer <- function(value, default = 0L) { encode_integer(value, default) }
setting_value_sheet <- function(value) { encode_sheet(value) }
setting_value_vector_text <- function(value, allow_null = FALSE) { encode_vector_text(value, allow_null) }
setting_value_vector_numeric <- function(value) { encode_vector_numeric(value) }

parse_num_robust <- function(x, decimal_mark = ".", grouping_mark = ",") {
  raw <- stringr::str_trim(as.character(x))
  
  suppressWarnings(
    readr::parse_number(
      raw,
      locale = readr::locale(decimal_mark = decimal_mark, grouping_mark = grouping_mark)
    )
  )
}
