# =============================================================================
# 04_metadata.R
# Metadata helpers
# =============================================================================

normalize_sex_value <- function(x) {
  s <- tolower(trimws(as.character(x)))
  s[s %in% c("", "na", "n/a", "nan", "unknown", "unk", "-", "null")] <- NA_character_

  dplyr::case_when(
    is.na(s) ~ NA_character_,
    s %in% c("f", "fem", "female") ~ "F",
    s %in% c("m", "mal", "male") ~ "M",
    TRUE ~ NA_character_
  )
}

clean_metadata <- function(metadata) {
  # Normalize column names: lowercase, trim, replace non-alphanumeric with underscores, remove leading/trailing underscores
  normalize_colname <- function(x) {
    x %>%
      tolower() %>%
      stringr::str_trim() %>%
      stringr::str_replace_all("[^a-z0-9]+", "_") %>%
      stringr::str_replace_all("^_+|_+$", "")
  }

  original_names <- names(metadata)
  clean_names <- normalize_colname(original_names)
  names(metadata) <- clean_names

  # Define expected columns and their possible aliases
  alias_map <- list(
    sample    = c("sample", "sample_id", "sample_name", "id_sample", "id", "name"),
    weight = c("weight_mg", "weight", "mass_mg", "mass", "mg", "sample_weight", "sample_mass", "weight_g", "mass_g", "sample_weight_g", "sample_mass_g"),
    group     = c("group", "treatment", "treat"),
    sex       = c("sex", "gender"),
    model     = c("model", "disease", "condition", "phenotype", "status"),
    injection_order = c("injection_order", "injectionorder", "injection", "injection_number", "injection_no", "run_order", "runorder", "order", "sequence_order", "acquisition_order")
  )

  # Function to standardize column names based on aliases
  standardize_column <- function(df, target, aliases) {
    found <- intersect(aliases, names(df))

    if (length(found) > 1 && !(target %in% found)) {
      stop(
        "Multiple possible columns found for '", target, "': ",
        paste(found, collapse = ", "),
        ". Please keep only one."
      )
    }

    if (length(found) >= 1) {
      chosen <- if (target %in% found) target else found[1]
      names(df)[names(df) == chosen] <- target
    }

    df
  }

  for (target in names(alias_map)) {
    metadata <- standardize_column(metadata, target, alias_map[[target]])
  }

  needed <- c("sample", "weight", "group", "sex", "model")
  miss <- setdiff(needed, names(metadata))

  if (length(miss) > 0) {
    stop("Metadata missing required columns: ", paste(miss, collapse = ", "))
  }

  metadata %>%
    dplyr::mutate(
      sample = stringr::str_trim(as.character(sample)),
      group = stringr::str_trim(as.character(group)),
      model = stringr::str_trim(as.character(model)),
      sex = normalize_sex_value(sex),
      weight = suppressWarnings(as.numeric(weight)),
      injection_order = if ("injection_order" %in% names(.)) {
        suppressWarnings(as.numeric(injection_order))
      } else {
        NA_real_
      },
      type = if (!("type" %in% names(.))) {
        dplyr::if_else(stringr::str_detect(sample, "^QC"), "QC", "Sample")
      } else {
        stringr::str_trim(as.character(type))
      }
    )
}

read_injection_order_file <- function(path, metadata_aligned) {
  path <- trimws(as.character(path)[1])
  if (!nzchar(path)) {
    return(NULL)
  }
  if (!file.exists(path)) {
    stop("Injection order file not found: ", path)
  }

  order_tbl <- safe_read_table(path)
  normalize_colname <- function(x) {
    x %>%
      tolower() %>%
      stringr::str_trim() %>%
      stringr::str_replace_all("[^a-z0-9]+", "_") %>%
      stringr::str_replace_all("^_+|_+$", "")
  }
  names(order_tbl) <- normalize_colname(names(order_tbl))

  normalize_sample_key <- function(x) {
    s <- trimws(as.character(x))
    s <- gsub("\\\\", "/", s)
    s <- sub("^.*/", "", s)
    s <- sub("(?i)\\.raw.*$", "", s, perl = TRUE)
    s <- sub("\\s*\\(.*\\)$", "", s)
    s <- trimws(s)
    s <- ifelse(grepl("^QC", s), s, sub("_.*$", "", s))
    tolower(s)
  }

  parse_creation_value <- function(x) {
    if (exists("parse_file_creation_value", mode = "function")) {
      return(parse_file_creation_value(x))
    }

    x <- trimws(as.character(x))
    if (length(x) == 0 || is.na(x) || !nzchar(x)) return(NA_real_)
    numeric_value <- suppressWarnings(as.numeric(x))
    if (is.finite(numeric_value)) return(numeric_value)

    formats <- c(
      "%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M",
      "%Y/%m/%d %H:%M:%S", "%Y/%m/%d %H:%M",
      "%d/%m/%Y %H:%M:%S", "%d/%m/%Y %H:%M",
      "%m/%d/%Y %H:%M:%S", "%m/%d/%Y %H:%M",
      "%Y-%m-%d", "%d/%m/%Y", "%m/%d/%Y"
    )
    for (fmt in formats) {
      parsed <- suppressWarnings(as.POSIXct(x, format = fmt, tz = "UTC"))
      if (!is.na(parsed)) return(as.numeric(parsed))
    }
    NA_real_
  }

  sample_aliases <- c(
    "sample", "samples", "sample_id", "sample_name", "id_sample", "id", "name",
    "file_name", "filename", "file", "raw_file", "input_file", "study_file"
  )
  order_aliases <- c(
    "injection_order", "injectionorder", "injection", "injection_number",
    "injection_no", "run_order", "runorder", "order", "sequence_order",
    "acquisition_order"
  )
  creation_aliases <- c(
    "creation_date", "creation_time", "created", "created_date", "created_time",
    "file_creation_date", "file_creation_time", "acquisition_date",
    "acquisition_time", "acquired", "acquired_date", "acquired_time",
    "injection_date", "injection_time", "run_date", "run_time"
  )

  sample_col <- intersect(sample_aliases, names(order_tbl))
  sample_col <- if (length(sample_col) > 0) sample_col[1] else NA_character_
  order_col <- intersect(order_aliases, names(order_tbl))
  order_col <- if (length(order_col) > 0) order_col[1] else NA_character_
  creation_col <- intersect(creation_aliases, names(order_tbl))
  creation_col <- if (length(creation_col) > 0) creation_col[1] else NA_character_

  if (is.na(order_col)) {
    if (!is.na(creation_col)) {
      if (is.na(sample_col)) {
        stop("Injection order file with creation dates must also contain a sample/file-name column.")
      }

      creation_values <- vapply(order_tbl[[creation_col]], parse_creation_value, numeric(1))
      sample_keys <- normalize_sample_key(order_tbl[[sample_col]])
      valid_rows <- is.finite(creation_values) & nzchar(sample_keys)

      if (sum(valid_rows) < 2) {
        stop("Injection order file creation-date column does not contain enough valid date/time values.")
      }

      ordered_rows <- which(valid_rows)[order(creation_values[valid_rows], which(valid_rows))]
      order_values <- rep(NA_real_, nrow(order_tbl))
      order_values[ordered_rows] <- seq_along(ordered_rows)
      order_col <- "__derived_from_creation_date__"
      order_tbl[[order_col]] <- order_values
    } else {
    numeric_cols <- names(order_tbl)[vapply(order_tbl, function(x) {
      any(is.finite(suppressWarnings(as.numeric(x))), na.rm = TRUE)
    }, logical(1))]

    if (length(numeric_cols) == 1) {
      order_col <- numeric_cols[1]
    } else {
      stop(
        "Injection order file must contain one order column such as ",
          "'injection_order', 'run_order', 'order', or 'acquisition_order', ",
          "or a creation-date column such as 'Creation Date'."
      )
    }
    }
  }

  order_values <- suppressWarnings(as.numeric(order_tbl[[order_col]]))
  fallback_order <- if ("injection_order" %in% names(metadata_aligned)) {
    suppressWarnings(as.numeric(metadata_aligned$injection_order))
  } else {
    seq_len(nrow(metadata_aligned))
  }
  fallback_order[!is.finite(fallback_order)] <- seq_len(nrow(metadata_aligned))[!is.finite(fallback_order)]

  if (!is.na(sample_col)) {
    sample_values <- as.character(order_tbl[[sample_col]])
    order_samples <- normalize_sample_key(sample_values)
    preferred_rows <- rep(TRUE, length(order_samples))
    non_ms2_rows <- !grepl("(?i)ddms2", sample_values, perl = TRUE)
    if (any(non_ms2_rows & is.finite(order_values) & nzchar(order_samples))) {
      preferred_rows <- non_ms2_rows
    }

    lookup_tbl <- data.frame(
      sample = order_samples[preferred_rows],
      order = order_values[preferred_rows],
      stringsAsFactors = FALSE
    )
    lookup_tbl <- lookup_tbl[is.finite(lookup_tbl$order) & nzchar(lookup_tbl$sample), , drop = FALSE]
    lookup_tbl <- lookup_tbl[order(lookup_tbl$order), , drop = FALSE]
    lookup_tbl <- lookup_tbl[!duplicated(lookup_tbl$sample), , drop = FALSE]

    metadata_samples <- normalize_sample_key(metadata_aligned$sample)
    matched_order <- lookup_tbl$order[match(metadata_samples, lookup_tbl$sample)]
    missing_order <- !is.finite(matched_order)

    if (any(missing_order)) {
      warning(
        "Injection order file is missing valid order values for samples: ",
        paste(head(metadata_aligned$sample[missing_order], 10), collapse = ", "),
        if (sum(missing_order) > 10) ", ..." else "",
        ". Falling back to metadata/data order for those samples."
      )
      matched_order[missing_order] <- fallback_order[missing_order]
    }

    return(matched_order)
  }

  order_values <- order_values[is.finite(order_values)]
  if (length(order_values) != nrow(metadata_aligned)) {
    stop(
      "Injection order file without a sample column must have one valid order ",
      "value per aligned sample. Expected ", nrow(metadata_aligned),
      " values, found ", length(order_values), "."
    )
  }

  order_values
}
