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
    model     = c("model", "disease", "condition", "phenotype", "status")
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
      type = if (!("type" %in% names(.))) {
        dplyr::if_else(stringr::str_detect(sample, "^QC"), "QC", "Sample")
      } else {
        stringr::str_trim(as.character(type))
      }
    )
}
