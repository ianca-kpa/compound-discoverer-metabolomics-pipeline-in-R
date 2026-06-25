validate_metadata_columns <- function(path,
                                      metadata_mapping = NULL,
                                      allowed_groups = c("WT", "TG"),
                                      model_allowed_groups_by_model = NULL) {
  alias_map <- list(
    sample = c("sample", "sample_id", "sample_name", "id_sample", "id", "name"),
    weight = c("weight", "weight_mg", "mass", "mass_mg", "mg", "sample_weight", "sample_mass", "weight_g", "mass_g", "sample_weight_g", "sample_mass_g"),
    group = c("group", "treatment", "treat"),
    sex = c("sex", "gender"),
    model = c("model", "disease", "condition", "phenotype", "status")
  )

  md <- safe_read_table(path)
  actual <- tolower(trimws(as.character(names(md))))

  missing <- character(0)
  resolved_cols <- list()

  for (target in names(alias_map)) {
    mapped_col <- ""
    if (!is.null(metadata_mapping) && !is.null(metadata_mapping[[target]])) {
      mapped_col <- normalize_name(metadata_mapping[[target]])
    }

    if (nzchar(mapped_col)) {
      if (!(mapped_col %in% actual)) {
        missing <- c(missing, target)
      } else {
        resolved_cols[[target]] <- mapped_col
      }
    } else {
      aliases <- normalize_name(alias_map[[target]])
      found <- intersect(aliases, actual)
      if (length(found) == 0) {
        missing <- c(missing, target)
      } else {
        resolved_cols[[target]] <- found[1]
      }
    }
  }

  if (length(missing) > 0) {
    return(list(
      ok = FALSE,
      message = paste(
        "Metadata missing required columns or mappings:",
        paste(missing, collapse = ", ")
      )
    ))
  }

  allowed_groups <- unique(trimws(as.character(allowed_groups)))
  allowed_groups <- allowed_groups[nzchar(allowed_groups)]

  if (length(allowed_groups) < 2 && !is.null(model_allowed_groups_by_model) && length(model_allowed_groups_by_model) > 0) {
    inferred_groups <- unique(trimws(unlist(strsplit(as.character(model_allowed_groups_by_model), ",", fixed = TRUE), use.names = FALSE)))
    inferred_groups <- inferred_groups[nzchar(inferred_groups)]
    if (length(inferred_groups) >= 2) {
      allowed_groups <- inferred_groups
    }
  }

  allowed_groups_norm <- toupper(allowed_groups)
  if (length(allowed_groups) < 2) {
    return(list(
      ok = FALSE,
      message = "Please provide at least two allowed group values in this order: control, test (e.g. WT, TG), or define them in model_allowed_groups_by_model."
    ))
  }

  group_col <- resolved_cols[["group"]]
  if (!is.null(group_col) && nzchar(group_col)) {
    col_idx <- which(actual == group_col)
    if (length(col_idx) == 0) {
      return(list(ok = FALSE, message = paste0("Group column '", group_col, "' not found in metadata.")))
    }
    groups_raw <- as.character(md[[col_idx[1]]])
    groups_raw <- trimws(groups_raw)
    qc_aliases <- metadata_qc_group_aliases()
    row_is_qc <- toupper(groups_raw) %in% qc_aliases

    sample_col <- resolved_cols[["sample"]]
    sample_idx <- which(actual == sample_col)
    if (length(sample_idx) > 0) {
      sample_raw <- trimws(as.character(md[[sample_idx[1]]]))
      row_is_qc <- row_is_qc | grepl("^QC", sample_raw, ignore.case = TRUE)
    }

    type_idx <- which(actual %in% c("type", "sample_type"))
    if (length(type_idx) > 0) {
      type_raw <- trimws(as.character(md[[type_idx[1]]]))
      row_is_qc <- row_is_qc | (toupper(type_raw) %in% qc_aliases)
    }

    if (!is.null(model_allowed_groups_by_model) && length(model_allowed_groups_by_model) > 0 && "model" %in% names(resolved_cols)) {
      model_idx <- which(actual == resolved_cols[["model"]])
      if (length(model_idx) > 0) {
        model_raw <- trimws(as.character(md[[model_idx[1]]]))
        valid_rows <- !is.na(groups_raw) & !is_missing_like(groups_raw) & !is.na(model_raw) & !is_missing_like(model_raw)

        if (any(valid_rows)) {
          groups_raw[valid_rows] <- normalize_model_group_pairs(
            groups_raw[valid_rows],
            model_raw[valid_rows],
            model_allowed_groups_by_model,
            allowed_groups[1],
            allowed_groups[2]
          )
        }
      }
    }

    groups_raw <- groups_raw[!row_is_qc]
    groups_raw <- groups_raw[!is.na(groups_raw)]
    groups_raw <- groups_raw[!is_missing_like(groups_raw)]

    groups_norm <- toupper(groups_raw)
    invalid_groups <- sort(unique(groups_raw[!(groups_norm %in% allowed_groups_norm)]))

    if (length(invalid_groups) > 0) {
      return(list(
        ok = FALSE,
        message = paste0(
          "Invalid values in metadata group column ('",
          group_col,
          "'). Allowed values: ",
          paste(allowed_groups, collapse = ", "),
          ". Found: ",
          paste(invalid_groups, collapse = ", "),
          ". Please review your metadata file and fix the group column."
        )
      ))
    }
  }

  list(ok = TRUE, message = "Metadata validation passed.")
}
