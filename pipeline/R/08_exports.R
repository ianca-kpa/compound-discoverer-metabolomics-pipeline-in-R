# =============================================================================
# 08_exports.R
# Export helpers (MetaboAnalyst)
# =============================================================================

# This module contains helper functions for exporting data in formats suitable for MetaboAnalyst, as specified in the export settings in config/settings.R.

# log2_transform(): Applies log2 transformation to the assay matrix with a specified offset to avoid log(0) issues.
log2_transform <- function(mat, offset) {
  log2(mat + offset)
}

# Rename feature columns in a data frame based on the feature_tbl mapping, while ensuring unique column names.
rename_feature_cols <- function(df, feature_tbl, sample_col = "sample") {
  map_vec <- feature_tbl$display_name
  names(map_vec) <- feature_tbl$featureID

  out <- df
  feat_cols <- setdiff(names(out), sample_col)
  hit <- feat_cols %in% names(map_vec)

  new_names <- feat_cols
  new_names[hit] <- unname(map_vec[feat_cols[hit]])
  new_names <- make.unique(new_names, sep = "_dup")

  names(out)[match(feat_cols, names(out))] <- new_names
  out
}

# This function creates two CSV files for a given model:
  # 1. One with only Sample type (no QC), named "MA_ACTIVE_log2_NO_QC_model_<model_name>.csv"
  # 2. One with both Sample and QC types, named "MA_ACTIVE_log2_WITH_QC_model_<model_name>.csv"

export_metaboanalyst_one_model <- function(log2_df,
                                           metadata_aligned,
                                           feature_tbl,
                                           model_name,
                                           export_dir,
                                           log_path = NULL) {
  
  df_named <- rename_feature_cols(log2_df, feature_tbl, sample_col = "sample")

  df <- df_named %>%
    dplyr::left_join(
      metadata_aligned %>% dplyr::select(sample, type, group, model),
      by = "sample"
    ) %>%
    dplyr::mutate(
      Class = dplyr::if_else(
        type == "QC",
        "QC",
        map_comparison_group_display_values(group, model_name = model_name)
      )
    )

  # Create export directory if it doesn't exist
  dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

  df_no_qc <- df %>%
    dplyr::filter(type == "Sample", model == model_name) %>%
    dplyr::select(sample, Class, dplyr::everything(), -type, -group, -model)

  out_no_qc <- file.path(export_dir, paste0("MA_ACTIVE_log2_NO_QC_model_", model_name, ".csv"))
  if (nrow(df_no_qc) >= 2) {
    write_csv_safe(df_no_qc, out_no_qc)
    if (!is.null(log_path)) {
      log_written_object(
        log_path,
        out_no_qc,
        df_no_qc,
        note = paste0("MetaboAnalyst export (ACTIVE, log2, NO_QC) | model=", model_name)
      )
    }
  }

  df_with_qc <- df %>%
    dplyr::filter(model == model_name | type == "QC") %>%
    dplyr::select(sample, Class, dplyr::everything(), -type, -group, -model)

  # With QC samples, but only if there are at least 2 rows (MetaboAnalyst requires at least 2 samples).
  out_with_qc <- file.path(export_dir, paste0("MA_ACTIVE_log2_WITH_QC_model_", model_name, ".csv"))
  if (nrow(df_with_qc) >= 2) {
    write_csv_safe(df_with_qc, out_with_qc)
    if (!is.null(log_path)) {
      log_written_object(
        log_path,
        out_with_qc,
        df_with_qc,
        note = paste0("MetaboAnalyst export (ACTIVE, log2, WITH_QC) | model=", model_name)
      )
    }
  }

  message("  ✓ MetaboAnalyst exports created for model ", model_name, ": ", export_dir)
}
