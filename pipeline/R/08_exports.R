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

# Build a MetaboAnalyst-ready export table from the raw assay and aligned metadata.
build_metaboanalyst_export_df <- function(raw_df,
                                          metadata_aligned,
                                          feature_tbl) {
  df_named <- rename_feature_cols(raw_df, feature_tbl, sample_col = "sample")

  df_named %>%
    dplyr::left_join(
      metadata_aligned %>%
        dplyr::transmute(
          sample,
          type,
          model,
          group,
          sex,
          model_group = dplyr::if_else(type == "QC", "QC", paste(model, group, sep = "-")),
          Class = dplyr::if_else(type == "QC", "QC", model_group)
        ),
      by = "sample"
    )
}

# Export a single global MetaboAnalyst table from the raw assay.
# Two CSV files are written:
# 1. One with only Sample type (no QC), named "MA_ACTIVE_raw_GLOBAL_NO_QC.csv"
# 2. One with both Sample and QC types, named "MA_ACTIVE_raw_GLOBAL_WITH_QC.csv"
export_metaboanalyst_global_raw <- function(raw_df,
                                            metadata_aligned,
                                            feature_tbl,
                                            export_dir,
                                            log_path = NULL,
                                            value_label = "raw",
                                            file_prefix = "MA_ACTIVE",
                                            include_with_qc = TRUE) {
  df <- build_metaboanalyst_export_df(raw_df, metadata_aligned, feature_tbl)

  # Create export directory if it doesn't exist.
  dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

  df_no_qc <- df %>%
    dplyr::filter(type == "Sample") %>%
    dplyr::select(-type) %>%
    dplyr::relocate(Class, model_group, model, group, sex, .after = sample)

  out_no_qc <- file.path(export_dir, paste0(file_prefix, "_", value_label, "_GLOBAL_NO_QC.csv"))
  if (nrow(df_no_qc) >= 2) {
    write_csv_safe(df_no_qc, out_no_qc)
    if (!is.null(log_path)) {
      log_written_object(
        log_path,
        out_no_qc,
        df_no_qc,
        note = paste0("MetaboAnalyst export (GLOBAL, ", value_label, ", NO_QC)")
      )
    }
  }

  if (isTRUE(include_with_qc)) {
    df_with_qc <- df %>%
      dplyr::select(-type) %>%
      dplyr::relocate(Class, model_group, model, group, sex, .after = sample)

    # With QC samples, but only if there are at least 2 rows (MetaboAnalyst requires at least 2 samples).
    out_with_qc <- file.path(export_dir, paste0(file_prefix, "_", value_label, "_GLOBAL_WITH_QC.csv"))
    if (nrow(df_with_qc) >= 2) {
      write_csv_safe(df_with_qc, out_with_qc)
      if (!is.null(log_path)) {
        log_written_object(
          log_path,
          out_with_qc,
          df_with_qc,
          note = paste0("MetaboAnalyst export (GLOBAL, ", value_label, ", WITH_QC)")
        )
      }
    }
  }

  message("  ✓ Global MetaboAnalyst exports created: ", export_dir)
}
