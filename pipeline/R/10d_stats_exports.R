# =============================================================================
# 10d_stats_exports.R
# Statistics exports
# =============================================================================

export_stats_excel_by_model <- function(stats_5sets_by_model, paths, p_value_cutoff, fdr_cutoff, fc_cutoff_log2,
                                        active_variant, log_path = NULL,
                                        statistical_test_type = "student",
                                        test_is_paired = FALSE,
                                        pvalue_correction_method = "FDR",
                                        comparison_mode = "pairwise",
                                        multigroup_test = "kruskal",
                                        multigroup_groups = character(0),
                                        multigroup_pairwise_mode = "selected",
                                        run_metrics = "FDR_and_p_value",
                                        include_multigroup_outputs = TRUE,
                                        output_level = "standard") {
  standard_comparisons <- get0("COMPARISON_NAMES", ifnotfound = character(0), inherits = TRUE)

  pretty_comparison_label <- function(comp_name, model_name, model_groups) {
    switch(
      comp_name,
      ALL_TGvsWT = paste0(model_groups$treatment, " vs ", model_groups$control, " | sex=ALL | model=", model_name),
      F_TGvsWT = paste0(model_groups$treatment, " vs ", model_groups$control, " | sex=F | model=", model_name),
      M_TGvsWT = paste0(model_groups$treatment, " vs ", model_groups$control, " | sex=M | model=", model_name),
      TG_FvsM = paste0("FvsM within ", model_groups$treatment, " | model=", model_name),
      WT_FvsM = paste0("FvsM within ", model_groups$control, " | model=", model_name),
      MULTIGROUP_GLOBAL = paste0("Multi-group global test | model=", model_name),
      comp_name
    )
  }

  safe_sheet_name <- function(name, used = character(0)) {
    name <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", as.character(name)[1])
    name <- substr(name, 1, 31)
    if (!nzchar(name)) {
      name <- "Sheet"
    }
    candidate <- name
    i <- 1
    while (candidate %in% used) {
      suffix <- paste0("_", i)
      candidate <- paste0(substr(name, 1, 31 - nchar(suffix)), suffix)
      i <- i + 1
    }
    candidate
  }

  col_idx <- function(df, colname) {
    idx <- match(colname, names(df))
    if (is.na(idx)) {
      return(NA_integer_)
    }
    idx
  }

  for (m in names(stats_5sets_by_model)) {
    mp <- get_model_paths(paths, m)
    out_dir <- mp$exports$stats
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    wb <- openxlsx::createWorkbook()

    # Get model-specific group values
    model_groups <- resolve_model_group_values(m)
    tabs <- stats_5sets_by_model[[m]]
    if (!isTRUE(include_multigroup_outputs)) {
      multigroup_tabs <- names(tabs) == "MULTIGROUP_GLOBAL" | grepl("^PAIR_", names(tabs))
      tabs <- tabs[!multigroup_tabs]
    }

    configured_multigroup_groups <- parse_multigroup_groups(multigroup_groups)
    groups_description <- if (length(configured_multigroup_groups) > 0) {
      paste(configured_multigroup_groups, collapse = " | ")
    } else if ("MULTIGROUP_GLOBAL" %in% names(tabs) &&
               !is.null(tabs[["MULTIGROUP_GLOBAL"]]) &&
               "groups_compared" %in% names(tabs[["MULTIGROUP_GLOBAL"]])) {
      detected <- unique(stats::na.omit(as.character(tabs[["MULTIGROUP_GLOBAL"]]$groups_compared)))
      if (length(detected) > 0) paste(detected, collapse = " | ") else "auto-detected per model"
    } else {
      "auto-detected per model"
    }

    stats_definition <- switch(
      tolower(as.character(comparison_mode)[1]),
      multigroup = paste0(
        "Global test=", multigroup_test,
        "; groups=", groups_description,
        "; pairwise follow-up=", multigroup_pairwise_mode,
        if (!identical(multigroup_pairwise_mode, "none")) paste0(" using ", statistical_test_type, " (paired=", test_is_paired, ")") else ""
      ),
      both = paste0(
        "Global test=", multigroup_test,
        "; groups=", groups_description,
        "; standard pairwise test=", statistical_test_type,
        " (paired=", test_is_paired, ")",
        "; multi-group pairwise follow-up=", multigroup_pairwise_mode
      ),
      paste0("Pairwise test=", statistical_test_type, " (paired=", test_is_paired, ")")
    )

    significance_logic <- switch(
      as.character(run_metrics)[1],
      p_value = paste0("significant when p_value < ", p_value_cutoff),
      FDR = paste0("significant when FDR < ", fdr_cutoff),
      paste0("significant when p_value < ", p_value_cutoff, " OR FDR < ", fdr_cutoff)
    )

    batch_adjusted_any <- any(vapply(tabs, function(tab) {
      !is.null(tab) && "batch_adjusted" %in% names(tab) && any(tab$batch_adjusted %in% TRUE, na.rm = TRUE)
    }, logical(1)))
    batch_levels_detected <- unique(unlist(lapply(tabs, function(tab) {
      if (is.null(tab) || !"batch_levels" %in% names(tab)) return(character(0))
      values <- as.character(tab$batch_levels)
      values[!is.na(values) & nzchar(values)]
    }), use.names = FALSE))

    readme <- tibble(
      field = c(
        "model", "Output_level", "Control group", "Treatment group", "Active_variant",
        "Analysis_mode", "Global_multigroup_test", "Multigroup_groups",
        "Pairwise_follow_up", "Pairwise_test", "Pairwise_test_paired",
        "Batch_adjusted", "Batch_levels",
        "P_value_correction", "Significance_metric", "p_value_cutoff",
        "fdr_cutoff", "fc_cutoff_log2", "log2FC_definition",
        "stats_definition", "significance_logic", "Multigroup_interpretation",
        "Volcano_policy"
      ),
      value = c(
        m, output_level, model_groups$control, model_groups$treatment, active_variant,
        comparison_mode,
        if (comparison_mode %in% c("multigroup", "both")) multigroup_test else "not run",
        if (comparison_mode %in% c("multigroup", "both")) groups_description else "not applicable",
        if (comparison_mode %in% c("multigroup", "both")) multigroup_pairwise_mode else "not applicable",
        statistical_test_type,
        as.character(test_is_paired),
        as.character(batch_adjusted_any),
        if (length(batch_levels_detected) > 0) paste(batch_levels_detected, collapse = "; ") else "not available",
        pvalue_correction_method,
        run_metrics,
        p_value_cutoff, fdr_cutoff, fc_cutoff_log2,
        "Pairwise only: log2FC = log2(mean(num_prelog)/mean(den_prelog)); for FvsM: log2(F/M). MULTIGROUP_GLOBAL keeps FC and log2FC as NA.",
        stats_definition,
        paste0(significance_logic, "; pairwise green rows additionally require |log2FC| >= ", fc_cutoff_log2),
        "MULTIGROUP_GLOBAL is exploratory/complementary. A significant global test shows that at least one group differs, but does not identify which group or provide a directional effect. Inspect group means and selected pairwise follow-up tests.",
        "Volcano plots and Up/Down direction are generated only for pairwise comparisons. MULTIGROUP_GLOBAL is always excluded."
      )
    )

    openxlsx::addWorksheet(wb, "README")
    openxlsx::writeData(wb, "README", readme)

    sig_rows <- list()
    comparisons <- unique(c(standard_comparisons[standard_comparisons %in% names(tabs)], setdiff(names(tabs), standard_comparisons)))
    used_sheets <- "README"

    for (nm in comparisons) {
      df <- tabs[[nm]]
      sheet_nm <- safe_sheet_name(nm, used = used_sheets)
      used_sheets <- c(used_sheets, sheet_nm)
      openxlsx::addWorksheet(wb, sheet_nm)

      if (is.null(df) || nrow(df) == 0) {
        tmp <- tibble(note = "Not enough samples or no valid tests")
        openxlsx::writeData(wb, sheet_nm, tmp)
        next
      }

      is_multigroup_global <- identical(nm, "MULTIGROUP_GLOBAL")
      if (is_multigroup_global) {
        # Preserve the explicit non-directional contract in the exported sheet.
        df_clean <- df %>%
          dplyr::mutate(
            FC_num_over_den = NA_real_,
            log2FC_num_over_den = NA_real_,
            row_sig_p = as.integer(!is.na(p_value) & p_value < p_value_cutoff),
            row_sig_fdr = as.integer(!is.na(FDR) & FDR < fdr_cutoff),
            row_sig_any = as.integer(row_sig_p == 1 | row_sig_fdr == 1),
            row_sig_and_fc = 0L
          )
      } else {
        df_clean <- df %>%
          dplyr::rename(
            FC = dplyr::any_of("FC_num_over_den"),
            log2FC = dplyr::any_of("log2FC_num_over_den")
          )
        if (!"FC" %in% names(df_clean)) {
          df_clean$FC <- NA_real_
        }
        if (!"log2FC" %in% names(df_clean)) {
          df_clean$log2FC <- NA_real_
        }

        df_clean <- df_clean %>%
          dplyr::mutate(
            up_down = dplyr::case_when(
              is.na(log2FC) ~ NA_character_,
              log2FC > 0 ~ "Up",
              log2FC < 0 ~ "Down",
              TRUE ~ "Flat"
            ),
            row_sig_p = as.integer(!is.na(p_value) & p_value < p_value_cutoff),
            row_sig_fdr = as.integer(!is.na(FDR) & FDR < fdr_cutoff),
            row_sig_any = as.integer(row_sig_p == 1 | row_sig_fdr == 1),
            row_sig_and_fc = as.integer(
              row_sig_any == 1 &
                !is.na(log2FC) &
                abs(log2FC) >= fc_cutoff_log2
            )
          )
      }

      openxlsx::writeData(wb, sheet_nm, df_clean)
      openxlsx::freezePane(wb, sheet_nm, firstRow = TRUE)
      openxlsx::addFilter(wb, sheet_nm, rows = 1, cols = 1:ncol(df_clean))

      style_red <- openxlsx::createStyle(fgFill = "#FF0000", fontColour = "#FFFFFF")
      style_yellow <- openxlsx::createStyle(fgFill = "#FFD966")
      style_green <- openxlsx::createStyle(fgFill = "#00B050", fontColour = "#FFFFFF")

      cols_full_row <- 1:ncol(df_clean)

      green_rows <- which(df_clean$row_sig_and_fc == 1) + 1
      yellow_rows <- which(df_clean$row_sig_any == 1 & df_clean$row_sig_and_fc == 0) + 1

      if (length(yellow_rows) > 0) {
        openxlsx::addStyle(
          wb, sheet_nm,
          style = style_yellow,
          rows = yellow_rows,
          cols = cols_full_row,
          gridExpand = TRUE,
          stack = FALSE
        )
      }

      if (length(green_rows) > 0) {
        openxlsx::addStyle(
          wb, sheet_nm,
          style = style_green,
          rows = green_rows,
          cols = cols_full_row,
          gridExpand = TRUE,
          stack = FALSE
        )
      }

      p_col <- col_idx(df_clean, "p_value")
      if (!is.na(p_col)) {
        p_rows <- which(!is.na(df_clean$p_value) & df_clean$p_value < p_value_cutoff) + 1
        if (length(p_rows) > 0) {
          openxlsx::addStyle(
          wb, sheet_nm,
            style = style_red,
            rows = p_rows,
            cols = p_col,
            gridExpand = TRUE,
            stack = FALSE
          )
        }
      }

      fdr_col <- col_idx(df_clean, "FDR")
      if (!is.na(fdr_col)) {
        fdr_rows <- which(!is.na(df_clean$FDR) & df_clean$FDR < fdr_cutoff) + 1
        if (length(fdr_rows) > 0) {
          openxlsx::addStyle(
          wb, sheet_nm,
            style = style_red,
            rows = fdr_rows,
            cols = fdr_col,
            gridExpand = TRUE,
            stack = FALSE
          )
        }
      }

      helper_cols <- c(
        col_idx(df_clean, "row_sig_p"),
        col_idx(df_clean, "row_sig_fdr"),
        col_idx(df_clean, "row_sig_any"),
        col_idx(df_clean, "row_sig_and_fc")
      )
      helper_cols <- helper_cols[!is.na(helper_cols)]

      visible_cols <- setdiff(seq_len(ncol(df_clean)), helper_cols)

      if (length(visible_cols) > 0) {
        openxlsx::setColWidths(wb, sheet_nm, cols = visible_cols, widths = "auto")
      }

      if (length(helper_cols) > 0) {
        openxlsx::setColWidths(wb, sheet_nm, cols = helper_cols, widths = 0)
      }
      # collect significant rows for aggregated sheet
      # The aggregate Significant sheet is directional/pairwise. Global
      # multi-group hits remain in their dedicated MULTIGROUP_GLOBAL sheet.
      sig_subset <- if (is_multigroup_global) df_clean[0, , drop = FALSE] else df_clean[df_clean$row_sig_any == 1, , drop = FALSE]
      if (nrow(sig_subset) > 0) {
        sig_subset$comparison <- pretty_comparison_label(
          comp_name = nm,
          model_name = m,
          model_groups = model_groups
        )
        sig_rows[[length(sig_rows) + 1]] <- sig_subset
      }
    }

    # aggregated sheet of significant features across comparisons
    openxlsx::addWorksheet(wb, "Significant")
    if (length(sig_rows) == 0) {
      openxlsx::writeData(wb, "Significant", tibble::tibble(note = "No significant features found for this model"))
    } else {
      agg <- tryCatch({ do.call(rbind, sig_rows) }, error = function(e) NULL)
      if (is.null(agg) || nrow(agg) == 0) {
        openxlsx::writeData(wb, "Significant", tibble::tibble(note = "No significant features found for this model"))
      } else {
        # move comparison and direction columns to front
        cols <- c("comparison", "up_down", setdiff(names(agg), c("comparison", "up_down")))
        agg <- agg[, intersect(cols, names(agg)), drop = FALSE]
        openxlsx::writeData(wb, "Significant", agg)
        openxlsx::freezePane(wb, "Significant", firstRow = TRUE)
      }
    }

    openxlsx::setColWidths(wb, "README", cols = 1:2, widths = "auto")

    if ("MULTIGROUP_GLOBAL" %in% names(tabs)) {
      multigroup_readme_path <- file.path(out_dir, "MULTIGROUP_README.txt")
      writeLines(
        c(
          "MULTIGROUP_GLOBAL — exploratory global analysis",
          "",
          paste0("Model: ", m),
          paste0("Global test: ", multigroup_test),
          paste0("Groups: ", groups_description),
          "",
          "A significant p-value or FDR indicates that at least one analyzed group differs from the others.",
          "The global ANOVA, Welch ANOVA, or Kruskal-Wallis test does not identify a single numerator/denominator and therefore has no directional fold change.",
          "FC_num_over_den and log2FC_num_over_den are intentionally NA; Up/Down direction and volcano plots are not produced.",
          "Use the per-group mean columns for exploration and selected, biologically interpretable pairwise comparisons for directional follow-up.",
          "Primary pairwise outputs remain ALL_TGvsWT, F_TGvsWT, M_TGvsWT, TG_FvsM, and WT_FvsM."
        ),
        con = multigroup_readme_path,
        useBytes = TRUE
      )
      if (!is.null(log_path)) {
        append_log_line(log_path, paste0("- MULTIGROUP_README.txt -> explanatory multi-group README | path: ", multigroup_readme_path))
      }
    }

    out_xlsx <- file.path(out_dir, paste0("STATS_model_", m, ".xlsx"))
    openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

    if (!is.null(log_path)) {
      append_log_line(
        log_path,
        paste0("- ", basename(out_xlsx), " -> (Excel workbook) | path: ", out_xlsx)
      )
    }
  }

  message("  ✓ Stats Excel per model saved in: ", out_dir)
}

# export_significant_metabolites_txt_by_model(): Writes one TXT list per
# comparison, with the significance metric encoded in the filename.
compact_stats_path_component <- function(value, max_chars = 72L) {
  value <- gsub("[^A-Za-z0-9_-]+", "_", as.character(value)[1])
  value <- gsub("_+", "_", value)
  value <- gsub("^_+|_+$", "", value)
  if (!nzchar(value)) value <- "comparison"

  max_chars <- max(20L, as.integer(max_chars)[1])
  if (nchar(value) <= max_chars) return(value)

  code_points <- utf8ToInt(value)
  checksum <- sum((as.numeric(code_points) * seq_along(code_points)) %% 10000019) %% 100000000
  suffix <- sprintf("_%08d", checksum)
  paste0(substr(value, 1, max_chars - nchar(suffix)), suffix)
}

make_significant_metabolite_path <- function(out_dir, metric, comparison_name) {
  metric_tag <- if (identical(metric, "FDR")) "fdr" else "p"
  target_dir <- file.path(out_dir, "significant")
  target_dir_abs <- normalizePath(target_dir, winslash = "/", mustWork = FALSE)

  # Keep comfortable headroom below the traditional Windows MAX_PATH limit.
  fixed_chars <- nchar(target_dir_abs) + nchar(paste0("/SIG__", metric_tag, ".txt"))
  comparison_budget <- max(20L, min(72L, 235L - fixed_chars))
  comparison_tag <- compact_stats_path_component(comparison_name, comparison_budget)

  file.path(target_dir, paste0("SIG_", comparison_tag, "_", metric_tag, ".txt"))
}

write_significant_metabolites_txt_safe <- function(lines, path) {
  if (is.null(path) || length(path) != 1 || is.na(path) || !nzchar(path)) {
    stop("write_significant_metabolites_txt_safe() received an invalid path.")
  }

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(dirname(path))) {
    stop("Cannot create significant metabolites output directory: ", dirname(path))
  }

  write_path <- path
  abs_path <- path

  path_dir_abs <- tryCatch(
    normalizePath(dirname(path), winslash = "/", mustWork = TRUE),
    error = function(e) NULL
  )
  cwd_abs <- tryCatch(
    normalizePath(getwd(), winslash = "/", mustWork = TRUE),
    error = function(e) NULL
  )

  if (!is.null(path_dir_abs)) {
    abs_path <- paste0(path_dir_abs, "/", basename(path))
  }

  if (!is.null(path_dir_abs) && !is.null(cwd_abs)) {
    cwd_prefix <- paste0(tolower(cwd_abs), "/")
    abs_lower <- tolower(abs_path)
    if (startsWith(abs_lower, cwd_prefix)) {
      write_path <- substr(abs_path, nchar(cwd_abs) + 2, nchar(abs_path))
    }
  }

  tryCatch(
    {
      con <- file(write_path, open = "w", encoding = "UTF-8")
      on.exit(close(con), add = TRUE)
      writeLines(lines, con = con, useBytes = TRUE)
    },
    error = function(e) {
      stop(
        "Cannot open significant metabolites TXT for writing.\n",
        "Requested path: ", path, "\n",
        "Write path: ", write_path, "\n",
        "Absolute path length: ", nchar(abs_path), "\n",
        "Original error: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )
}

export_significant_metabolites_txt_by_model <- function(stats_5sets_by_model, paths, p_value_cutoff, fdr_cutoff,
                                                        fc_cutoff_log2, active_variant,
                                                        log_path = NULL,
                                                        require_fc_cutoff = FALSE) {
  standard_comparisons <- get0("COMPARISON_NAMES", ifnotfound = character(0), inherits = TRUE)
  metrics <- c("p_value", "FDR")

  for (m in names(stats_5sets_by_model)) {
    mp <- get_model_paths(paths, m)
    out_dir <- mp$exports$stats
    dir.create(file.path(out_dir, "significant"), recursive = TRUE, showWarnings = FALSE)

    tabs <- stats_5sets_by_model[[m]]
    comparisons <- unique(c(standard_comparisons[standard_comparisons %in% names(tabs)], setdiff(names(tabs), standard_comparisons)))

    for (nm in comparisons) {
      df <- tabs[[nm]]

      if (identical(nm, "MULTIGROUP_GLOBAL")) {
        message("  - Significant-metabolite TXT skipped for MULTIGROUP_GLOBAL: the global test is non-directional; use its Excel sheet and top-feature heatmap.")
        next
      }

      for (metric in metrics) {
        out_path <- make_significant_metabolite_path(
          out_dir = out_dir,
          metric = metric,
          comparison_name = nm
        )
        dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

        out_names <- character(0)

        if (!is.null(df) && nrow(df) > 0) {
          if (!"log2FC_num_over_den" %in% names(df)) {
            df$log2FC_num_over_den <- NA_real_
          }
          out_names <- df %>%
            dplyr::mutate(
                cutoff = dplyr::if_else(metric == "FDR", fdr_cutoff, p_value_cutoff),
                passes_metric = !is.na(.data[[metric]]) & .data[[metric]] < cutoff,
              passes_fc_cutoff = !is.na(log2FC_num_over_den) & abs(log2FC_num_over_den) >= fc_cutoff_log2,
              primary_name = dplyr::na_if(trimws(as.character(Name)), "")
            ) %>%
            dplyr::filter(passes_metric) %>%
            dplyr::filter(!require_fc_cutoff | passes_fc_cutoff) %>%
            dplyr::arrange(.data[[metric]], dplyr::desc(abs(log2FC_num_over_den)), featureID) %>%
            dplyr::pull(primary_name) %>%
            unique() %>%
            stats::na.omit() %>%
            as.character()
        }

        if (length(out_names) > 0) {
          write_significant_metabolites_txt_safe(out_names, out_path)

          if (!is.null(log_path)) {
            log_written_object(
              log_path,
              out_path,
              out_names,
              note = paste0(
                "Significant metabolites list export | model=", m,
                " | comparison=", nm,
                " | metric=", metric,
                " | p_value_cutoff=", p_value_cutoff,
                " | fdr_cutoff=", fdr_cutoff,
                " | fc_cutoff_log2=", fc_cutoff_log2,
                " | require_fc_cutoff=", require_fc_cutoff,
                " | active_variant=", active_variant
              )
            )
          }
        }
      }
    }
  }

  message("  Significant metabolite TXT exports saved under: ", file.path(out_dir, "significant"))
}
