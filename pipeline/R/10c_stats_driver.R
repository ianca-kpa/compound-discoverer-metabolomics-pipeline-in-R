# =============================================================================
# 10c_stats_driver.R
# Statistics driver
# =============================================================================

run_all_stats_5sets_per_model <- function(mat_log2, mat_raw, metadata_aligned, feature_tbl, paths,
                                          p_value_cutoff = p_value_cutoff, fdr_cutoff = fdr_cutoff, fc_cutoff_log2 = fc_cutoff_log2,
                                          run_metrics = run_metrics,
                                          make_volcano_plots = TRUE,
                                          volcano_style = volcano_style,
                                          comparison_configs = COMPARISON_CONFIGS,
                                          statistical_test_type = "student",
                                          test_is_paired = FALSE,
                                          pvalue_correction_method = "FDR",
                                          comparison_mode = "pairwise",
                                          multigroup_groups = character(0),
                                          multigroup_test = "kruskal",
                                          multigroup_pairwise_mode = "selected",
                                          multigroup_pairwise_pairs = NULL,
                                          export_all_pairwise_multigroup = FALSE) {
  metrics <- normalize_metric_modes(run_metrics)
  styles_to_export <- resolve_volcano_styles(volcano_style)
  comparison_mode <- tolower(trimws(as.character(comparison_mode)[1]))
  if (!comparison_mode %in% c("pairwise", "multigroup", "both")) {
    comparison_mode <- "pairwise"
  }
  models <- sort(unique(metadata_aligned$model[metadata_aligned$type == "Sample"]))
  out <- list()

  get_volcano_dir <- function(prefix, mp) {
    # Accept both legacy lowercase prefixes and comparison-config uppercase prefixes
    if (prefix %in% c("ALL", "tg_vs_wt", "tg-vs-wt", "ALL_TGvsWT")) {
      return(mp$plots$volcano_tg_vs_wt)
    }
    if (prefix %in% c("sex", "tg-f_vs_wt-f", "tg_f_vs_wt_f", "tg-m_vs_wt-m", "tg_m_vs_wt_m", "F_TGvsWT", "M_TGvsWT")) {
      return(mp$plots$volcano_by_sex)
    }
    if (prefix %in% c("sex_FvsM", "sex_MvsF", "tg-f_vs_tg-m", "tg_f_vs_tg_m", "TG_FvsM", "TG_MvsF")) {
      return(mp$plots$volcano_tg_f_vs_tg_m)
    }
    if (prefix %in% c("wt-f_vs_wt-m", "wt_f_vs_wt_m", "WT_FvsM", "WT_MvsF")) {
      return(mp$plots$volcano_wt_f_vs_wt_m)
    }

    # fallback: try lowercase-normalized prefix
    np <- tolower(prefix)
    if (np %in% c("all", "tg_vs_wt", "tg-vs-wt")) return(mp$plots$volcano_tg_vs_wt)
    if (np %in% c("sex", "tg-f_vs_wt-f", "tg_f_vs_wt_f", "tg-m_vs_wt-m", "tg_m_vs_wt_m")) return(mp$plots$volcano_by_sex)
    if (np %in% c("sex_fvsm", "sex_mvsf", "tg-f_vs_tg-m", "tg_f_vs_tg_m", "tg_fvsm", "tg_mvsf")) return(mp$plots$volcano_tg_f_vs_tg_m)
    if (np %in% c("wt-f_vs_wt-m", "wt_f_vs_wt_m", "wt_fvsm", "wt_mvsf")) return(mp$plots$volcano_wt_f_vs_wt_m)
    if (np %in% c("multigroup_pairwise")) return(mp$plots$volcano_by_group)

    stop("Unknown comparison prefix: ", prefix)
  }

  for (m in models) {
    mp <- get_model_paths(paths, m)
    meta_m <- metadata_aligned %>% filter(type == "Sample", model == m)
    model_groups <- resolve_model_group_values(m)
    message("  - Stats model groups: model=", m, " | control=", model_groups$control, " | treatment=", model_groups$treatment)
    message("    - Test type: ", statistical_test_type, " | Paired: ", test_is_paired, " | P-value correction: ", pvalue_correction_method)
    out[[m]] <- list()

    groups_for_model <- parse_multigroup_groups(multigroup_groups)
    if (length(groups_for_model) == 0) {
      groups_for_model <- sort(unique(as.character(meta_m$group[!is.na(meta_m$group) & nzchar(as.character(meta_m$group))])))
    }

    # The five predefined pairwise comparisons are the primary analysis in every
    # mode. Multi-group analysis is exploratory and is added to, never substituted
    # for, these biologically defined comparisons.
    active_comparison_configs <- comparison_configs
    if (comparison_mode %in% c("multigroup", "both")) {
      active_comparison_configs <- c(
        active_comparison_configs,
        build_multigroup_pairwise_configs(
          groups = groups_for_model,
          pairwise_mode = multigroup_pairwise_mode,
          selected_pairs = multigroup_pairwise_pairs
        )
      )

      st_multi <- compute_multigroup_stats_general(
        mat_log2 = mat_log2,
        mat_raw = mat_raw,
        meta_sub = meta_m,
        feature_tbl = feature_tbl,
        groups = groups_for_model,
        multigroup_test = multigroup_test,
        pvalue_correction_method = pvalue_correction_method
      )
      out[[m]][["MULTIGROUP_GLOBAL"]] <- st_multi
      if (is.null(st_multi)) {
        message("    - MULTIGROUP_GLOBAL skipped: at least 3 groups with at least 2 samples each are required.")
      } else {
        message("    - Multi-group global test: ", multigroup_test, " | groups=", paste(groups_for_model, collapse = ", "))
      }
      if (!is.null(st_multi) && isTRUE(make_volcano_plots)) {
        message("    - MULTIGROUP_GLOBAL volcano skipped: the global ", multigroup_test,
                " test detects differences among groups but has no single numerator/denominator or directional log2FC. Volcano plots are generated only for pairwise comparisons.")
      }
    }

    for (comp_name in names(active_comparison_configs)) {
      cfg <- active_comparison_configs[[comp_name]]
      meta_sub <- cfg$meta_filter(meta_m, model_name = m)

      compare_den <- cfg$stats_den
      compare_num <- cfg$stats_num
      use_model_groups <- is.null(cfg$use_model_groups) || isTRUE(cfg$use_model_groups)
      if (identical(cfg$stats_compare_var, "group") && isTRUE(use_model_groups)) {
        compare_den <- model_groups$control
        compare_num <- model_groups$treatment
      }

      st <- compute_ttest_stats_general(
        mat_log2 = mat_log2,
        mat_raw = mat_raw,
        meta_sub = meta_sub,
        feature_tbl = feature_tbl,
        compare_var = cfg$stats_compare_var,
        num_level = compare_num,
        den_level = compare_den,
        statistical_test_type = statistical_test_type,
        test_is_paired = test_is_paired,
        pvalue_correction_method = pvalue_correction_method
      )

      out[[m]][[comp_name]] <- st

      if (!isTRUE(make_volcano_plots)) next
      is_primary_comparison <- comp_name %in% names(comparison_configs)
      if (!is_primary_comparison && !isTRUE(export_all_pairwise_multigroup)) {
        message("    - Additional multi-group pairwise volcano skipped at this output level: ", comp_name)
        next
      }

      for (met in metrics) {
        out_dir <- get_volcano_dir(cfg$prefix, mp)

          comp_label <- switch(
            comp_name,
            ALL_TGvsWT = paste0(model_groups$treatment, " vs ", model_groups$control, " | sex=ALL"),
            F_TGvsWT = paste0(model_groups$treatment, " vs ", model_groups$control, " | sex=F"),
            M_TGvsWT = paste0(model_groups$treatment, " vs ", model_groups$control, " | sex=M"),
            TG_FvsM = paste0("FvsM within ", model_groups$treatment),
            WT_FvsM = paste0("FvsM within ", model_groups$control),
            cfg$label
          )

        for (style_name in styles_to_export) {
          base_name <- paste0("volcano_ACTIVE_model_", m, "_", comp_name, "_metric_", met)
          file_name <- if (length(styles_to_export) > 1) {
            paste0(base_name, "_", style_name, ".png")
          } else {
            paste0(base_name, ".png")
          }

          out_path <- file.path(out_dir, file_name)

          title <- paste0(
            "Volcano (", style_name, ") | model=", m,
            " | ", comp_label,
            " | metric=", met,
            " | log2FC=log2(", compare_num, "/", compare_den, ")"
          )

          if (is.null(st) || nrow(st) == 0) {
            save_placeholder_volcano(
              out_path = out_path,
              title = title,
              reason = "Not enough samples for this comparison (need >=2 per group and >=4 total), or no valid tests after filtering."
            )
            next
          }

          tryCatch(
            {
              plot_volcano_metric(
                stats_df = st,
                title = title,
                out_path = out_path,
                metric = met,
                alpha = p_value_cutoff,
                fdr_alpha = fdr_cutoff,
                fc_cutoff_log2 = fc_cutoff_log2
              )
            },
            error = function(e) {
              message("  - Volcano failed: ", out_path)
              message("    Error: ", conditionMessage(e))

              save_placeholder_volcano(
                out_path = out_path,
                title = title,
                reason = paste0("Volcano failed with error: ", conditionMessage(e))
              )
            }
          )
        }
      }
    }
  }

  out
}

# export_stats_excel_by_model() exports the computed statistics for each model in an Excel workbook 
