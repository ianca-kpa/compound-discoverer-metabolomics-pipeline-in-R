# =============================================================================
# 10_stats_volcano.R
# Statistics + volcano + Excel export
# =============================================================================

# Summary:
  # This module contains functions for performing statistical tests (t-tests) between groups, 
  # generating volcano plots, and exporting results to Excel. The main function run_all_stats_
  # 5sets_per_model() computes statistics for predefined comparisons and creates volcano plots, 
  # while export_stats_excel_by_model() saves the results in an Excel workbook with conditional formatting.

resolve_volcano_styles <- function(volcano_style) {
  style <- as.character(volcano_style)[1]
  if (is.na(style) || !nzchar(style)) {
    style <- "classic"
  }

  if (identical(style, "both")) {
    return(c("classic", "gradual"))
  }
  if (style %in% c("classic", "gradual")) {
    return(style)
  }

  warning("Invalid volcano_style: ", style, ". Falling back to 'classic'.")
  "classic"
}

# normalize_metric_modes(): Normalizes the input for metrics to run, 
# expanding "FDR_and_p_value" into its components and validating the selection.
normalize_metric_modes <- function(metrics) {
  valid <- c("FDR", "p_value")

  if (is.null(metrics) || length(metrics) == 0) {
    return(valid)
  }

  expanded <- unlist(lapply(as.character(metrics), function(m) {
    if (identical(m, "FDR_and_p_value")) {
      return(c("FDR", "p_value"))
    }
    m
  }), use.names = FALSE)

  expanded <- unique(expanded)
  expanded <- expanded[expanded %in% valid]

  if (length(expanded) == 0) {
    stop("No valid metric selected. Use one of: FDR, p_value, FDR_and_p_value")
  }

  expanded
}

FLOOR_P <- 1e-300

resolve_statistical_test_type <- function(test_type) {
  test <- tolower(as.character(test_type)[1])

  if (is.na(test) || !nzchar(test)) {
    test <- "student"
  }

  if (!test %in% c("student", "welch", "wilcoxon", "limma")) {
    warning("Invalid statistical_test_type: ", test, ". Falling back to 'student'.")
    test <- "student"
  }

  test
}

resolve_pvalue_correction_method <- function(correction_method) {
  method <- tolower(as.character(correction_method)[1])

  if (is.na(method) || !nzchar(method)) {
    method <- "FDR"
  }

  valid_methods <- c("raw", "fdr", "bonferroni", "holm", "hochberg", "hommel", "by")
  if (!tolower(method) %in% valid_methods) {
    warning("Invalid pvalue_correction_method: ", method, ". Falling back to 'FDR'.")
    method <- "FDR"
  }

  # Map "raw" to appropriate method for p.adjust
  if (identical(tolower(method), "raw")) {
    return("none")
  }

  # Ensure proper capitalization for p.adjust
  switch(tolower(method),
    "fdr" = "BH",
    "bonferroni" = "bonferroni",
    "holm" = "holm",
    "hochberg" = "hochberg",
    "hommel" = "hommel",
    "by" = "BY",
    "BH"
  )
}

perform_statistical_test <- function(x, y, test_type = "student", paired = FALSE) {
  test_type <- tolower(test_type)

  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  if (length(x) < 2 || length(y) < 2) {
    return(list(p_value = NA_real_, test_type_used = test_type, error = "Insufficient data"))
  }

  if (identical(test_type, "student")) {
    result <- tryCatch(
      {
        t <- t.test(y, x, var.equal = TRUE, paired = paired)
        list(p_value = as.numeric(t$p.value), test_type_used = "student", error = NULL)
      },
      error = function(e) {
        list(p_value = NA_real_, test_type_used = "student", error = as.character(e))
      }
    )
    return(result)
  }

  if (identical(test_type, "welch")) {
    result <- tryCatch(
      {
        t <- t.test(y, x, var.equal = FALSE, paired = paired)
        list(p_value = as.numeric(t$p.value), test_type_used = "welch", error = NULL)
      },
      error = function(e) {
        list(p_value = NA_real_, test_type_used = "welch", error = as.character(e))
      }
    )
    return(result)
  }

  if (identical(test_type, "wilcoxon")) {
    result <- tryCatch(
      {
        w <- wilcox.test(y, x, paired = paired, exact = FALSE)
        list(p_value = as.numeric(w$p.value), test_type_used = "wilcoxon", error = NULL)
      },
      error = function(e) {
        list(p_value = NA_real_, test_type_used = "wilcoxon", error = as.character(e))
      }
    )
    return(result)
  }

  if (identical(test_type, "limma")) {
    # Limma requires a different approach - we'll use it only if available and data permits
    if (!requireNamespace("limma", quietly = TRUE)) {
      warning("Limma package not available. Falling back to Welch's t-test.")
      return(perform_statistical_test(x, y, test_type = "welch", paired = paired))
    }

    result <- tryCatch(
      {
        # Create design matrix for limma
        group <- factor(c(rep("den", length(x)), rep("num", length(y))))
        design <- model.matrix(~group)

        # Combine data
        y_combined <- c(x, y)

        # Fit linear model
        fit <- limma::lmFit(matrix(y_combined, ncol = 1), design)
        fit <- limma::eBayes(fit)

        # Extract p-value from the coefficient (second column, treatment effect)
        p_val <- fit$p.value[1, 2]

        list(p_value = as.numeric(p_val), test_type_used = "limma", error = NULL)
      },
      error = function(e) {
        # Fallback to Welch if limma fails
        list(p_value = NA_real_, test_type_used = "limma", error = as.character(e))
      }
    )
    return(result)
  }

  # Default fallback
  list(p_value = NA_real_, test_type_used = test_type, error = "Unknown test type")
}

# compute_ttest_stats_general(): Computes fold changes, log2 fold changes, p-values, FDR, 
# for a specified comparison between groups in the metadata. 
# It returns a data frame with these statistics along with feature information.
# Supports multiple statistical test types: Student's t-test, Welch's t-test, Wilcoxon rank-sum test, and limma.
compute_ttest_stats_general <- function(mat_log2, mat_prelog, meta_sub, feat_info,
                                        compare_var = c("group", "sex"),
                                        num_level, den_level,
                                        statistical_test_type = "student",
                                        test_is_paired = FALSE,
                                        pvalue_correction_method = "FDR") {
  compare_var <- match.arg(compare_var)
  statistical_test_type <- resolve_statistical_test_type(statistical_test_type)
  pvalue_correction_method <- resolve_pvalue_correction_method(pvalue_correction_method)

  meta_sub <- meta_sub %>%
    filter(sample %in% rownames(mat_log2), sample %in% rownames(mat_prelog)) %>%
    mutate(.ord = match(sample, rownames(mat_log2))) %>%
    arrange(.ord) %>%
    select(-.ord)

  if (nrow(meta_sub) < 4) {
    return(NULL)
  }

  v <- meta_sub[[compare_var]]
  if (sum(v == den_level, na.rm = TRUE) < 2 || sum(v == num_level, na.rm = TRUE) < 2) {
    return(NULL)
  }

  s <- meta_sub$sample
  sub_log2 <- mat_log2[s, , drop = FALSE]
  sub_pre <- mat_prelog[s, , drop = FALSE]
  v2 <- meta_sub[[compare_var]]

  mean_den_pre <- colMeans(sub_pre[v2 == den_level, , drop = FALSE], na.rm = TRUE)
  mean_num_pre <- colMeans(sub_pre[v2 == num_level, , drop = FALSE], na.rm = TRUE)

  fc <- rep(NA_real_, length(mean_den_pre))
  ok <- is.finite(mean_den_pre) & is.finite(mean_num_pre) & mean_den_pre > 0 & mean_num_pre > 0
  fc[ok] <- mean_num_pre[ok] / mean_den_pre[ok]
  log2FC <- log2(fc)

  pvals <- rep(NA_real_, ncol(sub_log2))

  for (j in seq_len(ncol(sub_log2))) {
    x <- sub_log2[v2 == den_level, j]
    y <- sub_log2[v2 == num_level, j]

    if (all(is.na(x)) || all(is.na(y))) next
    if (length(unique(na.omit(c(x, y)))) < 2) next

    test_result <- perform_statistical_test(
      x = x,
      y = y,
      test_type = statistical_test_type,
      paired = test_is_paired
    )

    pvals[j] <- test_result$p_value
  }

  # Apply p-value correction method
  adjusted_pvals <- if (identical(pvalue_correction_method, "none")) {
    pvals
  } else {
    p.adjust(pvals, method = pvalue_correction_method)
  }

  result_df <- tibble(
    featureID = colnames(sub_log2),
    FC_num_over_den = as.numeric(fc),
    log2FC_num_over_den = as.numeric(log2FC),
    p_value = pvals,
    FDR = adjusted_pvals
  )

  # If correction method is "none" (raw p-values), add FDR column for compatibility
  if (identical(pvalue_correction_method, "none")) {
    # Compute BH FDR for compatibility even when not used for significance
    result_df <- result_df %>%
      mutate(FDR_BH = p.adjust(p_value, method = "BH"))
  }

  result_df %>%
    left_join(
      feat_info %>%
        select(
          featureID,
          any_of(c(
            "display_name", "mz", "RT", "Name", "Name_canon",
            "Metabolika_pathways", "Formula"
          ))
        ),
      by = "featureID"
    )
}

# plot_volcano_metric(): Generates a volcano plot for a given metric (FDR or p-value) using ggplot2.
plot_volcano_metric <- function(stats_df, title, out_path,
                                metric = c("FDR", "p_value"),
                                alpha = p_value_cutoff,
                                fdr_alpha = fdr_cutoff,
                                fc_cutoff_log2 = fc_cutoff_log2,
                                xlab = "log2FC",
                                style = "classic") {
  metric <- match.arg(metric)
  style <- resolve_volcano_styles(style)[1]

  # Determine which alpha to use depending on the metric plotted
  alpha_used <- if (identical(metric, "FDR")) {
    # use the FDR cutoff when plotting FDR
    fdr_alpha
  } else {
    # default: p-value cutoff
    alpha
  }

  if (is.null(stats_df) || nrow(stats_df) == 0) {
    stop("stats_df is NULL or empty.")
  }

  df <- stats_df %>%
    dplyr::mutate(
      metric_val = .data[[metric]],
      metric_plot = dplyr::if_else(
        is.na(metric_val),
        NA_real_,
        pmax(metric_val, FLOOR_P)
      ),
      minus_log10_metric = dplyr::if_else(
        is.na(metric_plot),
        NA_real_,
        -log10(metric_plot)
      ),
      label = dplyr::case_when(
        !is.na(display_name) & trimws(display_name) != "" ~ display_name,
        !is.na(Name) & trimws(Name) != "" ~ Name,
        TRUE ~ featureID
      ),

      # ---------------------------------------------------------
      # STRICT LOGIC:
      # Up/Down must respect BOTH significance and FC cutoff
      # ---------------------------------------------------------
      regulation = dplyr::case_when(
        !is.na(metric_val) &
          metric_val < alpha_used &
          !is.na(log2FC_num_over_den) &
          log2FC_num_over_den >= fc_cutoff_log2 ~ "Up",
        !is.na(metric_val) &
          metric_val < alpha_used &
          !is.na(log2FC_num_over_den) &
          log2FC_num_over_den <= -fc_cutoff_log2 ~ "Down",
        TRUE ~ "Normal"
      ),
      regulation = factor(regulation, levels = c("Down", "Normal", "Up"))
    )

  # -------------------------------------------------------------
  # Labels only for points that pass BOTH significance and FC cutoff
  # -------------------------------------------------------------
  df_labels <- df %>%
      dplyr::filter(
      !is.na(metric_val),
      !is.na(log2FC_num_over_den),
      metric_val < alpha_used,
      abs(log2FC_num_over_den) >= fc_cutoff_log2
    ) %>%
    dplyr::arrange(metric_val, dplyr::desc(abs(log2FC_num_over_den)))

  if (is.null(volcano_custom_labels)) {
    if (nrow(df_labels) > volcano_label_number) {
      df_labels <- df_labels %>% dplyr::slice_head(n = volcano_label_number)
    }
  } else {
    custom_labels <- as.character(volcano_custom_labels)
    df_labels <- df %>%
      dplyr::filter(label %in% custom_labels)
  }

  if (!isTRUE(volcano_add_labels)) {
    df_labels <- df[0, , drop = FALSE]
  }

  # -------------------------------------------------------------
  # Axis limits
  # -------------------------------------------------------------
  x_vals <- df$log2FC_num_over_den[is.finite(df$log2FC_num_over_den)]
  y_vals <- df$minus_log10_metric[is.finite(df$minus_log10_metric)]

  max_abs_x <- if (length(x_vals) > 0) max(abs(x_vals), na.rm = TRUE) else 1
  max_abs_x <- max(max_abs_x, fc_cutoff_log2, 0.5)
  if (isTRUE(volcano_auto_axis)) {
    max_abs_x <- max_abs_x * (1 + volcano_axis_expand_mult)
  }

  sig_y <- -log10(alpha_used)
  max_y <- if (length(y_vals) > 0) max(y_vals, na.rm = TRUE) else 1
  max_y <- max(max_y, sig_y, 1)
  if (isTRUE(volcano_auto_axis)) {
    max_y <- max_y * (1 + volcano_axis_expand_mult)
  }

  # -------------------------------------------------------------
  # Dummy points only for legend
  # This forces Down / Normal / Up to always appear in legend
  # -------------------------------------------------------------
  legend_dummy <- tibble::tibble(
    log2FC_num_over_den = c(-999, -999, -999),
    minus_log10_metric = c(-999, -999, -999),
    regulation = factor(c("Down", "Normal", "Up"),
      levels = c("Down", "Normal", "Up")
    )
  )

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = log2FC_num_over_den,
      y = minus_log10_metric
    )
  )

  if (isTRUE(volcano_add_cutoff_lines)) {
    p <- p +
      ggplot2::geom_hline(
        yintercept = sig_y,
        linetype = "dashed",
        linewidth = 0.5,
        color = "black"
      ) +
      ggplot2::geom_vline(
        xintercept = c(-fc_cutoff_log2, fc_cutoff_log2),
        linetype = "dashed",
        linewidth = 0.5,
        color = "black"
      )
  }

  if (identical(style, "classic")) {
    classic_fill_values <- setNames(volcano_classic_fills, c("Down", "Normal", "Up"))
    classic_color_values <- setNames(volcano_classic_colors, c("Down", "Normal", "Up"))

    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(fill = regulation, color = regulation),
        shape = volcano_classic_point_shape,
        size = volcano_classic_point_size,
        stroke = 0.25,
        alpha = 0.90,
        na.rm = TRUE
      ) +
      ggplot2::geom_point(
        data = legend_dummy,
        ggplot2::aes(
          x = log2FC_num_over_den,
          y = minus_log10_metric,
          fill = regulation,
          color = regulation
        ),
        shape = volcano_classic_point_shape,
        size = volcano_classic_point_size,
        stroke = 0.25,
        alpha = 0.90,
        inherit.aes = FALSE,
        show.legend = TRUE
      ) +
      ggplot2::scale_fill_manual(
        values = classic_fill_values,
        drop = FALSE,
        name = volcano_classic_legend_title
      ) +
      ggplot2::scale_color_manual(
        values = classic_color_values,
        drop = FALSE,
        name = volcano_classic_legend_title
      ) +
      ggplot2::guides(
        color = "none",
        fill = ggplot2::guide_legend(
          override.aes = list(
            shape = volcano_classic_point_shape,
            size = volcano_classic_point_size,
            alpha = 1
          )
        )
      )
  }

  if (identical(style, "gradual")) {
    n_grad <- max(3, as.integer(volcano_gradual_brewer_n))
    grad_fill_values <- volcano_gradual_fills
    grad_color_values <- volcano_gradual_colors

    if (isTRUE(volcano_gradual_use_RColorBrewer) &&
      requireNamespace("RColorBrewer", quietly = TRUE)) {
      max_colors <- RColorBrewer::brewer.pal.info[volcano_gradual_brewer_palette, "maxcolors"]
      safe_n <- min(max_colors, n_grad)
      pal <- RColorBrewer::brewer.pal(safe_n, volcano_gradual_brewer_palette)
      if (isTRUE(volcano_gradual_reverse_brewer)) {
        pal <- rev(pal)
      }
      grad_fill_values <- pal
      grad_color_values <- pal
    }

    grad_fill_fun <- grDevices::colorRampPalette(grad_fill_values)
    grad_color_fun <- grDevices::colorRampPalette(grad_color_values)

    grad_breaks <- as.numeric(volcano_gradual_legend_breaks)
    grad_limits <- as.numeric(volcano_gradual_legend_limits)

    df <- df %>%
      dplyr::mutate(
        sig_score = pmin(pmax(minus_log10_metric, grad_limits[1]), grad_limits[2]),
        size_score = pmin(sig_score, grad_limits[2])
      )

    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = log2FC_num_over_den,
        y = minus_log10_metric
      )
    )

    if (isTRUE(volcano_add_cutoff_lines)) {
      p <- p +
        ggplot2::geom_hline(
          yintercept = sig_y,
          linetype = "dashed",
          linewidth = 0.5,
          color = "black"
        ) +
        ggplot2::geom_vline(
          xintercept = c(-fc_cutoff_log2, fc_cutoff_log2),
          linetype = "dashed",
          linewidth = 0.5,
          color = "black"
        )
    }

    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(fill = sig_score, color = sig_score, size = size_score),
        shape = volcano_gradual_point_shape,
        stroke = 0.2,
        alpha = 0.92,
        na.rm = TRUE
      ) +
      ggplot2::scale_fill_gradientn(
        colors = grad_fill_fun(max(3, n_grad)),
        limits = grad_limits,
        breaks = grad_breaks,
        name = volcano_gradual_legend_title,
        oob = scales::squish
      ) +
      ggplot2::scale_color_gradientn(
        colors = grad_color_fun(max(3, n_grad)),
        limits = grad_limits,
        breaks = grad_breaks,
        guide = "none",
        oob = scales::squish
      ) +
      ggplot2::scale_size_continuous(
        range = volcano_gradual_point_size_range,
        guide = "none"
      )
  }

  p <- p +
    ggrepel::geom_text_repel(
      data = df_labels,
      ggplot2::aes(label = label),
      size = 3.8,
      box.padding = 0.35,
      point.padding = 0.20,
      segment.color = "grey40",
      segment.size = 0.30,
      max.overlaps = Inf,
      min.segment.length = 0,
      show.legend = FALSE
    ) +
    ggplot2::coord_cartesian(
      xlim = c(-max_abs_x, max_abs_x),
      ylim = c(0, max_y),
      expand = TRUE
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste0(if (identical(metric, "FDR")) {
        paste0("Cutoff: FDR < ", fdr_alpha)
      } else {
        paste0("Cutoff: p-value < ", alpha)
      }),
      x = xlab,
      y = paste0("-log10(", metric, ")")
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10),
      axis.title = ggplot2::element_text(face = "bold", size = 13),
      legend.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.key = ggplot2::element_rect(fill = "white", color = NA)
    )

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

  ggplot2::ggsave(
    filename = out_path,
    plot = p,
    width = 10,
    height = 7,
    dpi = 300,
    bg = "white"
  )
}

save_placeholder_volcano <- function(out_path, title, reason) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

  p <- ggplot() +
    theme_void() +
    labs(title = title, subtitle = reason) +
    theme(
      plot.title = element_text(size = 14, face = "bold", color = "black"),
      plot.subtitle = element_text(size = 12, color = "black"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )

  ggsave(
    filename = out_path,
    plot = p,
    width = 10,
    height = 7,
    dpi = 300,
    bg = "white"
  )
}

# run_all_stats_5sets_per_model(): For each model, computes statistics for predefined comparisons and 
# creates volcano plots for specified metrics (FDR, p-value, FDR_and_p_value). 
# It handles cases with insufficient samples gracefully by saving placeholder volcano plots with explanations.
run_all_stats_5sets_per_model <- function(mat_log2, mat_prelog, metadata_aligned, feat_info, paths,
                                          p_value_cutoff = p_value_cutoff, fdr_cutoff = fdr_cutoff, fc_cutoff_log2 = fc_cutoff_log2,
                                          run_metrics = run_metrics,
                                          make_volcano_plots = TRUE,
                                          volcano_style = volcano_style,
                                          comparison_configs = COMPARISON_CONFIGS,
                                          statistical_test_type = "student",
                                          test_is_paired = FALSE,
                                          pvalue_correction_method = "FDR") {
  metrics <- normalize_metric_modes(run_metrics)
  styles_to_export <- resolve_volcano_styles(volcano_style)
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
    if (prefix %in% c("sex_FvsM", "tg-f_vs_tg-m", "tg_f_vs_tg_m", "TG_FvsM")) {
      return(mp$plots$volcano_tg_f_vs_tg_m)
    }
    if (prefix %in% c("wt-f_vs_wt-m", "wt_f_vs_wt_m", "WT_FvsM")) {
      return(mp$plots$volcano_wt_f_vs_wt_m)
    }

    # fallback: try lowercase-normalized prefix
    np <- tolower(prefix)
    if (np %in% c("all", "tg_vs_wt", "tg-vs-wt")) return(mp$plots$volcano_tg_vs_wt)
    if (np %in% c("sex", "tg-f_vs_wt-f", "tg_f_vs_wt_f", "tg-m_vs_wt-m", "tg_m_vs_wt_m")) return(mp$plots$volcano_by_sex)
    if (np %in% c("sex_fvsm", "tg-f_vs_tg-m", "tg_f_vs_tg_m")) return(mp$plots$volcano_tg_f_vs_tg_m)
    if (np %in% c("wt-f_vs_wt-m", "wt_f_vs_wt_m")) return(mp$plots$volcano_wt_f_vs_wt_m)

    stop("Unknown comparison prefix: ", prefix)
  }

  for (m in models) {
    mp <- get_model_paths(paths, m)
    meta_m <- metadata_aligned %>% filter(type == "Sample", model == m)
    model_groups <- resolve_model_group_values(m)
    message("  - Stats model groups: model=", m, " | control=", model_groups$control, " | treatment=", model_groups$treatment)
    message("    - Test type: ", statistical_test_type, " | Paired: ", test_is_paired, " | P-value correction: ", pvalue_correction_method)
    out[[m]] <- list()

    for (comp_name in names(comparison_configs)) {
      cfg <- comparison_configs[[comp_name]]
      meta_sub <- cfg$meta_filter(meta_m, model_name = m)

      compare_den <- cfg$stats_den
      compare_num <- cfg$stats_num
      if (identical(cfg$stats_compare_var, "group")) {
        compare_den <- model_groups$control
        compare_num <- model_groups$treatment
      }

      st <- compute_ttest_stats_general(
        mat_log2 = mat_log2,
        mat_prelog = mat_prelog,
        meta_sub = meta_sub,
        feat_info = feat_info,
        compare_var = cfg$stats_compare_var,
        num_level = compare_num,
        den_level = compare_den,
        statistical_test_type = statistical_test_type,
        test_is_paired = test_is_paired,
        pvalue_correction_method = pvalue_correction_method
      )

      out[[m]][[comp_name]] <- st

      if (!isTRUE(make_volcano_plots)) next

      for (met in metrics) {
        out_dir <- get_volcano_dir(cfg$prefix, mp)

        comp_label <- switch(
          comp_name,
          ALL_TGvsWT = paste0(model_groups$treatment, " vs ", model_groups$control, " | sex=ALL"),
          F_TGvsWT = paste0(model_groups$treatment, " vs ", model_groups$control, " | sex=F"),
          M_TGvsWT = paste0(model_groups$treatment, " vs ", model_groups$control, " | sex=M"),
          TG_FvsM = paste0("F vs M within ", model_groups$treatment),
          WT_FvsM = paste0("F vs M within ", model_groups$control),
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
                fc_cutoff_log2 = fc_cutoff_log2,
                style = style_name
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
# with conditional formatting to highlight significant results.
export_stats_excel_by_model <- function(stats_5sets_by_model, paths, p_value_cutoff, fdr_cutoff, fc_cutoff_log2,
                                        active_variant, log_path = NULL,
                                        statistical_test_type = "student",
                                        test_is_paired = FALSE,
                                        pvalue_correction_method = "FDR") {
  comparisons <- COMPARISON_NAMES

  pretty_comparison_label <- function(comp_name, model_name, model_groups) {
    switch(
      comp_name,
      ALL_TGvsWT = paste0(model_groups$treatment, " vs ", model_groups$control, " | sex=ALL | model=", model_name),
      F_TGvsWT = paste0(model_groups$treatment, " vs ", model_groups$control, " | sex=F | model=", model_name),
      M_TGvsWT = paste0(model_groups$treatment, " vs ", model_groups$control, " | sex=M | model=", model_name),
      TG_FvsM = paste0("F vs M within ", model_groups$treatment, " | model=", model_name),
      WT_FvsM = paste0("F vs M within ", model_groups$control, " | model=", model_name),
      comp_name
    )
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

    readme <- tibble(
      field = c(
        "model", "Control group (WT)", "Treatment group (TG)", "Active_variant", "p_value_cutoff", "fdr_cutoff", "fc_cutoff_log2",
        "log2FC_definition", "Statistical_test", "Test_paired", "P_value_correction", "stats_definition", "significance_logic"
      ),
      value = c(
        m, model_groups$control, model_groups$treatment, active_variant, p_value_cutoff, fdr_cutoff, fc_cutoff_log2,
        "log2FC = log2(mean(num_prelog)/mean(den_prelog)); for FvsM: log2(F/M)",
        statistical_test_type,
        as.character(test_is_paired),
        pvalue_correction_method,
        paste0(
          "Test=", statistical_test_type, " (paired=", test_is_paired, "); ",
          "Correction method=", pvalue_correction_method
        ),
        "yellow row = (p_value<alpha OR FDR<alpha); green row = yellow row AND |log2FC|>=cutoff; red cell = p_value or FDR < alpha"
      )
    )

    openxlsx::addWorksheet(wb, "README")
    openxlsx::writeData(wb, "README", readme)

    tabs <- stats_5sets_by_model[[m]]
    sig_rows <- list()

    for (nm in comparisons) {
      df <- tabs[[nm]]
      openxlsx::addWorksheet(wb, nm)

      if (is.null(df) || nrow(df) == 0) {
        tmp <- tibble(note = "Not enough samples or no valid tests")
        openxlsx::writeData(wb, nm, tmp)
        next
      }

      df_clean <- df %>%
        rename(
          FC = FC_num_over_den,
          log2FC = log2FC_num_over_den
        ) %>%
        mutate(
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

      openxlsx::writeData(wb, nm, df_clean)
      openxlsx::freezePane(wb, nm, firstRow = TRUE)
      openxlsx::addFilter(wb, nm, rows = 1, cols = 1:ncol(df_clean))

      style_red <- openxlsx::createStyle(fgFill = "#FF0000", fontColour = "#FFFFFF")
      style_yellow <- openxlsx::createStyle(fgFill = "#FFD966")
      style_green <- openxlsx::createStyle(fgFill = "#00B050", fontColour = "#FFFFFF")

      cols_full_row <- 1:ncol(df_clean)

      green_rows <- which(df_clean$row_sig_and_fc == 1) + 1
      yellow_rows <- which(df_clean$row_sig_any == 1 & df_clean$row_sig_and_fc == 0) + 1

      if (length(yellow_rows) > 0) {
        openxlsx::addStyle(
          wb, nm,
          style = style_yellow,
          rows = yellow_rows,
          cols = cols_full_row,
          gridExpand = TRUE,
          stack = FALSE
        )
      }

      if (length(green_rows) > 0) {
        openxlsx::addStyle(
          wb, nm,
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
            wb, nm,
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
            wb, nm,
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
        openxlsx::setColWidths(wb, nm, cols = visible_cols, widths = "auto")
      }

      if (length(helper_cols) > 0) {
        openxlsx::setColWidths(wb, nm, cols = helper_cols, widths = 0)
      }
      # collect significant rows for aggregated sheet
      sig_subset <- df_clean[df_clean$row_sig_any == 1, , drop = FALSE]
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

# export_significant_metabolites_txt_by_model(): Writes one TXT list per comparison
# in metric-specific folders instead of mixing p_value and FDR.
export_significant_metabolites_txt_by_model <- function(stats_5sets_by_model, paths, p_value_cutoff, fdr_cutoff,
                                                        fc_cutoff_log2, active_variant,
                                                        log_path = NULL,
                                                        require_fc_cutoff = FALSE) {
  comparisons <- COMPARISON_NAMES
  metric_dirs <- c(p_value = "p_value", FDR = "FDR")

  for (m in names(stats_5sets_by_model)) {
    mp <- get_model_paths(paths, m)
    out_dir <- mp$exports$stats

    tabs <- stats_5sets_by_model[[m]]

    for (nm in comparisons) {
      df <- tabs[[nm]]

      for (metric in names(metric_dirs)) {
        out_path <- file.path(
          out_dir,
          "significant_metabolites",
          metric_dirs[[metric]],
          paste0("SIGNIFICANT_METABOLITES_model_", m, "_", nm, "_", metric, ".txt")
        )

        out_names <- character(0)

        if (!is.null(df) && nrow(df) > 0) {
          out_names <- df %>%
            dplyr::mutate(
              passes_metric = !is.na(.data[[metric]]) & .data[[metric]] < p_value_cutoff,
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
          dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
          readr::write_lines(out_names, out_path)

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

  message("  Significant metabolite txt exports saved in metric-specific folders under: ", out_dir)
}
