# =============================================================================
# 10b_volcano_plots.R
# Volcano plot helpers
# =============================================================================

plot_volcano_metric <- function(stats_df, title, out_path,
                                metric = c("FDR", "p_value"),
                                alpha = p_value_cutoff,
                                fdr_alpha = fdr_cutoff,
                                fc_cutoff_log2 = fc_cutoff_log2,
                                xlab = "log2FC") {
  metric <- match.arg(metric)

  if (!is.null(stats_df) && "comparison_type" %in% names(stats_df) &&
      any(stats_df$comparison_type %in% "multigroup_global", na.rm = TRUE)) {
    message("  - MULTIGROUP_GLOBAL volcano skipped: a global multi-group test has no single numerator/denominator or directional log2FC.")
    return(invisible(FALSE))
  }

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
