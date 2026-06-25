# =============================================================================
# 06b_normalization_plots.R
# Normalization diagnostic plots
# =============================================================================

normalization_plot_long_df <- function(mat,
                                       method,
                                       sample_labels = NULL,
                                       method_group = NULL) {
  df <- as.data.frame(mat, check.names = FALSE)
  if (is.null(sample_labels)) {
    sample_labels <- rownames(mat)
  }
  df$sample <- sample_labels

  out <- df %>%
    tidyr::pivot_longer(
      cols = -sample,
      names_to = "feature",
      values_to = "value"
    ) %>%
    dplyr::mutate(method = method)

  if (!is.null(method_group)) {
    out <- dplyr::mutate(out, method_group = method_group)
  }

  out
}

plot_normalization_comparison <- function(assay_num_base,
                                         metadata_aligned,
                                         sample_idx,
                                         qc_idx,
                                         out_png,
                                         qc_loess_span = 0.75,
                                         qc_loess_min_qc_points = 4,
                                         injection_order = NULL,
                                         log2_offset = 1) {
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

  sample_labels <- rownames(assay_num_base)
  if (is.null(sample_labels) || any(!nzchar(sample_labels))) {
    sample_labels <- metadata_aligned$sample
  }

  add_long <- function(mat, method, method_group) {
    normalization_plot_long_df(
      mat = mat,
      method = method,
      sample_labels = sample_labels,
      method_group = method_group
    )
  }

  qc_loess <- normalize_qc_loess_ref(
    assay_num_base,
    qc_idx = qc_idx,
    qc_loess_span = qc_loess_span,
    min_qc_points = qc_loess_min_qc_points,
    injection_order = injection_order
  )$assay_num_qc_loess

  qcrsc_qc <- normalize_qcrsc_qc_ref(
    assay_num_base,
    qc_idx = qc_idx,
    min_qc_points = qc_loess_min_qc_points,
    injection_order = injection_order,
    batch = if ("batch" %in% names(metadata_aligned)) metadata_aligned$batch else NULL
  )$assay_num_qcrsc

  pqn_qc <- tryCatch(
    normalize_pqn_qc_ref(
      assay_num_base,
      qc_idx = qc_idx,
      min_qc_points = qc_loess_min_qc_points
    )$assay_num_pqn,
    error = function(e) NULL
  )

  plot_inputs <- list(
    add_long(log2_transform(assay_num_base, log2_offset), "No normalization", "Baseline"),
    add_long(log2_transform(normalize_by_weight(
      assay_num_base,
      metadata_aligned,
      sample_idx,
      use_weight_normalization = TRUE,
      stop_on_invalid_weight = FALSE,
      invalid_weight_to_NA = TRUE
    ), log2_offset), "Weight normalization", "Scale / reference"),
    if (!is.null(pqn_qc)) add_long(log2_transform(pqn_qc, log2_offset), "PQN (QC median)", "Scale / reference"),
    add_long(normalize_cyclic_loess_ref(assay_num_base, log2_offset = log2_offset), "Cyclic LOESS", "Drift correction"),
    add_long(log2_transform(qc_loess, log2_offset), "QC-LOESS", "Drift correction"),
    add_long(log2_transform(qcrsc_qc, log2_offset), "QC-RSC", "Drift correction")
  )

  plot_df <- dplyr::bind_rows(plot_inputs) %>%
    dplyr::mutate(
      sample = factor(sample, levels = sample_labels),
      method = factor(
        method,
        levels = c(
          "No normalization",
          "Weight normalization",
          "PQN (QC median)",
          "Cyclic LOESS",
          "QC-LOESS",
          "QC-RSC"
        )
      )
    )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = sample, y = value)) +
    ggplot2::geom_boxplot(fill = "#BFD7EA", color = "#1F4E79", outlier.size = 0.25, linewidth = 0.25) +
    ggplot2::facet_wrap(ggplot2::vars(method_group, method), ncol = 2) +
    ggplot2::labs(
      title = "Normalization comparison",
      subtitle = "Drift-correction methods are shown separately: cyclic LOESS, QC-LOESS, and QC-RSC.",
      x = NULL,
      y = "log2 intensity"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1, size = 6),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9)
    )

  ggplot2::ggsave(out_png, p, width = 14, height = 9, dpi = 300)

  invisible(list(path = out_png, plot = p, data = plot_df))
}

plot_tutorial_style_normalization_comparison <- function(assay_num_raw,
                                                        metadata_aligned,
                                                        sample_idx,
                                                        qc_idx,
                                                        out_png,
                                                        qc_loess_span = 0.75,
                                                        qc_loess_min_qc_points = 4,
                                                        injection_order = NULL,
                                                        log2_offset = 1,
                                                        audit_path = NULL,
                                                        sample_order_path = NULL) {
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

  sample_labels <- rownames(assay_num_raw)
  if (is.null(sample_labels) || any(!nzchar(sample_labels))) {
    sample_labels <- paste0("sample_", seq_len(nrow(assay_num_raw)))
  }

  add_long <- function(mat, method, method_group) {
    normalization_plot_long_df(
      mat = mat,
      method = method,
      sample_labels = rownames(mat),
      method_group = method_group
    )
  }

  weight_norm <- tryCatch(
    normalize_by_weight(
      assay_num_raw,
      metadata_aligned,
      sample_idx,
      use_weight_normalization = TRUE,
      stop_on_invalid_weight = FALSE,
      invalid_weight_to_NA = TRUE
    ),
    error = function(e) NULL
  )

  pqn_qc <- tryCatch(
    normalize_pqn_qc_ref(
      assay_num_raw,
      qc_idx = qc_idx,
      min_qc_points = qc_loess_min_qc_points
    )$assay_num_pqn,
    error = function(e) NULL
  )

  cyclic_loess <- tryCatch(
    normalize_cyclic_loess_ref(
      assay_num_raw,
      log2_offset = log2_offset
    ),
    error = function(e) NULL
  )

  qc_loess <- tryCatch(
    normalize_qc_loess_ref(
      assay_num_raw,
      qc_idx = qc_idx,
      qc_loess_span = qc_loess_span,
      min_qc_points = qc_loess_min_qc_points,
      injection_order = injection_order
    )$assay_num_qc_loess,
    error = function(e) NULL
  )

  qcrsc_qc <- tryCatch(
    normalize_qcrsc_qc_ref(
      assay_num_raw,
      qc_idx = qc_idx,
      min_qc_points = qc_loess_min_qc_points,
      injection_order = injection_order,
      batch = if ("batch" %in% names(metadata_aligned)) metadata_aligned$batch else NULL
    )$assay_num_qcrsc,
    error = function(e) NULL
  )

  plot_inputs <- list(
    add_long(assay_num_raw, "No normalization (raw AUC)", "Baseline"),
    if (!is.null(weight_norm)) add_long(log2_transform(weight_norm, log2_offset), "Weight normalization", "Scale / reference"),
    if (!is.null(pqn_qc)) add_long(log2_transform(pqn_qc, log2_offset), "PQN (QC median)", "Scale / reference"),
    if (!is.null(cyclic_loess)) add_long(cyclic_loess, "Cyclic LOESS", "Drift correction"),
    if (!is.null(qc_loess)) add_long(log2_transform(qc_loess, log2_offset), "QC-LOESS", "Drift correction"),
    if (!is.null(qcrsc_qc)) add_long(log2_transform(qcrsc_qc, log2_offset), "QC-RSC", "Drift correction")
  )

  plot_df <- dplyr::bind_rows(plot_inputs) %>%
    dplyr::mutate(
      sample = factor(sample, levels = sample_labels),
      method_group = factor(method_group, levels = c("Baseline", "Scale / reference", "Drift correction")),
      method = factor(
        method,
        levels = c(
          "No normalization (raw AUC)",
          "Weight normalization",
          "PQN (QC median)",
          "Cyclic LOESS",
          "QC-LOESS",
          "QC-RSC"
        )
      )
    )

  audit <- plot_df %>%
    dplyr::group_by(method, sample) %>%
    dplyr::summarise(
      n_values = sum(is.finite(value)),
      median_value = stats::median(value, na.rm = TRUE),
      iqr_value = stats::IQR(value, na.rm = TRUE),
      .groups = "drop"
    )

  if (!is.null(audit_path)) {
    write_csv_safe(audit, audit_path)
  }

  if (!is.null(sample_order_path)) {
    write_csv_safe(
      tibble::tibble(sample_order = seq_along(sample_labels), sample = sample_labels),
      sample_order_path
    )
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = sample, y = value)) +
    ggplot2::geom_boxplot(fill = "#D8E8D2", color = "#2F5D50", outlier.size = 0.25, linewidth = 0.25) +
    ggplot2::facet_wrap(~method, ncol = 3, scales = "free_y") +
    ggplot2::labs(
      title = "Tutorial-style normalization comparison",
      subtitle = "Raw panel stays in AUC; the normalized panels are shown in the cyclic LOESS log2 unit.",
      x = NULL,
      y = "Raw AUC / log2 intensity"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1, size = 6),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9)
    )

  ggplot2::ggsave(out_png, p, width = 14, height = 9, dpi = 300)

  invisible(list(path = out_png, plot = p, data = plot_df, audit = audit))
}

plot_qc_loess_weight_comparison <- function(assay_num_base,
                                            metadata_aligned,
                                            sample_idx,
                                            qc_idx,
                                            out_png,
                                            qc_loess_span = 0.75,
                                            qc_loess_min_qc_points = 4,
                                            injection_order = NULL,
                                            log2_offset = 1) {
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

  sample_labels <- rownames(assay_num_base)
  if (is.null(sample_labels) || any(!nzchar(sample_labels))) {
    sample_labels <- metadata_aligned$sample
  }

  add_long <- function(mat, method) {
    normalization_plot_long_df(
      mat = mat,
      method = method,
      sample_labels = sample_labels
    )
  }

  weight_norm <- tryCatch(
    normalize_by_weight(
      assay_num_base,
      metadata_aligned,
      sample_idx,
      use_weight_normalization = TRUE,
      stop_on_invalid_weight = FALSE,
      invalid_weight_to_NA = TRUE
    ),
    error = function(e) NULL
  )

  qc_loess_raw <- tryCatch(
    normalize_qc_loess_ref(
      assay_num_base,
      qc_idx = qc_idx,
      qc_loess_span = qc_loess_span,
      min_qc_points = qc_loess_min_qc_points,
      injection_order = injection_order
    )$assay_num_qc_loess,
    error = function(e) NULL
  )

  qc_loess_weight <- NULL
  if (!is.null(weight_norm)) {
    qc_loess_weight <- tryCatch(
      normalize_qc_loess_ref(
        weight_norm,
        qc_idx = qc_idx,
        qc_loess_span = qc_loess_span,
        min_qc_points = qc_loess_min_qc_points,
        injection_order = injection_order
      )$assay_num_qc_loess,
      error = function(e) NULL
    )
  }

  plot_inputs <- list()
  if (!is.null(qc_loess_raw)) plot_inputs <- c(plot_inputs, list(add_long(log2_transform(qc_loess_raw, log2_offset), "QC-LOESS (no weight)")))
  if (!is.null(qc_loess_weight)) plot_inputs <- c(plot_inputs, list(add_long(log2_transform(qc_loess_weight, log2_offset), "QC-LOESS (weight)")))

  if (length(plot_inputs) == 0) stop("No QC-LOESS results available to plot.")

  plot_df <- dplyr::bind_rows(plot_inputs) %>%
    dplyr::mutate(
      sample = factor(sample, levels = sample_labels),
      method = factor(method, levels = unique(method))
    )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = sample, y = value)) +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "#F7D6D0", color = "#8B1A1A", linewidth = 0.25) +
    ggplot2::geom_jitter(width = 0.15, height = 0, alpha = 0.35, size = 0.6, color = "#6E0B0B") +
    ggplot2::facet_wrap(~method, ncol = 2, scales = "free_y") +
    ggplot2::labs(
      title = "QC-LOESS: weight vs no-weight",
      subtitle = "Comparison of QC-LOESS applied to raw data and to weight-normalized data.",
      x = NULL,
      y = "log2 intensity"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1, size = 6),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9)
    )

  ggplot2::ggsave(out_png, p, width = 10, height = 5, dpi = 300)

  invisible(list(path = out_png, plot = p, data = plot_df))
}

plot_qc_qcrsc_weight_comparison <- function(assay_num_base,
                                             metadata_aligned,
                                             sample_idx,
                                             qc_idx,
                                             out_png,
                                             min_qc_points = 4,
                                             injection_order = NULL,
                                             log2_offset = 1) {
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

  sample_labels <- rownames(assay_num_base)
  if (is.null(sample_labels) || any(!nzchar(sample_labels))) {
    sample_labels <- metadata_aligned$sample
  }

  add_long <- function(mat, method) {
    normalization_plot_long_df(
      mat = mat,
      method = method,
      sample_labels = sample_labels
    )
  }

  weight_norm <- tryCatch(
    normalize_by_weight(
      assay_num_base,
      metadata_aligned,
      sample_idx,
      use_weight_normalization = TRUE,
      stop_on_invalid_weight = FALSE,
      invalid_weight_to_NA = TRUE
    ),
    error = function(e) NULL
  )

  qcrsc_raw <- tryCatch(
    normalize_qcrsc_qc_ref(
      assay_num_base,
      qc_idx = qc_idx,
      min_qc_points = min_qc_points,
      injection_order = injection_order,
      batch = if ("batch" %in% names(metadata_aligned)) metadata_aligned$batch else NULL
    )$assay_num_qcrsc,
    error = function(e) NULL
  )

  qcrsc_weight <- NULL
  if (!is.null(weight_norm)) {
    qcrsc_weight <- tryCatch(
      normalize_qcrsc_qc_ref(
        weight_norm,
        qc_idx = qc_idx,
        min_qc_points = min_qc_points,
        injection_order = injection_order,
        batch = if ("batch" %in% names(metadata_aligned)) metadata_aligned$batch else NULL
      )$assay_num_qcrsc,
      error = function(e) NULL
    )
  }

  plot_inputs <- list()
  if (!is.null(qcrsc_raw)) plot_inputs <- c(plot_inputs, list(add_long(log2_transform(qcrsc_raw, log2_offset), "QC-RSC (no weight)")))
  if (!is.null(qcrsc_weight)) plot_inputs <- c(plot_inputs, list(add_long(log2_transform(qcrsc_weight, log2_offset), "QC-RSC (weight)")))

  if (length(plot_inputs) == 0) stop("No QC-RSC results available to plot.")

  plot_df <- dplyr::bind_rows(plot_inputs) %>%
    dplyr::mutate(
      sample = factor(sample, levels = sample_labels),
      method = factor(method, levels = unique(method))
    )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = sample, y = value)) +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "#DCEBF7", color = "#1A567A", linewidth = 0.25) +
    ggplot2::geom_jitter(width = 0.15, height = 0, alpha = 0.35, size = 0.6, color = "#0B3B59") +
    ggplot2::facet_wrap(~method, ncol = 2, scales = "free_y") +
    ggplot2::labs(
      title = "QC-RSC: weight vs no-weight",
      subtitle = "Comparison of QC-RSC applied to raw data and to weight-normalized data.",
      x = NULL,
      y = "log2 intensity"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1, size = 6),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9)
    )

  ggplot2::ggsave(out_png, p, width = 10, height = 5, dpi = 300)

  invisible(list(path = out_png, plot = p, data = plot_df))
}
