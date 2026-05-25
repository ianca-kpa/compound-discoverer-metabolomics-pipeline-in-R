# =============================================================================
# 09_pca.R
# PCA functions
# =============================================================================

# Summary:
  # This module contains functions for performing PCA on the log2-transformed assay matrix, with different
  # scaling options (none, pareto, autoscale) as specified in the settings. It includes a main driver function 
  # that generates PCA plots for each model and subsets of samples, and saves them to the appropriate output directories.

# scale_uv(): Autoscaling (unit variance scaling) - centers to mean and scales to unit variance.
scale_uv <- function(mat) {
  out <- scale(mat, center = TRUE, scale = TRUE)
  out <- as.matrix(out)
  out[is.na(out)] <- 0
  out
}

# scale_pareto(): Pareto scaling - centers to mean and scales to square root of standard deviation.
scale_pareto <- function(mat) {
  mu <- colMeans(mat, na.rm = TRUE)
  sdv <- apply(mat, 2, stats::sd, na.rm = TRUE)
  denom <- sqrt(sdv)
  denom[is.na(denom) | denom == 0] <- 1

  out <- sweep(sweep(mat, 2, mu, "-"), 2, denom, "/")
  out <- as.matrix(out)
  out[is.na(out)] <- 0
  out
}

# apply_pca_scaling(): Apply the specified scaling method to the matrix before PCA, 
# and return both the scaled matrix and a label for the method used.
apply_pca_scaling <- function(mat, method = c("pareto", "autoscale", "none")) {
  method <- match.arg(method)

  if (method == "none") {
    out <- as.matrix(mat)
    out[is.na(out)] <- 0
    return(list(mat = out, label = "none"))
  }

  if (method == "pareto") {
    return(list(mat = scale_pareto(mat), label = "pareto"))
  }

  list(mat = scale_uv(mat), label = "autoscale")
}

# -----------------------------------------------------------------------------
# Remove invalid / constant columns before PCA
# -----------------------------------------------------------------------------
prepare_matrix_for_pca <- function(mat) {
  if (is.null(mat) || nrow(mat) < 3 || ncol(mat) < 2) {
    return(NULL)
  }

  keep_cols <- apply(mat, 2, function(v) {
    vv <- stats::na.omit(v)
    length(vv) > 0 && length(unique(vv)) > 1
  })

  mat2 <- mat[, keep_cols, drop = FALSE]

  if (ncol(mat2) < 2) {
    return(NULL)
  }

  mat2 <- as.matrix(mat2)
  mat2[is.na(mat2)] <- 0
  mat2
}

pca_legend_title <- function(var_name) {
  switch(
    var_name,
    group = "Group",
    sex = "Sex",
    tools::toTitleCase(gsub("_", " ", as.character(var_name)))
  )
}

pca_color_values <- function(var_name, levels, model_name = NULL) {
  levels <- as.character(levels)
  levels <- levels[!is.na(levels) & nzchar(levels)]

  if (length(levels) == 0) {
    return(c(ALL = "#6B7280"))
  }

  vals <- setNames(rep("#6B7280", length(levels)), levels)

  if (identical(var_name, "group")) {
    groups <- resolve_model_group_values(model_name)

    # Blue for control and orange for treatment, keeping the order stable across models.
    ctrl_idx <- which(tolower(levels) == tolower(groups$control))[1]
    trt_idx <- which(tolower(levels) == tolower(groups$treatment))[1]

    if (!is.na(ctrl_idx)) vals[ctrl_idx] <- "#4EEE94"
    if (!is.na(trt_idx)) vals[trt_idx] <- "#FFA54F"

    remaining <- which(vals == "#6B7280")
    fallback <- c("#14B8A6", "#A855F7", "#EF4444", "#22C55E")
    if (length(remaining) > 0) {
      vals[remaining] <- rep_len(fallback, length(remaining))
    }

    return(vals)
  }

  if (identical(var_name, "sex")) {
    female_idx <- which(tolower(levels) == "f")[1]
    male_idx <- which(tolower(levels) == "m")[1]

    if (!is.na(female_idx)) vals[female_idx] <- "#CD0000"
    if (!is.na(male_idx)) vals[male_idx] <- "#009ACD"

    remaining <- which(vals == "#6B7280")
    fallback <- c("#6366F1", "#EC4899", "#16A34A")
    if (length(remaining) > 0) {
      vals[remaining] <- rep_len(fallback, length(remaining))
    }

    return(vals)
  }

  fallback <- c("#2563EB", "#F97316", "#14B8A6", "#A855F7")
  vals[] <- rep_len(fallback, length(levels))
  vals
}

pca_shape_values <- function(var_name, levels) {
  levels <- as.character(levels)
  levels <- levels[!is.na(levels) & nzchar(levels)]

  if (length(levels) == 0) {
    return(c(ALL = 16))
  }

  vals <- setNames(rep(16, length(levels)), levels)

  if (identical(var_name, "sex")) {
    female_idx <- which(tolower(levels) == "f")[1]
    male_idx <- which(tolower(levels) == "m")[1]
    if (!is.na(female_idx)) vals[female_idx] <- 16
    if (!is.na(male_idx)) vals[male_idx] <- 17
  } else {
    shape_seq <- c(16, 17, 15, 18, 8, 3, 7)
    vals <- setNames(rep_len(shape_seq, length(levels)), levels)
  }

  vals
}

pca_apply_legend_theme <- function(p) {
  p + ggplot2::theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.title = ggplot2::element_text(size = 10),
    legend.text = ggplot2::element_text(size = 9)
  ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(order = 1),
      fill = ggplot2::guide_legend(order = 1),
      shape = ggplot2::guide_legend(order = 2)
    )
}

# -----------------------------------------------------------------------------
# Draw and save one PCA plot
# -----------------------------------------------------------------------------
plot_one_pca_subset <- function(mat_log2 = mat_log2,
                                meta = meta,
                                out_png = out_png,
                                title_main = title_main,
                                pca_scaling = pca_scaling,
                                color_var = "group",
                                shape_var = "sex",
                                ellipse_color_var = NULL,
                                draw_ellipse = TRUE,
                                ellipse_positive = TRUE,
                                log_path = NULL) {
  meta <- meta %>%
    dplyr::filter(sample %in% rownames(mat_log2)) %>%
    dplyr::distinct(sample, .keep_all = TRUE)

  if (nrow(meta) < 3) {
    message("  - PCA skipped: fewer than 3 samples for ", basename(out_png))
    return(FALSE)
  }

  mat_sub <- mat_log2[meta$sample, , drop = FALSE]
  mat_sub <- prepare_matrix_for_pca(mat_sub)

  if (is.null(mat_sub)) {
    message("  - PCA skipped: matrix too small or without variable features for ", basename(out_png))
    return(FALSE)
  }

  sc <- apply_pca_scaling(mat_sub, method = pca_scaling)

  pca <- tryCatch(
    stats::prcomp(sc$mat, center = FALSE, scale. = FALSE),
    error = function(e) NULL
  )

  if (is.null(pca)) {
    message("  - PCA failed for ", basename(out_png))
    return(FALSE)
  }

  if (is.null(pca$x) || ncol(pca$x) < 2) {
    message("  - PCA skipped: PC1/PC2 not available for ", basename(out_png))
    return(FALSE)
  }

  var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
  pc1 <- round(100 * var_exp[1], 1)
  pc2 <- round(100 * var_exp[2], 1)

  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE]) %>%
    tibble::rownames_to_column("sample") %>%
    dplyr::left_join(meta, by = "sample")

  if (!(color_var %in% names(scores))) {
    scores[[color_var]] <- "ALL"
  }
  if (!is.null(shape_var) && length(shape_var) > 0 && !(shape_var %in% names(scores))) {
    scores[[shape_var]] <- "ALL"
  }

  model_name <- if ("model" %in% names(meta) && length(unique(meta$model)) >= 1) unique(meta$model)[1] else NULL

  scores[[color_var]] <- as.factor(scores[[color_var]])
  if (!is.null(shape_var) && length(shape_var) > 0) {
    scores[[shape_var]] <- as.factor(scores[[shape_var]])
  }

  if (identical(color_var, "group")) {
    groups <- resolve_model_group_values(model_name)
    lvl <- levels(scores[[color_var]])
    ordered <- c(groups$control, groups$treatment)
    ordered <- ordered[ordered %in% lvl]
    ordered <- c(ordered, setdiff(lvl, ordered))
    scores[[color_var]] <- factor(scores[[color_var]], levels = ordered)
  } else if (identical(color_var, "sex")) {
    lvl <- levels(scores[[color_var]])
    ordered <- c("F", "M")
    ordered <- ordered[ordered %in% lvl]
    ordered <- c(ordered, setdiff(lvl, ordered))
    scores[[color_var]] <- factor(scores[[color_var]], levels = ordered)
  }

  color_breaks <- levels(scores[[color_var]])
  color_values <- pca_color_values(color_var, color_breaks, model_name = model_name)

  # determine which variable should provide ellipse colors (allow separate mapping)
  ellipse_var <- if (!is.null(ellipse_color_var) && nzchar(as.character(ellipse_color_var))) ellipse_color_var else color_var
  if (!(ellipse_var %in% names(scores))) {
    scores[[ellipse_var]] <- scores[[color_var]]
  }
  scores[[ellipse_var]] <- as.factor(scores[[ellipse_var]])

  if (identical(ellipse_var, "group")) {
    groups <- resolve_model_group_values(model_name)
    lvl <- levels(scores[[ellipse_var]])
    ordered <- c(groups$control, groups$treatment)
    ordered <- ordered[ordered %in% lvl]
    ordered <- c(ordered, setdiff(lvl, ordered))
    scores[[ellipse_var]] <- factor(scores[[ellipse_var]], levels = ordered)
  } else if (identical(ellipse_var, "sex")) {
    lvl <- levels(scores[[ellipse_var]])
    ordered <- c("F", "M")
    ordered <- ordered[ordered %in% lvl]
    ordered <- c(ordered, setdiff(lvl, ordered))
    scores[[ellipse_var]] <- factor(scores[[ellipse_var]], levels = ordered)
  }

  ellipse_breaks <- levels(scores[[ellipse_var]])
  ellipse_values <- pca_color_values(ellipse_var, ellipse_breaks, model_name = model_name)

  # If ellipses should only be drawn for 'positive' case and the flag is FALSE,
  # when coloring by sex we prefer black points and rely on shape to distinguish sex.
  if (identical(color_var, "sex") && isFALSE(ellipse_positive)) {
    color_values[] <- "#000000"
  }

  if (!is.null(shape_var) && length(shape_var) > 0) {
    if (identical(shape_var, "group")) {
      groups <- resolve_model_group_values(model_name)
      lvl <- levels(scores[[shape_var]])
      ordered <- c(groups$control, groups$treatment)
      ordered <- ordered[ordered %in% lvl]
      ordered <- c(ordered, setdiff(lvl, ordered))
      scores[[shape_var]] <- factor(scores[[shape_var]], levels = ordered)
    } else if (identical(shape_var, "sex")) {
      lvl <- levels(scores[[shape_var]])
      ordered <- c("F", "M")
      ordered <- ordered[ordered %in% lvl]
      ordered <- c(ordered, setdiff(lvl, ordered))
      scores[[shape_var]] <- factor(scores[[shape_var]], levels = ordered)
    }
  }

  if (is.null(shape_var) || length(shape_var) == 0) {
    p <- ggplot2::ggplot(
      scores,
      ggplot2::aes(
        x = PC1,
        y = PC2,
        color = .data[[color_var]]
      )
    ) +
      ggplot2::geom_point(size = 3, alpha = 0.9) +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = sample),
        show.legend = FALSE,
        max.overlaps = 20
      ) +
      ggplot2::labs(
        title = title_main,
        x = paste0("PC1 (", pc1, "%)"),
        y = paste0("PC2 (", pc2, "%)"),
        color = pca_legend_title(color_var)
      ) +
      ggplot2::theme_minimal()

    p <- p + ggplot2::scale_color_manual(
      values = color_values,
      breaks = color_breaks,
      drop = TRUE
    )
    # add filled group ellipses if requested and group sizes allow
    if (isTRUE(draw_ellipse) && isTRUE(ellipse_positive)) {
      grp_counts <- table(scores[[color_var]])
      valid_grps <- names(grp_counts[grp_counts >= 3])
      if (length(valid_grps) > 0) {
        ell_data <- scores[scores[[color_var]] %in% valid_grps, , drop = FALSE]
        p <- p + ggplot2::stat_ellipse(
          data = ell_data,
          mapping = ggplot2::aes(x = PC1, y = PC2, fill = .data[[ellipse_var]], group = .data[[ellipse_var]]),
          inherit.aes = FALSE,
          geom = "polygon",
          alpha = 0.15,
          level = 0.95,
          show.legend = TRUE,
          type = "t"
        )
        p <- p + ggplot2::scale_fill_manual(values = ellipse_values, breaks = ellipse_breaks, name = paste0(pca_legend_title(ellipse_var), " (ellipse)"))
      }
    }
  } else {
    p <- ggplot2::ggplot(
      scores,
      ggplot2::aes(
        x = PC1,
        y = PC2,
        color = .data[[color_var]],
        shape = .data[[shape_var]]
      )
    ) +
      ggplot2::geom_point(size = 3, alpha = 0.9) +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = sample),
        show.legend = FALSE,
        max.overlaps = 20
      ) +
      ggplot2::labs(
        title = title_main,
        x = paste0("PC1 (", pc1, "%)"),
        y = paste0("PC2 (", pc2, "%)"),
        color = pca_legend_title(color_var),
        shape = pca_legend_title(shape_var)
      ) +
      ggplot2::theme_minimal()

    shape_breaks <- levels(scores[[shape_var]])
    shape_values <- pca_shape_values(shape_var, shape_breaks)

    p <- p + ggplot2::scale_color_manual(
      values = color_values,
      breaks = color_breaks,
      drop = TRUE
    ) +
      ggplot2::scale_shape_manual(
        values = shape_values,
        breaks = shape_breaks,
        drop = TRUE
      )
    # add filled group ellipses if requested and group sizes allow
    if (isTRUE(draw_ellipse) && isTRUE(ellipse_positive)) {
      grp_counts <- table(scores[[color_var]])
      valid_grps <- names(grp_counts[grp_counts >= 3])
      if (length(valid_grps) > 0) {
        ell_data <- scores[scores[[color_var]] %in% valid_grps, , drop = FALSE]
        p <- p + ggplot2::stat_ellipse(
          data = ell_data,
          mapping = ggplot2::aes(x = PC1, y = PC2, fill = .data[[ellipse_var]], group = .data[[ellipse_var]]),
          inherit.aes = FALSE,
          geom = "polygon",
          alpha = 0.12,
          level = 0.95,
          show.legend = TRUE,
          type = "t"
        )
        p <- p + ggplot2::scale_fill_manual(values = ellipse_values, breaks = ellipse_breaks, name = paste0(pca_legend_title(ellipse_var), " (ellipse)"))
      }
    }
  }

  p <- pca_apply_legend_theme(p)

  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

  ggplot2::ggsave(
    filename = out_png,
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )

  if (!is.null(log_path)) {
    log_written_object(
      log_path = log_path,
      file_path = out_png,
      object = scores,
      note = paste0(
        "PCA plot | title=", title_main,
        " | scaling=", sc$label,
        " | color=", color_var,
        " | shape=", shape_var
      )
    )
  }

  message("  ✓ PCA saved: ", out_png)
  TRUE
}

# -----------------------------------------------------------------------------
# PCA driver
# - one PCA per model for each comparison defined in COMPARISON_CONFIGS
# -----------------------------------------------------------------------------
plot_pca_per_model <- function(mat_log2 = mat_log2,
                               metadata_aligned = metadata_aligned,
                               paths = paths,
                               pca_scaling = pca_scaling,
                               ellipse_positive = TRUE,
                               log_path = NULL) {
  
  models <- metadata_aligned %>%
    dplyr::filter(type == "Sample") %>%
    dplyr::pull(model) %>%
    unique() %>%
    sort()
  n_done <- 0

  get_pca_dir_for_comparison <- function(prefix, mp) {
    # accept both legacy lowercase-style prefixes and comparison config prefixes
    if (prefix %in% c("tg_vs_wt", "tg-vs-wt", "ALL_TGvsWT", "ALL_TGvsWT")) {
      return(mp$plots$pca_tg_vs_wt)
    }
    if (prefix %in% c("tg-f_vs_wt-f", "tg_f_vs_wt_f", "tg-m_vs_wt-m", "tg_m_vs_wt_m", "F_TGvsWT", "M_TGvsWT")) {
      return(mp$plots$pca_by_sex)
    }
    if (prefix %in% c("tg-f_vs_tg-m", "tg_f_vs_tg_m", "TG_FvsM")) {
      return(mp$plots$pca_tg_f_vs_tg_m)
    }
    if (prefix %in% c("wt-f_vs_wt-m", "wt_f_vs_wt_m", "WT_FvsM")) {
      return(mp$plots$pca_wt_f_vs_wt_m)
    }

    # try lowercase-normalized fallback
    np <- tolower(prefix)
    if (np %in% c("tg_vs_wt", "tg-vs-wt")) return(mp$plots$pca_tg_vs_wt)
    if (np %in% c("tg-f_vs_wt-f", "tg_f_vs_wt_f", "tg-m_vs_wt-m", "tg_m_vs_wt_m")) return(mp$plots$pca_by_sex)
    if (np %in% c("tg-f_vs_tg-m", "tg_f_vs_tg_m")) return(mp$plots$pca_tg_f_vs_tg_m)
    if (np %in% c("wt-f_vs_wt-m", "wt_f_vs_wt_m")) return(mp$plots$pca_wt_f_vs_wt_m)

    stop("Unknown comparison prefix for PCA: ", prefix)
  }

  for (m in models) {
    mp <- get_model_paths(paths, m)
    model_groups <- resolve_model_group_values(m)
    message("  - PCA model groups: model=", m, " | control=", model_groups$control, " | treatment=", model_groups$treatment)

    # ensure only the PCA directories required by the comparison set exist
    pca_dirs <- vapply(
      COMPARISON_CONFIGS,
      function(cfg) get_pca_dir_for_comparison(cfg$prefix, mp),
      character(1)
    )

    for (d in unique(pca_dirs)) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }

    meta_model <- metadata_aligned %>%
      dplyr::filter(type == "Sample", model == m)

    for (comp_name in names(COMPARISON_CONFIGS)) {
      cfg <- COMPARISON_CONFIGS[[comp_name]]
      model_groups <- resolve_model_group_values(m)
      meta_sub <- cfg$meta_filter(meta_model, model_name = m)
      out_dir <- get_pca_dir_for_comparison(cfg$prefix, mp)

      comp_label <- switch(
        comp_name,
        tg_vs_wt = paste0(model_groups$treatment, " vs ", model_groups$control, " | sex=ALL"),
        "tg-f_vs_wt-f" = paste0(model_groups$treatment, " vs ", model_groups$control, " | sex=F"),
        "tg-m_vs_wt-m" = paste0(model_groups$treatment, " vs ", model_groups$control, " | sex=M"),
        "tg-f_vs_tg-m" = paste0("F vs M within ", model_groups$treatment),
        "wt-f_vs_wt-m" = paste0("F vs M within ", model_groups$control),
        cfg$label
      )

      # For ALL_TGvsWT produce two PCA plots: ellipse colored by group, and ellipse colored by sex
      if (identical(comp_name, "ALL_TGvsWT") || identical(cfg$prefix, "ALL_TGvsWT")) {
        out_base_group <- file.path(
          out_dir,
          paste0("PCA_ACTIVE_model_", m, "_", comp_name, "_scaling_", pca_scaling, "_ellipse_group.png")
        )
        ok1 <- plot_one_pca_subset(
          mat_log2 = mat_log2,
          meta = meta_sub,
          out_png = out_base_group,
          title_main = paste0("PCA - model=", m, " | ", comp_label, " (", pca_scaling, ") | ellipse=group"),
          pca_scaling = pca_scaling,
          color_var = cfg$pca_color_var,
          shape_var = cfg$pca_shape_var,
          ellipse_color_var = "group",
          draw_ellipse = TRUE,
          ellipse_positive = ellipse_positive,
          log_path = log_path
        )
        if (isTRUE(ok1)) n_done <- n_done + 1

        out_base_sex <- file.path(
          out_dir,
          paste0("PCA_ACTIVE_model_", m, "_", comp_name, "_scaling_", pca_scaling, "_ellipse_sex.png")
        )
        ok2 <- plot_one_pca_subset(
          mat_log2 = mat_log2,
          meta = meta_sub,
          out_png = out_base_sex,
          title_main = paste0("PCA - model=", m, " | ", comp_label, " (", pca_scaling, ") | ellipse=sex"),
          pca_scaling = pca_scaling,
          color_var = cfg$pca_color_var,
          shape_var = cfg$pca_shape_var,
          ellipse_color_var = "sex",
          draw_ellipse = TRUE,
          ellipse_positive = ellipse_positive,
          log_path = log_path
        )
        if (isTRUE(ok2)) n_done <- n_done + 1
      } else {
        title_base <- paste0("PCA - model=", m, " | ", comp_label, " (", pca_scaling, ")")

        ok <- plot_one_pca_subset(
          mat_log2 = mat_log2,
          meta = meta_sub,
          out_png = file.path(
            out_dir,
            paste0("PCA_ACTIVE_model_", m, "_", comp_name, "_scaling_", pca_scaling, ".png")
          ),
          title_main = title_base,
          pca_scaling = pca_scaling,
          color_var = cfg$pca_color_var,
          shape_var = cfg$pca_shape_var,
          ellipse_color_var = cfg$pca_ellipse_color_var,
          draw_ellipse = TRUE,
          ellipse_positive = ellipse_positive,
          log_path = log_path
        )
        if (isTRUE(ok)) n_done <- n_done + 1
      }
    }
  }

  message("  ✓ PCA plots created: ", n_done)
}
