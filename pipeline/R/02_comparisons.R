# =============================================================================
# 02_comparisons.R
# Comparison definitions
# =============================================================================

comparison_group_control <- get0("comparison_group_control", ifnotfound = "WT", inherits = TRUE)
comparison_group_treatment <- get0("comparison_group_treatment", ifnotfound = "TG", inherits = TRUE)

comparison_group_levels <- c(comparison_group_control, comparison_group_treatment)

COMPARISON_CONFIGS <- list(
  ALL_TGvsWT = list(
    meta_filter = function(m, model_name = NULL) {
      groups <- resolve_model_group_values(model_name)
      dplyr::filter(m, group %in% c(groups$control, groups$treatment))
    },
    stats_compare_var = "group",
    stats_den = comparison_group_control,
    stats_num = comparison_group_treatment,
    pca_color_var = "group",
    pca_shape_var = "sex",
    pca_ellipse_color_var = "group",
    prefix = "ALL_TGvsWT",
    label = paste0(comparison_group_treatment, " vs ", comparison_group_control, " | sex=ALL")
  ),
  
  "F_TGvsWT" = list(
    meta_filter = function(m, model_name = NULL) {
      groups <- resolve_model_group_values(model_name)
      dplyr::filter(m, sex == "F", group %in% c(groups$control, groups$treatment))
    },
    stats_compare_var = "group",
    stats_den = comparison_group_control,
    stats_num = comparison_group_treatment,
    pca_color_var = "group",
    pca_shape_var = NULL,
    prefix = "F_TGvsWT",
    label = paste0(comparison_group_treatment, " vs ", comparison_group_control, " | sex=F")
  ),
  
  "M_TGvsWT" = list(
    meta_filter = function(m, model_name = NULL) {
      groups <- resolve_model_group_values(model_name)
      dplyr::filter(m, sex == "M", group %in% c(groups$control, groups$treatment))
    },
    stats_compare_var = "group",
    stats_den = comparison_group_control,
    stats_num = comparison_group_treatment,
    pca_color_var = "group",
    pca_shape_var = NULL,
    prefix = "M_TGvsWT",
    label = paste0(comparison_group_treatment, " vs ", comparison_group_control, " | sex=M")
  ),
  
  "TG_FvsM" = list(
    meta_filter = function(m, model_name = NULL) {
      groups <- resolve_model_group_values(model_name)
      dplyr::filter(m, group == groups$treatment, sex %in% c("F", "M"))
    },
    stats_compare_var = "sex",
    stats_den = "M",
    stats_num = "F",
    pca_color_var = "sex",
    pca_shape_var = NULL,
    prefix = "TG_FvsM",
    label = paste0("F vs M within ", comparison_group_treatment)
  ),
  
  "WT_FvsM" = list(
    meta_filter = function(m, model_name = NULL) {
      groups <- resolve_model_group_values(model_name)
      dplyr::filter(m, group == groups$control, sex %in% c("F", "M"))
    },
    stats_compare_var = "sex",
    stats_den = "M",
    stats_num = "F",
    pca_color_var = "sex",
    pca_shape_var = NULL,
    prefix = "WT_FvsM",
    label = paste0("F vs M within ", comparison_group_control)
  )
)

COMPARISON_NAMES <- names(COMPARISON_CONFIGS)
