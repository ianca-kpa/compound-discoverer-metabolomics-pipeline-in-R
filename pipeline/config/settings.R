# Generated from app/global.R defaults.
# This is the app startup state; settings.R is only updated when you save or run.

# Normalization
normalization_mode <- "qc_loess"
make_qc_diagnostics <- TRUE

# Statistics thresholds
p_value_cutoff <- 0.05
fdr_cutoff <- 0.05
fc_cutoff_log2 <- 0

# Statistical analysis
comparison_mode <- "pairwise"
statistical_test_type <- "student"
test_is_paired <- FALSE
run_metrics <- "FDR_and_p_value"

# Multi-group statistics
multigroup_groups <- c()
multigroup_test <- "anova"
multigroup_pairwise_mode <- "selected"
multigroup_pairwise_pairs <- NULL

# PCA
pca_scaling <- "pareto"
pca_label_samples <- TRUE

# Heatmap
heatmap_top_n <- 50
make_heatmap_by_model <- TRUE
make_heatmap_by_model_sex <- TRUE
heatmap_scale_method <- "zscore"

# Output controls
output_level <- "full_debug"

# Volcano
make_volcano_plots <- TRUE
volcano_add_labels <- TRUE
volcano_add_cutoff_lines <- TRUE

# Feature filters
active_variant <- "QC_RSD"
rsd_thresholds <- c(20)
low_variance_filter_method <- "iqr"
low_variance_filter_fraction <- 0.2
use_only_known <- TRUE
duplicate_name_strategy <- "collapse_best_qc_rsd"
output_dir <- "output/TESTE"
use_reference_file <- TRUE
rsd_filter_type <- "QC_RSD"
use_weight_normalization <- TRUE
heatmap_rank_metrics <- "FDR_and_p_value"
pvalue_correction_method <- "FDR"
reference_col_metabolite <- ""
reference_col_ref_ion <- ""
reference_col_mz <- ""
reference_col_rt <- ""
model_allowed_groups_by_model <- c("AETA-m" = "WT, AETA-m", "APPdelETA" = "WT, APPdelETA")
comparison_group_control <- "WT"
comparison_group_treatment <- "AETA-m"
alpha_sig <- fdr_cutoff
cd_file_path <- "data/251114_HM_OP_untargeted_POS_new_appr.xlsx"
metadata_path <- "data/AETA-m_APPdelETA_weight_new.xlsx"
reference_path <- "data/MetabolitesList2023.xlsx"
