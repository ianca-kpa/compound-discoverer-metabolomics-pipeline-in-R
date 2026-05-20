# =============================================================================
# 05_features_assay.R
# Feature table + assay build
# =============================================================================

get_mz_col <- function(df) {
  nms <- names(df)
  if ("m/z" %in% nms) return("m/z")
  if ("mz" %in% nms) return("mz")
  
  cand <- nms[str_detect(tolower(nms), "m\\s*/\\s*z|mass\\s*to\\s*charge")]
  if (length(cand) > 0) return(cand[1])
  
  NA_character_
}

get_rt_col <- function(df) {
  nms <- names(df)
  if ("RT [min]" %in% nms) return("RT [min]")
  if ("RT" %in% nms) return("RT")
  
  cand <- nms[str_detect(tolower(nms), "^rt\\b|retention")]
  if (length(cand) > 0) return(cand[1])
  
  NA_character_
}

get_metabolika_col <- function(df) {
  nms <- names(df)
  nms_low <- tolower(nms)
  
  exact_idx <- which(nms_low == "metabolika pathways")
  if (length(exact_idx) > 0) return(nms[exact_idx[1]])
  
  cand <- nms[
    str_detect(nms_low, "metabolika") &
      str_detect(nms_low, "pathway")
  ]
  if (length(cand) > 0) return(cand[1])
  
  NA_character_
}

get_ref_ion_col <- function(df) {
  nms <- names(df)
  nms_low <- tolower(nms)

  exact <- c("ref ion", "reference ion", "ref_ion", "reference_ion")
  exact_idx <- which(nms_low %in% exact)
  if (length(exact_idx) > 0) return(nms[exact_idx[1]])

  cand <- nms[str_detect(nms_low, "ref\\s*ion|reference\\s*ion")]
  if (length(cand) > 0) return(cand[1])

  NA_character_
}

build_feature_table <- function(cd_raw,
                                sanitize_names_for_exports = TRUE,
                                sanitize_mode = "greek_latin_ascii",
                                mz_digits = 4,
                                rt_digits = 2) {
  mz_col <- get_mz_col(cd_raw)
  rt_col <- get_rt_col(cd_raw)
  
  if (is.na(mz_col) || is.na(rt_col)) {
    stop("Could not detect mz/RT columns.")
  }
  
  name_col <- if ("Name" %in% names(cd_raw)) "Name" else NA_character_
  formula_col <- if ("Formula" %in% names(cd_raw)) "Formula" else NA_character_
  metabolika_col <- get_metabolika_col(cd_raw)
  ref_ion_col <- get_ref_ion_col(cd_raw)
  
  Name_clean <- if (!is.na(name_col)) clean_text(cd_raw[[name_col]]) else rep(NA_character_, nrow(cd_raw))
  Formula_clean <- if (!is.na(formula_col)) clean_text(cd_raw[[formula_col]]) else rep(NA_character_, nrow(cd_raw))
  Metabolika_clean <- if (!is.na(metabolika_col)) clean_text(cd_raw[[metabolika_col]]) else rep(NA_character_, nrow(cd_raw))
  Ref_ion_clean <- if (!is.na(ref_ion_col)) clean_text(cd_raw[[ref_ion_col]]) else rep(NA_character_, nrow(cd_raw))
  
  Formula_clean <- str_replace_all(Formula_clean, "\\s+", "")
  Name_canon <- strip_v_suffix_end(Name_clean)
  
  mz_num <- parse_num_robust(cd_raw[[mz_col]])
  rt_num <- parse_num_robust(cd_raw[[rt_col]])
  
  mz_txt <- ifelse(is.finite(mz_num), format(round(mz_num, mz_digits), nsmall = mz_digits, trim = TRUE), NA_character_)
  rt_txt <- ifelse(is.finite(rt_num), format(round(rt_num, rt_digits), nsmall = rt_digits, trim = TRUE), NA_character_)
  
  display_raw <- case_when(
    !is.na(Name_canon) ~ Name_canon,
    is.na(Name_canon) & !is.na(Formula_clean) & !is.na(mz_txt) & !is.na(rt_txt) ~ paste0(Formula_clean, "_mz", mz_txt, "_rt", rt_txt),
    is.na(Name_canon) & !is.na(Formula_clean) & !is.na(mz_txt) ~ paste0(Formula_clean, "_mz", mz_txt),
    TRUE ~ NA_character_
  )
  
  feature_tbl <- tibble(
    Name = Name_clean,
    Name_canon = Name_canon,
    `Ref ion` = Ref_ion_clean,
    Formula = Formula_clean,
    Metabolika_pathways = Metabolika_clean,
    mz = mz_num,
    RT = rt_num
  ) %>%
    mutate(
      featureID_raw = case_when(
        !is.na(mz) & !is.na(RT) ~ paste0("mz", sprintf("%.4f", mz), "_rt", sprintf("%.2f", RT)),
        TRUE ~ paste0("row", row_number())
      ),
      featureID = make.unique(featureID_raw, sep = "_r"),
      display_name = display_raw
    ) %>%
    select(-featureID_raw)
  
  if (isTRUE(sanitize_names_for_exports)) {
    feature_tbl <- feature_tbl %>%
      mutate(
        Name = sanitize_text_for_exports(Name, mode = sanitize_mode),
        Name_canon = sanitize_text_for_exports(Name_canon, mode = sanitize_mode),
        `Ref ion` = sanitize_text_for_exports(`Ref ion`, mode = sanitize_mode),
        Formula = sanitize_text_for_exports(Formula, mode = sanitize_mode),
        Metabolika_pathways = sanitize_text_for_exports(Metabolika_pathways, mode = sanitize_mode),
        display_name = sanitize_text_for_exports(display_name, mode = sanitize_mode)
      )
  }
  
  feature_tbl <- feature_tbl %>%
    mutate(display_name = make.unique(display_name, sep = "__dup"))
  
  feature_tbl
}

clean_sample_from_area_col <- function(x) {
  s <- str_remove(x, "^Area\\s*:\\s*") %>% str_trim()
  s <- str_replace(s, "(?i)\\.raw.*$", "")
  s <- str_replace(s, "\\s*\\(.*\\)$", "")
  s <- str_trim(s)
  
  if (str_detect(s, "^QC")) return(s)
  if (str_detect(s, "_")) s <- str_split_fixed(s, "_", 2)[, 1]
  
  s
}

build_assay_from_cd <- function(cd_raw, feature_tbl, metadata, paths) {
  area_cols <- names(cd_raw)[str_detect(names(cd_raw), "^Area\\s*:")]
  
  if (length(area_cols) == 0) {
    stop("No 'Area:' columns found.")
  }
  
  clean_names <- make.unique(map_chr(area_cols, clean_sample_from_area_col), sep = "__rep")
  sample_map <- tibble(area_col = area_cols, sample = clean_names)
  write_csv_safe(sample_map, file.path(paths$global$exports, "01_area_column_to_sample_map.csv"))
  
  cd_area <- cd_raw %>%
    select(all_of(area_cols)) %>%
    mutate(across(everything(), ~parse_num_robust(.x))) %>%
    rlang::set_names(clean_names)
  
  if ("ExtBlk" %in% names(cd_area)) {
    cd_area <- cd_area %>% select(-ExtBlk)
  }
  
  area_mat <- as.matrix(cd_area)
  rownames(area_mat) <- feature_tbl$featureID
  
  assay_mat <- t(area_mat)
  assay_df <- as.data.frame(assay_mat, check.names = FALSE) %>%
    tibble::rownames_to_column("sample") %>%
    mutate(sample = str_trim(sample))
  
  metadata$sample <- str_trim(metadata$sample)
  common <- intersect(assay_df$sample, metadata$sample)
  
  assay_df <- assay_df %>% filter(sample %in% common)
  metadata_aligned <- metadata %>%
    filter(sample %in% common) %>%
    slice(match(assay_df$sample, sample))
  
  write_csv_safe(metadata_aligned, file.path(paths$global$exports, "04_sampleData_aligned.csv"))
  
  qc_idx <- which(metadata_aligned$type == "QC")
  sample_idx <- which(metadata_aligned$type == "Sample")
  
  if (length(qc_idx) < 2) stop("Need at least 2 QCs.")
  if (length(sample_idx) < 2) stop("Need at least 2 biological samples.")
  
  assay_num <- as.matrix(assay_df %>% select(-sample))
  mode(assay_num) <- "numeric"
  rownames(assay_num) <- assay_df$sample
  
  list(
    assay_num_raw = assay_num,
    metadata_aligned = metadata_aligned,
    qc_idx = qc_idx,
    sample_idx = sample_idx
  )
}