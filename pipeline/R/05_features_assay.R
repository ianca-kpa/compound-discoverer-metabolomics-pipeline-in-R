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
                                rt_digits = 2,
                                remove_greek_letters_from_names = TRUE,
                                remove_isomer_descriptors_from_names = TRUE) {
  mz_col <- get_mz_col(cd_raw)
  rt_col <- get_rt_col(cd_raw)
  
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
  Name_canon <- standardize_metabolite_name(
    Name_canon,
    remove_greek_letters = remove_greek_letters_from_names,
    remove_isomer_descriptors = remove_isomer_descriptors_from_names
  )
  
  mz_num <- if (!is.na(mz_col)) parse_num_robust(cd_raw[[mz_col]]) else rep(NA_real_, nrow(cd_raw))
  rt_num <- if (!is.na(rt_col)) parse_num_robust(cd_raw[[rt_col]]) else rep(NA_real_, nrow(cd_raw))
  
  mz_txt <- ifelse(is.finite(mz_num), format(round(mz_num, mz_digits), nsmall = mz_digits, trim = TRUE), NA_character_)
  rt_txt <- ifelse(is.finite(rt_num), format(round(rt_num, rt_digits), nsmall = rt_digits, trim = TRUE), NA_character_)
  row_id <- paste0("row", seq_len(nrow(cd_raw)))
  
  display_raw <- case_when(
    !is.na(Name_canon) ~ Name_canon,
    is.na(Name_canon) & !is.na(Formula_clean) & !is.na(mz_txt) & !is.na(rt_txt) ~ paste0(Formula_clean, "_mz", mz_txt, "_rt", rt_txt),
    is.na(Name_canon) & !is.na(Formula_clean) & !is.na(mz_txt) ~ paste0(Formula_clean, "_mz", mz_txt),
    TRUE ~ row_id
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

is_primary_ms_area_col <- function(x) {
  s <- str_remove(x, "^Area\\s*:\\s*") %>% str_trim()
  s <- str_replace(s, "(?i)\\.raw.*$", "")
  s <- str_replace(s, "\\s*\\(.*\\)$", "")
  s <- str_trim(s)

  str_detect(s, "(^|_)MS(_|$)") && !str_detect(s, regex("ddMS|MS2", ignore_case = TRUE))
}

build_assay_from_cd <- function(cd_raw,
                                feature_tbl,
                                metadata,
                                paths,
                                require_qc = TRUE,
                                require_injection_order = FALSE,
                                export_intermediate_tables = FALSE) {
  area_cols_all <- names(cd_raw)[str_detect(names(cd_raw), "^Area\\s*:")]
  
  if (length(area_cols_all) == 0) {
    stop("No 'Area:' columns found.")
  }

  area_clean_all <- map_chr(area_cols_all, clean_sample_from_area_col)
  area_is_primary_ms <- map_lgl(area_cols_all, is_primary_ms_area_col)
  area_has_primary_ms <- area_clean_all %in% names(which(tapply(area_is_primary_ms, area_clean_all, any)))
  area_is_duplicate_sample <- duplicated(area_clean_all) | duplicated(area_clean_all, fromLast = TRUE)
  area_keep <- !(area_is_duplicate_sample & area_has_primary_ms & !area_is_primary_ms)
  area_cols <- area_cols_all[area_keep]
  
  clean_names <- make.unique(map_chr(area_cols, clean_sample_from_area_col), sep = "__rep")
  sample_map <- tibble(area_col = area_cols, sample = clean_names)
  if (isTRUE(export_intermediate_tables)) {
    write_csv_safe(sample_map, file.path(paths$global$exports, "01_area_column_to_sample_map.csv"))
  }
  
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

  injection_order_config_path <- if (exists("injection_order_path", inherits = TRUE)) {
    trimws(as.character(get("injection_order_path", inherits = TRUE))[1])
  } else {
    ""
  }

  input_order_candidate_paths <- unique(c(
    if (nzchar(injection_order_config_path)) injection_order_config_path else character(0),
    file.path("data", "Input Files.xlsx"),
    list.files("data", pattern = "^.*inputfiles.*\\.xlsx$", full.names = TRUE, ignore.case = TRUE),
    list.files("data", pattern = "^input_order.*\\.xlsx$", full.names = TRUE, ignore.case = TRUE)
  ))
  input_order_candidate_paths <- input_order_candidate_paths[file.exists(input_order_candidate_paths)]
  input_order_candidates <- lapply(input_order_candidate_paths, function(path) {
    ref <- read_input_files_reference(path)
    if (is.null(ref) || nrow(ref) == 0) {
      return(NULL)
    }
    ref$input_order_path <- path
    ref
  })
  input_order_candidates <- Filter(Negate(is.null), input_order_candidates)

  input_files_ref <- NULL
  if (length(input_order_candidates) > 0) {
    match_scores <- vapply(
      input_order_candidates,
      function(ref) sum(metadata_aligned$sample %in% ref$sample),
      numeric(1)
    )
    best_idx <- which.max(match_scores)
    if (length(best_idx) > 0 && is.finite(match_scores[best_idx]) && match_scores[best_idx] > 0) {
      input_files_ref <- input_order_candidates[[best_idx]]
      message(
        "  - Injection order reference selected: ",
        unique(input_files_ref$input_order_path)[1],
        " (matched ",
        match_scores[best_idx],
        " aligned sample(s))."
      )
    }
  }

  drift_injection_order <- rep(NA_real_, nrow(metadata_aligned))
  drift_injection_order_source <- if (isTRUE(require_injection_order)) {
    "required_real_injection_order_not_available"
  } else {
    "fallback_aligned_row_order"
  }
  missing_injection_order_samples <- metadata_aligned$sample
  if (!is.null(input_files_ref)) {
    order_map <- input_files_ref %>%
      mutate(sample_key = clean_input_file_sample_name(sample)) %>%
      distinct(sample_key, .keep_all = TRUE)

    drift_injection_order <- order_map$input_order[match(metadata_aligned$sample, order_map$sample_key)]
    missing_injection_order_samples <- metadata_aligned$sample[!is.finite(drift_injection_order)]
    if (length(missing_injection_order_samples) == 0) {
      drift_injection_order_source <- paste0(
        "input_files_reference:",
        unique(input_files_ref$input_order_source)[1]
      )
    }
  }

  if (any(!is.finite(drift_injection_order))) {
    if (isTRUE(require_injection_order)) {
      details <- if (is.null(input_files_ref)) {
          "No usable injection order file was found."
      } else {
        paste0(
          "Missing injection order for ",
          length(missing_injection_order_samples),
          " aligned sample(s): ",
          paste(head(missing_injection_order_samples, 12), collapse = ", "),
          if (length(missing_injection_order_samples) > 12) ", ..." else "",
          "."
        )
      }

        stop(
          "Real injection order is required and must be complete for normalization_mode = 'qcrsc' or 'qc_loess'. ",
          details,
          " Please provide data/Input Files.xlsx or data/input_order*.xlsx with all aligned sample names and a real injection order ",
          "(either an Order column or file creation date/time that can be sorted). Fallback to aligned row order is only allowed ",
          "for modes without drift correction."
        )
      }

    drift_injection_order <- seq_len(nrow(metadata_aligned))
    message("  - Drift-correction input-file reference not found or incomplete; falling back to current aligned row order.")
  }

  drift_order_debug <- tibble(
    sample = metadata_aligned$sample,
    type = metadata_aligned$type,
    injection_order = drift_injection_order
  )

  if (isTRUE(export_intermediate_tables)) {
    write_csv_safe(metadata_aligned, file.path(paths$global$exports, "04_sampleData_aligned.csv"))
    write_csv_safe(drift_order_debug, file.path(paths$global$exports, "00_drift_injection_order_used.csv"))
  }
  
  qc_idx <- which(metadata_aligned$type == "QC")
  sample_idx <- which(metadata_aligned$type == "Sample")

  if (isTRUE(require_qc) && length(qc_idx) < 2) stop("Need at least 2 QCs.")
  if (length(sample_idx) < 2) stop("Need at least 2 biological samples.")
  
  assay_num <- as.matrix(assay_df %>% select(-sample))
  mode(assay_num) <- "numeric"
  rownames(assay_num) <- assay_df$sample
  
  list(
    assay_num_raw = assay_num,
    metadata_aligned = metadata_aligned,
    qc_idx = qc_idx,
    sample_idx = sample_idx,
    drift_injection_order = drift_injection_order,
    drift_injection_order_source = drift_injection_order_source,
    missing_injection_order_samples = missing_injection_order_samples
  )
}
