# =============================================================================
# 07_duplicates.R
# Duplicate handling
# =============================================================================

# Functions for handling duplicate features in the feature table, based on the
# settings specified in config/settings.R.
#
# The main function is collapse_duplicate_names(), which takes the assay matrix
# and feature table, identifies duplicates based on the "Name_canon" column,
# and handles them according to the specified strategy:
# - reference_or_best_qc_rsd
# - keep_separate
# - collapse_mean
# - collapse_sum
# - collapse_best_qc_rsd
#
# An audit table is also generated to show how duplicates were handled.

utils::globalVariables(c(
  "Name_canon",
  "display_name",
  "featureID",
  ".name_key",
  "is_named",
  "dup_key",
  "n_in_group",
  ".name_norm",
  ".name_keys",
  ".refion_norm",
  ".mz_round",
  ".rt_round",
  ".rt_raw"
))

collapse_duplicate_names <- function(mat, feature_tbl,
                                     strategy = c("reference_or_best_qc_rsd", "keep_separate", "collapse_mean", "collapse_sum", "collapse_best_qc_rsd"),
                                     reference_tbl = NULL,
                                     reference_col_overrides = NULL,
                                     sanitize_mode = "greek_latin_ascii",
                                     qc_rsd = NULL,
                                     audit_path = NULL) {
  strategy <- match.arg(strategy)
  orig_feature_cols <- colnames(feature_tbl)

  # ---------------------------------------------------------------------------
  # Basic validation
  # ---------------------------------------------------------------------------
  required_cols <- c("featureID", "Name_canon", "display_name")
  if (strategy == "reference_or_best_qc_rsd") {
    required_cols <- c(required_cols, "mz", "RT")
  }
  missing_cols <- setdiff(required_cols, names(feature_tbl))
  if (length(missing_cols) > 0) {
    stop(
      "collapse_duplicate_names(): feature_tbl is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (is.null(colnames(mat))) {
    stop("collapse_duplicate_names(): mat must have column names matching featureID.")
  }

  if (anyDuplicated(feature_tbl$featureID) > 0) {
    dup_ids <- unique(feature_tbl$featureID[duplicated(feature_tbl$featureID)])
    stop(
      "collapse_duplicate_names(): feature_tbl has duplicated featureID values: ",
      paste(utils::head(dup_ids, 10), collapse = ", "),
      if (length(dup_ids) > 10) " ..." else ""
    )
  }

  if (anyDuplicated(colnames(mat)) > 0) {
    dup_cols <- unique(colnames(mat)[duplicated(colnames(mat))])
    stop(
      "collapse_duplicate_names(): mat has duplicated column names: ",
      paste(utils::head(dup_cols, 10), collapse = ", "),
      if (length(dup_cols) > 10) " ..." else ""
    )
  }

  missing_in_mat <- setdiff(feature_tbl$featureID, colnames(mat))
  if (length(missing_in_mat) > 0) {
    stop(
      "collapse_duplicate_names(): some featureIDs from feature_tbl are missing in mat: ",
      paste(utils::head(missing_in_mat, 10), collapse = ", "),
      if (length(missing_in_mat) > 10) " ..." else ""
    )
  }

  extra_in_mat <- setdiff(colnames(mat), feature_tbl$featureID)
  if (length(extra_in_mat) > 0) {
    stop(
      "collapse_duplicate_names(): mat has columns not present in feature_tbl$featureID: ",
      paste(utils::head(extra_in_mat, 10), collapse = ", "),
      if (length(extra_in_mat) > 10) " ..." else ""
    )
  }

  make_unique_id <- function(base_id, used_ids) {
    candidate <- base_id
    i <- 1
    while (candidate %in% used_ids) {
      candidate <- paste0(base_id, "_r", i)
      i <- i + 1
    }
    candidate
  }

  normalize_name_key <- function(x) {
    if (exists("sanitize_text_for_exports", mode = "function")) {
      x <- sanitize_text_for_exports(x, mode = sanitize_mode)
    }
    x <- as.character(x)
    x <- trimws(tolower(x))
    x[x %in% c("", "na", "nan", "null")] <- NA_character_
    x
  }

  ft <- feature_tbl |>
    dplyr::mutate(
      .name_key = normalize_name_key(Name_canon),
      is_named = !is.na(.name_key) & .name_key != "",
      dup_key = dplyr::if_else(is_named, .name_key, featureID)
    )

  audit <- ft |>
    dplyr::group_by(dup_key) |>
    dplyr::mutate(n_in_group = dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      will_collapse = (is_named & n_in_group > 1),
      duplicate_action = dplyr::case_when(
        strategy == "keep_separate" ~ "kept_separate",
        !is_named ~ "kept_unnamed",
        is_named & n_in_group == 1 ~ "kept_single_named",
        strategy == "reference_or_best_qc_rsd" & is_named & n_in_group > 1 ~ "selected_reference_or_best_qc_rsd",
        strategy == "collapse_best_qc_rsd" & is_named & n_in_group > 1 ~ "selected_best_qc_rsd",
        strategy %in% c("collapse_mean", "collapse_sum") & is_named & n_in_group > 1 ~ "collapsed",
        TRUE ~ "unknown"
      )
    )

  # ---------------------------------------------------------------------------
  # Keep all features separate
  # ---------------------------------------------------------------------------
  if (strategy == "keep_separate") {
    mat2 <- mat[, feature_tbl$featureID, drop = FALSE]

    ft2 <- ft |>
      dplyr::mutate(display_name = make.unique(display_name, sep = "_dup")) |>
      dplyr::select(dplyr::all_of(orig_feature_cols))

    if (!is.null(audit_path)) write_csv_safe(audit, audit_path)

    return(list(mat = mat2, feature = ft2, audit = audit))
  }

  named_groups <- split(audit$featureID[audit$is_named], audit$dup_key[audit$is_named])
  unnamed_ids <- audit$featureID[!audit$is_named]

  keep_named_ids <- character(0)
  collapsed_vecs <- list()
  collapsed_feats <- list()
  used_ids <- feature_tbl$featureID

  # ---------------------------------------------------------------------------
  # Strategy: use reference file or best QC RSD to choose representative
  # ---------------------------------------------------------------------------
if (strategy == "reference_or_best_qc_rsd") {
  normalize_text <- function(x) {
    if (exists("sanitize_text_for_exports", mode = "function")) {
      x <- sanitize_text_for_exports(x, mode = sanitize_mode)
    }
    x <- as.character(x)
    x <- trimws(tolower(x))
    x[x %in% c("", "na", "nan", "null")] <- NA_character_
    x
  }

  normalize_compact_text <- function(x) {
    x <- normalize_text(x)
    if (length(x) == 0 || all(is.na(x))) {
      return(character(0))
    }

    x <- gsub("[’'`]+", "", x, perl = TRUE)
    x <- gsub("\\((?:\\+|-)\\)", "", x, perl = TRUE)
    x <- gsub("[\\[\\]{}]", "", x, perl = TRUE)
    x <- gsub("^\\d+[a-z]?\\s*[-_:.\\s]+", "", x, perl = TRUE)
    x <- gsub("(^|[\\s\\-_/,])(dl|d|l)(?=([\\s\\-_/,]|$))", "\\1", x, perl = TRUE)
    x <- gsub("\\b([+-])\\b", "", x, perl = TRUE)
    x <- gsub("\\s+", " ", x, perl = TRUE)
    x <- trimws(x)
    x <- gsub("[^a-z0-9]+", "", x, perl = TRUE)
    x <- x[!is.na(x) & nzchar(x)]

    if (length(x) == 0) {
      return(character(0))
    }

    x
  }

  expand_compact_synonyms <- function(keys) {
    if (length(keys) == 0) {
      return(character(0))
    }

    keys <- unique(keys[!is.na(keys) & nzchar(keys)])
    if (length(keys) == 0) {
      return(character(0))
    }

    synonym_map <- c(
      "adenosine5monophosphate" = "amp",
      "adenosine5diphosphate" = "adp",
      "adenosine5triphosphate" = "atp",
      "inosine5monophosphate" = "imp",
      "nicotinamideadeninedinucleotide" = "nad",
      "nicotinamideadeninedinucleotideh" = "nadh",
      "acetyllcarnitine" = "acetylcarnitine",
      "glutamicacid" = "glutamate",
      "asparticacid" = "aspartate"
    )

    mapped <- unname(synonym_map[keys])
    mapped <- mapped[!is.na(mapped) & nzchar(mapped)]

    acid_to_ate <- character(0)
    acid_keys <- keys[grepl("icacid$", keys, perl = TRUE)]
    if (length(acid_keys) > 0) {
      acid_to_ate <- paste0(sub("icacid$", "", acid_keys, perl = TRUE), "ate")
      acid_to_ate <- acid_to_ate[nzchar(acid_to_ate)]
    }

    unique(c(keys, mapped, acid_to_ate))
  }

  # Build comparable name variants to improve matching against reference names.
  # Includes: base normalized form, without leading stereo/number prefixes,
  # and acronym/full-name variants when present in parentheses.
  build_name_keys <- function(x) {
    base <- normalize_text(x)
    if (is.na(base) || !nzchar(base)) {
      return(character(0))
    }

    keys <- base

    stripped <- base
    stripped <- gsub("^(?:l|d|dl)\\s*[-_\\s]+", "", stripped, perl = TRUE)
    stripped <- gsub("^\\d+\\s*[-_:.\\s]+", "", stripped, perl = TRUE)
    if (nzchar(stripped)) {
      keys <- c(keys, stripped)
    }

    paren <- regmatches(base, gregexpr("\\([^()]+\\)", base, perl = TRUE))[[1]]
    if (length(paren) > 0) {
      paren_clean <- trimws(gsub("^\\(|\\)$", "", paren))
      paren_clean <- normalize_text(paren_clean)
      paren_clean <- paren_clean[!is.na(paren_clean) & nzchar(paren_clean)]
      if (length(paren_clean) > 0) {
        keys <- c(keys, paren_clean)
      }
    }

    no_paren <- trimws(gsub("\\s*\\([^)]*\\)", "", base, perl = TRUE))
    no_paren <- normalize_text(no_paren)
    if (!is.na(no_paren) && nzchar(no_paren)) {
      keys <- c(keys, no_paren)
    }

    stripped_all <- gsub("^(?:l|d|dl)\\s*[-_\\s]+", "", keys, perl = TRUE)
    stripped_all <- gsub("^\\d+\\s*[-_:.\\s]+", "", stripped_all, perl = TRUE)
    keys <- c(keys, stripped_all)

    stereo_free <- keys
    stereo_free <- gsub("\\((?:\\+|-)\\)", "", stereo_free, perl = TRUE)
    stereo_free <- gsub("[\\[\\]{}]", "", stereo_free, perl = TRUE)
    stereo_free <- gsub("(^|[\\s\\-_/,])(dl|d|l)(?=([\\s\\-_/,]|$))", "\\1", stereo_free, perl = TRUE)
    stereo_free <- gsub("\\s+", " ", stereo_free, perl = TRUE)
    stereo_free <- trimws(stereo_free)
    stereo_free <- stereo_free[nzchar(stereo_free)]
    keys <- c(keys, stereo_free)

    compact_keys <- unlist(lapply(keys, normalize_compact_text), use.names = FALSE)
    compact_keys <- compact_keys[!is.na(compact_keys) & nzchar(compact_keys)]
    if (length(compact_keys) > 0) {
      keys <- c(keys, compact_keys)
      keys <- c(keys, expand_compact_synonyms(compact_keys))
    }

    keys <- unique(trimws(keys))
    keys[!is.na(keys) & nzchar(keys)]
  }

  parse_num_local <- function(x) {
    x <- as.character(x)
    x <- gsub(",", ".", x, fixed = TRUE)
    suppressWarnings(as.numeric(x))
  }

  normalize_header <- function(x) {
    x <- as.character(x)
    if (requireNamespace("stringi", quietly = TRUE)) {
      x <- stringi::stri_trans_general(x, "Any-Latin; Latin-ASCII")
    } else {
      x <- iconv(x, from = "", to = "ASCII//TRANSLIT", sub = "")
    }
    tolower(gsub("[^a-z0-9]+", "", x))
  }

  find_col_by_pattern <- function(df, patterns, fallback_first = FALSE) {

    nms <- names(df)
    nms_norm <- normalize_header(nms)

    # First try exact normalized aliases (useful when headers carry punctuation/spacing variants).
    pat_norm <- normalize_header(patterns)
    idx_exact <- which(nms_norm %in% pat_norm[nzchar(pat_norm)])
    if (length(idx_exact) > 0) {
      return(nms[idx_exact[1]])
    }

    for (p in patterns) {
      idx <- which(grepl(p, nms_norm, perl = TRUE))
      if (length(idx) > 0) return(nms[idx[1]])
    }
    if (fallback_first && length(nms) > 0) return(nms[1])
    NA_character_
  }

  resolve_override_col <- function(df, override_value) {
    if (is.null(override_value) || length(override_value) != 1 || is.na(override_value)) {
      return(NA_character_)
    }

    override_value <- trimws(as.character(override_value))
    if (!nzchar(override_value)) {
      return(NA_character_)
    }

    nms <- names(df)
    if (override_value %in% nms) {
      return(override_value)
    }

    nms_norm <- normalize_header(nms)
    override_norm <- normalize_header(override_value)
    idx <- which(nms_norm == override_norm)
    if (length(idx) > 0) {
      return(nms[idx[1]])
    }

    NA_character_
  }

  choose_best_qc <- function(ids, qc_rsd_vec) {
    if (!is.null(qc_rsd_vec) && !is.null(names(qc_rsd_vec))) {
      rsd <- unname(qc_rsd_vec[ids])
      rsd_sort <- ifelse(is.finite(rsd), rsd, Inf)
      return(ids[order(rsd_sort, ids)][1])
    }
    ids[order(ids)][1]
  }

  mz_digits_local <- if (exists("dup_mz_digits", inherits = TRUE)) {
    as.integer(get("dup_mz_digits", inherits = TRUE))
  } else {
    4L
  }

  rt_digits_local <- if (exists("dup_rt_digits", inherits = TRUE)) {
    as.integer(get("dup_rt_digits", inherits = TRUE))
  } else {
    2L
  }

  if (is.null(reference_tbl)) {
    reference_path_local <- if (exists("reference_path", inherits = TRUE)) {
      get("reference_path", inherits = TRUE)
    } else {
      NULL
    }

    reference_sheet_local <- if (exists("reference_sheet", inherits = TRUE)) {
      get("reference_sheet", inherits = TRUE)
    } else {
      1
    }

    if (!is.null(reference_path_local) &&
        nzchar(as.character(reference_path_local)) &&
        exists("read_any_table", mode = "function")) {
      reference_tbl <- tryCatch(
        read_any_table(reference_path_local, reference_sheet_local),
        error = function(e) NULL
      )
    }
  }

  if (!is.null(qc_rsd) && is.null(names(qc_rsd))) {
    stop("reference_or_best_qc_rsd requires qc_rsd to be a named vector by featureID when provided.")
  }

  if (is.null(reference_tbl) || nrow(reference_tbl) == 0) {
    has_qc_rsd <- !is.null(qc_rsd) && !is.null(names(qc_rsd))
    fallback_msg <- if (has_qc_rsd) {
      "reference_or_best_qc_rsd: no reference table found. Falling back to best QC RSD for all duplicate groups."
    } else {
      "reference_or_best_qc_rsd: no reference table found and qc_rsd is unavailable. Falling back to deterministic featureID ordering for duplicate groups."
    }
    warning(
      fallback_msg,
      call. = FALSE
    )
  }

  ft_local <- feature_tbl
  ft_local$.name_norm <- normalize_text(ft_local$Name_canon)
  ft_local$.name_keys <- lapply(ft_local$Name_canon, build_name_keys)
  ft_local$.mz_round <- round(ft_local$mz, digits = mz_digits_local)
  ft_local$.rt_round <- round(ft_local$RT, digits = rt_digits_local)
  ft_local$.rt_raw <- parse_num_local(ft_local$RT)

  ref_proc <- NULL
  if (!is.null(reference_tbl) && nrow(reference_tbl) > 0) {
    ref_nms <- names(reference_tbl)
    if (requireNamespace("stringi", quietly = TRUE)) {
      ref_nms_norm <- stringi::stri_trans_general(ref_nms, "Any-Latin; Latin-ASCII")
    } else {
      ref_nms_norm <- iconv(ref_nms, from = "", to = "ASCII//TRANSLIT", sub = "")
    }
    ref_nms_norm <- tolower(gsub("[^a-z0-9]+", "", ref_nms_norm))

      metab_col <- resolve_override_col(reference_tbl, reference_col_overrides$metabolite)
      if (is.na(metab_col)) {
        metab_col <- find_col_by_pattern(
          reference_tbl,
          c("^metabolite$", "^compound$", "^name$", "metabolite", "compound", "Name"),
          fallback_first = TRUE
        )
      }

    refion_col <- resolve_override_col(reference_tbl, reference_col_overrides$ref_ion)
    if (is.na(refion_col)) {
      refion_col <- find_col_by_pattern(
        reference_tbl,
        c("^refion$", "^referenceion$", "ref-ion", "referenceion", "^ion$", "Ref Ion", "Reference Ion", "ref_ion")
      )
    }

    if (is.na(refion_col)) {
      idx_refion <- which(ref_nms_norm %in% c("refion", "referenceion", "refions", "referenceions", "efion", "eferenceion"))
      if (length(idx_refion) > 0) {
        refion_col <- ref_nms[idx_refion[1]]
      }
    }

    mz_col <- resolve_override_col(reference_tbl, reference_col_overrides$mz)
    if (is.na(mz_col)) {
      mz_col <- find_col_by_pattern(
        reference_tbl,
        c("^mz$", "mz", "masstocharge", "moverz", "masscharge", "m/z")
      )
    }

    rt_col <- resolve_override_col(reference_tbl, reference_col_overrides$rt)
    if (is.na(rt_col)) {
      rt_col <- find_col_by_pattern(
        reference_tbl,
        c("^rt", "rt", "^rtmin$", "retentiontime", "retention", "Retention Time", "RT [min]")
      )
    }

    if (is.na(rt_col)) {
      idx_rt <- which(
        grepl("^rt", ref_nms_norm) |
          grepl("retentiontime", ref_nms_norm) |
          grepl("retention", ref_nms_norm) |
          grepl("min", ref_nms_norm)
      )
      if (length(idx_rt) > 0) {
        rt_col <- ref_nms[idx_rt[1]]
      }
    }

    step_info("Reference column detected for metabolite: ", metab_col)
    step_info("Reference column detected for Ref ion: ", refion_col)
    step_info("Reference column detected for m/z: ", mz_col)
    step_info("Reference column detected for RT: ", rt_col)

    if (is.na(refion_col) || is.na(rt_col)) {
      step_info("Reference normalized headers: ", paste(ref_nms_norm, collapse = " | "))
    }

    if (is.na(metab_col) || is.na(rt_col)) {
      warning(
        "reference_or_best_qc_rsd: could not detect required columns in reference table (metabolite/name and RT). Falling back to best QC RSD.",
        call. = FALSE
      )
    } else {
      refion_vec <- if (!is.na(refion_col)) reference_tbl[[refion_col]] else rep(NA_character_, nrow(reference_tbl))
      mz_vec <- if (!is.na(mz_col)) reference_tbl[[mz_col]] else rep(NA_real_, nrow(reference_tbl))

      metab_vec <- reference_tbl[[metab_col]]
      rt_vec <- reference_tbl[[rt_col]]

      ref_proc <- reference_tbl |>
        dplyr::mutate(
          .name_norm = normalize_text(metab_vec),
          .name_keys = lapply(metab_vec, build_name_keys),
          .refion_norm = normalize_text(refion_vec),
          .mz_round = round(parse_num_local(mz_vec), digits = mz_digits_local),
          .rt_round = round(parse_num_local(rt_vec), digits = rt_digits_local),
          .rt_raw = parse_num_local(rt_vec)
        ) |>
        dplyr::filter(!is.na(.name_norm), is.finite(.rt_round))
    }
  }

    refion_candidates <- c("Ref ion", "Ref Ion", "Reference ion", "Reference Ion", "ref_ion")
    refion_found <- refion_candidates[refion_candidates %in% names(ft_local)]
    if (length(refion_found) > 0) {
      ft_local$.refion_norm <- normalize_text(ft_local[[refion_found[1]]])
    } else {
      ft_local$.refion_norm <- NA_character_
    }

    group_keys <- names(named_groups)
    ref_option_count <- setNames(integer(length(group_keys)), group_keys)
    full_match_option_count <- setNames(integer(length(group_keys)), group_keys)
    selection_source <- setNames(rep("", length(group_keys)), group_keys)
    selected_feature_rt <- setNames(rep(NA_real_, length(group_keys)), group_keys)
    selected_reference_rt <- setNames(rep(NA_real_, length(group_keys)), group_keys)
    selected_rt_abs_diff <- setNames(rep(NA_real_, length(group_keys)), group_keys)

    get_nearest_reference_metrics <- function(candidate_row, refs_for_name) {
      if (is.null(refs_for_name) || nrow(refs_for_name) == 0) {
        return(c(ref_rt = NA_real_, abs_diff = Inf))
      }

      rt_candidate <- candidate_row$.rt_raw[[1]]
      if (!is.finite(rt_candidate)) {
        return(c(ref_rt = NA_real_, abs_diff = Inf))
      }

      refs_pool <- refs_for_name

      if (!is.na(candidate_row$.refion_norm[[1]]) && nzchar(candidate_row$.refion_norm[[1]])) {
        refs_by_refion <- refs_pool |>
          dplyr::filter(.refion_norm == candidate_row$.refion_norm[[1]])
        if (nrow(refs_by_refion) > 0) {
          refs_pool <- refs_by_refion
        }
      }

      if (is.finite(candidate_row$.mz_round[[1]])) {
        refs_by_mz <- refs_pool |>
          dplyr::filter(.mz_round == candidate_row$.mz_round[[1]])
        if (nrow(refs_by_mz) > 0) {
          refs_pool <- refs_by_mz
        }
      }

      if (!(".rt_raw" %in% names(refs_pool))) {
        return(c(ref_rt = NA_real_, abs_diff = Inf))
      }

      ref_rt <- refs_pool$.rt_raw
      if (!any(is.finite(ref_rt))) {
        return(c(ref_rt = NA_real_, abs_diff = Inf))
      }

      ref_rt_valid <- ref_rt[is.finite(ref_rt)]
      diffs <- abs(rt_candidate - ref_rt_valid)
      idx <- which.min(diffs)

      c(ref_rt = ref_rt_valid[idx], abs_diff = diffs[idx])
    }

    rt_distance_to_reference <- function(candidate_row, refs_for_name) {
      # Distance by subtraction: choose candidate with smallest |RT_candidate - RT_reference|.
      get_nearest_reference_metrics(candidate_row, refs_for_name)[["abs_diff"]]
    }

    for (k in names(named_groups)) {
      ids <- named_groups[[k]]
      name_keys_k <- build_name_keys(k)

      refs_k <- NULL
      if (!is.null(ref_proc) && nrow(ref_proc) > 0) {
        refs_k <- ref_proc |>
          dplyr::filter(vapply(.name_keys, function(keys) {
            length(intersect(keys, name_keys_k)) > 0
          }, logical(1)))
      }

      n_ref <- if (is.null(refs_k)) 0L else nrow(refs_k)
      ref_option_count[k] <- n_ref

      candidates <- ft_local |>
        dplyr::filter(featureID %in% ids)
      matched_ids <- character(0)

      if (n_ref > 0) {
        has_match <- vapply(seq_len(nrow(candidates)), function(i) {
          rr <- candidates[i, , drop = FALSE]
          if (is.na(rr$.refion_norm) || !is.finite(rr$.mz_round) || !is.finite(rr$.rt_round)) {
            return(FALSE)
          }

          refs_full <- refs_k |>
            dplyr::filter(!is.na(.refion_norm), is.finite(.mz_round), is.finite(.rt_round))
          if (nrow(refs_full) == 0) {
            return(FALSE)
          }

          any(
            refs_full$.refion_norm == rr$.refion_norm &
              refs_full$.mz_round == rr$.mz_round &
              refs_full$.rt_round == rr$.rt_round,
            na.rm = TRUE
          )
        }, logical(1))
        matched_ids <- candidates$featureID[has_match]
      }

      full_match_option_count[k] <- length(matched_ids)

      nearest_ids <- character(0)
      if (n_ref > 0) {
        dists <- vapply(seq_len(nrow(candidates)), function(i) {
          rt_distance_to_reference(candidates[i, , drop = FALSE], refs_k)
        }, numeric(1))

        if (any(is.finite(dists))) {
          min_dist <- min(dists[is.finite(dists)], na.rm = TRUE)
          nearest_ids <- candidates$featureID[is.finite(dists) & dists == min_dist]
        }
      }

      chosen <- NA_character_
      # Rule priority:
      # 1) If a reference exists for this metabolite group, choose the candidate with
      #    the smallest absolute RT difference to the reference (tie -> best QC RSD).
      # 2) If RT cannot be evaluated but full-match candidates exist, use them.
      # 3) Otherwise, fall back to best QC RSD.
      if (length(ids) == 1) {
        chosen <- ids[1]
        selection_source[k] <- "kept_single_named"
      } else if (length(nearest_ids) > 0) {
        chosen <- choose_best_qc(nearest_ids, qc_rsd)
        selection_source[k] <- "reference_nearest_rt"
      } else if (length(matched_ids) > 0) {
        chosen <- choose_best_qc(matched_ids, qc_rsd)
        selection_source[k] <- if (length(matched_ids) == 1) {
          "reference_full_match"
        } else {
          "reference_tie_best_qc_rsd"
        }
      } else {
        chosen <- choose_best_qc(ids, qc_rsd)
        selection_source[k] <- "best_qc_rsd_fallback"
      }

      chosen_row <- candidates |>
        dplyr::filter(featureID == chosen) |>
        dplyr::slice(1)
      if (nrow(chosen_row) == 1) {
        if (is.finite(chosen_row$.rt_raw[[1]])) {
          selected_feature_rt[k] <- chosen_row$.rt_raw[[1]]
        }

        nearest_metrics <- get_nearest_reference_metrics(chosen_row, refs_k)
        if (is.finite(nearest_metrics[["ref_rt"]])) {
          selected_reference_rt[k] <- nearest_metrics[["ref_rt"]]
        }
        if (is.finite(nearest_metrics[["abs_diff"]])) {
          selected_rt_abs_diff[k] <- nearest_metrics[["abs_diff"]]
        }
      }

      keep_named_ids <- c(keep_named_ids, chosen)
    }

    group_stats <- tibble::tibble(
      dup_key = names(ref_option_count),
      reference_has_metabolite = as.logical(ref_option_count > 0),
      reference_option_count = as.integer(ref_option_count),
      full_match_option_count = as.integer(full_match_option_count),
      selection_source = as.character(selection_source),
      selected_feature_rt = as.numeric(selected_feature_rt),
      selected_reference_rt = as.numeric(selected_reference_rt),
      selected_rt_abs_diff = as.numeric(selected_rt_abs_diff)
    )

    audit <- audit |>
      dplyr::left_join(group_stats, by = "dup_key")

    keep_ids <- c(keep_named_ids, unnamed_ids)
    mat2 <- mat[, keep_ids, drop = FALSE]

    ft2 <- feature_tbl |>
      dplyr::filter(featureID %in% keep_ids) |>
      dplyr::slice(match(keep_ids, featureID)) |>
      dplyr::mutate(
        display_name = dplyr::if_else(
          !is.na(Name_canon) & Name_canon != "",
          Name_canon,
          display_name
        ),
        display_name = make.unique(display_name, sep = "__dup")
      ) |>
      dplyr::select(dplyr::all_of(orig_feature_cols))

    if (!is.null(audit_path)) {
      write_csv_safe(audit, audit_path)

      # Export unresolved metabolite names (not found in reference by name keys).
      unmatched_idx <- with(
        audit,
        is_named & !is.na(reference_has_metabolite) & !reference_has_metabolite
      )

      unmatched_audit <- audit[unmatched_idx, , drop = FALSE]
      if (nrow(unmatched_audit) > 0) {
        keep_cols <- intersect(
          c(
            "featureID", "Name", "Name_canon", "display_name", "Ref ion",
            "mz", "RT", "dup_key", "n_in_group", "reference_has_metabolite",
            "reference_option_count", "full_match_option_count", "selection_source"
          ),
          names(unmatched_audit)
        )

        unmatched_path <- file.path(
          dirname(audit_path),
          "unrecognized_metabolite_name_audit.csv"
        )
        write_csv_safe(unmatched_audit[, keep_cols, drop = FALSE], unmatched_path)
      }
    }

    return(list(mat = mat2, feature = ft2, audit = audit))
  }
  
  
  # ---------------------------------------------------------------------------
  # Strategy: choose the best representative by QC RSD
  # ---------------------------------------------------------------------------
  if (strategy == "collapse_best_qc_rsd") {
    if (is.null(qc_rsd)) {
      stop("collapse_best_qc_rsd requires qc_rsd (named by featureID).")
    }

    if (is.null(names(qc_rsd))) {
      stop("collapse_best_qc_rsd requires qc_rsd to be a named vector by featureID.")
    }

    named_ids <- unique(audit$featureID[audit$is_named])
    missing_qc <- setdiff(named_ids, names(qc_rsd))

    if (length(named_ids) > 0 && length(missing_qc) == length(named_ids)) {
      stop("collapse_best_qc_rsd found no qc_rsd values for named features.")
    }

    if (length(missing_qc) > 0) {
      warning(
        paste0(
          "collapse_best_qc_rsd: ",
          length(missing_qc),
          " named featureIDs are missing qc_rsd values. ",
          "Missing/non-finite values will be treated as Inf, ",
          "and ties will be resolved by featureID order."
        ),
        call. = FALSE
      )
    }

    for (k in names(named_groups)) {
      ids <- named_groups[[k]]

      if (length(ids) == 1) {
        keep_named_ids <- c(keep_named_ids, ids)
      } else {
        rsd <- unname(qc_rsd[ids])
        rsd_sort <- ifelse(is.finite(rsd), rsd, Inf)
        chosen <- ids[order(rsd_sort, ids)][1]
        keep_named_ids <- c(keep_named_ids, chosen)
      }
    }

    keep_ids <- c(keep_named_ids, unnamed_ids)
    mat2 <- mat[, keep_ids, drop = FALSE]

    ft2 <- feature_tbl |>
      dplyr::filter(featureID %in% keep_ids) |>
      dplyr::slice(match(keep_ids, featureID)) |>
      dplyr::mutate(
        display_name = dplyr::if_else(
          !is.na(Name_canon) & Name_canon != "",
          Name_canon,
          display_name
        ),
        display_name = make.unique(display_name, sep = "_dup")
      ) |>
      dplyr::select(dplyr::all_of(orig_feature_cols))

    if (!is.null(audit_path)) write_csv_safe(audit, audit_path)

    return(list(mat = mat2, feature = ft2, audit = audit))
  }

  # ---------------------------------------------------------------------------
  # Strategies: collapse_mean / collapse_sum
  # ---------------------------------------------------------------------------
  for (k in names(named_groups)) {
    ids <- named_groups[[k]]

    if (length(ids) == 1) {
      keep_named_ids <- c(keep_named_ids, ids)
    } else {
      sub <- mat[, ids, drop = FALSE]
      agg <- if (strategy == "collapse_sum") {
        rowSums(sub, na.rm = TRUE)
      } else {
        rowMeans(sub, na.rm = TRUE)
      }

      # Preserve NA when all duplicate values are NA for a sample row
      all_na_rows <- apply(sub, 1, function(x) all(is.na(x)))
      agg[all_na_rows] <- NA

      base_id <- paste0("COLLAPSED_", make.names(k))
      new_id <- make_unique_id(base_id, used_ids)
      used_ids <- c(used_ids, new_id)
      collapsed_vecs[[new_id]] <- agg

      base_row <- feature_tbl |>
        dplyr::filter(featureID %in% ids) |>
        dplyr::slice(1)

      base_row$featureID <- new_id
      if ("Name" %in% names(base_row)) {
        base_row$Name <- k
      }
      base_row$Name_canon <- k
      base_row$display_name <- k

      collapsed_feats[[new_id]] <- base_row
    }
  }

  mat_keep <- mat[, c(keep_named_ids, unnamed_ids), drop = FALSE]

  if (length(collapsed_vecs) > 0) {
    mat_coll <- do.call(cbind, collapsed_vecs)
    mat_out <- cbind(mat_keep, mat_coll)

    ft_keep <- feature_tbl |>
      dplyr::filter(featureID %in% colnames(mat_keep)) |>
      dplyr::slice(match(colnames(mat_keep), featureID))

    ft_coll <- dplyr::bind_rows(collapsed_feats)
    ft_out <- dplyr::bind_rows(ft_keep, ft_coll)
  } else {
    mat_out <- mat_keep
    ft_out <- feature_tbl |>
      dplyr::filter(featureID %in% colnames(mat_out)) |>
      dplyr::slice(match(colnames(mat_out), featureID))
  }

  ft_out <- ft_out |>
    dplyr::filter(featureID %in% colnames(mat_out)) |>
    dplyr::slice(match(colnames(mat_out), featureID)) |>
    dplyr::mutate(
      display_name = dplyr::if_else(
        !is.na(Name_canon) & Name_canon != "",
        Name_canon,
        display_name
      ),
      display_name = make.unique(display_name, sep = "__dup")
    ) |>
    dplyr::select(dplyr::all_of(orig_feature_cols))

  if (!is.null(audit_path)) write_csv_safe(audit, audit_path)

  list(mat = mat_out, feature = ft_out, audit = audit)
}