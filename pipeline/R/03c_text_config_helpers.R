# =============================================================================
# 03c_text_config_helpers.R
# Text and config helpers
# =============================================================================

# -----------------------------------------------------------------------------
# Text cleaning / sanitization helpers
# -----------------------------------------------------------------------------
clean_text <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_trim(x)
  
  x_low <- tolower(x)
  
  x[x_low %in% c("not named","notnamed","unnamed","no name","noname")] <- NA_character_
  x[x %in% c("", "NA","N/A","n/a","-","Unknown","unknown","No results","no results")] <- NA_character_
  
  x
}

strip_v_suffix_end <- function(x) {
  x <- as.character(x)
  stringr::str_replace(x, "_v\\d+$", "")
}

standardize_metabolite_name <- function(x,
                                        remove_greek_letters = TRUE,
                                        remove_isomer_descriptors = TRUE) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_squish(x)
  
  missing <- is.na(x) | tolower(x) %in% c("", "na", "n/a", "nan", "null")
  out <- x
  
  if (isTRUE(remove_greek_letters)) {
    greek_map <- c(
      "\u03b1" = " alpha ", "\u0391" = " alpha ",
      "\u03b2" = " beta ",  "\u0392" = " beta ",
      "\u03b3" = " gamma ", "\u0393" = " gamma ",
      "\u03b4" = " delta ", "\u0394" = " delta ",
      "\u03b5" = " epsilon ", "\u0395" = " epsilon ",
      "\u03b6" = " zeta ", "\u0396" = " zeta ",
      "\u03b7" = " eta ", "\u0397" = " eta ",
      "\u03b8" = " theta ", "\u0398" = " theta ",
      "\u03b9" = " iota ", "\u0399" = " iota ",
      "\u03ba" = " kappa ", "\u039a" = " kappa ",
      "\u03bb" = " lambda ", "\u039b" = " lambda ",
      "\u03bc" = " mu ", "\u039c" = " mu ",
      "\u03bd" = " nu ", "\u039d" = " nu ",
      "\u03be" = " xi ", "\u039e" = " xi ",
      "\u03bf" = " omicron ", "\u039f" = " omicron ",
      "\u03c0" = " pi ", "\u03a0" = " pi ",
      "\u03c1" = " rho ", "\u03a1" = " rho ",
      "\u03c3" = " sigma ", "\u03a3" = " sigma ",
      "\u03c2" = " sigma ",
      "\u03c4" = " tau ", "\u03a4" = " tau ",
      "\u03c5" = " upsilon ", "\u03a5" = " upsilon ",
      "\u03c6" = " phi ", "\u03a6" = " phi ",
      "\u03c7" = " chi ", "\u03a7" = " chi ",
      "\u03c8" = " psi ", "\u03a8" = " psi ",
      "\u03c9" = " omega ", "\u03a9" = " omega "
    )
    out <- stringr::str_replace_all(out, greek_map)
    greek_words <- paste(
      c(
        "alpha", "alfa", "beta", "gamma", "gama", "delta", "epsilon",
        "zeta", "eta", "theta", "teta", "iota", "kappa", "lambda",
        "mu", "nu", "xi", "omicron", "pi", "rho", "sigma", "tau",
        "upsilon", "phi", "chi", "psi", "omega"
      ),
      collapse = "|"
    )
    out <- stringr::str_replace_all(
      out,
      stringr::regex(paste0("(^|[\\s,;_\\-/\\(\\[]+)(", greek_words, ")(?=([\\s,;_\\-/\\)\\]]+|$))"), ignore_case = TRUE),
      "\\1"
    )
  }
  
  if (isTRUE(remove_isomer_descriptors)) {
    out <- stringr::str_replace_all(out, stringr::regex("\\((?:\\+|-|\\+/-|\\u00b1|rac|r|s|e|z|cis|trans)\\)", ignore_case = TRUE), " ")
    out <- stringr::str_replace_all(out, stringr::regex("\\[(?:\\+|-|\\+/-|\\u00b1|rac|r|s|e|z|cis|trans)\\]", ignore_case = TRUE), " ")
    out <- stringr::str_replace_all(out, stringr::regex("(^|[\\s,;_\\-/]+)(?:d|l|dl|ld|cis|trans|rac|endo|exo)(?=([\\s,;_\\-/]+|$))", ignore_case = TRUE), "\\1")
    out <- stringr::str_replace_all(out, stringr::regex("\\b(?:isomer|isomeric)\\s*[A-Za-z0-9._-]*", ignore_case = TRUE), " ")
  }
  
  out <- stringr::str_replace_all(out, "[\\[\\]{}()]+", " ")
  out <- stringr::str_replace_all(out, "\\s*[;/]+\\s*", " ")
  out <- stringr::str_replace_all(out, "\\s*[-_]+\\s*", "-")
  out <- stringr::str_replace_all(out, "^-+|-+$", "")
  out <- stringr::str_squish(out)
  out[missing | out == ""] <- NA_character_
  out
}

sanitize_text_for_exports <- function(x, mode = c("greek_latin_ascii","ascii_translit")) {
  mode <- match.arg(mode)
  
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_squish(x)
  
  if (all(is.na(x))) return(x)
  
  if (mode == "greek_latin_ascii") {
    x <- stringi::stri_trans_general(x, "Greek-Latin; Latin-ASCII")
  } else {
    x <- iconv(x, from = "", to = "ASCII//TRANSLIT", sub = "")
  }
  
  x <- stringr::str_replace_all(x, "[[:cntrl:]]+", "")
  x
}

# Normalize simple names (lowercase trim)
normalize_name <- function(x) {
  tolower(trimws(as.character(x)))
}

# Normalize model group pairs: map per-model alias values to control/treatment labels
normalize_model_group_pairs <- function(groups_vec, model_vec, pair_map, control_label, treatment_label) {
  if (is.null(pair_map) || length(pair_map) == 0) {
    return(groups_vec)
  }
  
  if (is.null(model_vec) || length(model_vec) == 0) {
    return(groups_vec)
  }
  
  out <- as.character(groups_vec)
  model_vec <- trimws(as.character(model_vec))
  model_keys <- names(pair_map)
  
  if (is.null(model_keys) || length(model_keys) == 0) {
    return(groups_vec)
  }
  
  for (model_name in model_keys) {
    pair_raw <- as.character(pair_map[[model_name]])
    pair_vals <- unlist(strsplit(pair_raw, ",", fixed = TRUE), use.names = FALSE)
    pair_vals <- trimws(pair_vals)
    pair_vals <- pair_vals[nzchar(pair_vals)]
    
    if (length(pair_vals) == 0) {
      next
    }
    
    model_idx <- !is.na(model_vec) & trimws(model_vec) == trimws(model_name)
    if (!any(model_idx)) {
      next
    }
    
    group_vals <- trimws(out[model_idx])
    group_norm <- toupper(group_vals)
    if (length(pair_vals) >= 1) {
      control_idx <- which(model_idx)[group_norm %in% toupper(pair_vals[1])]
      if (length(control_idx) > 0) {
        out[control_idx] <- control_label
      }
    }
    
    if (length(pair_vals) >= 2) {
      treatment_idx <- which(model_idx)[group_norm %in% toupper(pair_vals[2])]
      if (length(treatment_idx) > 0) {
        out[treatment_idx] <- treatment_label
      }
    } else {
      treatment_idx <- which(model_idx)[group_norm %in% toupper(pair_vals[1])]
      if (length(treatment_idx) > 0) {
        out[treatment_idx] <- treatment_label
      }
    }
  }
  
  out
}

normalize_config_text <- function(text) {
  text <- gsub("\u201C|\u201D", '"', text, perl = TRUE)
  text <- gsub("\u2018|\u2019", "'", text, perl = TRUE)
  text
}

normalize_normalization_mode <- function(value, default = "qcrsc") {
  if (is.null(value) || length(value) == 0 || all(is.na(value))) {
    value <- default
  }
  if (is.null(value) || length(value) == 0 || all(is.na(value))) {
    return("")
  }

  mode_key <- gsub("[-_ ]", "", tolower(trimws(as.character(value)[1])))
  if (!nzchar(mode_key) && !is.null(default)) {
    mode_key <- gsub("[-_ ]", "", tolower(trimws(as.character(default)[1])))
  }

  if (identical(mode_key, "pqnqc")) {
    return("pqn_qc")
  }
  if (mode_key %in% c("pqnsample", "pqnnoqc", "pqnnqc")) {
    return("pqn_sample")
  }
  if (identical(mode_key, "cyclicloess")) {
    return("cyclic_loess")
  }
  if (identical(mode_key, "qcloess")) {
    return("qc_loess")
  }
  if (identical(mode_key, "qcrsc")) {
    return("qcrsc")
  }

  valid_modes <- c("none", "weight", "qc_loess", "cyclic_loess", "qcrsc", "pqn_qc", "pqn_sample")
  if (mode_key %in% valid_modes) {
    return(mode_key)
  }

  if (!is.null(default)) {
    return(normalize_normalization_mode(default, default = NULL))
  }

  mode_key
}

