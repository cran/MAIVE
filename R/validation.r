#' Validate MAIVE Data Structure and Scientific Constraints
#'
#' Internal function that validates the input data frame for MAIVE analysis.
#' Checks data structure, numeric types, minimum sample size, and other
#' scientific requirements.
#'
#' @param dat Data frame with required columns: bs, sebs, Ns
#' @return Cleaned data frame (invisible) with empty rows removed
#' @keywords internal
#' @noRd
validate_maive_data <- function(dat,
                                estimate = NULL,
                                se = NULL,
                                n = NULL,
                                study_id = NULL) {
  # Data frame type check
  if (!is.data.frame(dat)) {
    cli::cli_abort("Input 'dat' must be a data frame.", call. = FALSE)
  }

  # Required columns (check before row count so we can clean properly)
  required_cols <- c("bs", "sebs", "Ns")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    cli::cli_abort(
      sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  # Remove empty rows
  original_nrow <- nrow(dat)
  dat <- dat[rowSums(is.na(dat)) != ncol(dat), ]
  if (nrow(dat) < original_nrow) {
    cli::cli_alert_info(
      sprintf(
        "Removed %d completely empty row%s",
        original_nrow - nrow(dat),
        if (original_nrow - nrow(dat) == 1) "" else "s"
      )
    )
  }

  # Numeric types (strict validation, no silent coercion) after removing empty rows
  for (col in required_cols) {
    if (!is.numeric(dat[[col]])) {
      cli::cli_abort(
        sprintf("Column '%s' must be numeric. Convert your data before calling MAIVE.", col),
        call. = FALSE
      )
    }
    if (any(is.na(dat[[col]]))) {
      cli::cli_abort(
        sprintf("Column '%s' must not contain missing values.", col),
        call. = FALSE
      )
    }
    if (any(!is.finite(dat[[col]]))) {
      cli::cli_abort(
        sprintf("Column '%s' must contain finite values only.", col),
        call. = FALSE
      )
    }
  }

  # Minimum rows after cleaning (scientific requirement for reliable estimation)
  if (nrow(dat) < 4) {
    if (nrow(dat) < original_nrow) {
      # Some rows were removed, use specific message
      cli::cli_abort(
        "After removing empty rows, insufficient data remains (< 4 observations).",
        call. = FALSE
      )
    } else {
      # No rows were removed, use general message
      cli::cli_abort(
        "Insufficient data: MAIVE requires at least 4 observations for reliable estimation.",
        call. = FALSE
      )
    }
  }

  # Positive standard errors (on non-missing values)
  invalid_se <- !is.na(dat$sebs) & dat$sebs <= 0
  if (any(invalid_se)) {
    n_invalid <- sum(invalid_se)
    cli::cli_abort(
      sprintf(
        "Standard errors ('sebs') must be positive. Found %d non-positive value%s.",
        n_invalid,
        if (n_invalid == 1) "" else "s"
      ),
      call. = FALSE
    )
  }

  # Positive sample sizes (on non-missing values)
  invalid_n <- !is.na(dat$Ns) & dat$Ns <= 0
  if (any(invalid_n)) {
    n_invalid <- sum(invalid_n)
    cli::cli_abort(
      sprintf(
        "Sample sizes ('Ns') must be positive. Found %d non-positive value%s.",
        n_invalid,
        if (n_invalid == 1) "" else "s"
      ),
      call. = FALSE
    )
  }

  # Study ID degrees of freedom check
  if ("study_id" %in% names(dat)) {
    n_studies <- length(unique(dat$study_id))
    min_required <- n_studies + 3
    if (nrow(dat) < min_required) {
      cli::cli_abort(
        sprintf(
          "Insufficient degrees of freedom: %d observations with %d unique studies requires at least %d rows.",
          nrow(dat),
          n_studies,
          min_required
        ),
        call. = FALSE
      )
    }
  }

  invisible(dat)
}

#' Resolve column mappings for MAIVE inputs
#'
#' Allows users to specify custom column names for estimate, SE, Ns, and study_id.
#' Falls back to defaults (bs, sebs, Ns, study_id positional) when mappings are
#' not provided. Performs presence, numeric, non-missing, and finite validation.
#'
#' @keywords internal
#' @noRd
resolve_maive_columns <- function(dat, estimate = NULL, se = NULL, n = NULL, study_id = NULL) { # nolint: object_name_linter.
  col_for <- function(arg, default_name) {
    if (is.null(arg) || is.na(arg) || identical(arg, "")) {
      return(default_name)
    }
    as.character(arg)
  }

  est_col <- col_for(estimate, "bs")
  se_col <- col_for(se, "sebs")
  n_col <- col_for(n, "Ns")

  required <- c(est_col, se_col, n_col)
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0) {
    cli::cli_abort(
      sprintf("Missing required columns: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }

  check_vector <- function(vec, name) {
    if (!is.numeric(vec)) {
      cli::cli_abort(sprintf("Column '%s' must be numeric.", name), call. = FALSE)
    }
    if (any(is.na(vec))) {
      cli::cli_abort(sprintf("Column '%s' must not contain missing values.", name), call. = FALSE)
    }
    if (any(!is.finite(vec))) {
      cli::cli_abort(sprintf("Column '%s' must contain finite values only.", name), call. = FALSE)
    }
  }

  bs <- dat[[est_col]]
  sebs <- dat[[se_col]]
  Ns <- dat[[n_col]]

  check_vector(bs, est_col)
  check_vector(sebs, se_col)
  check_vector(Ns, n_col)

  # Study ID handling: use mapped name if provided; otherwise fallback to positional 4th col if present
  if (!is.null(study_id) && !identical(study_id, "")) {
    study_col <- as.character(study_id)
    if (!study_col %in% names(dat)) {
      cli::cli_abort(sprintf("Column '%s' (study_id) is missing.", study_col), call. = FALSE)
    }
    studyid <- dat[[study_col]]
  } else if (ncol(dat) >= 4) {
    studyid <- dat[[4]]
  } else {
    studyid <- NULL
  }

  # Validate study_id if present
  if (!is.null(studyid)) {
    if (any(is.na(studyid))) {
      cli::cli_abort("Column 'study_id' must not contain missing values.", call. = FALSE)
    }
    if (!is.atomic(studyid)) {
      cli::cli_abort("Column 'study_id' must be a vector.", call. = FALSE)
    }
    if (length(studyid) != length(bs)) {
      cli::cli_abort("Column 'study_id' must have the same length as the estimates.", call. = FALSE)
    }
  }

  # Rebuild a clean data.frame with standardized names for downstream use
  cleaned <- data.frame(
    bs = bs,
    sebs = sebs,
    Ns = Ns,
    stringsAsFactors = FALSE
  )
  if (!is.null(studyid)) {
    cleaned$study_id <- studyid
  }

  list(dat = cleaned)
}

#' Validate MAIVE Parameter Values
#'
#' Internal function that validates parameter values are within expected ranges.
#'
#' @param method Method choice (1=PET, 2=PEESE, 3=PET-PEESE, 4=EK)
#' @param weight Weighting scheme (0=equal, 1=standard, 2=adjusted, 3=study)
#' @param instrument Variance instrumentation (0=disabled, 1=enabled)
#' @param studylevel Study-level effects (0-3)
#' @param SE Standard error treatment (0-3)
#' @param AR Anderson-Rubin confidence interval (0=no, 1=yes)
#' @return TRUE invisibly if all validations pass
#' @keywords internal
#' @noRd
validate_maive_parameters <- function(method, weight, instrument, studylevel, SE, AR) {
  # Method: 1=PET, 2=PEESE, 3=PET-PEESE, 4=EK
  if (!method %in% c(1, 2, 3, 4)) {
    cli::cli_abort(
      "Parameter 'method' must be 1 (PET), 2 (PEESE), 3 (PET-PEESE), or 4 (EK).",
      call. = FALSE
    )
  }

  # Weight: 0=equal, 1=standard, 2=adjusted, 3=study
  if (!weight %in% c(0, 1, 2, 3)) {
    cli::cli_abort(
      "Parameter 'weight' must be 0 (equal), 1 (standard), 2 (adjusted), or 3 (study).",
      call. = FALSE
    )
  }

  # Binary/categorical parameters
  if (!instrument %in% c(0, 1)) {
    cli::cli_abort("Parameter 'instrument' must be 0 or 1.", call. = FALSE)
  }

  if (!studylevel %in% c(0, 1, 2, 3)) {
    cli::cli_abort("Parameter 'studylevel' must be 0, 1, 2, or 3.", call. = FALSE)
  }

  if (!SE %in% c(0, 1, 2, 3)) {
    cli::cli_abort("Parameter 'SE' must be 0, 1, 2, or 3.", call. = FALSE)
  }

  if (!AR %in% c(0, 1)) {
    cli::cli_abort("Parameter 'AR' must be 0 or 1.", call. = FALSE)
  }

  invisible(TRUE)
}

#' Normalize MAIVE options and apply cross-parameter guardrails
#'
#' Coerces parameter inputs to integers, maps string first-stage names, applies
#' AR/instrument guardrails, and returns a normalized option list for downstream
#' pipeline use. Relies on existing data/parameter validators.
#'
#' @param dat Data frame with required columns: bs, sebs, Ns (or mapped), study_id optional
#' @param estimate Optional column name to use instead of 'bs'
#' @param se Optional column name to use instead of 'sebs'
#' @param n Optional column name to use instead of 'Ns'
#' @param study_id Optional column name for study identifiers (defaults to 4th col if absent)
#' @param method Method choice (1=PET, 2=PEESE, 3=PET-PEESE, 4=EK)
#' @param weight Weighting scheme (0=equal, 1=standard, 2=adjusted, 3=study)
#' @param instrument Variance instrumentation (0=disabled, 1=enabled)
#' @param studylevel Study-level effects (0-3)
#' @param SE Standard error treatment (0-3)
#' @param AR Anderson-Rubin confidence interval (0=no, 1=yes)
#' @param first_stage First-stage specification: 0/\"levels\" or 1/\"log\"
#' @return List of normalized options and derived metadata
#' @keywords internal
#' @noRd
normalize_maive_options <- function(dat,
                                    method,
                                    weight,
                                    instrument,
                                    studylevel,
                                    SE,
                                    AR,
                                    first_stage = 0L,
                                    estimate = NULL,
                                    se = NULL,
                                    n = NULL,
                                    study_id = NULL) {
  dat <- as.data.frame(dat)
  dat <- validate_maive_data(dat)

  scalar_int <- function(value, name) {
    if (length(value) != 1L || is.na(value)) {
      cli::cli_abort(sprintf("Parameter '%s' must be a single non-missing value.", name), call. = FALSE)
    }
    as.integer(value)
  }

  normalize_first_stage <- function(first_stage) {
    if (missing(first_stage) || is.null(first_stage)) {
      return(0L)
    }
    if (is.character(first_stage)) {
      match_idx <- match(tolower(first_stage), c("levels", "log"))
      if (is.na(match_idx)) {
        cli::cli_abort("Parameter 'first_stage' must be one of 'levels' or 'log'.", call. = FALSE)
      }
      return(as.integer(match_idx - 1L))
    }
    scalar_int(first_stage, "first_stage")
  }

  apply_guardrails <- function(method, weight, instrument, AR, dat) {
    # Disable AR when incompatible with method/weights/instrument choice
    if (method == 4L || weight == 1L || instrument == 0L) {
      AR <- 0L
    }

    # Disable instrumentation (and AR) if Ns has no variation
    if (instrument == 1L) {
      # Check minimum sample size for IV estimation
      min_n_iv <- 10L
      if (nrow(dat) < min_n_iv) {
        cli::cli_warn(
          "Sample size ({nrow(dat)}) is small for IV estimation. Results may be unreliable. Consider using instrument=0 for small samples.",
          call. = FALSE
        )
      }

      Ns <- dat$Ns
      unique_Ns <- unique(stats::na.omit(Ns))
      if (length(unique_Ns) < 2L) {
        cli::cli_warn("Ns has no variation; disabling variance instrumentation (instrument=0).")
        instrument <- 0L
        AR <- 0L
      }
    }
    list(instrument = instrument, AR = AR)
  }

  method <- scalar_int(method, "method")
  weight <- scalar_int(weight, "weight")
  instrument <- scalar_int(instrument, "instrument")
  studylevel <- scalar_int(studylevel, "studylevel")
  SE <- scalar_int(SE, "SE")
  first_stage <- normalize_first_stage(first_stage)

  guardrails <- apply_guardrails(method, weight, instrument, AR, dat)
  instrument <- guardrails$instrument
  AR <- guardrails$AR

  if (!first_stage %in% c(0L, 1L)) {
    cli::cli_abort("Parameter 'first_stage' must be 0 (levels) or 1 (log).", call. = FALSE)
  }

  type_map <- c("CR0", "CR1", "CR2")
  type_choice <- if (SE == 3L) "CR0" else type_map[SE + 1L]
  first_stage_type <- c("levels", "log")[first_stage + 1L]

  validate_maive_parameters(method, weight, instrument, studylevel, SE, AR)

  resolved <- resolve_maive_columns(dat, estimate = estimate, se = se, n = n, study_id = study_id)
  dat <- resolved$dat

  list(
    dat = dat,
    method = method,
    weight = weight,
    instrument = instrument,
    studylevel = studylevel,
    SE = SE,
    AR = AR,
    type_choice = type_choice,
    alpha_s = 0.05,
    first_stage = first_stage,
    first_stage_type = first_stage_type
  )
}

#' @keywords internal
maive_prepare_data <- function(dat, studylevel) {
  bs <- dat$bs
  sebs <- dat$sebs
  Ns <- dat$Ns

  M <- length(bs)

  cluster <- studylevel %/% 2L
  dummy <- studylevel %% 2L

  if (ncol(dat) >= 4) {
    studyid <- dat[[4]]
  } else {
    studyid <- seq_len(M)
    dummy <- 0L
    cluster <- 0L
  }

  D <- maive_center_dummy_matrix(studyid)
  g <- if (cluster == 0L) seq_len(M) else studyid
  dat$g <- g

  list(
    dat = dat,
    bs = bs,
    sebs = sebs,
    Ns = Ns,
    M = M,
    studyid = studyid,
    dummy = dummy,
    cluster = cluster,
    g = g,
    D = D
  )
}
