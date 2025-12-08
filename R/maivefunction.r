#' @keywords internal
maive_validate_inputs <- function(dat, method, weight, instrument, studylevel, SE, AR, first_stage) {
  dat <- as.data.frame(dat)
  if (ncol(dat) < 3) {
    stop("dat must contain at least three columns: bs, sebs, and Ns.")
  }

  scalar_int <- function(value, name) {
    if (length(value) != 1L || is.na(value)) {
      stop(sprintf("%s must be a single non-missing value.", name))
    }
    as.integer(value)
  }

  method <- scalar_int(method, "method")
  weight <- scalar_int(weight, "weight")
  instrument <- scalar_int(instrument, "instrument")
  studylevel <- scalar_int(studylevel, "studylevel")
  SE <- scalar_int(SE, "SE")
  AR <- scalar_int(AR, "AR")
  if (missing(first_stage)) {
    first_stage <- 0L
  }

  if (is.character(first_stage)) {
    match_idx <- match(tolower(first_stage), c("levels", "log"))
    if (is.na(match_idx)) {
      stop("first_stage must be one of 'levels' or 'log'.")
    }
    first_stage <- match_idx - 1L
  }
  first_stage <- scalar_int(first_stage, "first_stage")

  if (!method %in% 1:4) stop("method must be between 1 and 4.")
  if (!weight %in% 0:3) stop("weight must be 0, 1, 2, or 3.")
  if (!instrument %in% 0:1) stop("instrument must be 0 or 1.")
  if (!studylevel %in% 0:3) stop("studylevel must be between 0 and 3.")
  if (!SE %in% 0:3) stop("SE must be between 0 and 3.")
  if (!AR %in% 0:1) stop("AR must be 0 or 1.")
  if (!first_stage %in% 0:1) stop("first_stage must be 0 (levels) or 1 (log).")

  type_map <- c("CR0", "CR1", "CR2")
  type_choice <- if (SE == 3L) "CR0" else type_map[SE + 1L]
  first_stage_type <- c("levels", "log")[first_stage + 1L]

  if (method == 4L || weight == 1L || instrument == 0L) {
    AR <- 0L
  }

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
maive_build_dummy_matrix <- function(values) {
  f <- factor(values)
  mm <- stats::model.matrix(~ f - 1)
  colnames(mm) <- paste0("studyid.", levels(f))
  rownames(mm) <- NULL
  attr(mm, "assign") <- NULL
  attr(mm, "contrasts") <- NULL
  mm
}

#' @keywords internal
maive_center_dummy_matrix <- function(values) {
  D <- maive_build_dummy_matrix(values)
  if (is.null(dim(D)) || ncol(D) == 0L) {
    return(matrix(0, nrow = length(values), ncol = 0))
  }
  centered <- sweep(D, 2, colMeans(D), "-")
  if (ncol(centered) <= 1L) {
    centered[, 0, drop = FALSE]
  } else {
    centered[, seq_len(ncol(centered) - 1L), drop = FALSE]
  }
}

#' @keywords internal
maive_prepare_data <- function(dat, studylevel) {
  bs <- dat[[1]]
  sebs <- dat[[2]]
  Ns <- dat[[3]]

  if (!is.numeric(bs) || !is.numeric(sebs) || !is.numeric(Ns)) {
    stop("bs, sebs, and Ns must be numeric.")
  }
  if (any(!is.finite(bs)) || any(!is.finite(sebs)) || any(!is.finite(Ns))) {
    stop("bs, sebs, and Ns must be finite.")
  }
  if (any(sebs <= 0)) {
    stop("sebs must be strictly positive.")
  }
  if (any(Ns <= 0)) {
    stop("Ns must be strictly positive.")
  }

  M <- length(bs)
  if (length(sebs) != M || length(Ns) != M) {
    stop("bs, sebs, and Ns must have the same length.")
  }

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

#' @keywords internal
maive_compute_variance_instrumentation <- function(sebs, Ns, g, type_choice, instrument, first_stage_type) {
  invNs <- 1 / Ns
  sebs2 <- sebs^2
  if (first_stage_type == "log") {
    log_sebs2 <- log(sebs2)
    log_Ns <- log(Ns)
    Xiv <- cbind(1, log_Ns)
    varreg1 <- lm(log_sebs2 ~ 0 + Xiv)
    resid_varreg1 <- stats::residuals(varreg1)
    weights_varreg1 <- stats::weights(varreg1)
    if (is.null(weights_varreg1)) {
      weights_varreg1 <- rep(1, length(resid_varreg1))
    }
    smearing <- sum(weights_varreg1 * exp(resid_varreg1)) / sum(weights_varreg1)
    sebs2fit1 <- exp(stats::fitted(varreg1)) * smearing
    slope_index <- 2L
    # Use log(Ns) as instrument for AR test (matches first-stage)
    instrument_for_ar <- log_Ns
  } else {
    Xiv <- cbind(1, invNs)
    varreg1 <- lm(sebs2 ~ 0 + Xiv)
    dimiv <- 2L
    if (varreg1$coefficients[1] < 0) {
      Xiv <- as.matrix(invNs)
      varreg1 <- lm(sebs2 ~ 0 + Xiv)
      dimiv <- 1L
    }

    sebs2fit1 <- stats::fitted(varreg1)
    slope_index <- dimiv
    # Use 1/Ns as instrument for AR test (matches first-stage)
    instrument_for_ar <- invNs
  }
  if (instrument == 0L) {
    F_hac <- "NA"
  } else {
    V <- clubSandwich::vcovCR(varreg1, cluster = g, type = type_choice)
    F_hac <- unname(round(varreg1$coefficients[slope_index]^2 / V[slope_index, slope_index], 3))
  }

  list(
    invNs = invNs,
    instrument_for_ar = instrument_for_ar,
    sebs2fit1 = sebs2fit1,
    F_hac = F_hac,
    first_stage_model = varreg1
  )
}

#' @keywords internal
maive_compute_weights <- function(weight, sebs, sebs2fit1, studyid = NULL) {
  if (weight == 0L) {
    rep(1, length(sebs))
  } else if (weight == 1L) {
    sebs
  } else if (weight == 2L) {
    sqrt(sebs2fit1)
  } else if (weight == 3L) {
    if (is.null(studyid)) {
      studyid <- seq_along(sebs)
    }
    if (length(studyid) != length(sebs)) {
      stop("studyid must align with sebs when using study weights.")
    }
    counts <- ave(rep(1, length(studyid)), studyid, FUN = length)
    1 / counts
  } else {
    stop("Invalid weight option.")
  }
}

#' Compute exponential-decay weights from first-stage residuals
#'
#' @param first_stage_model Fitted lm object from first-stage regression
#' @return Exponential-decay weights normalized to mean 1
#' @keywords internal
#' @noRd
maive_compute_waive_weights <- function(first_stage_model) {
  if (is.null(first_stage_model)) {
    stop("first_stage_model must be supplied for weighting.")
  }

  nu <- stats::residuals(first_stage_model)
  if (length(nu) == 0L) {
    stop("first_stage_model must provide residuals.")
  }

  sigma <- 1.4826 * stats::mad(nu, constant = 1, na.rm = TRUE)

  if (!is.finite(sigma) || sigma <= 0) {
    sigma <- stats::sd(nu, na.rm = TRUE)
    if (!is.finite(sigma) || sigma <= 0) {
      sigma <- 1e-12
    } else {
      sigma <- sigma + 1e-12
    }
  }

  z <- nu / sigma
  z_neg <- pmax(-z, 0)
  z_out <- pmax(abs(z) - 2, 0)

  w <- exp(-1.0 * z_neg - 0.25 * z_out^2)
  w <- pmax(w, 0.05)

  mean_w <- mean(w)
  if (!is.finite(mean_w) || mean_w <= 0) {
    stop("Failed to compute valid weights.")
  }

  w / mean_w
}

#' Run the shared analysis pipeline
#'
#' @param opts Validated options
#' @param prepared Prepared data
#' @param instrumentation First-stage results
#' @param w Weights for second-stage regression
#' @return List of analysis results
#' @keywords internal
#' @noRd
maive_run_pipeline <- function(opts, prepared, instrumentation, w) {
  if (!is.numeric(w) || length(w) != prepared$M) {
    stop("w must be a numeric vector aligned with the input data.")
  }

  x <- if (opts$instrument == 0L) prepared$sebs else sqrt(instrumentation$sebs2fit1)
  x2 <- if (opts$instrument == 0L) prepared$sebs^2 else instrumentation$sebs2fit1

  design <- maive_build_design_matrices(prepared$bs, prepared$sebs, w, x, x2, prepared$D, prepared$dummy)
  fits <- maive_fit_models(design)
  selection <- maive_select_petpeese(fits, design, opts$alpha_s, opts$SE, prepared$dat, opts$type_choice)
  sighats <- maive_compute_sigma_h(fits, design$w, design$sebs)
  ek <- maive_fit_ek(selection, design, sighats, opts$method)

  slope_info <- maive_slope_information(opts$method, fits, selection, ek)
  slope_summary <- maive_quadratic_summary(opts$method, selection, slope_info)

  # Determine PET-PEESE selection (method == 3 only)
  petpeese_selected <- if (opts$method == 3L) {
    if (isTRUE(selection$quadratic_decision)) "PEESE" else "PET"
  } else {
    NA_character_
  }

  # Compute PEESE SE^2 coefficient and standard error when PEESE is final model
  peese_is_final <- (opts$method == 2L) || (opts$method == 3L && isTRUE(selection$quadratic_decision))
  if (peese_is_final) {
    peese_se2_inf <- maive_infer_coef(fits$peese, 2L, opts$SE, prepared$dat, "g", opts$type_choice)
    peese_se2_coef <- round(peese_se2_inf$b, 3)
    peese_se2_se <- round(peese_se2_inf$se, 3)
  } else {
    peese_se2_coef <- NA_real_
    peese_se2_se <- NA_real_
  }

  egger_inf <- maive_infer_coef(fits$fatpet, 2L, opts$SE, prepared$dat, "g", opts$type_choice)
  egger_boot_ci <- round(egger_inf$ci, 3)
  egger_ar_ci <- maive_compute_egger_ar_ci(
    opts,
    fits,
    prepared,
    instrumentation$instrument_for_ar,
    adjusted_variance = instrumentation$sebs2fit1,
    f_stat = instrumentation$F_hac
  )
  cfg <- maive_get_config(opts$method, fits, selection, ek)
  hausman_cfg <- maive_get_hausman_models(opts$method, cfg, selection, design)
  if (is.null(cfg$maive) || is.null(cfg$std)) {
    stop("Failed to identify models for the selected method.")
  }

  beta <- unname(coef(hausman_cfg$maive)[1])
  beta0 <- unname(coef(hausman_cfg$std)[1])

  se_ma <- maive_get_intercept_se(cfg$maive, opts$SE, prepared$dat, "g", opts$type_choice)
  se_std <- maive_get_intercept_se(cfg$std, opts$SE, prepared$dat, "g", opts$type_choice)

  hausman <- maive_compute_hausman(beta, beta0, hausman_cfg$maive, hausman_cfg$std, prepared$g, opts$type_choice)
  chi2 <- qchisq(p = 0.05, df = 1, lower.tail = FALSE)

  ar_ci_res <- maive_compute_ar_ci(
    opts,
    fits,
    selection,
    prepared,
    instrumentation$instrument_for_ar,
    opts$type_choice,
    adjusted_variance = instrumentation$sebs2fit1
  )

  list(
    "beta" = round(beta, 3),
    "SE" = round(as.numeric(se_ma$se), 3),
    "F-test" = instrumentation$F_hac,
    "beta_standard" = round(beta0, 3),
    "SE_standard" = round(as.numeric(se_std$se), 3),
    "Hausman" = round(hausman, 3),
    "Chi2" = round(chi2, 3),
    "SE_instrumented" = sqrt(instrumentation$sebs2fit1),
    "AR_CI" = ar_ci_res$b0_CI,
    "pub bias p-value" = round(egger_inf$p, 3),
    "egger_coef" = round(egger_inf$b, 3),
    "egger_se" = round(egger_inf$se, 3),
    "egger_boot_ci" = egger_boot_ci,
    "egger_ar_ci" = egger_ar_ci,
    "is_quadratic_fit" = slope_summary,
    "boot_result" = se_ma$boot_result,
    "slope_coef" = slope_info$coefficient,
    "petpeese_selected" = petpeese_selected,
    "peese_se2_coef" = peese_se2_coef,
    "peese_se2_se" = peese_se2_se,
    "weights" = w
  )
}

#' @keywords internal
maive_build_design_matrices <- function(bs, sebs, w, x, x2, D, dummy) {
  y <- bs / w
  X <- cbind(1, x) / w
  X2 <- cbind(1, x2) / w

  y0 <- bs / sebs
  X0 <- cbind(1, sebs) / sebs
  X20 <- cbind(1, sebs^2) / sebs

  if (dummy == 1L) {
    X <- cbind(X, D / w)
    X2 <- cbind(X2, D / w)
    cD <- cbind(1, D) / w
    X0 <- cbind(X0, D / sebs)
    X20 <- cbind(X20, D / sebs)
    cD0 <- cbind(1, D) / sebs
  } else {
    cD <- matrix(1 / w, ncol = 1)
    cD0 <- matrix(1 / sebs, ncol = 1)
  }

  list(
    y = y,
    X = X,
    X2 = X2,
    y0 = y0,
    X0 = X0,
    X20 = X20,
    cD = cD,
    cD0 = cD0,
    w = w,
    sebs = sebs,
    x = x,
    D = D,
    dummy = dummy
  )
}

#' @keywords internal
maive_build_auxiliary_petpeese_matrix <- function(term, design) {
  base <- cbind(1, term) / design$w
  if (design$dummy == 1L && !is.null(design$D) && ncol(design$D) > 0L) {
    base <- cbind(base, design$D / design$w)
  }
  base
}

#' @keywords internal
maive_fit_auxiliary_petpeese <- function(selection, design) {
  if (isTRUE(selection$quadratic_decision)) {
    X_aux <- maive_build_auxiliary_petpeese_matrix(design$sebs^2, design)
  } else {
    X_aux <- maive_build_auxiliary_petpeese_matrix(design$sebs, design)
  }
  lm(design$y ~ 0 + X_aux)
}

#' @keywords internal
maive_fit_models <- function(design) {
  y <- design$y
  X <- design$X
  X2 <- design$X2
  y0 <- design$y0
  X0 <- design$X0
  X20 <- design$X20
  cD <- design$cD
  cD0 <- design$cD0

  list(
    fatpet = lm(y ~ 0 + X),
    peese = lm(y ~ 0 + X2),
    fatpet0 = lm(y0 ~ 0 + X0),
    peese0 = lm(y0 ~ 0 + X20),
    wlsreg = lm(y ~ 0 + cD),
    wlsreg0 = lm(y0 ~ 0 + cD0)
  )
}

#' @keywords internal
maive_select_petpeese <- function(fits, design, alpha_s, SE = NULL, data = NULL, type_choice = NULL) {
  M <- length(design$y)

  # Use proper SE method for PET/PEESE selection to respect clustering/bootstrapping
  if (!is.null(SE) && !is.null(data) && !is.null(type_choice)) {
    # Get robust SE for the intercept using the user's chosen method
    intercept_inf <- maive_infer_coef(fits$fatpet, 1L, SE, data, "g", type_choice)
    quad_stat <- abs(intercept_inf$b / intercept_inf$se)
  } else {
    # Fallback to unclustered variance (backward compatibility)
    quad_stat <- abs(coef(fits$fatpet)[1] / sqrt(vcov(fits$fatpet)[1, 1]))
  }
  quad_cutoff <- qt(1 - alpha_s / 2, M - ncol(design$X))
  quadratic_decision <- quad_stat > quad_cutoff
  petpeese <- if (quadratic_decision) fits$peese else fits$fatpet

  # For the standard model (fatpet0), also use proper SE if available
  if (!is.null(SE) && !is.null(data) && !is.null(type_choice)) {
    intercept_inf0 <- maive_infer_coef(fits$fatpet0, 1L, SE, data, "g", type_choice)
    quad_stat0 <- abs(intercept_inf0$b / intercept_inf0$se)
  } else {
    quad_stat0 <- abs(coef(fits$fatpet0)[1] / sqrt(vcov(fits$fatpet0)[1, 1]))
  }
  quad_cutoff0 <- qt(1 - alpha_s / 2, M - ncol(design$X0))
  quadratic_decision0 <- quad_stat0 > quad_cutoff0
  petpeese0 <- if (quadratic_decision0) fits$peese0 else fits$fatpet0

  list(
    quadratic_decision = quadratic_decision,
    petpeese = petpeese,
    quadratic_decision0 = quadratic_decision0,
    petpeese0 = petpeese0
  )
}

#' @keywords internal
maive_compute_sigma_h <- function(fits, w, sebs) {
  M <- length(w)
  wis0 <- 1 / (w^2)
  Qfe0 <- sum(residuals(fits$wlsreg)^2)
  denom0 <- M - ncol(model.matrix(fits$wlsreg)) - 1
  sigh2hat0 <- max(0, M * ((Qfe0 / denom0) - 1) / sum(wis0))
  sighhat0 <- sqrt(sigh2hat0)

  wis00 <- 1 / (sebs^2)
  Qfe00 <- sum(residuals(fits$wlsreg0)^2)
  denom00 <- M - ncol(model.matrix(fits$wlsreg0)) - 1
  sigh2hat00 <- max(0, M * ((Qfe00 / denom00) - 1) / sum(wis00))
  sighhat00 <- sqrt(sigh2hat00)

  list(sighhat0 = sighhat0, sighhat00 = sighhat00)
}

#' @keywords internal
maive_fit_ek <- function(selection, design, sighats, method) {
  if (method != 4L) {
    return(list(ekreg = NULL, ekreg0 = NULL, structure = "linear", a0 = NA_real_, a00 = NA_real_))
  }

  intercept <- coef(selection$petpeese)[1]
  threshold <- 1.96 * sighats$sighhat0
  if (intercept > threshold) {
    a0 <- (intercept - threshold) * (intercept + threshold) / (2 * 1.96 * intercept)
  } else {
    a0 <- 0
  }

  intercept0 <- coef(selection$petpeese0)[1]
  threshold0 <- 1.96 * sighats$sighhat00
  if (intercept0 > threshold0) {
    a00 <- (intercept0 - threshold0) * (intercept0 + threshold0) / (2 * 1.96 * intercept0)
  } else {
    a00 <- 0
  }

  y <- design$y
  cD <- design$cD
  y0 <- design$y0
  cD0 <- design$cD0

  if (!is.na(a0) && a0 > min(design$x) && a0 < max(design$x)) {
    xx_w <- (design$x - a0) * (design$x > a0) / design$w
    ekreg <- lm(y ~ 0 + cD + xx_w)
    ek_structure <- "kink"
  } else if (!is.na(a0) && a0 < min(design$x)) {
    x_w <- design$x / design$w
    ekreg <- lm(y ~ 0 + cD + x_w)
    ek_structure <- "linear"
  } else {
    ekreg <- lm(y ~ 0 + cD)
    ek_structure <- "intercept"
  }

  if (a00 > min(design$sebs) && a00 < max(design$sebs)) {
    xx0_w <- (design$sebs - a00) * (design$sebs > a00) / design$sebs
    ekreg0 <- lm(y0 ~ 0 + cD0 + xx0_w)
  } else if (a00 < min(design$sebs)) {
    x0_w <- design$sebs / design$sebs
    ekreg0 <- lm(y0 ~ 0 + cD0 + x0_w)
  } else {
    ekreg0 <- lm(y0 ~ 0 + cD0)
  }

  list(ekreg = ekreg, ekreg0 = ekreg0, structure = ek_structure, a0 = a0, a00 = a00)
}

#' @keywords internal
maive_infer_coef <- function(model, coef_index, SE, data, cluster_var, type_choice, alpha = 0.05) {
  if (SE == 3L) {
    boot <- manual_wild_cluster_boot_se(
      model = model,
      data = data,
      cluster_var = cluster_var,
      B = 999
    )
    if (is.null(boot$boot_se)) {
      stop("Bootstrap helper must return boot_se.")
    }
    se <- unname(boot$boot_se[coef_index])
    boot_result <- boot
  } else {
    V <- clubSandwich::vcovCR(model, cluster = data[[cluster_var]], type = type_choice)
    se <- unname(sqrt(V[coef_index, coef_index]))
    boot_result <- NULL
  }
  b <- unname(coef(model)[coef_index])
  t <- as.numeric(b / se)
  p <- 2 * pnorm(-abs(t))
  ci <- maive_prepare_confidence_interval(
    model = model,
    coef_index = coef_index,
    estimate = b,
    se = se,
    boot_result = boot_result,
    alpha = alpha
  )
  list(b = b, se = se, p = p, ci = ci, boot_result = boot_result)
}

#' @keywords internal
maive_get_intercept_se <- function(model, SE, data, cluster_var, type_choice) {
  inf <- maive_infer_coef(model, 1L, SE, data, cluster_var, type_choice)
  list(se = inf$se, ci = inf$ci, boot_result = inf$boot_result)
}

#' @keywords internal
maive_prepare_confidence_interval <- function(model, coef_index, estimate, se, boot_result, alpha) {
  if (!is.null(boot_result) && !is.null(boot_result$boot_ci)) {
    boot_ci <- boot_result$boot_ci
    coef_names <- names(coef(model))
    ci_row <- NULL
    if (!is.null(coef_names) && length(coef_names) >= coef_index) {
      coef_name <- coef_names[coef_index]
      if (!is.null(rownames(boot_ci)) && coef_name %in% rownames(boot_ci)) {
        ci_row <- boot_ci[coef_name, ]
      }
    }
    if (is.null(ci_row)) {
      ci_row <- boot_ci[coef_index, ]
    }
    ci_vals <- maive_normalize_ci_bounds(ci_row)
    ci_vals
  } else {
    crit <- stats::qnorm(1 - alpha / 2)
    ci_vals <- c(estimate - crit * se, estimate + crit * se)
    names(ci_vals) <- c("lower", "upper")
    ci_vals
  }
}

maive_normalize_ci_bounds <- function(ci_row) {
  ci_vals <- as.numeric(ci_row)
  if (length(ci_vals) == 0) {
    ci_vals <- rep(NA_real_, 2L)
  } else if (length(ci_vals) == 1) {
    ci_vals <- rep(ci_vals[1], 2L)
  } else if (length(ci_vals) > 2) {
    ci_vals <- ci_vals[seq_len(2L)]
  }
  if (length(ci_vals) < 2) {
    ci_vals <- c(ci_vals, rep(NA_real_, 2L - length(ci_vals)))
  }
  names(ci_vals) <- c("lower", "upper")
  ci_vals
}

#' @keywords internal
maive_compute_egger_ar_ci <- function(opts, fits, prepared, invNs, adjusted_variance = NULL, f_stat = NULL) {
  if (opts$AR != 1L || opts$weight == 1L || opts$instrument == 0L || prepared$dummy == 1L) {
    return("NA")
  }
  if (is.null(fits$fatpet)) {
    return("NA")
  }
  ar_weights <- NULL
  if (opts$weight == 2L) {
    if (is.null(adjusted_variance)) {
      stop("Adjusted variance estimates are required when computing weighted AR intervals.")
    }
    ar_weights <- 1 / adjusted_variance
  }

  # Determine which SE to use for AR test

  # When instrument=1, the fatpet model uses instrumented SE (sqrt(adjusted_variance))
  # The AR test must use the same SE for consistency
  if (opts$instrument == 1L && !is.null(adjusted_variance)) {
    sebs_for_ar <- sqrt(adjusted_variance)
  } else {
    sebs_for_ar <- prepared$sebs
  }

  # Always use subset AR for Egger slope CI

  # The joint method can produce spuriously narrow CIs due to banana-projection
  # even when F-stat is strong. Subset AR is more robust for slope inference.
  ar_method <- "slope_only"

  ar_result <- compute_AR_CI_optimized(
    model = fits$fatpet,
    adjust_fun = PET_adjust,
    bs = prepared$bs,
    sebs = sebs_for_ar,
    invNs = invNs,
    g = prepared$g,
    type_choice = opts$type_choice,
    weights = ar_weights,
    method = ar_method
  )
  if (is.null(ar_result$b1_CI) || identical(ar_result$b1_CI, "NA")) {
    ci_vals <- c(NA_real_, NA_real_)
    names(ci_vals) <- c("lower", "upper")
    return(ci_vals)
  }
  if (all(is.na(ar_result$b1_CI)) || any(is.infinite(ar_result$b1_CI))) {
    ci_vals <- c(NA_real_, NA_real_)
    names(ci_vals) <- c("lower", "upper")
    return(ci_vals)
  }
  ci_vals <- round(ar_result$b1_CI, 3)
  names(ci_vals) <- c("lower", "upper")
  ci_vals
}

#' @keywords internal
maive_slope_information <- function(method, fits, selection, ek) {
  method_str <- as.character(method)
  if (method_str == "1") {
    return(list(type = "linear", coefficient = round(as.numeric(coef(fits$fatpet)[2]), 3), detail = NULL))
  }
  if (method_str == "2") {
    return(list(type = "quadratic", coefficient = round(as.numeric(coef(fits$peese)[2]), 3), detail = NULL))
  }
  if (method_str == "3") {
    if (identical(selection$petpeese, fits$peese)) {
      return(list(type = "quadratic", coefficient = round(as.numeric(coef(fits$peese)[2]), 3), detail = NULL))
    }
    return(list(type = "linear", coefficient = round(as.numeric(coef(fits$fatpet)[2]), 3), detail = NULL))
  }
  if (method_str == "4") {
    if (is.null(ek$ekreg)) {
      return(list(type = "linear", coefficient = 0, detail = NULL))
    }
    if (ek$structure == "kink") {
      kink_effect <- round(as.numeric(tail(coef(ek$ekreg), 1)), 3)
      kink_location <- as.numeric(round(ek$a0, 3))
      detail <- list(kink_location = kink_location, kink_effect = kink_effect)
      return(list(
        type = "kinked",
        coefficient = list(kink_effect = kink_effect, kink_location = kink_location),
        detail = detail
      ))
    }
    if (ek$structure == "linear") {
      slope <- round(as.numeric(tail(coef(ek$ekreg), 1)), 3)
      return(list(type = "linear", coefficient = slope, detail = NULL))
    }
    return(list(type = "linear", coefficient = 0, detail = NULL))
  }
  stop("Invalid method")
}

#' @keywords internal
maive_quadratic_summary <- function(method, selection, slope_info) {
  quadratic_flag <- switch(as.character(method),
    "1" = FALSE,
    "2" = TRUE,
    "3" = isTRUE(selection$quadratic_decision),
    "4" = FALSE,
    stop("Invalid method")
  )
  list(
    quadratic = quadratic_flag,
    slope_type = slope_info$type,
    slope_detail = slope_info$detail
  )
}

#' @keywords internal
maive_get_config <- function(method, fits, selection, ek) {
  switch(as.character(method),
    "1" = list(maive = fits$fatpet, std = fits$fatpet0),
    "2" = list(maive = fits$peese, std = fits$peese0),
    "3" = list(maive = selection$petpeese, std = selection$petpeese0),
    "4" = list(maive = ek$ekreg, std = ek$ekreg0),
    stop("Invalid method")
  )
}

#' @keywords internal
maive_get_hausman_models <- function(method, cfg, selection, design) {
  if (as.character(method) != "3") {
    return(cfg)
  }

  aux_std <- maive_fit_auxiliary_petpeese(selection, design)
  list(maive = cfg$maive, std = aux_std)
}

#' @keywords internal
maive_compute_hausman <- function(beta_iv, beta_ols, model_iv, model_ols, g, type_choice) {
  V_iv <- clubSandwich::vcovCR(model_iv, cluster = g, type = type_choice)
  V_ols <- clubSandwich::vcovCR(model_ols, cluster = g, type = type_choice)
  var_diff <- V_iv[1, 1] - V_ols[1, 1]
  if (!is.finite(var_diff) || var_diff <= 0) {
    return(NA_real_)
  }
  unname((beta_iv - beta_ols)^2 / var_diff)
}

#' @keywords internal
maive_compute_ar_ci <- function(opts, fits, selection, prepared, invNs, type_choice, adjusted_variance = NULL) {
  if (opts$AR != 1L || opts$method == 4L || opts$weight == 1L || prepared$dummy == 1L) {
    return(list(b0_CI = "NA", b1_CI = "NA"))
  }

  cfg_ar <- switch(as.character(opts$method),
    "1" = list(model = fits$fatpet, adjust_fun = PET_adjust),
    "2" = list(model = fits$peese, adjust_fun = PEESE_adjust),
    "3" = if (identical(selection$petpeese, fits$peese)) {
      list(model = fits$peese, adjust_fun = PEESE_adjust)
    } else {
      list(model = fits$fatpet, adjust_fun = PET_adjust)
    },
    stop("Invalid method")
  )

  ar_weights <- NULL
  if (opts$weight == 2L) {
    if (is.null(adjusted_variance)) {
      stop("Adjusted variance estimates are required when computing weighted AR intervals.")
    }
    ar_weights <- 1 / adjusted_variance
  }

  # When instrument=1, use instrumented SE for consistency with fitted model
  if (opts$instrument == 1L && !is.null(adjusted_variance)) {
    sebs_for_ar <- sqrt(adjusted_variance)
  } else {
    sebs_for_ar <- prepared$sebs
  }

  do.call(
    compute_AR_CI_optimized,
    c(cfg_ar, list(
      bs = prepared$bs,
      sebs = sebs_for_ar,
      invNs = invNs,
      g = prepared$g,
      type_choice = type_choice,
      weights = ar_weights
    ))
  )
}

#' R code for MAIVE
#'
#' R package for MAIVE: "Spurious Precision in Meta-Analysis of Observational Research" by
#' Zuzana Irsova, Pedro Bom, Tomas Havranek, Petr Cala, and Heiko Rachinger.
#'
#' @param dat Data frame with columns bs, sebs, Ns, study_id (optional).
#' @param method 1 FAT-PET, 2 PEESE, 3 PET-PEESE, 4 EK.
#' @param weight 0 no weights, 1 standard weights, 2 MAIVE adjusted weights, 3 study weights.
#' @param instrument 1 yes, 0 no.
#' @param studylevel Correlation at study level: 0 none, 1 fixed effects, 2 cluster.
#' @param SE SE estimator: 0 CR0 (Huber-White), 1 CR1 (Standard empirical correction),
#' 2 CR2 (Bias-reduced estimator), 3 wild bootstrap.
#' @param AR Anderson Rubin corrected CI for weak instruments (available for unweighted and MAIVE-adjusted weight versions of
#' PET, PEESE, PET-PEESE, not available for fixed effects): 0 no, 1 yes.
#' @param first_stage First-stage specification for the variance model: 0 levels, 1 log.
#'
#' @details Data \code{dat} can be imported from an Excel file via:
#' \code{dat <- read_excel("inputdata.xlsx")} or from a csv file via: \code{dat <- read.csv("inputdata.csv")}
#' It should contain:
#' \itemize{
#'   \item Estimates: bs
#'   \item Standard errors: sebs
#'   \item Number of observations: Ns
#'   \item Optional: study_id
#' }
#' Default option for MAIVE: MAIVE-PET-PEESE, unweighted, instrumented, cluster SE, wild bootstrap, AR.
#'
#' @return \itemize{
#'   \item beta: MAIVE meta-estimate
#'   \item SE: MAIVE standard error
#'   \item F-test: heteroskedastic robust F-test of the first step instrumented SEs
#'   \item beta_standard: point estimate from the method chosen
#'   \item SE_standard: standard error from the method chosen
#'   \item Hausman: Hausman type test: comparison between MAIVE and standard version
#'   \item Chi2: 5% critical value for Hausman test
#'   \item SE_instrumented: instrumented standard errors
#'   \item AR_CI: Anderson-Rubin confidence interval for weak instruments
#'   \item pub bias p-value: p-value of test for publication bias / p-hacking based on instrumented FAT
#'   \item egger_coef: Egger Coefficient (PET estimate)
#'   \item egger_se: Egger Standard Error (PET standard error)
#'   \item egger_boot_ci: Confidence interval for the Egger coefficient using the selected resampling scheme
#'   \item egger_ar_ci: Anderson-Rubin confidence interval for the Egger coefficient (when available)
#'   \item is_quadratic_fit: Details on quadratic selection and slope behaviour
#'   \item boot_result: Boot result
#'   \item slope_coef: Slope coefficient
#'   \item petpeese_selected: Which model (PET or PEESE) was selected when method=3 (NA otherwise)
#'   \item peese_se2_coef: Coefficient on SE^2 when PEESE is the final model (NA otherwise)
#'   \item peese_se2_se: Standard error of the PEESE SE^2 coefficient (NA otherwise)
#' }
#'
#' @examples
#' dat <- data.frame(
#'   bs = c(0.5, 0.45, 0.55, 0.6),
#'   sebs = c(0.25, 0.2, 0.22, 0.27),
#'   Ns = c(50, 80, 65, 90)
#' )
#'
#' result <- maive(dat,
#'   method = 3, weight = 0, instrument = 1,
#'   studylevel = 0, SE = 0, AR = 0, first_stage = 0
#' )
#'
#' @export
maive <- function(dat, method, weight, instrument, studylevel, SE, AR, first_stage = 0L) {
  opts <- maive_validate_inputs(dat, method, weight, instrument, studylevel, SE, AR, first_stage)
  prepared <- maive_prepare_data(opts$dat, opts$studylevel)
  instrumentation <- maive_compute_variance_instrumentation(prepared$sebs, prepared$Ns, prepared$g, opts$type_choice, opts$instrument, opts$first_stage_type)

  w <- maive_compute_weights(opts$weight, prepared$sebs, instrumentation$sebs2fit1, prepared$studyid)

  maive_run_pipeline(opts, prepared, instrumentation, w)
}

#' WAIVE: Weighted Adjusted Instrumental Variable Estimator
#'
#' WAIVE extends MAIVE by applying exponential-decay weights that downweight
#' studies with spurious precision or extreme outlier behavior.
#'
#' @inheritParams maive
#' @return List with the same structure as \code{maive()}. See \code{?maive} for details.
#'
#' @details
#' Computes robust downweighting based on first-stage residuals. Studies with
#' negative residuals (spurious precision) or extreme residuals (outliers) receive
#' reduced influence in the meta-analytic estimate.
#'
#' @examples
#' dat <- data.frame(
#'   bs = c(0.5, 0.45, 0.55, 0.6),
#'   sebs = c(0.25, 0.2, 0.22, 0.27),
#'   Ns = c(50, 80, 65, 90)
#' )
#'
#' result <- waive(dat,
#'   method = 3, weight = 0, instrument = 1,
#'   studylevel = 0, SE = 0, AR = 0, first_stage = 0
#' )
#'
#' @seealso \code{\link{maive}}
#' @export
waive <- function(dat, method, weight, instrument, studylevel, SE, AR, first_stage = 0L) {
  opts <- maive_validate_inputs(dat, method, weight, instrument, studylevel, SE, AR, first_stage)
  prepared <- maive_prepare_data(opts$dat, opts$studylevel)
  instrumentation <- maive_compute_variance_instrumentation(prepared$sebs, prepared$Ns, prepared$g, opts$type_choice, opts$instrument, opts$first_stage_type)

  base_w <- maive_compute_weights(opts$weight, prepared$sebs, instrumentation$sebs2fit1, prepared$studyid)
  decay_weights <- maive_compute_waive_weights(instrumentation$first_stage_model)

  if (opts$weight == 0L) {
    w <- sqrt(decay_weights)
  } else {
    w <- base_w * sqrt(decay_weights)
  }

  maive_run_pipeline(opts, prepared, instrumentation, w)
}
