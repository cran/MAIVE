#' @keywords internal
PET_adjust <- function(bs, b0, b1, sebs) bs - b0 - b1 * sebs

#' @keywords internal
PEESE_adjust <- function(bs, b0, b1, sebs) bs - b0 - b1 * sebs^2

#' Joint Anderson-Rubin CI computation
#'
#' Computes AR confidence intervals using a 2D grid search over (b0, b1)
#' with chi^2_2 critical value. For the slope (Egger) coefficient under
#' weak instruments, use method = "slope_only" instead.
#'
#' @param model Fitted lm object from second-stage regression
#' @param adjust_fun PET_adjust or PEESE_adjust function
#' @param bs Effect estimates
#' @param sebs Standard errors
#' @param invNs Inverse sample sizes (instrument)
#' @param g Cluster variable
#' @param type_choice CR variance type ("CR0", "CR1", "CR2")
#' @param weights Optional weights for weighted AR
#' @param method "joint" for 2D grid, "slope_only" for subset AR (robust under weak IV)
#'
#' @keywords internal
compute_AR_CI_optimized <- function(model, adjust_fun, bs, sebs, invNs, g, type_choice, weights = NULL, method = "joint") {
  # Dispatch to subset AR for slope-only inference (recommended under weak instruments)
  if (method == "slope_only") {
    return(compute_AR_CI_slope_only(model, adjust_fun, bs, sebs, invNs, g, type_choice, weights))
  }

  # Joint 2D grid method for intercept (corrected mean) inference
  beta0 <- model$coefficients[1]
  beta1 <- model$coefficients[2]
  vc <- clubSandwich::vcovCR(model, cluster = g, type = type_choice)
  beta0se <- sqrt(vc[1, 1])
  beta1se <- sqrt(vc[2, 2])

  M <- length(bs)

  # Validate and process weights
  if (!is.null(weights)) {
    if (length(weights) != M) {
      stop("weights must have the same length as the data.")
    }
    if (any(!is.finite(weights) | weights <= 0)) {
      stop("weights must be positive and finite.")
    }
    weight_ratio <- max(weights) / min(weights)
    if (is.finite(weight_ratio) && weight_ratio > 1000) {
      warning("Adjusted weights vary substantially (ratio > 1000); AR CI may be unstable.")
    }
    sqrt_weights <- sqrt(weights)
  } else {
    sqrt_weights <- rep(1, M)
  }

  # Transform variables
  bs_t <- bs * sqrt_weights
  sebs_t <- sebs * sqrt_weights
  invNs_t <- invNs * sqrt_weights
  ones_vec <- sqrt_weights

  # Large dataset check
  if (M > 5000) {
    warning("Dataset too large for AR computation (", M, " observations). AR CI disabled.")
    return(list(b0_CI = c(NA_real_, NA_real_), b1_CI = c(NA_real_, NA_real_)))
  }

  # Grid resolution based on sample size
  base_resolution <- if (M > 1000) 30 else min(60, max(25, round(sqrt(M))))

  # Pre-compute projection matrices
  Z <- cbind(ones_vec, invNs_t)
  ZtZ_inv <- tryCatch(solve(t(Z) %*% Z), error = function(e) NULL)
  if (is.null(ZtZ_inv)) {
    warning("Instrument matrix is singular; AR CI unavailable.")
    return(list(b0_CI = c(NA_real_, NA_real_), b1_CI = c(NA_real_, NA_real_)))
  }
  PZ <- Z %*% ZtZ_inv %*% t(Z)
  MZ <- diag(M) - PZ

  # Pre-compute sebs term for adjustment
  sebs_term <- if (identical(adjust_fun, PET_adjust)) sebs_t else sebs_t^2

  # Vectorized AR statistic computation
  compute_AR_stats <- function(b0_vals, b1_vals) {
    n_b0 <- length(b0_vals)
    n_b1 <- length(b1_vals)
    stats_mat <- matrix(NA_real_, nrow = n_b0, ncol = n_b1)

    for (i in seq_len(n_b0)) {
      for (j in seq_len(n_b1)) {
        bs_star <- bs_t - b0_vals[i] - b1_vals[j] * sebs_term
        PZ_bs_star <- as.vector(PZ %*% bs_star)
        MZ_bs_star <- as.vector(MZ %*% bs_star)
        num <- sum(bs_star * PZ_bs_star)
        denom <- sum(bs_star * MZ_bs_star)
        if (abs(denom) > 1e-12) {
          stats_mat[i, j] <- (M - 2) * num / denom
        }
      }
    }
    stats_mat
  }

  # Chi^2_2(0.95) critical value
  crit_value <- qchisq(0.95, df = 2)

  # Adaptive grid search
  grid_mult <- 5
  max_mult <- 50
  accepted <- FALSE

  while (!accepted && grid_mult <= max_mult) {
    range0 <- if (is.finite(beta0se) && beta0se > 0) grid_mult * beta0se else grid_mult * max(1, abs(beta0))
    range1 <- if (is.finite(beta1se) && beta1se > 0) grid_mult * beta1se else grid_mult * max(1, abs(beta1))

    b0_grid <- seq(beta0 - range0, beta0 + range0, length.out = base_resolution)
    b1_grid <- seq(beta1 - range1, beta1 + range1, length.out = base_resolution)

    AR_stats <- compute_AR_stats(b0_grid, b1_grid)
    AR_accept <- !is.na(AR_stats) & AR_stats < crit_value

    accepted <- any(AR_accept)
    if (!accepted) grid_mult <- grid_mult * 2
  }

  if (!accepted) {
    warning("AR grid search failed to locate acceptance region.")
    return(list(b0_CI = c(NA_real_, NA_real_), b1_CI = c(NA_real_, NA_real_)))
  }

  # Extract CI bounds
  b0_accept_idx <- which(rowSums(AR_accept) > 0)
  b1_accept_idx <- which(colSums(AR_accept) > 0)

  if (length(b0_accept_idx) == 0 || length(b1_accept_idx) == 0) {
    return(list(b0_CI = c(NA_real_, NA_real_), b1_CI = c(NA_real_, NA_real_)))
  }

  # Check for disjoint regions
  has_gaps <- function(idx) length(idx) > 1 && any(diff(idx) > 1)
  if (has_gaps(b0_accept_idx) || has_gaps(b1_accept_idx)) {
    warning("AR acceptance region is disjoint; returning conservative interval.")
  }

  b0_CI <- round(c(min(b0_grid[b0_accept_idx]), max(b0_grid[b0_accept_idx])), 3)
  b1_CI <- round(c(min(b1_grid[b1_accept_idx]), max(b1_grid[b1_accept_idx])), 3)

  list(b0_CI = b0_CI, b1_CI = b1_CI)
}

#' Subset Anderson-Rubin CI for the slope coefficient
#'
#' Computes weak-instrument-robust AR confidence interval for the slope (Egger)
#' coefficient by treating the intercept as a nuisance parameter.
#'
#' This method is robust under weak instruments and avoids the "banana projection"
#' artifact that can produce spuriously narrow CIs when projecting from 2D
#' joint AR regions.
#'
#' Algorithm (per tjhavranek, Nov 2025):
#' 1. For each candidate slope b1, form residuals: r_i = y_i - b1 * x_i
#' 2. Run auxiliary regression: r ~ 1 + z (intercept absorbs b0)
#' 3. Use cluster-robust (CR2) variance for the z coefficient
#' 4. Test statistic: t_z^2 from the clustered test
#' 5. Accept b1 if t_z^2 <= qchisq(0.95, df = 1) = 3.84
#'
#' References:
#' - Guggenberger, P., Kleibergen, F. (2012). Econometrica.
#' - Andrews, D. W. K., & Mikusheva, A. (2016). Econometrica.
#'
#' @keywords internal
compute_AR_CI_slope_only <- function(model, adjust_fun, bs, sebs, invNs, g, type_choice, weights = NULL) {
  # Extract beta estimates and robust SEs for grid initialization

  beta1 <- model$coefficients[2]
  beta1se <- sqrt(clubSandwich::vcovCR(model, cluster = g, type = type_choice)[2, 2])

  M <- length(bs)

  # Handle weights
  if (!is.null(weights)) {
    if (length(weights) != M) {
      stop("weights must have the same length as the data used in the AR test.")
    }
    if (any(!is.finite(weights) | weights <= 0)) {
      stop("weights supplied to the AR test must be positive and finite.")
    }
    weight_ratio <- max(weights) / min(weights)
    if (is.finite(weight_ratio) && weight_ratio > 1000) {
      warning(
        "Adjusted weights vary substantially (ratio > 1000); Anderson-Rubin CI may be unstable."
      )
    }
    sqrt_weights <- sqrt(weights)
  } else {
    sqrt_weights <- rep(1, M)
  }

  # For extremely large datasets, disable AR computation
  if (M > 5000) {
    warning("Dataset too large for AR computation (", M, " observations). AR CI disabled.")
    return(list(b0_CI = c(NA_real_, NA_real_), b1_CI = c(NA_real_, NA_real_)))
  }

  # Transform variables using weight scaling
  # y = bs / w, x = sebs / w (for PET) or sebs^2 / w (for PEESE), z = invNs / w
  y <- bs * sqrt_weights
  if (identical(adjust_fun, PET_adjust)) {
    x <- sebs * sqrt_weights
  } else {
    x <- (sebs^2) * sqrt_weights
  }
  z <- invNs * sqrt_weights

  # Subset AR test function: for a given b1, test if z coefficient is significant
  # when regressing residuals r = y - b1*x on intercept + z
  compute_subset_ar_stat <- function(b1_val) {
    # Form residuals under the null slope
    r <- y - b1_val * x

    # Auxiliary regression: r ~ 1 + z
    # The intercept absorbs b0 (nuisance parameter)
    fit <- tryCatch(
      lm(r ~ z),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      return(Inf)
    }

    # Get cluster-robust variance for z coefficient
    vc <- tryCatch(
      clubSandwich::vcovCR(fit, cluster = g, type = type_choice),
      error = function(e) NULL
    )

    if (is.null(vc)) {
      return(Inf)
    }

    # Extract z coefficient and its SE
    coef_z <- coef(fit)["z"]
    se_z_sq <- vc["z", "z"]

    if (!is.finite(coef_z) || !is.finite(se_z_sq) || se_z_sq <= 0) {
      return(Inf)
    }

    # Test statistic: t_z^2
    t_z_sq <- (coef_z^2) / se_z_sq

    if (!is.finite(t_z_sq)) {
      return(Inf)
    }

    t_z_sq
  }

  # Critical value: chi^2_1(0.95) = 3.84

  crit_value <- qchisq(0.95, df = 1)

  # Create grid over b1 (slope) values
  # Start with +-10 SE, expand if needed
  grid_multiplier <- 10
  max_multiplier <- 100
  accepted <- FALSE
  b1_grid <- NULL
  ar_stats <- NULL

  # Use finer grid for better precision
  base_resolution <- min(200, max(50, round(sqrt(M) * 2)))

  while (!accepted && grid_multiplier <= max_multiplier) {
    range1 <- if (is.finite(beta1se) && beta1se > 0) {
      grid_multiplier * beta1se
    } else {
      grid_multiplier * max(1, abs(beta1))
    }

    l1 <- beta1 - range1
    u1 <- beta1 + range1

    b1_grid <- seq(l1, u1, length.out = base_resolution)

    # Compute subset AR statistic for each b1
    ar_stats <- sapply(b1_grid, compute_subset_ar_stat)

    # Check for accepted values
    accepted <- any(is.finite(ar_stats) & ar_stats < crit_value)

    if (!accepted) {
      grid_multiplier <- grid_multiplier * 2
    }
  }

  if (!accepted) {
    warning("Subset AR search failed to locate an acceptance region; returning NA interval.")
    return(list(b0_CI = c(NA_real_, NA_real_), b1_CI = c(NA_real_, NA_real_)))
  }

  # Find accepted b1 values
  b1_accept <- is.finite(ar_stats) & ar_stats < crit_value
  b1_accept_idx <- which(b1_accept)

  if (length(b1_accept_idx) == 0) {
    return(list(b0_CI = c(NA_real_, NA_real_), b1_CI = c(NA_real_, NA_real_)))
  }

  # Check for disjoint acceptance region
  detect_disjoint <- function(indices) {
    if (length(indices) <= 1) {
      return(FALSE)
    }
    any(diff(indices) > 1)
  }

  if (detect_disjoint(b1_accept_idx)) {
    warning("Subset AR acceptance region is disjoint; returning conservative interval spanning all segments.")
  }

  # Return slope CI
  b1_CI <- c(min(b1_grid[b1_accept_idx]), max(b1_grid[b1_accept_idx]))

  # b0_CI is not computed by this method (nuisance parameter)
  b0_CI <- c(NA_real_, NA_real_)

  list(
    b0_CI = b0_CI,
    b1_CI = round(b1_CI, 3)
  )
}
