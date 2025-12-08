#' @keywords internal
manual_wild_cluster_boot_se <- function(model, data, cluster_var, B = 500, seed = 123) {
  set.seed(seed)

  # Extract residuals and fitted values
  resids <- residuals(model)
  fitted_vals <- fitted(model)

  # Extract weights (if any)
  orig_weights <- model$weights
  has_weights <- !is.null(orig_weights) && length(orig_weights) > 0

  # Get cluster IDs
  clusters <- unique(as.character(data[[cluster_var]]))
  G <- length(clusters)

  # Coefficient names and original estimates
  coef_names <- names(coef(model))
  k <- length(coef_names)
  base_coefs <- coef(model)

  # Compute cluster-robust SE (CR1) for the base model
  # This is used as the denominator in t-statistic computation
  V_cr1 <- clubSandwich::vcovCR(model, cluster = data[[cluster_var]], type = "CR1")
  base_se_cr1 <- sqrt(diag(V_cr1))

  # Matrix to store bootstrap coefficients, bootstrap SEs, and t-statistics
  boot_coefs <- matrix(NA, nrow = B, ncol = k)
  boot_se_cr1 <- matrix(NA, nrow = B, ncol = k)
  boot_t_stats <- matrix(NA, nrow = B, ncol = k)
  colnames(boot_coefs) <- coef_names
  colnames(boot_se_cr1) <- coef_names
  colnames(boot_t_stats) <- coef_names

  # Loop over bootstrap replications
  for (b in 1:B) {
    # Draw Rademacher multipliers per cluster (one per study)
    u_g <- sample(c(-1, 1), size = G, replace = TRUE)
    names(u_g) <- as.character(clusters)

    # Create bootstrap outcome
    data$y_boot <- fitted_vals + resids * u_g[as.character(data[[cluster_var]])]

    # Refit the same model formula on bootstrap sample
    # IMPORTANT: preserve weights from original model
    form_boot <- update(formula(model), y_boot ~ .)
    if (has_weights) {
      fit_boot <- lm(form_boot, data = data, weights = orig_weights)
    } else {
      fit_boot <- lm(form_boot, data = data)
    }

    # Store bootstrap coefficients
    coef_boot <- coef(fit_boot)
    if (is.null(names(coef_boot))) {
      names(coef_boot) <- coef_names
    }
    boot_coefs[b, ] <- coef_boot[coef_names]

    # Compute bootstrap clustered SEs (CR1) for studentized statistics
    V_boot <- tryCatch(
      clubSandwich::vcovCR(fit_boot, cluster = data[[cluster_var]], type = "CR1"),
      error = function(e) NULL
    )

    if (!is.null(V_boot)) {
      se_boot <- sqrt(diag(V_boot))
      se_boot_named <- rep(NA_real_, k)
      names(se_boot_named) <- coef_names
      if (!is.null(names(se_boot))) {
        se_boot_named[names(se_boot)] <- se_boot
      } else if (length(se_boot) == k) {
        se_boot_named <- se_boot
        names(se_boot_named) <- coef_names
      }
      boot_se_cr1[b, ] <- se_boot_named
    }
  }

  # Compute t-statistics for each bootstrap replication
  # t_b = (beta_b - beta_base) / SE*_b
  for (j in 1:k) {
    diffs <- boot_coefs[, j] - base_coefs[j]
    denom <- boot_se_cr1[, j]
    valid <- !is.na(denom) & denom != 0
    boot_t_stats[, j] <- NA_real_
    boot_t_stats[valid, j] <- diffs[valid] / denom[valid]
  }

  # Compute t-percentile CI
  # CI = [beta_base - t_upper * SE_CR1, beta_base - t_lower * SE_CR1]
  alpha <- 0.05
  boot_ci <- matrix(NA, nrow = k, ncol = 2)
  rownames(boot_ci) <- coef_names
  colnames(boot_ci) <- c("lower", "upper")

  for (j in 1:k) {
    if (all(is.na(boot_t_stats[, j]))) {
      t_quantiles <- c(NA_real_, NA_real_)
    } else {
      t_quantiles <- stats::quantile(
        boot_t_stats[, j],
        probs = c(1 - alpha / 2, alpha / 2),
        na.rm = TRUE,
        names = FALSE
      )
    }
    boot_ci[j, "lower"] <- base_coefs[j] - t_quantiles[1] * base_se_cr1[j]
    boot_ci[j, "upper"] <- base_coefs[j] - t_quantiles[2] * base_se_cr1[j]
  }

  # Return CR1 SE (not bootstrap SD) as the SE
  # This ensures consistency across SE methods
  boot_se <- base_se_cr1

  # Return list
  return(list(
    boot_se = boot_se,
    boot_ci = boot_ci,
    boot_coefs = boot_coefs,
    boot_t_stats = boot_t_stats,
    boot_rep_se = boot_se_cr1
  ))
}
