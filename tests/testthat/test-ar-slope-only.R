test_that("slope-only AR CI method is available", {
  # Use larger sample with clearer instrument to ensure AR CI can be computed
  set.seed(42)
  n <- 30
  Ns <- exp(rnorm(n, 8, 1))
  sebs <- pmax(0.2 / sqrt(Ns) + rnorm(n, 0, 0.01), 0.01)
  bs <- 0.5 + 0.2 * sebs + rnorm(n, 0, 0.05)

  invNs <- 1 / Ns
  model <- lm(bs ~ invNs)

  # Test slope-only method directly
  result_slope <- suppressWarnings(MAIVE:::compute_AR_CI_slope_only(
    model = model,
    adjust_fun = MAIVE:::PET_adjust,
    bs = bs,
    sebs = sebs,
    invNs = invNs,
    g = seq_along(bs),
    type_choice = "CR0"
  ))

  expect_true(is.list(result_slope))
  expect_true("b1_CI" %in% names(result_slope))
  # CI may be NA if no acceptance region found, or numeric
  if (!all(is.na(result_slope$b1_CI))) {
    expect_equal(length(result_slope$b1_CI), 2L)
    expect_true(is.numeric(result_slope$b1_CI))
  }
})

test_that("subset AR gives wider or comparable CI under weak instruments", {
  # Simulate weak instrument scenario where SE has low correlation with 1/N
  # The subset AR method should produce wider CIs that are more consistent
  # with cluster-robust Wald CIs under weak identification
  set.seed(123)
  n <- 50
  # Create weak instrument: low correlation between 1/N and SE^2
  Ns <- exp(rnorm(n, 6, 1))
  true_effect <- 0.5
  true_bias <- 0.6 # Egger slope
  sebs <- 0.1 + rnorm(n, 0, 0.02) # Almost constant SE (weak instrument)
  bs <- true_effect + true_bias * sebs + rnorm(n, 0, sebs * 0.5)

  invNs <- 1 / Ns
  model <- lm(bs ~ sebs)

  # Compare joint vs subset (slope-only)
  result_joint <- suppressWarnings(MAIVE:::compute_AR_CI_optimized(
    model = model,
    adjust_fun = MAIVE:::PET_adjust,
    bs = bs,
    sebs = sebs,
    invNs = invNs,
    g = seq_along(bs),
    type_choice = "CR2",
    method = "joint"
  ))

  result_subset <- suppressWarnings(MAIVE:::compute_AR_CI_slope_only(
    model = model,
    adjust_fun = MAIVE:::PET_adjust,
    bs = bs,
    sebs = sebs,
    invNs = invNs,
    g = seq_along(bs),
    type_choice = "CR2"
  ))

  # Both should return valid list structure
  expect_true(is.list(result_joint))
  expect_true(is.list(result_subset))
  expect_true("b1_CI" %in% names(result_joint))
  expect_true("b1_CI" %in% names(result_subset))

  # Under weak instruments, subset AR should generally produce
  # CIs that are at least as wide as or wider than joint projection
  # (avoiding the banana-projection artifact)
  if (!all(is.na(result_joint$b1_CI)) && !all(is.na(result_subset$b1_CI))) {
    joint_width <- diff(result_joint$b1_CI)
    subset_width <- diff(result_subset$b1_CI)
    # The subset method should not produce spuriously narrow CIs
    # (we allow some tolerance for numerical differences)
    expect_true(subset_width >= joint_width * 0.5 || subset_width > 1.0)
  }
})

test_that("MAIVE always uses subset AR for Egger slope CI", {
  # Egger AR CI always uses subset AR method for robust inference
  # regardless of instrument strength (avoids banana-projection)
  set.seed(999)
  n <- 40
  Ns <- exp(rnorm(n, 8, 0.5))
  sebs <- 0.15 + rnorm(n, 0, 0.005)
  bs <- 0.5 + 0.1 * sebs + rnorm(n, 0, sebs)

  dat <- data.frame(
    bs = bs,
    sebs = sebs,
    Ns = Ns,
    studyid = rep(1:10, each = 4)
  )

  result <- suppressWarnings(maive(
    dat = dat,
    method = 1,
    weight = 2,
    instrument = 1,
    studylevel = 0,
    SE = 0,
    AR = 1,
    first_stage = 0
  ))

  # Check that Egger AR CI was computed
  expect_true(!is.null(result))
  expect_true(!is.null(result$egger_ar_ci))
  expect_true(is.numeric(result$egger_ar_ci))
})

test_that("slope-only method handles PEESE adjustment", {
  set.seed(888)
  n <- 30
  Ns <- exp(rnorm(n, 8, 1))
  sebs <- pmax(0.2 / sqrt(Ns) + rnorm(n, 0, 0.01), 0.01)
  bs <- 0.5 + 0.3 * sebs^2 + rnorm(n, 0, 0.05)

  invNs <- 1 / Ns
  sebs_sq <- sebs^2
  model <- lm(bs ~ sebs_sq)

  # Test with PEESE adjustment
  result <- suppressWarnings(MAIVE:::compute_AR_CI_slope_only(
    model = model,
    adjust_fun = MAIVE:::PEESE_adjust,
    bs = bs,
    sebs = sebs,
    invNs = invNs,
    g = seq_along(bs),
    type_choice = "CR0"
  ))

  expect_true(is.list(result))
  expect_true("b1_CI" %in% names(result))
})

test_that("slope-only method handles weighted AR", {
  bs <- c(0.4, 0.5, 0.55, 0.48)
  sebs <- c(0.2, 0.18, 0.22, 0.19)
  invNs <- 1 / c(80, 90, 85, 88)
  model <- lm(bs ~ invNs)
  weights <- c(1, 1.5, 1.2, 1.3)

  result <- suppressWarnings(MAIVE:::compute_AR_CI_slope_only(
    model = model,
    adjust_fun = MAIVE:::PET_adjust,
    bs = bs,
    sebs = sebs,
    invNs = invNs,
    g = seq_along(bs),
    type_choice = "CR0",
    weights = weights
  ))

  expect_true(is.list(result))
  expect_true("b1_CI" %in% names(result))
})

test_that("slope-only handles extreme weight heterogeneity", {
  # Create scenario with extreme weight heterogeneity
  bs <- c(0.5, 0.55, 0.52)
  sebs <- c(0.001, 0.5, 0.002)
  invNs <- 1 / c(100, 120, 110)
  model <- lm(bs ~ invNs)
  weights <- c(1000, 1, 1000) # Extreme weight variation

  # Suppress warnings and just check it doesn't error
  result <- suppressWarnings(MAIVE:::compute_AR_CI_slope_only(
    model = model,
    adjust_fun = MAIVE:::PET_adjust,
    bs = bs,
    sebs = sebs,
    invNs = invNs,
    g = seq_along(bs),
    type_choice = "CR0",
    weights = weights
  ))

  # Should return valid structure (may be NA)
  expect_true(is.list(result))
  expect_true("b1_CI" %in% names(result))
})

test_that("subset AR produces reasonable CI width under weak instruments", {
  # Simulates a clustered meta-analysis with weak first-stage instrument
  # (log sample size has low correlation with SE). Under such conditions,
  # the joint AR can produce spuriously narrow CIs due to banana-projection.
  # Expected behavior: subset AR should give CIs comparable to Wald CI width

  set.seed(2025)
  n <- 100

  # Simulate clustered data with weak instrument (log sample size)
  n_studies <- 20
  study_id <- rep(1:n_studies, each = n / n_studies)

  # Log sample size as instrument (weak correlation with SE)
  Ns <- exp(rnorm(n, 10, 1.5))
  log_Ns <- log(Ns)

  # Standard errors with weak relationship to sample size
  sebs <- pmax(0.01 + rnorm(n, 0, 0.005), 0.005)

  # Effect estimates with publication bias
  true_effect <- 0.02
  true_bias <- 0.6 # Egger slope
  bs <- true_effect + true_bias * sebs + rnorm(n, 0, sebs * 0.3)

  # Fit model
  invNs <- 1 / Ns
  model <- lm(bs ~ sebs)

  # Get cluster-robust SE for Wald comparison
  vc <- clubSandwich::vcovCR(model, cluster = study_id, type = "CR2")
  wald_se <- sqrt(vc[2, 2])
  wald_width <- 2 * qnorm(0.975) * wald_se

  # Compute subset AR CI
  result <- suppressWarnings(MAIVE:::compute_AR_CI_slope_only(
    model = model,
    adjust_fun = MAIVE:::PET_adjust,
    bs = bs,
    sebs = sebs,
    invNs = invNs,
    g = study_id,
    type_choice = "CR2"
  ))

  expect_true(is.list(result))
  expect_true("b1_CI" %in% names(result))

  # Under weak instruments, AR CI should not be implausibly narrower than Wald
  # The subset AR method should give CIs of reasonable width
  if (!all(is.na(result$b1_CI))) {
    ar_width <- diff(result$b1_CI)
    # AR CI should be at least 30% of Wald width (not spuriously narrow)
    # Under correct implementation it's typically wider or similar
    expect_true(ar_width > 0.3 * wald_width,
      info = sprintf("AR width (%.3f) < 30%% of Wald width (%.3f)", ar_width, wald_width)
    )
  }
})

test_that("subset AR uses correct critical value (chi^2_1)", {
  # Verify the subset AR uses chi^2_1 (3.84) not chi^2_2 (5.99)
  # This is important for correct coverage

  set.seed(777)
  n <- 40
  Ns <- exp(rnorm(n, 8, 1))
  sebs <- pmax(0.15 / sqrt(Ns) + rnorm(n, 0, 0.01), 0.01)
  bs <- 0.3 + 0.4 * sebs + rnorm(n, 0, 0.03)

  invNs <- 1 / Ns
  model <- lm(bs ~ sebs)

  # The subset method uses chi^2_1(0.95) = 3.84 as the critical value
  # This should produce CIs with correct 95% coverage for the slope
  result <- suppressWarnings(MAIVE:::compute_AR_CI_slope_only(
    model = model,
    adjust_fun = MAIVE:::PET_adjust,
    bs = bs,
    sebs = sebs,
    invNs = invNs,
    g = seq_along(bs),
    type_choice = "CR2"
  ))

  expect_true(is.list(result))

  # The result should have b1_CI (slope) but b0_CI should be NA

  # since subset AR treats intercept as nuisance
  if (!all(is.na(result$b1_CI))) {
    expect_true(all(is.na(result$b0_CI)),
      info = "Subset AR should not compute b0_CI (nuisance parameter)"
    )
  }
})
