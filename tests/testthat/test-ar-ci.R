test_that("AR CI is available with MAIVE-adjusted weights", {
  dat <- data.frame(
    bs = c(0.2, 0.25, 0.22, 0.3, 0.27),
    sebs = c(0.1, 0.12, 0.11, 0.13, 0.12),
    Ns = c(100, 120, 110, 130, 115)
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

  expect_true(is.numeric(result$AR_CI))
  expect_equal(length(result$AR_CI), 2L)
  expect_true(all(is.finite(result$AR_CI)))
})

test_that("standard weights continue to disable AR CI", {
  dat <- data.frame(
    bs = c(0.2, 0.25, 0.22, 0.3),
    sebs = c(0.1, 0.12, 0.11, 0.13),
    Ns = c(100, 120, 110, 130)
  )

  result <- suppressWarnings(maive(
    dat = dat,
    method = 1,
    weight = 1,
    instrument = 1,
    studylevel = 0,
    SE = 0,
    AR = 1,
    first_stage = 0
  ))

  expect_identical(result$AR_CI, "NA")
})

test_that("weighted AR computation handles weight diagnostics", {
  bs <- c(0.4, 0.5, 0.55, 0.48)
  sebs <- c(0.2, 0.18, 0.22, 0.19)
  invNs <- 1 / c(80, 90, 85, 88)
  model <- lm(bs ~ invNs)

  base <- MAIVE:::compute_AR_CI_optimized(
    model = model,
    adjust_fun = MAIVE:::PET_adjust,
    bs = bs,
    sebs = sebs,
    invNs = invNs,
    g = seq_along(bs),
    type_choice = "CR0"
  )

  identical_weights <- MAIVE:::compute_AR_CI_optimized(
    model = model,
    adjust_fun = MAIVE:::PET_adjust,
    bs = bs,
    sebs = sebs,
    invNs = invNs,
    g = seq_along(bs),
    type_choice = "CR0",
    weights = rep(1, length(bs))
  )

  expect_equal(base, identical_weights)

  expect_warning(
    MAIVE:::compute_AR_CI_optimized(
      model = model,
      adjust_fun = MAIVE:::PET_adjust,
      bs = bs,
      sebs = sebs,
      invNs = invNs,
      g = seq_along(bs),
      type_choice = "CR0",
      weights = c(1e-6, 1, 5, 10)
    ),
    "Adjusted weights vary substantially"
  )
})

test_that("Egger AR CI uses instrumented SE when instrument=1", {
  # When instrument=1, the AR CI should use instrumented SE (sqrt(sebs2fit1))

  # not the original SE. This ensures consistency with the fitted model.
  set.seed(456)
  n <- 30
  Ns <- exp(rnorm(n, 7, 1))
  sebs <- 0.15 / sqrt(Ns) + abs(rnorm(n, 0, 0.02))
  bs <- 0.3 + 0.5 * sebs + rnorm(n, 0, 0.05)

  dat <- data.frame(bs = bs, sebs = sebs, Ns = Ns)

  # Run with instrumentation on
  result <- suppressWarnings(maive(
    dat = dat,
    method = 1,
    weight = 0,
    instrument = 1,
    studylevel = 0,
    SE = 0,
    AR = 1,
    first_stage = 0
  ))

  # Egger AR CI should be computed and not truncated at point estimate
  expect_true(!is.null(result$egger_ar_ci))
  expect_true(is.numeric(result$egger_ar_ci))

  if (all(is.finite(result$egger_ar_ci))) {
    # The CI upper bound should NOT equal the point estimate
    # (which was the bug when using wrong SE)
    expect_true(
      abs(result$egger_ar_ci["upper"] - result$egger_coef) > 0.01 ||
        result$egger_ar_ci["upper"] > result$egger_coef,
      info = "AR CI upper bound should not be truncated at point estimate"
    )
  }
})

test_that("AR CI uses correct instrument based on first_stage", {
  # When first_stage=1 (log), the AR test should use log(Ns) as instrument
  # When first_stage=0 (levels), the AR test should use 1/Ns as instrument
  set.seed(789)
  n <- 40
  Ns <- exp(rnorm(n, 8, 1))
  sebs <- 0.1 / sqrt(Ns) + abs(rnorm(n, 0, 0.01))
  bs <- 0.2 + 0.4 * sebs + rnorm(n, 0, 0.03)

  dat <- data.frame(bs = bs, sebs = sebs, Ns = Ns)

  # Run with levels first-stage
  result_levels <- suppressWarnings(maive(
    dat = dat,
    method = 1,
    weight = 0,
    instrument = 1,
    studylevel = 0,
    SE = 0,
    AR = 1,
    first_stage = 0
  ))

  # Run with log first-stage
  result_log <- suppressWarnings(maive(
    dat = dat,
    method = 1,
    weight = 0,
    instrument = 1,
    studylevel = 0,
    SE = 0,
    AR = 1,
    first_stage = 1
  ))

  # Both should produce valid AR CIs
  expect_true(is.numeric(result_levels$egger_ar_ci))
  expect_true(is.numeric(result_log$egger_ar_ci))

  # The CIs may differ due to different instruments and instrumented SE
  # Both should be reasonable (not NA, not truncated)
  if (all(is.finite(result_levels$egger_ar_ci)) && all(is.finite(result_log$egger_ar_ci))) {
    levels_width <- diff(result_levels$egger_ar_ci)
    log_width <- diff(result_log$egger_ar_ci)
    # Both should produce CIs with positive width
    expect_true(levels_width > 0, info = "Levels first-stage AR CI should have positive width")
    expect_true(log_width > 0, info = "Log first-stage AR CI should have positive width")
  }
})

test_that("Egger AR CI always uses subset method for robust inference", {
  # The egger_ar_ci should always use subset AR method (not joint)
  # This avoids banana-projection artifacts even with strong instruments
  set.seed(321)
  n <- 50
  Ns <- exp(rnorm(n, 9, 0.5)) # Strong instrument (clear relationship with SE)
  sebs <- 0.2 / sqrt(Ns) + abs(rnorm(n, 0, 0.005))
  bs <- 0.5 + 0.3 * sebs + rnorm(n, 0, 0.02)

  dat <- data.frame(bs = bs, sebs = sebs, Ns = Ns)

  result <- suppressWarnings(maive(
    dat = dat,
    method = 1,
    weight = 0,
    instrument = 1,
    studylevel = 0,
    SE = 2, # CR2
    AR = 1,
    first_stage = 0
  ))

  # Should have valid AR CI
  expect_true(is.numeric(result$egger_ar_ci))

  if (all(is.finite(result$egger_ar_ci))) {
    # Compute Wald CI for comparison
    wald_width <- 2 * qnorm(0.975) * result$egger_se
    ar_width <- diff(result$egger_ar_ci)

    # Subset AR should produce reasonable CI width (not spuriously narrow)
    # It should be at least 50% of Wald width
    expect_true(
      ar_width >= 0.5 * wald_width || ar_width > 0.5,
      info = sprintf(
        "Egger AR CI width (%.3f) is suspiciously narrow compared to Wald (%.3f)",
        ar_width, wald_width
      )
    )
  }
})

test_that("instrumentation returns correct instrument_for_ar", {
  # Verify that maive_compute_variance_instrumentation returns the correct
# instrument based on first_stage_type
  set.seed(555)
  n <- 20
  Ns <- exp(rnorm(n, 7, 1))
  sebs <- 0.1 + abs(rnorm(n, 0, 0.02))
  g <- seq_len(n)

  # Levels first-stage should return 1/Ns
  result_levels <- MAIVE:::maive_compute_variance_instrumentation(
    sebs = sebs,
    Ns = Ns,
    g = g,
    type_choice = "CR0",
    instrument = 1L,
    first_stage_type = "levels"
  )
  expect_equal(result_levels$instrument_for_ar, 1 / Ns)

  # Log first-stage should return log(Ns)
  result_log <- MAIVE:::maive_compute_variance_instrumentation(
    sebs = sebs,
    Ns = Ns,
    g = g,
    type_choice = "CR0",
    instrument = 1L,
    first_stage_type = "log"
  )
  expect_equal(result_log$instrument_for_ar, log(Ns))
})
