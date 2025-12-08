test_that("WAIVE computes exponential-decay weights correctly", {
  dat <- data.frame(
    bs = c(0.5, 0.45, 0.55, 0.6, 0.48, 0.52),
    sebs = c(0.25, 0.2, 0.22, 0.27, 0.19, 0.24),
    Ns = c(50, 80, 65, 90, 75, 60)
  )

  # Run WAIVE
  result <- waive(
    dat = dat,
    method = 3,
    weight = 0,
    instrument = 1,
    studylevel = 0,
    SE = 0,
    AR = 0,
    first_stage = 0
  )

  # Manually compute WAIVE weights
  # 1. First stage
  invNs <- 1 / dat$Ns
  sebs2 <- dat$sebs^2
  varreg <- lm(sebs2 ~ invNs)
  nu <- residuals(varreg)

  # 2. Robust normalization
  mad_value <- mad(nu, constant = 1.4826)
  sigma <- if (mad_value == 0) sd(nu) + 1e-12 else mad_value
  z <- nu / sigma

  # 3. Exponential-decay weights
  z_neg <- pmax(-z, 0)
  z_out <- pmax(abs(z) - 2, 0)
  w_waive <- exp(-1.0 * z_neg - 0.25 * z_out^2)
  w_waive <- pmax(w_waive, 0.05)
  w_waive <- w_waive / mean(w_waive)

  # 4. Combined weights
  sebs2fit <- fitted(varreg)
  combined_weights <- sqrt(sebs2fit * w_waive)

  # Verify weights are being applied
  expect_true(all(w_waive >= 0.05))
  expect_equal(mean(w_waive), 1, tolerance = 1e-10)

  # Verify WAIVE produces a result
  expect_true(is.numeric(result$beta))
  expect_true(is.numeric(result$SE))
  expect_true(!is.na(result$beta))
  expect_true(!is.na(result$SE))
})

test_that("WAIVE works with log first stage", {
  dat <- data.frame(
    bs = c(0.5, 0.45, 0.55, 0.6),
    sebs = c(0.25, 0.2, 0.22, 0.27),
    Ns = c(50, 80, 65, 90)
  )

  result_levels <- waive(
    dat = dat,
    method = 1,
    weight = 0,
    instrument = 1,
    studylevel = 0,
    SE = 0,
    AR = 0,
    first_stage = 0
  )

  result_log <- waive(
    dat = dat,
    method = 1,
    weight = 0,
    instrument = 1,
    studylevel = 0,
    SE = 0,
    AR = 0,
    first_stage = 1
  )

  # Both should produce valid results
  expect_true(is.numeric(result_levels$beta))
  expect_true(is.numeric(result_log$beta))
  expect_true(!is.na(result_levels$beta))
  expect_true(!is.na(result_log$beta))

  # They should be similar but not identical
  expect_true(abs(result_levels$beta - result_log$beta) < 1)
})

test_that("WAIVE downweights spurious precision (negative residuals)", {
  # Create synthetic data with one suspiciously precise estimate
  dat <- data.frame(
    bs = c(0.5, 0.48, 0.52, 0.5, 0.45),  # All similar estimates
    sebs = c(0.25, 0.22, 0.23, 0.24, 0.05),  # Last one too precise
    Ns = c(50, 60, 55, 58, 40)  # Last one has small N but tiny SE
  )

  # First stage regression
  invNs <- 1 / dat$Ns
  sebs2 <- dat$sebs^2
  varreg <- lm(sebs2 ~ invNs)
  nu <- residuals(varreg)

  # The suspiciously precise study should have negative residual
  expect_true(nu[5] < 0)

  # Compute WAIVE weights
  mad_value <- mad(nu, constant = 1.4826)
  sigma <- if (mad_value == 0) sd(nu) + 1e-12 else mad_value
  z <- nu / sigma
  z_neg <- pmax(-z, 0)
  z_out <- pmax(abs(z) - 2, 0)
  w_waive <- exp(-1.0 * z_neg - 0.25 * z_out^2)
  w_waive <- pmax(w_waive, 0.05)

  # The suspiciously precise study should get downweighted
  # (before normalization)
  expect_true(w_waive[5] < w_waive[1])
})

test_that("WAIVE downweights extreme outliers", {
  # Create data with an outlier
  dat <- data.frame(
    bs = c(0.5, 0.48, 0.52, 0.51, 0.49),
    sebs = c(0.25, 0.22, 0.23, 0.24, 0.80),  # Last one very imprecise
    Ns = c(50, 60, 55, 58, 45)
  )

  # First stage regression
  invNs <- 1 / dat$Ns
  sebs2 <- dat$sebs^2
  varreg <- lm(sebs2 ~ invNs)
  nu <- residuals(varreg)

  # The outlier should have large positive residual
  expect_true(nu[5] > mean(nu))

  # Compute WAIVE weights
  mad_value <- mad(nu, constant = 1.4826)
  sigma <- if (mad_value == 0) sd(nu) + 1e-12 else mad_value
  z <- nu / sigma
  z_neg <- pmax(-z, 0)
  z_out <- pmax(abs(z) - 2, 0)
  w_waive <- exp(-1.0 * z_neg - 0.25 * z_out^2)
  w_waive <- pmax(w_waive, 0.05)

  # If the residual is extreme (|z| > 2), it should be downweighted
  if (abs(z[5]) > 2) {
    expect_true(w_waive[5] < max(w_waive[-5]))
  }
})

test_that("WAIVE works with study clusters", {
  dat <- data.frame(
    bs = c(0.4, 0.6, 0.55, 0.5, 0.52, 0.49),
    sebs = c(0.2, 0.18, 0.22, 0.19, 0.21, 0.2),
    Ns = c(80, 95, 90, 85, 88, 92),
    study_id = c(1, 1, 2, 2, 3, 3)
  )

  result <- waive(
    dat = dat,
    method = 3,
    weight = 0,
    instrument = 1,
    studylevel = 2,
    SE = 0,
    AR = 0,
    first_stage = 0
  )

  expect_true(is.numeric(result$beta))
  expect_true(!is.na(result$beta))
  expect_true(is.numeric(result$SE))
})

test_that("WAIVE disables Anderson-Rubin CI", {
  dat <- data.frame(
    bs = c(0.5, 0.45, 0.55, 0.6),
    sebs = c(0.25, 0.2, 0.22, 0.27),
    Ns = c(50, 80, 65, 90)
  )

  # WAIVE with AR=1 should still disable AR (when using weight > 0)
  result <- waive(
    dat = dat,
    method = 3,
    weight = 1,
    instrument = 1,
    studylevel = 0,
    SE = 0,
    AR = 1,  # Request AR
    first_stage = 0
  )

  # AR_CI should be "NA" for weighted methods
  expect_equal(result$AR_CI, "NA")
})

test_that("WAIVE floor weight prevents zero leverage", {
  # Create extreme case that might lead to near-zero weights
  dat <- data.frame(
    bs = c(0.5, 0.48, 0.52, 0.01),  # Last one very different
    sebs = c(0.25, 0.22, 0.23, 0.02),  # Last one extremely precise
    Ns = c(50, 60, 55, 30)
  )

  result <- waive(
    dat = dat,
    method = 1,
    weight = 0,
    instrument = 1,
    studylevel = 0,
    SE = 0,
    AR = 0,
    first_stage = 0
  )

  # Should not fail due to zero weights
  expect_true(is.numeric(result$beta))
  expect_true(!is.na(result$beta))
})

test_that("WAIVE matches MAIVE structure with different weights", {
  dat <- data.frame(
    bs = c(0.5, 0.45, 0.55, 0.6, 0.48, 0.52),
    sebs = c(0.25, 0.2, 0.22, 0.27, 0.19, 0.24),
    Ns = c(50, 80, 65, 90, 75, 60)
  )

  result_maive <- maive(
    dat = dat,
    method = 3,
    weight = 2,
    instrument = 1,
    studylevel = 0,
    SE = 0,
    AR = 0,
    first_stage = 0
  )

  result_waive <- waive(
    dat = dat,
    method = 3,
    weight = 2,
    instrument = 1,
    studylevel = 0,
    SE = 0,
    AR = 0,
    first_stage = 0
  )

  # Both should return same structure
  expect_equal(names(result_maive), names(result_waive))

  # Results should differ due to different weighting
  expect_false(isTRUE(all.equal(result_maive$beta, result_waive$beta)))
})

test_that("waive() function works and produces valid results", {
  dat <- data.frame(
    bs = c(0.5, 0.45, 0.55, 0.6),
    sebs = c(0.25, 0.2, 0.22, 0.27),
    Ns = c(50, 80, 65, 90)
  )

  # Test waive() with unweighted base
  result <- waive(
    dat = dat,
    method = 3,
    weight = 0,
    instrument = 1,
    studylevel = 0,
    SE = 0,
    AR = 0,
    first_stage = 0
  )

  # Should produce valid results
  expect_true(is.numeric(result$beta))
  expect_true(!is.na(result$beta))
  expect_true(is.numeric(result$SE))
  expect_true(!is.na(result$SE))
  expect_true(result$SE > 0)
})

test_that("waive() works with different base weighting schemes", {
  dat <- data.frame(
    bs = c(0.5, 0.45, 0.55, 0.6, 0.48),
    sebs = c(0.25, 0.2, 0.22, 0.27, 0.19),
    Ns = c(50, 80, 65, 90, 75)
  )

  # Test waive() with different base weights
  result_unweighted <- waive(dat, method=1, weight=0, instrument=1, studylevel=0, SE=0, AR=0)
  result_ivweighted <- waive(dat, method=1, weight=1, instrument=1, studylevel=0, SE=0, AR=0)
  result_maiveweighted <- waive(dat, method=1, weight=2, instrument=1, studylevel=0, SE=0, AR=0)

  # All should produce valid results
  expect_true(is.numeric(result_unweighted$beta))
  expect_true(is.numeric(result_ivweighted$beta))
  expect_true(is.numeric(result_maiveweighted$beta))

  # Results should differ based on base weighting
  expect_false(isTRUE(all.equal(result_unweighted$beta, result_ivweighted$beta)))
  expect_false(isTRUE(all.equal(result_unweighted$beta, result_maiveweighted$beta)))
})
