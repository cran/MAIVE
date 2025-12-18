test_that("log first stage applies smearing retransformation", {
  dat <- data.frame(
    bs = c(0.5, 0.45, 0.55, 0.6),
    sebs = c(0.25, 0.2, 0.22, 0.27),
    Ns = c(50, 80, 65, 90)
  )

  result <- maive(
    dat = dat,
    method = 1,
    weight = 0,
    instrument = 1,
    studylevel = 0,
    SE = 0,
    AR = 0,
    first_stage = 1
  )

  log_model <- lm(log(dat$sebs^2) ~ log(dat$Ns))
  smearing <- mean(exp(residuals(log_model)))
  sehat_manual <- exp(predict(log_model)) * smearing

  expect_equal(result$SE_instrumented, sqrt(sehat_manual), tolerance = 1e-10)

  manual_vcov <- clubSandwich::vcovCR(log_model, cluster = seq_len(nrow(dat)), type = "CR0")
  slope <- coef(log_model)[2]
  manual_F <- unname(round(slope^2 / manual_vcov[2, 2], 3))
  expect_equal(result$`F-test`, manual_F)
})

test_that("first-stage F-test does not crash when Ns is constant (rank deficient)", {
  dat <- data.frame(
    bs = c(0.4, 0.6, 0.55, 0.5, 0.52, 0.49),
    sebs = c(0.2, 0.18, 0.22, 0.19, 0.21, 0.2),
    Ns = rep(100, 6)
  )

  # With constant Ns, the instrument (1/Ns) has no variation so IV is not identified.
  # MAIVE should not error; it should fall back to instrument=0 and report F-test as "NA".
  result <- maive(
    dat = dat,
    method = 1,
    weight = 0,
    instrument = 1,
    studylevel = 0,
    SE = 0,
    AR = 0,
    first_stage = 0
  )

  expect_identical(result$`F-test`, "NA")
  expect_true(is.finite(as.numeric(result$beta)))
  expect_true(is.finite(as.numeric(result$SE)))
})

test_that("Hausman statistic uses difference-in-estimators variance", {
  dat <- data.frame(
    bs = c(0.4, 0.6, 0.55, 0.5, 0.52, 0.49),
    sebs = c(0.2, 0.18, 0.22, 0.19, 0.21, 0.2),
    Ns = c(80, 95, 90, 85, 88, 92),
    study_id = c(1, 1, 2, 2, 3, 3)
  )

  result <- maive(
    dat = dat,
    method = 1,
    weight = 0,
    instrument = 1,
    studylevel = 2,
    SE = 0,
    AR = 1,
    first_stage = 0
  )

  opts <- MAIVE:::maive_validate_inputs(dat, 1, 0, 1, 2, 0, 1, 0)
  prepared <- MAIVE:::maive_prepare_data(opts$dat, opts$studylevel)
  instrumentation <- MAIVE:::maive_compute_variance_instrumentation(
    prepared$sebs,
    prepared$Ns,
    prepared$g,
    opts$type_choice,
    opts$instrument,
    opts$first_stage_type
  )
  w <- MAIVE:::maive_compute_weights(opts$weight, prepared$sebs, instrumentation$sebs2fit1, prepared$studyid)
  x <- if (opts$instrument == 0L) prepared$sebs else sqrt(instrumentation$sebs2fit1)
  x2 <- if (opts$instrument == 0L) prepared$sebs^2 else instrumentation$sebs2fit1
  design <- MAIVE:::maive_build_design_matrices(prepared$bs, prepared$sebs, w, x, x2, prepared$D, prepared$dummy)
  fits <- MAIVE:::maive_fit_models(design)
  selection <- MAIVE:::maive_select_petpeese(fits, design, opts$alpha_s)
  sighats <- MAIVE:::maive_compute_sigma_h(fits, design$w, design$sebs)
  ek <- MAIVE:::maive_fit_ek(selection, design, sighats, opts$method)
  cfg <- MAIVE:::maive_get_config(opts$method, fits, selection, ek)

  V_iv <- clubSandwich::vcovCR(cfg$maive, cluster = prepared$g, type = opts$type_choice)
  V_ols <- clubSandwich::vcovCR(cfg$std, cluster = prepared$g, type = opts$type_choice)
  var_diff <- V_iv[1, 1] - V_ols[1, 1]
  expected <- (coef(cfg$maive)[1] - coef(cfg$std)[1])^2 / var_diff
  expected_value <- as.numeric(round(expected, 3))

  expect_equal(as.numeric(result$Hausman), expected_value)
})

test_that("Hausman PET-PEESE uses MAIVE weights", {
  dat <- data.frame(
    bs = c(0.42, 0.5, 0.48, 0.55, 0.46, 0.53),
    sebs = c(0.18, 0.2, 0.19, 0.21, 0.22, 0.2),
    Ns = c(80, 95, 90, 85, 88, 92)
  )

  result <- maive(
    dat = dat,
    method = 3,
    weight = 2,
    instrument = 1,
    studylevel = 0,
    SE = 0,
    AR = 0,
    first_stage = 0
  )

  opts <- MAIVE:::maive_validate_inputs(dat, 3, 2, 1, 0, 0, 0, 0)
  prepared <- MAIVE:::maive_prepare_data(opts$dat, opts$studylevel)
  instrumentation <- MAIVE:::maive_compute_variance_instrumentation(
    prepared$sebs,
    prepared$Ns,
    prepared$g,
    opts$type_choice,
    opts$instrument,
    opts$first_stage_type
  )
  w <- MAIVE:::maive_compute_weights(opts$weight, prepared$sebs, instrumentation$sebs2fit1, prepared$studyid)
  x <- if (opts$instrument == 0L) prepared$sebs else sqrt(instrumentation$sebs2fit1)
  x2 <- if (opts$instrument == 0L) prepared$sebs^2 else instrumentation$sebs2fit1
  design <- MAIVE:::maive_build_design_matrices(prepared$bs, prepared$sebs, w, x, x2, prepared$D, prepared$dummy)
  fits <- MAIVE:::maive_fit_models(design)
  selection <- MAIVE:::maive_select_petpeese(fits, design, opts$alpha_s)
  sighats <- MAIVE:::maive_compute_sigma_h(fits, design$w, design$sebs)
  ek <- MAIVE:::maive_fit_ek(selection, design, sighats, opts$method)
  cfg <- MAIVE:::maive_get_config(opts$method, fits, selection, ek)
  hausman_cfg <- MAIVE:::maive_get_hausman_models(opts$method, cfg, selection, design)

  V_iv <- clubSandwich::vcovCR(hausman_cfg$maive, cluster = prepared$g, type = opts$type_choice)
  V_ols <- clubSandwich::vcovCR(hausman_cfg$std, cluster = prepared$g, type = opts$type_choice)
  var_diff <- V_iv[1, 1] - V_ols[1, 1]
  expected <- (coef(hausman_cfg$maive)[1] - coef(hausman_cfg$std)[1])^2 / var_diff
  expected_value <- as.numeric(round(expected, 3))

  expect_equal(as.numeric(result$Hausman), expected_value)
})
