test_that("Study weights give each study equal total weight", {
  dat <- data.frame(
    bs = c(0.4, 0.5, 0.45, 0.55, 0.47, 0.52),
    sebs = c(0.2, 0.19, 0.21, 0.18, 0.22, 0.2),
    Ns = c(80, 85, 90, 95, 88, 92),
    study_id = c("A", "A", "B", "B", "B", "C")
  )

  opts <- MAIVE:::maive_validate_inputs(dat, method = 1, weight = 3, instrument = 1, studylevel = 2, SE = 0, AR = 0, first_stage = 0)
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

  expect_equal(w, c(rep(1 / 2, 2), rep(1 / 3, 3), 1))

  totals <- tapply(w, prepared$studyid, sum)
  expect_true(all(abs(totals - 1) < 1e-12))
})
