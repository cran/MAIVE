test_that("maive reports Egger bootstrap and AR intervals when eligible", {
  dat <- data.frame(
    bs = c(0.4, 0.6, 0.55, 0.5, 0.52, 0.49),
    sebs = c(0.2, 0.18, 0.22, 0.19, 0.21, 0.2),
    Ns = c(80, 95, 90, 85, 88, 92),
    study_id = c(1, 1, 2, 2, 3, 3)
  )

  res <- maive(
    dat = dat,
    method = 1,
    weight = 0,
    instrument = 1,
    studylevel = 2,
    SE = 0,
    AR = 1
  )

  expect_true(is.numeric(res$egger_boot_ci))
  expect_length(res$egger_boot_ci, 2)
  expect_equal(sort(names(res$egger_boot_ci)), c("lower", "upper"))

  expect_true(is.numeric(res$egger_ar_ci))
  expect_length(res$egger_ar_ci, 2)
  expect_equal(sort(names(res$egger_ar_ci)), c("lower", "upper"))
})

test_that("maive returns NA Egger AR interval when AR is disabled", {
  dat <- data.frame(
    bs = c(0.4, 0.6, 0.55, 0.5, 0.52, 0.49),
    sebs = c(0.2, 0.18, 0.22, 0.19, 0.21, 0.2),
    Ns = c(80, 95, 90, 85, 88, 92),
    study_id = c(1, 1, 2, 2, 3, 3)
  )

  res <- maive(
    dat = dat,
    method = 1,
    weight = 0,
    instrument = 1,
    studylevel = 2,
    SE = 0,
    AR = 0
  )

  expect_true(is.numeric(res$egger_boot_ci))
  expect_length(res$egger_boot_ci, 2)
  expect_equal(sort(names(res$egger_boot_ci)), c("lower", "upper"))

  expect_identical(res$egger_ar_ci, "NA")
})

test_that("maive returns NA Egger AR interval when AR grid has no admissible points", {
  path <- testthat::test_path("fixtures", "egger_ar_no_acceptance.csv")
  dat <- read.csv(path, stringsAsFactors = FALSE)
  colnames(dat) <- c("bs", "sebs", "Ns", "study_id")

  res <- suppressWarnings(maive(
    dat = dat,
    method = 3,
    weight = 0,
    instrument = 1,
    studylevel = 3,
    SE = 0,
    AR = 1
  ))

  # When AR grid has no admissible points, function returns "NA" (string)
  expect_identical(res$egger_ar_ci, "NA")
})
