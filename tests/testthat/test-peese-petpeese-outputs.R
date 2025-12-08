run_maive <- function(data, method) {
  maive(
    dat = data,
    method = method,
    weight = 0,
    instrument = 1,
    studylevel = 0,
    SE = 0,
    AR = 0
  )
}

test_that("method=2 (PEESE) returns PEESE SE^2 coefficient and SE", {
  test_dat <- data.frame(
    bs = c(2.218944, 3.470763, 2.522442, 3.707544, 3.851168, 1.613891),
    sebs = c(0.2584316, 0.3677257, 0.2654305, 0.2369844, 0.38705, 0.2360002),
    Ns = c(74, 118, 106, 58, 75, 56)
  )

  res <- run_maive(test_dat, method = 2)

  expect_true(is.numeric(res$peese_se2_coef))
  expect_false(is.na(res$peese_se2_coef))
  expect_true(is.numeric(res$peese_se2_se))
  expect_false(is.na(res$peese_se2_se))
  expect_true(is.na(res$petpeese_selected))
})

test_that("method=1 (PET) returns NA for PEESE coefficients", {
  test_dat <- data.frame(
    bs = c(0.4, 0.45, 0.5, 0.55, 0.52),
    sebs = c(0.2, 0.22, 0.21, 0.24, 0.23),
    Ns = c(80, 90, 85, 95, 100)
  )

  res <- run_maive(test_dat, method = 1)

  expect_true(is.na(res$peese_se2_coef))
  expect_true(is.na(res$peese_se2_se))
  expect_true(is.na(res$petpeese_selected))
})

test_that("method=3 (PET-PEESE) returns selection indicator", {
  test_dat <- data.frame(
    bs = c(2.218944, 3.470763, 2.522442, 3.707544, 3.851168, 1.613891),
    sebs = c(0.2584316, 0.3677257, 0.2654305, 0.2369844, 0.38705, 0.2360002),
    Ns = c(74, 118, 106, 58, 75, 56)
  )

  res <- run_maive(test_dat, method = 3)

  expect_true(res$petpeese_selected %in% c("PET", "PEESE"))
  expect_type(res$petpeese_selected, "character")
})

test_that("method=3 returns PEESE coefficients only when PEESE is selected", {
  test_dat <- data.frame(
    bs = c(2.218944, 3.470763, 2.522442, 3.707544, 3.851168, 1.613891),
    sebs = c(0.2584316, 0.3677257, 0.2654305, 0.2369844, 0.38705, 0.2360002),
    Ns = c(74, 118, 106, 58, 75, 56)
  )

  res <- run_maive(test_dat, method = 3)

  if (res$petpeese_selected == "PEESE") {
    expect_true(is.numeric(res$peese_se2_coef))
    expect_false(is.na(res$peese_se2_coef))
    expect_true(is.numeric(res$peese_se2_se))
    expect_false(is.na(res$peese_se2_se))
  } else {
    expect_true(is.na(res$peese_se2_coef))
    expect_true(is.na(res$peese_se2_se))
  }
})

test_that("method=4 (EK) returns NA for PEESE coefficients", {
  test_dat <- data.frame(
    bs = c(0.4, 0.45, 0.5, 0.55, 0.52),
    sebs = c(0.2, 0.22, 0.21, 0.24, 0.23),
    Ns = c(80, 90, 85, 95, 100)
  )

  res <- suppressWarnings(run_maive(test_dat, method = 4))

  expect_true(is.na(res$peese_se2_coef))
  expect_true(is.na(res$peese_se2_se))
  expect_true(is.na(res$petpeese_selected))
})

