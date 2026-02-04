test_that("validate_maive_data rejects non-data frame input", {
  expect_error(
    MAIVE:::validate_maive_data(list(a = 1, b = 2)),
    "must be a data frame"
  )
})

test_that("validate_maive_data rejects < 4 rows", {
  dat <- data.frame(bs = c(1, 2, 3), sebs = c(0.1, 0.2, 0.3), Ns = c(100, 200, 300))
  expect_error(
    MAIVE:::validate_maive_data(dat),
    "at least 4 observations"
  )
})

test_that("validate_maive_data accepts exactly 4 rows", {
  dat <- data.frame(bs = c(1, 2, 3, 4), sebs = c(0.1, 0.2, 0.3, 0.4), Ns = c(100, 200, 300, 400))
  expect_silent(MAIVE:::validate_maive_data(dat))
})

test_that("validate_maive_data requires bs column", {
  dat <- data.frame(sebs = c(0.1, 0.2, 0.3, 0.4), Ns = c(100, 200, 300, 400))
  expect_error(
    MAIVE:::validate_maive_data(dat),
    "Missing required columns.*bs"
  )
})

test_that("validate_maive_data requires sebs column", {
  dat <- data.frame(bs = c(1, 2, 3, 4), Ns = c(100, 200, 300, 400))
  expect_error(
    MAIVE:::validate_maive_data(dat),
    "Missing required columns.*sebs"
  )
})

test_that("validate_maive_data requires Ns column", {
  dat <- data.frame(bs = c(1, 2, 3, 4), sebs = c(0.1, 0.2, 0.3, 0.4))
  expect_error(
    MAIVE:::validate_maive_data(dat),
    "Missing required columns.*Ns"
  )
})

test_that("validate_maive_data requires numeric bs", {
  dat <- data.frame(
    bs = c("1", "2", "3", "4"),
    sebs = c(0.1, 0.2, 0.3, 0.4),
    Ns = c(100, 200, 300, 400)
  )
  expect_error(
    MAIVE:::validate_maive_data(dat),
    "Column 'bs' must be numeric"
  )
})

test_that("validate_maive_data requires numeric sebs", {
  dat <- data.frame(
    bs = c(1, 2, 3, 4),
    sebs = c("0.1", "0.2", "0.3", "0.4"),
    Ns = c(100, 200, 300, 400)
  )
  expect_error(
    MAIVE:::validate_maive_data(dat),
    "Column 'sebs' must be numeric"
  )
})

test_that("validate_maive_data requires numeric Ns", {
  dat <- data.frame(
    bs = c(1, 2, 3, 4),
    sebs = c(0.1, 0.2, 0.3, 0.4),
    Ns = c("100", "200", "300", "400")
  )
  expect_error(
    MAIVE:::validate_maive_data(dat),
    "Column 'Ns' must be numeric"
  )
})

test_that("validate_maive_data removes empty rows", {
  dat <- data.frame(
    bs = c(1, NA, 3, 4, 5),
    sebs = c(0.1, NA, 0.3, 0.4, 0.5),
    Ns = c(100, NA, 300, 400, 500)
  )
  expect_message(MAIVE:::validate_maive_data(dat), "Removed 1 completely empty row")
})

test_that("validate_maive_data fails if too few rows after cleaning", {
  dat <- data.frame(
    bs = c(1, NA, 3),
    sebs = c(0.1, NA, 0.3),
    Ns = c(100, NA, 300)
  )
  # Remove empty rows (row 2 fully empty) leaves 2 rows, triggering size check
  expect_error(
    MAIVE:::validate_maive_data(dat),
    "insufficient data remains"
  )
})

test_that("validate_maive_data rejects non-positive standard errors", {
  dat <- data.frame(bs = 1:4, sebs = c(0.1, 0, 0.3, 0.4), Ns = c(100, 200, 300, 400))
  expect_error(
    MAIVE:::validate_maive_data(dat),
    "must be positive.*Found 1 non-positive value"
  )
})

test_that("validate_maive_data rejects multiple non-positive standard errors", {
  dat <- data.frame(bs = 1:4, sebs = c(0, -0.1, 0.3, 0.4), Ns = c(100, 200, 300, 400))
  expect_error(
    MAIVE:::validate_maive_data(dat),
    "must be positive.*Found 2 non-positive values"
  )
})

test_that("validate_maive_data rejects NA standard errors", {
  dat <- data.frame(bs = 1:4, sebs = c(0.1, NA, 0.3, 0.4), Ns = c(100, 200, 300, 400))
  expect_error(
    MAIVE:::validate_maive_data(dat),
    "must not contain missing values"
  )
})

test_that("validate_maive_data rejects non-positive sample sizes", {
  dat <- data.frame(bs = 1:4, sebs = c(0.1, 0.2, 0.3, 0.4), Ns = c(100, 0, 300, 400))
  expect_error(
    MAIVE:::validate_maive_data(dat),
    "Sample sizes.*must be positive.*Found 1 non-positive value"
  )
})

test_that("validate_maive_data checks study_id degrees of freedom", {
  dat <- data.frame(
    bs = 1:5,
    sebs = rep(0.1, 5),
    Ns = rep(100, 5),
    study_id = c(1, 1, 2, 3, 4) # 4 unique studies, need 7+ rows
  )
  expect_error(
    MAIVE:::validate_maive_data(dat),
    "Insufficient degrees of freedom.*requires at least 7 rows"
  )
})

test_that("validate_maive_data passes with sufficient study_id rows", {
  dat <- data.frame(
    bs = 1:8,
    sebs = rep(0.1, 8),
    Ns = rep(100, 8),
    study_id = c(1, 1, 2, 2, 3, 3, 4, 4) # 4 unique studies, 8 rows (>= 4+3)
  )
  expect_silent(MAIVE:::validate_maive_data(dat))
})

test_that("validate_maive_parameters rejects invalid method", {
  expect_error(
    MAIVE:::validate_maive_parameters(
      method = 99, weight = 1, instrument = 1,
      studylevel = 0, SE = 0, AR = 0
    ),
    "method.*must be.*PET.*PEESE.*EK"
  )
})

test_that("validate_maive_parameters accepts valid method values", {
  for (method in 1:4) {
    expect_silent(
      MAIVE:::validate_maive_parameters(
        method = method, weight = 1, instrument = 1,
        studylevel = 0, SE = 0, AR = 0
      )
    )
  }
})

test_that("validate_maive_parameters rejects invalid weight", {
  expect_error(
    MAIVE:::validate_maive_parameters(
      method = 1, weight = 99, instrument = 1,
      studylevel = 0, SE = 0, AR = 0
    ),
    "weight.*must be.*equal.*standard.*adjusted.*study"
  )
})

test_that("validate_maive_parameters accepts valid weight values", {
  for (weight in 0:3) {
    expect_silent(
      MAIVE:::validate_maive_parameters(
        method = 1, weight = weight, instrument = 1,
        studylevel = 0, SE = 0, AR = 0
      )
    )
  }
})

test_that("validate_maive_parameters rejects invalid instrument", {
  expect_error(
    MAIVE:::validate_maive_parameters(
      method = 1, weight = 1, instrument = 2,
      studylevel = 0, SE = 0, AR = 0
    ),
    "instrument.*must be 0 or 1"
  )
})

test_that("validate_maive_parameters rejects invalid studylevel", {
  expect_error(
    MAIVE:::validate_maive_parameters(
      method = 1, weight = 1, instrument = 1,
      studylevel = 99, SE = 0, AR = 0
    ),
    "studylevel.*must be.*0.*1.*2.*3"
  )
})

test_that("validate_maive_parameters rejects invalid SE", {
  expect_error(
    MAIVE:::validate_maive_parameters(
      method = 1, weight = 1, instrument = 1,
      studylevel = 0, SE = 99, AR = 0
    ),
    "SE.*must be.*0.*1.*2.*3"
  )
})

test_that("validate_maive_parameters rejects invalid AR", {
  expect_error(
    MAIVE:::validate_maive_parameters(
      method = 1, weight = 1, instrument = 1,
      studylevel = 0, SE = 0, AR = 2
    ),
    "AR.*must be 0 or 1"
  )
})

test_that("maive runs with valid data", {
  set.seed(123)
  dat <- data.frame(
    bs = rnorm(20, 0.3, 0.2),
    sebs = runif(20, 0.05, 0.15),
    Ns = round(runif(20, 100, 500))
  )

  result <- maive(dat,
    method = 3, weight = 1, instrument = 1,
    studylevel = 0, SE = 0, AR = 0
  )
  expect_true(!is.null(result))
  expect_true(is.list(result))
})

test_that("maive fails with insufficient data", {
  dat <- data.frame(bs = c(1, 2, 3), sebs = c(0.1, 0.2, 0.3), Ns = c(100, 200, 300))
  expect_error(
    maive(dat,
      method = 3, weight = 1, instrument = 1,
      studylevel = 0, SE = 0, AR = 0
    ),
    "at least 4 observations"
  )
})

test_that("waive runs with valid data", {
  set.seed(456)
  dat <- data.frame(
    bs = rnorm(20, 0.3, 0.2),
    sebs = runif(20, 0.05, 0.15),
    Ns = round(runif(20, 100, 500))
  )

  result <- waive(dat,
    method = 3, weight = 0, instrument = 1,
    studylevel = 0, SE = 0, AR = 0
  )
  expect_true(!is.null(result))
  expect_true(is.list(result))
})

test_that("waive fails with insufficient data", {
  dat <- data.frame(bs = c(1, 2), sebs = c(0.1, 0.2), Ns = c(100, 200))
  expect_error(
    waive(dat,
      method = 3, weight = 0, instrument = 1,
      studylevel = 0, SE = 0, AR = 0
    ),
    "at least 4 observations"
  )
})
