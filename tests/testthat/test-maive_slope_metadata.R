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

test_that("PET slope metadata remains linear", {
  pet_dat <- data.frame(
    bs = c(0.4, 0.45, 0.5, 0.55, 0.52),
    sebs = c(0.2, 0.22, 0.21, 0.24, 0.23),
    Ns = c(80, 90, 85, 95, 100)
  )

  pet_res <- run_maive(pet_dat, method = 1)

  expect_equal(pet_res$is_quadratic_fit$slope_type, "linear")
  expect_true(is.numeric(pet_res$slope_coef))
})

test_that("EK without kink stays linear", {
  ek_linear_dat <- data.frame(
    bs = c(0.1, 0.12, 0.11, 0.09, 0.13),
    sebs = c(0.4, 0.42, 0.41, 0.39, 0.43),
    Ns = c(150, 160, 170, 180, 190)
  )

  ek_linear_res <- suppressWarnings(run_maive(ek_linear_dat, method = 4))

  expect_equal(ek_linear_res$is_quadratic_fit$slope_type, "linear")
  expect_null(ek_linear_res$is_quadratic_fit$slope_detail)
  expect_true(is.numeric(ek_linear_res$slope_coef))
})

test_that("EK with kink exposes kink metadata", {
  ek_kink_dat <- data.frame(
    bs = c(2.218944, 3.470763, 2.522442, 3.707544, 3.851168, 1.613891),
    sebs = c(0.2584316, 0.3677257, 0.2654305, 0.2369844, 0.38705, 0.2360002),
    Ns = c(74, 118, 106, 58, 75, 56)
  )

  ek_kink_res <- suppressWarnings(run_maive(ek_kink_dat, method = 4))

  expect_equal(ek_kink_res$is_quadratic_fit$slope_type, "kinked")
  expect_true(is.list(ek_kink_res$slope_coef))
  expect_true(all(c("kink_effect", "kink_location") %in% names(ek_kink_res$slope_coef)))
  expect_false(is.null(ek_kink_res$is_quadratic_fit$slope_detail))
})
