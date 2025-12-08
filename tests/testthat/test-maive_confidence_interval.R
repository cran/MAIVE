test_that("maive_prepare_confidence_interval handles degenerate bootstrap ci", {
  model <- stats::lm(effect ~ 1, data = data.frame(effect = 1:3))

  base_args <- list(
    model = model,
    coef_index = 1,
    estimate = unname(coef(model)[1]),
    se = 1,
    alpha = 0.05
  )

  expect_silent({
    ci <- do.call(
      MAIVE:::maive_prepare_confidence_interval,
      c(base_args, list(boot_result = list(boot_ci = matrix(NA_real_, nrow = 1, ncol = 1))))
    )
    expect_identical(names(ci), c("lower", "upper"))
    expect_true(all(is.na(ci)))
  })

  expect_silent({
    ci <- do.call(
      MAIVE:::maive_prepare_confidence_interval,
      c(base_args, list(boot_result = list(boot_ci = matrix(0.5, nrow = 1, ncol = 1))))
    )
    expect_identical(names(ci), c("lower", "upper"))
    expect_equal(ci, setNames(rep(0.5, 2), c("lower", "upper")))
  })

  expect_silent({
    boot_ci <- matrix(c(0.2, 0.8), nrow = 1)
    dimnames(boot_ci) <- list(NULL, NULL)
    ci <- do.call(
      MAIVE:::maive_prepare_confidence_interval,
      c(base_args, list(boot_result = list(boot_ci = boot_ci)))
    )
    expect_identical(names(ci), c("lower", "upper"))
    expect_equal(ci, setNames(c(0.2, 0.8), c("lower", "upper")))
  })
})
