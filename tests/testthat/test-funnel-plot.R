test_that("get_funnel_plot draws to a graphics device (smoke test)", {
  dat <- data.frame(
    bs = c(0.5, 0.45, 0.55, 0.6),
    sebs = c(0.25, 0.2, 0.22, 0.27),
    Ns = c(50, 80, 65, 90)
  )

  result <- maive(dat,
    method = 3, weight = 0, instrument = 1,
    studylevel = 0, SE = 0, AR = 0, first_stage = 0
  )

  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp, width = 600, height = 600)

  expect_invisible(
    get_funnel_plot(dat = dat, result = result, instrument = 1, model_type = "MAIVE")
  )

  grDevices::dev.off()

  expect_true(file.exists(tmp))
  expect_true(file.info(tmp)$size > 0)
})

test_that("get_funnel_plot validates required inputs", {
  dat_ok <- data.frame(bs = c(0.1, 0.2, 0.15, 0.12), sebs = c(0.1, 0.2, 0.18, 0.16), Ns = c(10, 20, 15, 18))
  result_ok <- maive(dat_ok,
    method = 1, weight = 0, instrument = 0,
    studylevel = 0, SE = 0, AR = 0, first_stage = 0
  )

  expect_error(
    get_funnel_plot(dat = data.frame(bs = c(1, 2)), result = result_ok),
    "bs.*sebs"
  )
  expect_error(
    get_funnel_plot(dat = dat_ok, result = list(beta = 1)),
    "beta.*SE"
  )
  expect_error(
    get_funnel_plot(dat = data.frame(bs = c(0.1, 0.2), sebs = c(0, 0.2)), result = result_ok),
    "strictly positive"
  )
})
