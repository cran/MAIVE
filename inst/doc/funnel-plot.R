## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----synth-data, eval = FALSE-------------------------------------------------
# library(MAIVE)
# 
# set.seed(123)
# n <- 40
# dat <- data.frame(
#   bs = rnorm(n, mean = 0.2, sd = 0.25),
#   sebs = runif(n, min = 0.05, max = 0.3),
#   Ns = sample(80:800, n, replace = TRUE),
#   study_id = rep(1:10, each = 4)
# )
# 
# # Fit MAIVE (instrumented PET-PEESE, no weights, cluster-robust, wild bootstrap)
# res <- maive(
#   dat = dat,
#   method = 3,
#   weight = 0,
#   instrument = 1,
#   studylevel = 2,
#   SE = 3,
#   AR = 1
# )
# 
# # Draw to the current device
# get_funnel_plot(dat = dat, result = res, model_type = "MAIVE")

## ----save-png, eval = FALSE---------------------------------------------------
# png("maive-funnel.png", width = 1800, height = 1400, res = 200)
# get_funnel_plot(dat = dat, result = res, model_type = "MAIVE")
# dev.off()

