## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
# # Install from CRAN (once published)
# install.packages("MAIVE")
# 
# # Or install development version from GitHub
# install.packages("devtools")
# devtools::install_github("meta-analysis-es/maive")

## ----setup--------------------------------------------------------------------
library(MAIVE)

## ----example-data-------------------------------------------------------------
# Simulated meta-analysis data
set.seed(123)
n_studies <- 50

data <- data.frame(
  bs = rnorm(n_studies, mean = 0.3, sd = 0.2),
  sebs = runif(n_studies, min = 0.05, max = 0.3),
  Ns = sample(100:1000, n_studies, replace = TRUE),
  study_id = rep(1:10, each = 5)
)

head(data)

## ----default-maive, eval = FALSE----------------------------------------------
# # Run MAIVE with defaults
# result <- maive(
#   dat = data,
#   method = 3,      # PET-PEESE (default)
#   weight = 0,      # No weights (default)
#   instrument = 1,  # Instrument SEs (default)
#   studylevel = 2,  # Cluster-robust (default)
#   SE = 3,          # Wild bootstrap (default)
#   AR = 1           # Anderson-Rubin CI (default)
# )
# 
# # View key results
# cat("MAIVE Estimate:", round(result$Estimate, 3), "\n")
# cat("MAIVE SE:", round(result$SE, 3), "\n")
# cat("Standard Estimate:", round(result$StdEstimate, 3), "\n")
# cat("Hausman Test:", round(result$Hausman, 3), "\n")
# cat("First-stage F-test:", round(result$`F-test`, 3), "\n")

## ----pet, eval = FALSE--------------------------------------------------------
# result_pet <- maive(
#   dat = data,
#   method = 1,  # FAT-PET
#   weight = 0,
#   instrument = 1,
#   studylevel = 2,
#   SE = 3,
#   AR = 1
# )
# 
# cat("PET Estimate:", round(result_pet$Estimate, 3), "\n")

## ----peese, eval = FALSE------------------------------------------------------
# result_peese <- maive(
#   dat = data,
#   method = 2,  # PEESE
#   weight = 0,
#   instrument = 1,
#   studylevel = 2,
#   SE = 3,
#   AR = 1
# )
# 
# cat("PEESE Estimate:", round(result_peese$Estimate, 3), "\n")

## ----petpeese, eval = FALSE---------------------------------------------------
# result_petpeese <- maive(
#   dat = data,
#   method = 3,  # PET-PEESE (default)
#   weight = 0,
#   instrument = 1,
#   studylevel = 2,
#   SE = 3,
#   AR = 1
# )
# 
# cat("PET-PEESE Estimate:", round(result_petpeese$Estimate, 3), "\n")

## ----ek, eval = FALSE---------------------------------------------------------
# result_ek <- maive(
#   dat = data,
#   method = 4,  # EK
#   weight = 0,
#   instrument = 1,
#   studylevel = 2,
#   SE = 3,
#   AR = 0  # AR not available for EK
# )
# 
# cat("EK Estimate:", round(result_ek$Estimate, 3), "\n")

## ----no-weights, eval = FALSE-------------------------------------------------
# result_noweight <- maive(
#   dat = data,
#   method = 3,
#   weight = 0,  # No weights
#   instrument = 1,
#   studylevel = 2,
#   SE = 3,
#   AR = 1
# )

## ----iv-weights, eval = FALSE-------------------------------------------------
# result_ivweight <- maive(
#   dat = data,
#   method = 3,
#   weight = 1,  # Inverse-variance weights
#   instrument = 1,
#   studylevel = 2,
#   SE = 3,
#   AR = 0  # AR not available with weights
# )

## ----maive-weights, eval = FALSE----------------------------------------------
# result_maiveweight <- maive(
#   dat = data,
#   method = 3,
#   weight = 2,  # MAIVE-adjusted weights
#   instrument = 1,
#   studylevel = 2,
#   SE = 3,
#   AR = 1
# )

## ----studylevel, eval = FALSE-------------------------------------------------
# # No study-level adjustment
# result_none <- maive(data, method = 3, weight = 0, instrument = 1,
#                      studylevel = 0, SE = 0, AR = 1)
# 
# # Study fixed effects (demeaned)
# result_fe <- maive(data, method = 3, weight = 0, instrument = 1,
#                    studylevel = 1, SE = 1, AR = 0)  # AR not available with FE
# 
# # Cluster-robust standard errors
# result_cluster <- maive(data, method = 3, weight = 0, instrument = 1,
#                         studylevel = 2, SE = 3, AR = 1)
# 
# # Both fixed effects and clustering
# result_both <- maive(data, method = 3, weight = 0, instrument = 1,
#                      studylevel = 3, SE = 3, AR = 0)

## ----se-options, eval = FALSE-------------------------------------------------
# # CR0 (Huber-White)
# result_cr0 <- maive(data, method = 3, weight = 0, instrument = 1,
#                     studylevel = 2, SE = 0, AR = 1)
# 
# # CR1 (Standard empirical correction)
# result_cr1 <- maive(data, method = 3, weight = 0, instrument = 1,
#                     studylevel = 2, SE = 1, AR = 1)
# 
# # CR2 (Bias-reduced estimator)
# result_cr2 <- maive(data, method = 3, weight = 0, instrument = 1,
#                     studylevel = 2, SE = 2, AR = 1)
# 
# # Wild bootstrap (recommended, default)
# result_boot <- maive(data, method = 3, weight = 0, instrument = 1,
#                      studylevel = 2, SE = 3, AR = 1)

## ----first-stage-levels, eval = FALSE-----------------------------------------
# result_levels <- maive(data, method = 3, weight = 0, instrument = 1,
#                        studylevel = 2, SE = 3, AR = 1, first_stage = 0)
# 
# cat("First-stage (levels) F-test:", round(result_levels$`F-test`, 3), "\n")

## ----first-stage-log, eval = FALSE--------------------------------------------
# result_log <- maive(data, method = 3, weight = 0, instrument = 1,
#                     studylevel = 2, SE = 3, AR = 1, first_stage = 1)
# 
# cat("First-stage (log) F-test:", round(result_log$`F-test`, 3), "\n")

## ----waive, eval = FALSE------------------------------------------------------
# result_waive <- waive(
#   dat = data,
#   method = 3,
#   instrument = 1,
#   studylevel = 2,
#   SE = 3,
#   AR = 1
# )
# 
# cat("WAIVE Estimate:", round(result_waive$Estimate, 3), "\n")
# cat("WAIVE SE:", round(result_waive$SE, 3), "\n")

