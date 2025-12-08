test_that("wild cluster bootstrap uses bootstrap SEs for t statistics", {
  skip_if_not_installed("clubSandwich")

  set.seed(123)
  n_clusters <- 8
  cluster <- rep(seq_len(n_clusters), each = 5)
  x <- rnorm(length(cluster))
  cluster_effect <- rnorm(n_clusters, sd = 0.3)
  y <- 1 + 0.5 * x + cluster_effect[cluster] + rnorm(length(cluster), sd = 0.5)

  data <- data.frame(y = y, x = x, cluster = cluster)
  model <- lm(y ~ x, data = data)

  boot <- manual_wild_cluster_boot_se(
    model = model,
    data = data,
    cluster_var = "cluster",
    B = 50,
    seed = 2024
  )

  expect_true(!is.null(boot$boot_rep_se))
  expect_equal(dim(boot$boot_rep_se), dim(boot$boot_coefs))

  base_coefs <- coef(model)
  for (j in seq_along(base_coefs)) {
    denom <- boot$boot_rep_se[, j]
    numer <- boot$boot_coefs[, j] - base_coefs[j]
    manual_t <- rep(NA_real_, length(numer))
    ok <- !is.na(denom) & denom != 0
    manual_t[ok] <- numer[ok] / denom[ok]
    expect_equal(boot$boot_t_stats[, j], manual_t)
  }
})
