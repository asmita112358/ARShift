test_that("core test returns expected structure", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.random")

  set.seed(11)
  data <- simulate_arshift_pair(model = "thomas", dependant = FALSE, seed = 11)

  res <- test_spatial_association_pmclapply(
    data = data,
    n_perm = 9,
    workers = 1,
    include_RL = FALSE,
    type = "inhom"
  )

  expect_named(
    res,
    c(
      "pval_Kcross_tor", "pval_Kcross_vc", "pval_Kcross_uncorr",
      "pval_Kcross_area", "pval_RL", "n_valid", "n_valid_total",
      "n_attempts", "k"
    )
  )
  expect_true(res$n_attempts >= 9)
})

test_that("missing taxa returns NA p-values", {
  skip_if_not_installed("spatstat.geom")

  win <- spatstat.geom::square(1)
  x <- runif(20)
  y <- runif(20)
  marks <- factor(rep(1, 20), levels = c(1, 2))
  data <- spatstat.geom::ppp(x, y, window = win, marks = marks)

  res <- test_spatial_association_pmclapply(data = data, n_perm = 9, workers = 1)

  expect_true(is.na(res$pval_Kcross_tor))
  expect_equal(res$n_valid, 0)
})

test_that("invalid inputs are rejected", {
  expect_error(
    test_spatial_association_pmclapply(data = mtcars),
    "must be a spatstat point pattern"
  )
})
