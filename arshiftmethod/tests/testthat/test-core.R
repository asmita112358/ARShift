test_that("ARShift returns area p-value structure", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.random")

  set.seed(11)
  data <- simulate_arshift_pair(model = "thomas", dependant = FALSE, seed = 11)

  res <- ARShift(
    data = data,
    n_perm = 9,
    workers = 1,
    type = "inhom"
  )

  expect_named(
    res,
    c(
      "pval_Kcross_area", "n_valid", "n_valid_total",
      "n_attempts", "k", "r"
    )
  )
  expect_true(is.numeric(res$pval_Kcross_area) || is.na(res$pval_Kcross_area))
  expect_true(res$n_attempts >= 9)
})

test_that("missing taxa returns NA area p-value", {
  skip_if_not_installed("spatstat.geom")

  win <- spatstat.geom::square(1)
  x <- runif(20)
  y <- runif(20)
  marks <- factor(rep(1, 20), levels = c(1, 2))
  data <- spatstat.geom::ppp(x, y, window = win, marks = marks)

  res <- ARShift(data = data, n_perm = 9, workers = 1)

  expect_true(is.na(res$pval_Kcross_area))
  expect_equal(res$n_valid, 0)
})

test_that("invalid inputs are rejected", {
  expect_error(
    ARShift(data = mtcars),
    "must be a spatstat point pattern"
  )
})

test_that("area-only extraction matches wrapper area p-value", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.random")

  set.seed(23)
  data <- simulate_arshift_pair(model = "cauchy", dependant = FALSE, seed = 23)

  set.seed(101)
  area_only <- ARShift(data = data, n_perm = 9, workers = 1, type = "inhom")
  set.seed(101)
  wrapper <- test_spatial_association_pmclapply(
    data = data,
    n_perm = 9,
    workers = 1,
    include_RL = FALSE,
    type = "inhom"
  )

  expect_equal(area_only$pval_Kcross_area, wrapper$pval_Kcross_area, tolerance = 1e-10)
})
