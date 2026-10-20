make_test_ppp <- function(seed = 1, win = spatstat.geom::square(1)) {
  set.seed(seed)
  a <- spatstat.random::rThomas(kappa = 8, scale = 0.05, mu = 15, win = win)
  b <- spatstat.random::rThomas(kappa = 7, scale = 0.04, mu = 14, win = win)
  spatstat.geom::marks(a) <- factor(rep(1, spatstat.geom::npoints(a)), levels = c(1, 2))
  spatstat.geom::marks(b) <- factor(rep(2, spatstat.geom::npoints(b)), levels = c(1, 2))
  out <- spatstat.geom::superimpose(a, b, W = win)
  spatstat.geom::marks(out) <- factor(spatstat.geom::marks(out), levels = c(1, 2))
  out
}

test_that("ARShift returns expected area-output structure", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.random")

  data <- make_test_ppp(seed = 11)

  res <- ARShift(
    data = data,
    n_perm = 9,
    workers = 1,
    type = "inhom"
  )

  expect_named(res, c("pval_Kcross_area", "n_valid", "n_valid_total", "n_attempts", "k", "r"))
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
