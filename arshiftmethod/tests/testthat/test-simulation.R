test_that("simulation generator returns marked point pattern", {
  skip_if_not_installed("spatstat.geom")

  set.seed(123)
  X <- simulate_arshift_pair(model = "cauchy", dependant = TRUE, seed = 123)

  expect_s3_class(X, "ppp")
  expect_true(length(levels(spatstat.geom::marks(X))) == 2)
})

test_that("simulation runner produces expected output", {
  skip_if_not_installed("spatstat.geom")

  res <- run_arshift_simulation(
    n_sim = 3,
    seed = 1,
    test_args = list(n_perm = 9, workers = 1, type = "inhom")
  )

  expect_equal(nrow(res), 3)
  expect_true(all(c("reject_tor", "reject_vc", "reject_area", "reject_uncorr") %in% names(res)))
  expect_true(!is.null(attr(res, "summary")))
})
