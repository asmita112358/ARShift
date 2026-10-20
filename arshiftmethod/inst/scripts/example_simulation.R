library(arshiftmethod)

set.seed(2026)

X <- simulate_arshift_pair(model = "cauchy", dependant = FALSE)
res <- test_spatial_association_pmclapply(
  data = X,
  n_perm = 39,
  workers = 1,
  type = "inhom"
)
print(res)

sim_res <- run_arshift_simulation(
  n_sim = 10,
  seed = 2026,
  test_args = list(n_perm = 39, workers = 1, type = "inhom")
)
print(attr(sim_res, "summary"))
