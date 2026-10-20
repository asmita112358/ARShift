# arshiftmethod

Standalone R package extraction of ARShift's fast wrapper method from:

- `/home/runner/work/ARShift/ARShift/main_fun_fast.R`

## What is included

- `test_spatial_association_pmclapply()` (core extracted method)
- Required helper utilities (variance correction, envelope p-values, bandwidth selection)
- Simulation helpers:
  - `simulate_arshift_pair()`
  - `run_arshift_simulation()`
- Unit and integration tests under `tests/testthat`
- Example script at `inst/scripts/example_simulation.R`

## Install locally

```r
# from repository root
install.packages("/home/runner/work/ARShift/ARShift/arshiftmethod", repos = NULL, type = "source")
```

## Quick start

```r
library(arshiftmethod)

X <- simulate_arshift_pair(model = "cauchy", dependant = FALSE, seed = 42)
out <- test_spatial_association_pmclapply(X, n_perm = 39, workers = 1)
out

sim <- run_arshift_simulation(n_sim = 5, test_args = list(n_perm = 39, workers = 1))
attr(sim, "summary")
```

## Migration note

- Old direct wrapper usage:
  - `source(".../main_fun_fast.R")`
  - `test_spatial_association_pmclapply(...)`
- New package usage:
  - `library(arshiftmethod)`
  - `test_spatial_association_pmclapply(...)`
