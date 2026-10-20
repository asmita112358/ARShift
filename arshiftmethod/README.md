# ARShift

Standalone R package extraction of the `pval_Kcross_area` method path from:

- `/home/runner/work/ARShift/ARShift/main_fun_fast.R::test_spatial_association_pmclapply`

## Core extracted API

- `ARShift()`
- `arshift_pval_Kcross_area()`

Both compute the area-corrected ARShift p-value (`pval_Kcross_area`) and return permutation diagnostics.

## Included helper chain

The package includes all helper logic needed by the extracted method:

- point-process input validation
- K-cross computation helpers
- reduced-window shift permutations
- area-scaling transform (`T_0`, `T_scaled`)
- global envelope p-value computation

## Install locally

```r
install.packages("/home/runner/work/ARShift/ARShift/arshiftmethod", repos = NULL, type = "source")
```

## Quick start

```r
library(ARShift)

X <- simulate_arshift_pair(model = "cauchy", dependant = FALSE, seed = 42)
out <- ARShift(X, n_perm = 39, workers = 1)
out$pval_Kcross_area
```

## Migration note

- Old direct wrapper usage:
  - `source(".../main_fun_fast.R")`
  - `test_spatial_association_pmclapply(... )$pval_Kcross_area`
- New package usage:
  - `library(ARShift)`
  - `ARShift(... )$pval_Kcross_area`
