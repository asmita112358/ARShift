# ARShift

This package extracts only the area-corrected ARShift method from:

- `/home/runner/work/ARShift/ARShift/main_fun_fast.R::test_spatial_association_pmclapply`

## Public function

- `ARShift()`

`ARShift()` computes `pval_Kcross_area` and returns permutation diagnostics.

## Install locally

```r
install.packages("/home/runner/work/ARShift/ARShift/arshiftmethod", repos = NULL, type = "source")
```

## Quick start

```r
library(ARShift)

win <- spatstat.geom::square(1)
a <- spatstat.random::rThomas(8, 0.05, 15, win)
b <- spatstat.random::rThomas(7, 0.04, 14, win)
spatstat.geom::marks(a) <- factor(rep(1, spatstat.geom::npoints(a)), levels = c(1, 2))
spatstat.geom::marks(b) <- factor(rep(2, spatstat.geom::npoints(b)), levels = c(1, 2))
X <- spatstat.geom::superimpose(a, b, W = win)
spatstat.geom::marks(X) <- factor(spatstat.geom::marks(X), levels = c(1, 2))

out <- ARShift(X, n_perm = 39, workers = 1)
out$pval_Kcross_area
```

## Note on original code

The original repository scripts and wrapper files remain untouched. This extraction lives in a separate package directory: `/home/runner/work/ARShift/ARShift/arshiftmethod`.
