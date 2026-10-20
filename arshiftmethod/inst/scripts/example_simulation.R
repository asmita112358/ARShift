library(ARShift)

set.seed(2026)
win <- spatstat.geom::square(1)
a <- spatstat.random::rThomas(8, 0.05, 15, win)
b <- spatstat.random::rThomas(7, 0.04, 14, win)
spatstat.geom::marks(a) <- factor(rep(1, spatstat.geom::npoints(a)), levels = c(1, 2))
spatstat.geom::marks(b) <- factor(rep(2, spatstat.geom::npoints(b)), levels = c(1, 2))
X <- spatstat.geom::superimpose(a, b, W = win)
spatstat.geom::marks(X) <- factor(spatstat.geom::marks(X), levels = c(1, 2))

res <- ARShift(data = X, n_perm = 39, workers = 1, type = "inhom")
print(res)
