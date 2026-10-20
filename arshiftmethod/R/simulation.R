simulate_arshift_pair <- function(win = spatstat.geom::square(1),
                                  model = c("thomas", "cauchy"),
                                  dependant = FALSE,
                                  shift_radius = 0.2,
                                  kappa = c(10, 10),
                                  mu = c(20, 20),
                                  scale = c(0.05, 0.05),
                                  seed = NULL) {
  model <- match.arg(model)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (!inherits(win, "owin")) {
    stop("`win` must be a spatstat window (`owin`).", call. = FALSE)
  }

  if (length(kappa) != 2 || length(mu) != 2 || length(scale) != 2) {
    stop("`kappa`, `mu`, and `scale` must each have length 2.", call. = FALSE)
  }

  if (model == "thomas") {
    data1 <- spatstat.random::rThomas(kappa[1], scale[1], mu[1], win)
    data2 <- spatstat.random::rThomas(kappa[2], scale[2], mu[2], win)
  } else {
    data1 <- spatstat.random::rCauchy(kappa[1], scale[1], mu[1], win)
    data2 <- spatstat.random::rCauchy(kappa[2], scale[2], mu[2], win)
  }

  if (dependant) {
    data2 <- data1
    runif_sh <- spatstat.random::runifdisc(data1$n, radius = shift_radius)
    data2$x <- data1$x + runif_sh$x
    data2$y <- data1$y + runif_sh$y

    keep <- spatstat.geom::inside.owin(data2$x, data2$y, win)
    data2 <- data2[keep]
  }

  spatstat.geom::marks(data1) <- factor(rep(1, spatstat.geom::npoints(data1)), levels = c(1, 2))
  spatstat.geom::marks(data2) <- factor(rep(2, spatstat.geom::npoints(data2)), levels = c(1, 2))

  data_all <- spatstat.geom::superimpose(data1, data2, W = win)
  spatstat.geom::marks(data_all) <- factor(spatstat.geom::marks(data_all), levels = c(1, 2))
  data_all
}

run_arshift_simulation <- function(n_sim = 20,
                                   generator = NULL,
                                   alpha = 0.05,
                                   seed = 123,
                                   test_args = list(
                                     base_taxa = 1,
                                     shift_taxa = 2,
                                     n_perm = 39,
                                     type = "inhom",
                                     workers = 1
                                   )) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (is.null(generator)) {
    generator <- function(i) {
      simulate_arshift_pair(
        model = "cauchy",
        dependant = FALSE
      )
    }
  }

  out <- vector("list", n_sim)

  for (i in seq_len(n_sim)) {
    data_i <- generator(i)
    args_i <- c(list(data = data_i), test_args)
    result <- do.call(ARShift, args_i)

    out[[i]] <- data.frame(
      iter = i,
      pval_Kcross_area = result$pval_Kcross_area,
      stringsAsFactors = FALSE
    )
  }

  res <- do.call(rbind, out)
  res$reject_area <- res$pval_Kcross_area <= alpha

  attr(res, "summary") <- list(
    alpha = alpha,
    n_sim = n_sim,
    area_rate = mean(res$reject_area, na.rm = TRUE)
  )

  res
}
