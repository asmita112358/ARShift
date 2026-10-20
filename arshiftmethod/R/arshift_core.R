validate_arshift_inputs <- function(data,
                                    base_taxa,
                                    shift_taxa,
                                    r,
                                    n_perm,
                                    k,
                                    bw,
                                    type,
                                    lite,
                                    thin_prob,
                                    include_RL,
                                    jump_radius,
                                    correction,
                                    workers) {
  if (!inherits(data, "ppp")) {
    stop("`data` must be a spatstat point pattern (`ppp`).", call. = FALSE)
  }

  if (length(base_taxa) != 1 || length(shift_taxa) != 1) {
    stop("`base_taxa` and `shift_taxa` must be scalar mark values.", call. = FALSE)
  }

  if (identical(as.character(base_taxa), as.character(shift_taxa))) {
    stop("`base_taxa` and `shift_taxa` must be different.", call. = FALSE)
  }

  if (!is.null(r)) {
    if (!is.numeric(r) || length(r) < 2 || any(!is.finite(r)) || any(r < 0)) {
      stop("`r` must be a numeric vector of at least two finite non-negative values.", call. = FALSE)
    }
  }

  if (!is.numeric(n_perm) || length(n_perm) != 1 || !is.finite(n_perm) || n_perm < 1) {
    stop("`n_perm` must be a positive number.", call. = FALSE)
  }

  if (!is.numeric(k) || length(k) != 1 || !is.finite(k) || k < 1) {
    stop("`k` must be >= 1.", call. = FALSE)
  }

  if (!is.character(bw) && !is.numeric(bw)) {
    stop("`bw` must be one of 'silverman'/'scott' or a numeric bandwidth.", call. = FALSE)
  }

  if (!(type %in% c("inhom", "hom"))) {
    stop("`type` must be 'inhom' or 'hom'.", call. = FALSE)
  }

  if (!is.logical(lite) || length(lite) != 1) {
    stop("`lite` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.numeric(thin_prob) || length(thin_prob) != 1 || thin_prob <= 0 || thin_prob >= 1) {
    stop("`thin_prob` must be in (0, 1).", call. = FALSE)
  }

  if (!is.logical(include_RL) || length(include_RL) != 1) {
    stop("`include_RL` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.null(jump_radius) && (!is.numeric(jump_radius) || length(jump_radius) != 1 || jump_radius <= 0)) {
    stop("`jump_radius` must be a positive number when provided.", call. = FALSE)
  }

  if (!(correction %in% c("isotropic", "translate", "translation", "border"))) {
    stop("Unsupported correction. Use one of isotropic/translate/translation/border.", call. = FALSE)
  }

  if (!is.numeric(workers) || length(workers) != 1 || workers < 1) {
    stop("`workers` must be a positive integer.", call. = FALSE)
  }

  if (is.null(spatstat.geom::marks(data))) {
    stop("`data` must contain marks for taxa labels.", call. = FALSE)
  }

  invisible(TRUE)
}

safe_pmclapply <- function(X, FUN, workers, ...) {
  workers <- as.integer(workers)

  if (workers <= 1 || .Platform$OS.type == "windows") {
    return(lapply(X, FUN, ...))
  }

  mcprogress::pmclapply(
    X = X,
    FUN = FUN,
    mc.cores = workers,
    ...
  )
}

one_area_permutation_attempt <- function(attempt_id,
                                         points_shift_original,
                                         points_other_original,
                                         original_window,
                                         base_taxa,
                                         shift_taxa,
                                         r,
                                         type,
                                         correction,
                                         sigma_val,
                                         jump_radius) {
  base_chr <- as.character(base_taxa)
  shift_chr <- as.character(shift_taxa)
  out <- list(valid = FALSE)

  shift_vector <- spatstat.random::runifdisc(1, radius = jump_radius)
  shift_x <- shift_vector$x
  shift_y <- shift_vector$y

  points_shifted <- tryCatch(
    {
      spatstat.geom::shift(points_shift_original, vec = c(shift_x, shift_y))
    },
    error = function(e) NULL
  )
  if (is.null(points_shifted)) {
    return(out)
  }

  window_shifted <- tryCatch(
    {
      spatstat.geom::shift(original_window, vec = c(shift_x, shift_y))
    },
    error = function(e) NULL
  )
  if (is.null(window_shifted)) {
    return(out)
  }

  window_reduced <- tryCatch(
    {
      spatstat.geom::intersect.owin(original_window, window_shifted)
    },
    error = function(e) NULL
  )
  if (is.null(window_reduced)) {
    return(out)
  }

  reduced_area <- spatstat.geom::area.owin(window_reduced)
  if (!is.finite(reduced_area) || reduced_area <= 0) {
    return(out)
  }

  pp_combined <- tryCatch(
    {
      spatstat.geom::superimpose(points_other_original, points_shifted, check = FALSE)
    },
    error = function(e) NULL
  )
  if (is.null(pp_combined)) {
    return(out)
  }

  pp_reduced <- tryCatch(
    {
      pp_combined[window_reduced]
    },
    error = function(e) NULL
  )
  if (is.null(pp_reduced) || spatstat.geom::npoints(pp_reduced) == 0) {
    return(out)
  }

  freq_vc <- table(spatstat.geom::marks(pp_reduced))
  n_base_vc <- as.numeric(freq_vc[base_chr])
  n_shift_vc <- as.numeric(freq_vc[shift_chr])
  if (is.na(n_base_vc) || is.na(n_shift_vc) || n_base_vc == 0 || n_shift_vc == 0) {
    return(out)
  }

  K_uncorr <- tryCatch(
    {
      compute_Kcross_fast(
        X = pp_reduced,
        base_taxa = base_taxa,
        shift_taxa = shift_taxa,
        r = r,
        type = type,
        correction = correction,
        sigma_val = sigma_val
      )
    },
    error = function(e) NULL
  )
  if (is.null(K_uncorr) || length(K_uncorr) != length(r) || any(!is.finite(K_uncorr))) {
    return(out)
  }

  list(
    valid = TRUE,
    K_uncorr = K_uncorr,
    area_shift = reduced_area,
    shift_vector = c(shift_x, shift_y)
  )
}

ARShift <- function(data,
                    base_taxa = 1,
                    shift_taxa = 2,
                    r = NULL,
                    n_perm = 499,
                    k = 1.2,
                    type = "inhom",
                    lite = FALSE,
                    thin_prob = 0.2,
                    jump_radius = NULL,
                    correction = "isotropic",
                    workers = max(1L, parallel::detectCores() - 1L)) {
  validate_arshift_inputs(
    data = data,
    base_taxa = base_taxa,
    shift_taxa = shift_taxa,
    r = r,
    n_perm = n_perm,
    k = k,
    bw = "silverman",
    type = type,
    lite = lite,
    thin_prob = thin_prob,
    include_RL = FALSE,
    jump_radius = jump_radius,
    correction = correction,
    workers = workers
  )

  n_perm <- as.integer(n_perm)

  if (lite) {
    data <- spatstat.explore::rthin(data, thin_prob)
  }

  original_window <- spatstat.geom::Window(data)
  win_area <- spatstat.geom::area.owin(original_window)

  if (is.null(r)) {
    rmax <- 0.3 * spatstat.geom::incircle(original_window)$r
    r <- seq(0, rmax, length.out = 50)
  }

  base_chr <- as.character(base_taxa)
  shift_chr <- as.character(shift_taxa)
  freq_marks <- table(spatstat.geom::marks(data))
  n_base <- as.numeric(freq_marks[base_chr])
  n_shift <- as.numeric(freq_marks[shift_chr])

  if (is.na(n_base) || is.na(n_shift) || n_base == 0 || n_shift == 0) {
    return(list(
      pval_Kcross_area = NA_real_,
      n_valid = 0L,
      n_valid_total = 0L,
      n_attempts = 0L,
      k = k
    ))
  }

  if (is.null(jump_radius)) {
    jump_radius <- spatstat.geom::incircle(original_window)$r
  }

  sigma_val <- NULL
  if (type == "inhom") {
    sigma_val <- spatstat.explore::bw.CvL(data)
  }

  obs_Kcross <- compute_Kcross_fast(
    X = data,
    base_taxa = base_taxa,
    shift_taxa = shift_taxa,
    r = r,
    type = type,
    correction = correction,
    sigma_val = sigma_val
  )

  points_shift_original <- subset(data, marks == shift_taxa)
  points_other_original <- subset(data, marks != shift_taxa)

  n_attempts <- as.integer(ceiling(k * n_perm))
  attempt_results <- safe_pmclapply(
    X = seq_len(n_attempts),
    FUN = one_area_permutation_attempt,
    workers = workers,
    points_shift_original = points_shift_original,
    points_other_original = points_other_original,
    original_window = original_window,
    base_taxa = base_taxa,
    shift_taxa = shift_taxa,
    r = r,
    type = type,
    correction = correction,
    sigma_val = sigma_val,
    jump_radius = jump_radius
  )

  valid_results <- Filter(function(z) isTRUE(z$valid), attempt_results)
  n_valid_total <- length(valid_results)

  if (n_valid_total < n_perm) {
    warning(
      "Only ", n_valid_total, " valid permutations obtained out of requested ",
      n_perm, ". Consider increasing k. Current k = ", k, ".",
      call. = FALSE
    )
  }

  n_keep <- min(n_perm, n_valid_total)
  if (n_keep == 0L) {
    return(list(
      pval_Kcross_area = NA_real_,
      n_valid = 0L,
      n_valid_total = n_valid_total,
      n_attempts = n_attempts,
      k = k
    ))
  }

  valid_results <- valid_results[seq_len(n_keep)]
  Kcross_uncorr <- matrix(NA_real_, nrow = length(r), ncol = n_keep)
  area_shift <- numeric(n_keep)

  for (b in seq_len(n_keep)) {
    Kcross_uncorr[, b] <- valid_results[[b]]$K_uncorr
    area_shift[b] <- valid_results[[b]]$area_shift
  }

  meanK <- rowMeans(Kcross_uncorr, na.rm = TRUE)
  T_0 <- sqrt(win_area) * (obs_Kcross - meanK)
  T_sim <- Kcross_uncorr - meanK
  T_scaled <- sweep(T_sim, 2, sqrt(area_shift), "*")

  pval_Kcross_area <- compute_envelope_pval_fast(
    r = r,
    obs = T_0,
    sim_m = T_scaled
  )

  list(
    pval_Kcross_area = pval_Kcross_area,
    n_valid = n_keep,
    n_valid_total = n_valid_total,
    n_attempts = n_attempts,
    k = k,
    r = r
  )
}
