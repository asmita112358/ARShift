library(spatstat.geom)
library(spatstat.explore)
library(GET)
library(mcprogress)

test_spatial_association_pmclapply <- function(data,
                                               base_taxa = 1,
                                               shift_taxa = 2,
                                               r = NULL,
                                               n_perm = 499,
                                               k = 1.2,
                                               bw = "silverman",
                                               type = "inhom",
                                               lite = FALSE,
                                               thin_prob = 0.2,
                                               include_RL = FALSE,
                                               jump_radius = NULL,
                                               correction = "isotropic",
                                               workers = max(1, parallel::detectCores() - 1)) {
  
  if (lite) {
    data <- rthin(data, thin_prob)
  }
  
  original_window <- Window(data)
  win_area <- area.owin(original_window)
  
  if (is.null(r)) {
    rmax <- 0.3 * incircle(original_window)$r
    r <- seq(0, rmax, length.out = 50)
  }
  
  n_dist <- length(r)
  base_chr <- as.character(base_taxa)
  shift_chr <- as.character(shift_taxa)
  
  freq_marks <- table(marks(data))
  n_base <- as.numeric(freq_marks[base_chr])
  n_shift <- as.numeric(freq_marks[shift_chr])
  
  if (is.na(n_base) || is.na(n_shift) || n_base == 0 || n_shift == 0) {
    return(list(
      pval_Kcross_tor = NA_real_,
      pval_Kcross_vc = NA_real_,
      pval_Kcross_uncorr = NA_real_,
      pval_Kcross_area = NA_real_,
      pval_RL = NA_real_,
      n_valid = 0,
      n_valid_total = 0,
      n_attempts = 0,
      k = k
    ))
  }
  
  if (is.null(jump_radius)) {
    jump_radius <- incircle(original_window)$r
  }
  
  sigma_val <- NULL
  if (type == "inhom") {
    sigma_val <- bw.CvL(data)
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
  
  data_toroidal <- data
  
  if (!is.rectangle(original_window)) {
    Window(data_toroidal) <- boundingbox(original_window)
  }
  
  points_shift_original <- subset(data, marks == shift_taxa)
  points_other_original <- subset(data, marks != shift_taxa)
  
  n_attempts <- ceiling(k * n_perm)
  
  attempt_results <- mcprogress::pmclapply(
    X = seq_len(n_attempts),
    FUN = one_permutation_attempt_pmclapply,
    data = data,
    data_toroidal = data_toroidal,
    original_window = original_window,
    points_shift_original = points_shift_original,
    points_other_original = points_other_original,
    base_taxa = base_taxa,
    shift_taxa = shift_taxa,
    r = r,
    type = type,
    correction = correction,
    sigma_val = sigma_val,
    jump_radius = jump_radius,
    include_RL = include_RL,
    mc.cores = workers
  )
  
  valid_results <- Filter(function(z) isTRUE(z$valid), attempt_results)
  n_valid_total <- length(valid_results)
  
  if (n_valid_total < n_perm) {
    warning(
      "Only ", n_valid_total, " valid permutations obtained out of requested ",
      n_perm, ". Consider increasing k. Current k = ", k, "."
    )
  }
  
  n_keep <- min(n_perm, n_valid_total)
  
  if (n_keep == 0) {
    return(list(
      pval_Kcross_tor = NA_real_,
      pval_Kcross_vc = NA_real_,
      pval_Kcross_uncorr = NA_real_,
      pval_Kcross_area = NA_real_,
      pval_RL = NA_real_,
      n_valid = 0,
      n_valid_total = n_valid_total,
      n_attempts = n_attempts,
      k = k
    ))
  }
  
  valid_results <- valid_results[seq_len(n_keep)]
  
  Kcross_toroidal <- matrix(NA_real_, nrow = n_dist, ncol = n_keep)
  Kcross_uncorr <- matrix(NA_real_, nrow = n_dist, ncol = n_keep)
  shift_vectors <- matrix(NA_real_, nrow = 2, ncol = n_keep)
  area_shift <- numeric(n_keep)
  
  if (include_RL) {
    Kcross_RL <- matrix(NA_real_, nrow = n_dist, ncol = n_keep)
  } else {
    Kcross_RL <- NULL
  }
  
  for (b in seq_len(n_keep)) {
    Kcross_toroidal[, b] <- valid_results[[b]]$K_tor
    Kcross_uncorr[, b] <- valid_results[[b]]$K_uncorr
    shift_vectors[, b] <- valid_results[[b]]$shift_vector
    area_shift[b] <- valid_results[[b]]$area_shift
    
    if (include_RL) {
      Kcross_RL[, b] <- valid_results[[b]]$K_RL
    }
  }
  
  meanK <- rowMeans(Kcross_uncorr, na.rm = TRUE)
  
  T_0 <- sqrt(win_area) * (obs_Kcross - meanK)
  T_sim <- Kcross_uncorr - meanK
  T_scaled <- sweep(T_sim, 2, sqrt(area_shift), "*")
  
  pval_Kcross_tor <- compute_envelope_pval_fast(
    r = r,
    obs = obs_Kcross,
    sim_m = Kcross_toroidal
  )
  
  pval_Kcross_uncorr <- compute_envelope_pval_fast(
    r = r,
    obs = obs_Kcross,
    sim_m = Kcross_uncorr
  )
  
  pval_Kcross_area <- compute_envelope_pval_fast(
    r = r,
    obs = T_0,
    sim_m = T_scaled
  )
  
  if (include_RL) {
    pval_RL <- compute_envelope_pval_fast(
      r = r,
      obs = obs_Kcross,
      sim_m = Kcross_RL
    )
  } else {
    pval_RL <- NA_real_
  }
  
  shift_vectors_full <- cbind(shift_vectors, c(0, 0))
  bw_shift <- select_bandwidth(shift_vectors_full, method = bw)
  
  vc_results <- compute_variance_corrected_tests_fast(
    K_vc_matrix = Kcross_uncorr,
    obs_K = obs_Kcross,
    shift_vectors = shift_vectors_full,
    bw_shift = bw_shift,
    r = r
  )
  
  list(
    pval_Kcross_tor = pval_Kcross_tor,
    pval_Kcross_vc = vc_results$pval_vc_shift,
    pval_Kcross_uncorr = pval_Kcross_uncorr,
    pval_Kcross_area = pval_Kcross_area,
    pval_RL = pval_RL,
    n_valid = n_keep,
    n_valid_total = n_valid_total,
    n_attempts = n_attempts,
    k = k
  )
}

one_permutation_attempt_pmclapply <- function(attempt_id,
                                              data,
                                              data_toroidal,
                                              original_window,
                                              points_shift_original,
                                              points_other_original,
                                              base_taxa,
                                              shift_taxa,
                                              r,
                                              type,
                                              correction,
                                              sigma_val,
                                              jump_radius,
                                              include_RL = FALSE) {
  
  base_chr <- as.character(base_taxa)
  shift_chr <- as.character(shift_taxa)
  
  out <- list(valid = FALSE)
  
  data_shifted_tor <- tryCatch(
    {
      rshift(
        X = data_toroidal,
        which = shift_chr,
        edge = "torus"
      )
    },
    error = function(e) NULL
  )
  
  if (is.null(data_shifted_tor)) {
    return(out)
  }
  
  K_tor <- tryCatch(
    {
      compute_Kcross_fast(
        X = data_shifted_tor,
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
  
  if (is.null(K_tor) || length(K_tor) != length(r) || any(!is.finite(K_tor))) {
    return(out)
  }
  
  shift_vector <- runifdisc(1, radius = jump_radius)
  shift_x <- shift_vector$x
  shift_y <- shift_vector$y
  
  points_shifted <- tryCatch(
    {
      spatstat.geom::shift(
        X = points_shift_original,
        vec = c(shift_x, shift_y)
      )
    },
    error = function(e) NULL
  )
  
  if (is.null(points_shifted)) {
    return(out)
  }
  
  window_shifted <- tryCatch(
    {
      spatstat.geom::shift(
        X = original_window,
        vec = c(shift_x, shift_y)
      )
    },
    error = function(e) NULL
  )
  
  if (is.null(window_shifted)) {
    return(out)
  }
  
  window_reduced <- tryCatch(
    {
      intersect.owin(original_window, window_shifted)
    },
    error = function(e) NULL
  )
  
  if (is.null(window_reduced)) {
    return(out)
  }
  
  reduced_area <- area.owin(window_reduced)
  
  if (!is.finite(reduced_area) || reduced_area <= 0) {
    return(out)
  }
  
  pp_combined <- tryCatch(
    {
      superimpose(
        points_other_original,
        points_shifted,
        check = FALSE
      )
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
  
  if (is.null(pp_reduced) || npoints(pp_reduced) == 0) {
    return(out)
  }
  
  freq_vc <- table(marks(pp_reduced))
  n_base_vc <- as.numeric(freq_vc[base_chr])
  n_shift_vc <- as.numeric(freq_vc[shift_chr])
  
  if (
    is.na(n_base_vc) || is.na(n_shift_vc) ||
    n_base_vc == 0 || n_shift_vc == 0
  ) {
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
  
  if (include_RL) {
    permuted_marks <- sample(marks(data))
    
    pp_RL <- tryCatch(
      {
        ppp(
          x = data$x,
          y = data$y,
          window = original_window,
          marks = permuted_marks
        )
      },
      error = function(e) NULL
    )
    
    if (is.null(pp_RL)) {
      return(out)
    }
    
    K_RL <- tryCatch(
      {
        compute_Kcross_fast(
          X = pp_RL,
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
    
    if (is.null(K_RL) || length(K_RL) != length(r) || any(!is.finite(K_RL))) {
      return(out)
    }
  } else {
    K_RL <- NULL
  }
  
  list(
    valid = TRUE,
    K_tor = K_tor,
    K_uncorr = K_uncorr,
    K_RL = K_RL,
    shift_vector = c(shift_x, shift_y),
    area_shift = reduced_area
  )
}

compute_Kcross_fast <- function(X,
                                base_taxa,
                                shift_taxa,
                                r,
                                type = "inhom",
                                correction = "isotropic",
                                sigma_val = NULL) {
  
  i <- as.character(base_taxa)
  j <- as.character(shift_taxa)
  
  if (type == "inhom") {
    Kobj <- Kcross.inhom(
      X = X,
      i = i,
      j = j,
      r = r,
      correction = correction,
      sigma = sigma_val
    )
  } else {
    Kobj <- Kcross(
      X = X,
      i = i,
      j = j,
      r = r,
      correction = correction
    )
  }
  
  if (correction == "isotropic") {
    return(Kobj$iso)
  }
  
  if (correction %in% c("translate", "translation")) {
    return(Kobj$trans)
  }
  
  if (correction == "border") {
    return(Kobj$border)
  }
  
  stop("Unsupported correction: ", correction)
}

compute_envelope_pval_fast <- function(r, obs, sim_m) {
  
  if (is.null(sim_m)) {
    return(NA_real_)
  }
  
  if (!is.matrix(sim_m)) {
    sim_m <- as.matrix(sim_m)
  }
  
  keep <- is.finite(r) &
    is.finite(obs) &
    apply(sim_m, 1, function(z) all(is.finite(z)))
  
  if (!any(keep)) {
    return(NA_real_)
  }
  
  cset <- create_curve_set(list(
    r = r[keep],
    obs = obs[keep],
    sim_m = sim_m[keep, , drop = FALSE]
  ))
  
  result <- rank_envelope(
    cset,
    type = "erl",
    alternative = "greater"
  )
  
  attr(result, "p")
}

compute_variance_corrected_tests_fast <- function(K_vc_matrix,
                                                  obs_K,
                                                  shift_vectors,
                                                  bw_shift,
                                                  r) {
  
  n_perm <- ncol(K_vc_matrix)
  
  K_full <- cbind(K_vc_matrix, obs_K)
  K_mean <- rowMeans(K_full, na.rm = TRUE)
  K_centered <- K_full - K_mean
  
  var_shift <- t(nadaraya_watson(
    X = t(shift_vectors),
    Y = t(K_centered^2),
    X_pred = t(shift_vectors),
    bw = bw_shift
  ))
  
  var_shift[var_shift <= 1e-10] <- 1e-10
  
  K_std_shift <- K_centered / sqrt(var_shift)
  K_std_shift[!is.finite(K_std_shift)] <- 0
  
  CS_shift <- create_curve_set(list(
    r = r,
    obs = K_std_shift[, n_perm + 1],
    sim_m = K_std_shift[, seq_len(n_perm), drop = FALSE]
  ))
  
  rank_env_obj <- rank_envelope(
    CS_shift,
    type = "erl",
    alternative = "greater"
  )
  
  list(
    pval_vc_shift = attr(rank_env_obj, "p")
  )
}
nadaraya_watson <- function(X, Y, X_pred = X, bw) {
  
  n <- nrow(X)
  m <- nrow(X_pred)
  
  if (is.vector(Y)) {
    Y <- matrix(Y, ncol = 1)
  }
  
  n_response <- ncol(Y)
  predictions <- matrix(0, nrow = m, ncol = n_response)
  
  gaussian_kernel <- function(dist_sq, bandwidth) {
    exp(-0.5 * dist_sq / bandwidth^2)
  }
  
  for (i in seq_len(m)) {
    dist_sq <- rowSums(
      (
        X -
          matrix(
            X_pred[i, ],
            nrow = n,
            ncol = ncol(X),
            byrow = TRUE
          )
      )^2
    )
    
    weights <- gaussian_kernel(dist_sq, bw)
    sum_weights <- sum(weights)
    
    if (sum_weights > 1e-10) {
      predictions[i, ] <- as.numeric(crossprod(weights, Y) / sum_weights)
    } else {
      predictions[i, ] <- NA_real_
    }
  }
  
  if (n_response == 1) {
    predictions <- as.vector(predictions)
  }
  
  predictions
}

select_bandwidth <- function(X, method = "silverman") {
  
  n_sim <- ncol(X)
  d <- nrow(X)
  
  sd_X <- apply(X, 1, sd)
  sd_X[sd_X < 1e-10] <- 1e-10
  
  sigma_X <- exp(mean(log(sd_X)))
  
  if (method == "silverman") {
    bw <- sigma_X * (n_sim * (d + 2) / 4)^(-1 / (d + 4))
  } else if (method == "scott") {
    bw <- sigma_X * n_sim^(-1 / (d + 4))
  } else {
    bw <- as.numeric(method)
  }
  
  bw <- max(bw, 1e-10)
  
  bw
}
