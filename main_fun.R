library(GET)
test_spatial_association <- function(data, base_taxa = 1, shift_taxa = 2, r, 
                                     n_perm = 199, bw = "silverman", type = "inhom") {
  
  n_dist <- length(r)
  freq_marks <- table(data$marks)
  original_window <- data$window
  win_area <- area.owin(original_window)
  # Add safety checks
  n_base <- as.numeric(freq_marks[as.character(base_taxa)])
  n_shift <- as.numeric(freq_marks[as.character(shift_taxa)])
  
  if (is.na(n_base) || is.na(n_shift) || n_base == 0 || n_shift == 0) {
    stop("Base or shift taxa not found or has zero count")
  }
  
  # Compute observed statistics
  #obs_K <- compute_K(data, base.taxa = base_taxa, shift.taxa = shift_taxa, 
  #   lambda1 = n_base, lambda2 = n_shift, r = r)
  
  
  if(type == "inhom"){
    #obs_Kcross <- Kcross.inhom(data, i = as.character(base_taxa), j = as.character(shift_taxa), 
    # r = r, correction = "isotropic")[[3]]
    obs_Kcross <- Kcross.inhom(data, i = as.character(base_taxa), j = as.character(shift_taxa), r = r, correction = "isotropic")[[3]]
  }else{
    obs_Kcross <- Kcross(data, i = as.character(base_taxa), j = as.character(shift_taxa), 
                         r = r, correction = "isotropic")[[3]]
  }
  
  
  
  
  # Initialize storage for permutation results
  Kcross_toroidal <- matrix(NA, nrow = n_perm, ncol = n_dist)
  Kcross_uncorr <- matrix(NA, nrow = n_perm, ncol = n_dist)
  Kcross_areacorr <- matrix(NA, nrow = n_perm, ncol = n_dist)
  
  
  
  shift_vectors <- matrix(nrow = 2, ncol = n_perm)
  area_shift <- c()
  
  # Prepare toroidal shift data
  data_toroidal <- data
  if (!is.rectangle(Window(data))) {
    Window(data_toroidal) <- boundingbox(Window(data))
  }
  
  perm_idx <- 1
  
  while (perm_idx <= n_perm) {
    # Toroidal shift
    data_shifted_tor <- rshift(data_toroidal, which = as.character(shift_taxa), edge = "torus")
    freq_tor <- table(data_shifted_tor$marks)
    n_base_tor <- as.numeric(freq_tor[as.character(base_taxa)])
    n_shift_tor <- as.numeric(freq_tor[as.character(shift_taxa)])
    
    if (is.na(n_base_tor) || is.na(n_shift_tor) || n_base_tor == 0 || n_shift_tor == 0) next
    
    if(type == "inhom"){
      Kcross_toroidal[perm_idx, ] <- Kcross.inhom(data_shifted_tor, i = as.character(base_taxa), j = as.character(shift_taxa),
                                                  r = r, correction = "isotropic")[[3]]
    }else{
      Kcross_toroidal[perm_idx, ] <- Kcross(data_shifted_tor, i = as.character(base_taxa), j = as.character(shift_taxa),
                                            r = r, correction = "isotropic")[[3]]
    }
    
    
    
    # Non-toroidal shift (variance correction)
    jump_radius <- incircle(data$window)$r
    shift_vector <- runifdisc(1, radius = jump_radius)
    shift_x <- shift_vector$x
    shift_y <- shift_vector$y
    
    points_shift <- subset(data, marks == shift_taxa)
    points_other <- subset(data, marks != shift_taxa)
    points_shifted <- spatstat.geom::shift(points_shift, vec = c(shift_x, shift_y))
    
    window_shifted <- spatstat.geom::shift(original_window, vec = c(shift_x, shift_y))
    window_reduced <- intersect.owin(original_window, window_shifted)
    
    if (is.null(window_reduced) || area.owin(window_reduced) == 0) next
    
    pp_combined <- superimpose(points_other, points_shifted)
    pp_reduced <- pp_combined[window_reduced]
    pp_og <- data[window_reduced]
    
    if(type == "inhom"){
      Kcross_uncorr[perm_idx,] <- Kcross.inhom(pp_reduced, i = as.character(base_taxa), j = as.character(shift_taxa), 
                                               r = r, correction = "isotropic")[[3]]
    
    }else{
      Kcross_uncorr[perm_idx,] <- Kcross(pp_reduced, i = as.character(base_taxa), j = as.character(shift_taxa), 
                                         r = r, correction = "isotropic")[[3]]
   
    }
    
    area_shift[perm_idx] <- area.owin(window_reduced)
    
    
    freq_vc <- table(pp_reduced$marks)
    n_base_vc <- as.numeric(freq_vc[as.character(base_taxa)])
    n_shift_vc <- as.numeric(freq_vc[as.character(shift_taxa)])
    
    if (is.na(n_base_vc) || is.na(n_shift_vc) || n_base_vc == 0 || n_shift_vc == 0) next
    
    shift_vectors[,perm_idx] = c(shift_x, shift_y)
    
    
    
    
    perm_idx <- perm_idx + 1
  }
 
  #variance correction by area of window
  meanK <- colMeans(Kcross_uncorr, na.rm = TRUE)
  T_0 <- sqrt(win_area)*(obs_Kcross - meanK)
  T_sim <-sweep(t(Kcross_uncorr),1, meanK, "-")
  T_scaled <- T_sim*sqrt(area_shift)
  
  CS <- curve_set(r = r, obs = T_0, sim = T_scaled)
  
  # Global envelope tests for toroidal shift
  pval_Kcross_tor <- compute_envelope_pval(r, obs_Kcross, Kcross_toroidal)
  
  # Variance correction with bandwidth selection
  shift_vectors_full <- cbind(shift_vectors, c(0, 0))
  bw_shift <- select_bandwidth(shift_vectors_full, method = bw)
  
  # Variance-corrected tests 
  vc_results <- compute_variance_corrected_tests(
    Kcross_uncorr, obs_Kcross, shift_vectors_full, 
    bw_shift, r, n_perm
  )
  
  # Use compute_envelope_pval for Kcross_uncorr
  pval_Kcross_uncorr <- compute_envelope_pval(r, obs_Kcross, Kcross_uncorr)
  
  pval_Kcross_area <- compute_envelope_pval(r, T_0, t(T_scaled))
  
  
  return(list(
    pval_Kcross_tor = pval_Kcross_tor,
    pval_Kcross_vc = vc_results$pval_vc_shift,
    pval_Kcross_uncorr = pval_Kcross_uncorr, 
    pval_Kcross_area = pval_Kcross_area
  ))
}


## Helper functions


compute_variance_corrected_tests <- function(K_vc_matrix, obs_K, shift_vectors, 
                                             bw_shift, r, n_perm) {
  
  K_full <- cbind(t(K_vc_matrix), obs_K)
  K_mean <- rowMeans(K_full, na.rm = TRUE)
  K_centered <- sweep(K_full, 1, K_mean, "-")
  
  # Shift-based variance correction
  var_shift <- t(nadaraya_watson(X = t(shift_vectors), Y = t(K_centered^2),
                                 X_pred = t(shift_vectors), bw = bw_shift))
  var_shift[var_shift <= 1e-10] <- 1e-10
  K_std_shift <- K_centered / sqrt(var_shift)
  K_std_shift[!is.finite(K_std_shift)] <- 0
  
  CS_shift <- create_curve_set(list(r = r, obs = K_std_shift[, n_perm + 1],
                                    sim_m = K_std_shift[, 1:n_perm]))
  rank_env_obj <- rank_envelope(CS_shift, type = "erl", alternative = "greater")
  
  return(list(
    pval_vc_shift = attr(rank_env_obj, "p")
  ))
}

nadaraya_watson <- function(X, Y, X_pred = X, bw) {
  n <- nrow(X)
  m <- nrow(X_pred)
  
  if (is.vector(Y)) {
    Y <- matrix(Y, ncol = 1)
  }
  
  n_response <- ncol(Y)
  predictions <- matrix(0, nrow = m, ncol = n_response)
  
  # Gaussian kernel
  gaussian_kernel <- function(dist_sq, bandwidth) {
    exp(-0.5 * dist_sq / bandwidth^2)
  }
  
  for (i in seq_len(m)) {
    # Squared Euclidean distances
    dist_sq <- rowSums((X - matrix(X_pred[i, ], nrow = n, ncol = ncol(X), byrow = TRUE))^2)
    weights <- gaussian_kernel(dist_sq, bw)
    sum_weights <- sum(weights)
    
    if (sum_weights > 1e-10) {
      predictions[i, ] <- crossprod(weights, Y) / sum_weights
    } else {
      predictions[i, ] <- NA
    }
  }
  
  if (n_response == 1) {
    predictions <- as.vector(predictions)
  }
  
  return(predictions)
}

compute_envelope_pval <- function(r, obs, sim_matrix) {
  df <- data.frame(r = r, obs = obs, t(sim_matrix))
  clean_idx <- apply(df, 1, function(x) all(is.finite(x)))
  df <- df[clean_idx, ]
  
  if (nrow(df) == 0) {
    return(NA)
  }
  
  cset <- create_curve_set(list(r = df$r, obs = df$obs, sim_m = as.matrix(df[, -c(1, 2)])))
  result <- rank_envelope(cset, type = "erl", alternative = "greater")
  return(attr(result, "p"))
}

select_bandwidth <- function(X, method = "silverman") {
  n_sim <- ncol(X)
  d <- nrow(X)
  
  sd_X <- apply(X, 1, sd)
  sd_X[sd_X < 1e-10] <- 1e-10  # Handle zero variance
  
  sigma_X <- exp(mean(log(sd_X)))
  
  if (method == "silverman") {
    bw <- sigma_X * (n_sim * (d + 2) / 4)^(-1 / (d + 4))
  } else if (method == "scott") {
    bw <- sigma_X * n_sim^(-1 / (d + 4))
  } else {
    bw <- as.numeric(method)
  }
  
  bw <- max(bw, 1e-10)  # Ensure minimum bandwidth
  
  return(bw)
}