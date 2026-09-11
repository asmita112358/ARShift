
# ============================================================
# Stress Test: Type I Error of Uncorrected Random Shift
# Testing boundary conditions and extreme scenarios
# ============================================================
#
# Stress dimensions:
#   1. small random shifts (jump_radius = incircle(win)$r/2)
#   2. small random shifts + large cluster scales
#   3. large rmax = 0.5*incircle(win)$r
#   4. tiny rmax = 0.1*incircle(win)$r

# ============================================================

library(spatstat)
library(spatstat.geom)

set.seed(123)

# ============================================================
# FIXED PARAMETERS
# ============================================================

W       <- square(1)
alpha   <- 0.05
N_sims  <- 500   # per scenario; reduce to 200 for quick runs


ARShift_wrapper <- function(win = square(1), stress_case = 0) {
  
  if(stress_case == 0){   #base case
    rmax = 0.3*incircle(win)$r
    scale = c(0.05, 0.04)
    jump_radius = incircle(win)$r
  }else if(stress_case == 1){
    rmax = 0.3*incircle(win)$r
    scale = c(0.05, 0.04)
    jump_radius = incircle(win)$r/2
  }else if(stress_case == 2){
    rmax = 0.3*incircle(win)$r
    scale = c(0.25, 0.25)
    jump_radius = incircle(win)$r/2
  }else if(stress_case == 3){
    rmax = 0.5*incircle(win)$r
    scale = c(0.05, 0.04)
    jump_radius = incircle(win)$r
  }else if(stress_case == 4){
    rmax = 0.1*incircle(win)$r
    scale = c(0.05, 0.04)
    jump_radius = incircle(win)$r
  }
  data <- generate_Thomas(win, kappa = c(12,10), mu = c(30, 20), scale = scale, dependant = FALSE)
  r <- seq(0, rmax, length.out = 50)
 obj <- test_spatial_association(data, base_taxa = 1, shift_taxa = 2,
                                 r = r, n_perm = 499, bw = "silverman", type = "hom", include_RL = FALSE, 
                                 jump_radius = jump_radius)
 return(obj$pval_Kcross_area)
 
}

results = data.frame()

results = rbind(results, mean(replicate(N_sims, ARShift_wrapper(win = W, stress_case = 0)) <= 0.05))
results = rbind(results, mean(replicate(N_sims, ARShift_wrapper(win = W, stress_case = 1)) <= 0.05))
results = rbind(results, mean(replicate(N_sims, ARShift_wrapper(win = W, stress_case = 2)) <= 0.05))
results = rbind(results,mean(replicate(N_sims, ARShift_wrapper(win = W, stress_case = 3)) <= 0.05))
results = rbind(results, mean(replicate(N_sims, ARShift_wrapper(win = W, stress_case = 4)) <= 0.05))
results$stress_event = c("base_case", "small_shift", "small_shift_large_cluster", "large_rmax", "tiny_rmax")
beepr::beep(4)
colnames(results)[1] <- "Type_I_Error"
setwd("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_assoc_test/ARShift")
saveRDS(results, file = "stress_test_results.rds")
