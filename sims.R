##All simulations in manuscript
library(ggplot2)
library(mcprogress)
rm(list = ls())

source("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_assoc_test/ARShift/generate_data.R", echo = FALSE)
source("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_assoc_test/ARShift/main_fun.R", echo = FALSE)


#Load Irregular window if needed
path <- "2023_10_18_hsdm_slide_IIL_fov_01_centroid_sciname.csv"
load(paste0("/Users/asmitaroy/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome4/spatial_microbiome4/windows/", path, ".RData"))

W_new <- W
W_new$bdry[[1]]$x <- scale(W$bdry[[1]]$x, center = 0, scale = 6000)
W_new$bdry[[1]]$y <- scale(W$bdry[[1]]$y, center = 0, scale = 6000)
W_new$xrange <- range(W_new$bdry[[1]]$x)
W_new$yrange <- range(W_new$bdry[[1]]$y)


#win = square(1)
win = W_new
rmax = 0.3*incircle(win)$r
r = seq(0, rmax, length.out = 50)
n.sim = 1000
n_perm = 499

Kcross_tor_rej <- Kcross_vc_rej <- Kcross_area_rej <- Kcross_uncorr_rej <- Kcross_og <- c()
Kcross_RL <- c()
type = "inhom"

##NULL scenario: dependent = FALSE
##Alternative scenario: dependent = TRUE
i = 1

while(i <= n.sim){

  #data <- generate_LGCP(win, corrfun = "gauss", mu = 6, scale = 0.5, dependant = FALSE)
  
  #data <- generate_Thomas(win, kappa = c(12,10), mu = c(30, 20), scale = c(0.05, 0.04), dependant = FALSE)
  
  #data <- generate_inhom_LGCP(win, mu = 6, corrfun = "matern", dependant = TRUE)
  
  data <- generate_inhom_thomas(win,  scale = c(0.05, 0.04), dependant = FALSE)
  
  
  
  pvals <- tryCatch({
    test_spatial_association(data, base_taxa = 1, shift_taxa = 2,
                             r = r, n_perm = n_perm, bw = "silverman", type = type, include_RL = TRUE)
  }, error = function(e) {
    message(paste("Error in iteration", i, ":", e$message))
    return(NULL)
  }, warning = function(w) {
    message(paste("Warning in iteration", i, ":", w$message))
  })
  
  if(is.null(pvals)){
    next
  }
  
  Kcross_tor_rej[i] <- pvals$pval_Kcross_tor
  Kcross_vc_rej[i] <- pvals$pval_Kcross_vc
  Kcross_area_rej[i] <- pvals$pval_Kcross_area
  Kcross_uncorr_rej[i] <- pvals$pval_Kcross_uncorr
  Kcross_RL[i] <- pvals$pval_RL
  
  
  print(i)
  print(Kcross_area_rej[i])
  i = i+1
}

mean(Kcross_tor_rej <= 0.05, na.rm = TRUE)
mean(Kcross_vc_rej <= 0.05, na.rm = TRUE)
mean(Kcross_area_rej <= 0.05, na.rm = TRUE)
mean(Kcross_uncorr_rej<= 0.05, na.rm = TRUE)
mean(Kcross_RL <= 0.05, na.rm = TRUE)
beepr::beep(4)
df <- data.frame(data)
ggplot(df, aes(x = x, y = y, color = marks)) +
  geom_point(size = 1) +
  theme_minimal() +
  labs(title = "Homogeneous Neyman-scott, cauchy kernel") +
  scale_color_manual(values = c("red", "blue"))



##parallel version for inhomogeneous K
run_one_sim <- function(i) {
  
  tryCatch({
    
    withCallingHandlers({
      
      # data <- generate_inhom_thomas(
      #   win,
      #   scale = c(0.05, 0.04),
      #   dependant = TRUE
      # )
      #data <- generate_inhom_LGCP(win, mu = 6, corrfun = "matern", dependant = TRUE)
      
      # data <- generate_inhom_cauchy(
      #   win,
      #     scale = c(0.05, 0.04),
      #     dependant = TRUE
      # )
      
      data <- generate_rCauchy(win = win, kappa = c(20, 15), mu = c(30, 20), scale = c(0.15, 0.1), dependant = FALSE)
      #data <- generate_Thomas(win, kappa = c(10,10), mu = c(20, 20), scale = c(0.05, 0.04), dependant = FALSE)
      pvals <- test_spatial_association(
        data,
        base_taxa = 1,
        shift_taxa = 2,
        r = r,
        n_perm = n_perm,
        bw = "silverman",
        type = "hom", #type = "inhom"
        include_RL = TRUE
      )
      
      data.frame(
        iter = i,
        Kcross_tor_rej    = pvals$pval_Kcross_tor,
        Kcross_vc_rej     = pvals$pval_Kcross_vc,
        Kcross_area_rej   = pvals$pval_Kcross_area,
        Kcross_uncorr_rej = pvals$pval_Kcross_uncorr,
        Kcross_RL         = pvals$pval_RL,
        error             = NA_character_,
        stringsAsFactors  = FALSE
      )
      
    }, warning = function(w) {
      message(sprintf("Warning in iteration %s: %s", i, conditionMessage(w)))
      invokeRestart("muffleWarning")
    })
    
  }, error = function(e) {
    
    message(sprintf("Error in iteration %s: %s", i, conditionMessage(e)))
    
    data.frame(
      iter = i,
      Kcross_tor_rej    = NA_real_,
      Kcross_vc_rej     = NA_real_,
      Kcross_area_rej   = NA_real_,
      Kcross_uncorr_rej = NA_real_,
      Kcross_RL         = NA_real_,
      error             = conditionMessage(e),
      stringsAsFactors  = FALSE
    )
  })
}



n.cores <- max(1, parallel::detectCores() - 1)

res_list <- mcprogress::pmclapply(
  X = seq_len(n.sim),
  FUN = run_one_sim,
  mc.cores = n.cores,
  mc.set.seed = TRUE
)

res <- do.call(rbind, res_list)
Kcross_tor_rej    <- res$Kcross_tor_rej
Kcross_vc_rej     <- res$Kcross_vc_rej
Kcross_area_rej   <- res$Kcross_area_rej
Kcross_uncorr_rej <- res$Kcross_uncorr_rej
Kcross_RL         <- res$Kcross_RL
mean(Kcross_tor_rej <= 0.05, na.rm = TRUE)
mean(Kcross_vc_rej <= 0.05, na.rm = TRUE)
mean(Kcross_area_rej <= 0.05, na.rm = TRUE)
mean(Kcross_uncorr_rej<= 0.05, na.rm = TRUE)
mean(Kcross_RL <= 0.05, na.rm = TRUE)
