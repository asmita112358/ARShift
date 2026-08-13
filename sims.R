##All simulations in manuscript
library(ggplot2)
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


win = square(1)
#win = W_new
rmax = 0.3*incircle(win)$r
r = seq(0, rmax, length.out = 50)
n.sim = 100
n_perm = 199

Kcross_tor_rej <- Kcross_vc_rej <- Kcross_area_rej <- Kcross_uncorr_rej <- Kcross_og <- c()
type = "hom"

##NULL scenario: dependent = FALSE
##Alternative scenario: dependent = TRUE
i = 1

while(i <= n.sim){

  data <- generate_LGCP(win, corrfun = "gauss", mu = 6, scale = 0.5, dependant = FALSE)
  
  #data <- generate_Thomas(win, kappa = c(12,10), mu = c(30, 20), scale = c(0.05, 0.04), dependant = FALSE)
  
  
  
  pvals <- tryCatch({
    test_spatial_association(data, base_taxa = 1, shift_taxa = 2,
                             r = r, n_perm = n_perm, bw = "silverman", type = type)
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
  
  
  print(i)
  print(Kcross_area_rej[i])
  i = i+1
}

mean(Kcross_tor_rej <= 0.05, na.rm = TRUE)
mean(Kcross_vc_rej <= 0.05, na.rm = TRUE)
mean(Kcross_area_rej <= 0.05, na.rm = TRUE)
mean(Kcross_uncorr_rej<= 0.05, na.rm = TRUE)

df <- data.frame(data)
ggplot(df, aes(x = x, y = y, color = marks)) +
  geom_point(size = 1) +
  theme_minimal() +
  labs(title = "dependant Thomas, sigma1 = 0.02") +
  scale_color_manual(values = c("red", "blue"))
