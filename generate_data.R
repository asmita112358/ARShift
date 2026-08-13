library(spatstat)

generate_LGCP <- function(win, corrfun = "exponential",mu = 5, scale = 0.1, dependant = FALSE){
  
  data1 <- rLGCP(corrfun, mu = mu, win = win, scale = scale, var = 1)
  
  if(dependant == FALSE){
    data2 <- rLGCP(corrfun, mu = mu, win = win, scale = scale, var = 1)
  } else {
    data2 <- data1
    runif_sh <- runifdisc(data1$n, radius = 0.2)
    data2$x <- data1$x + runif_sh$x
    data2$y <- data1$y + runif_sh$y
    
    keep <- inside.owin(data2$x, data2$y, win)
    data2 <- data2[keep]
  }
  
  marks(data1) <- factor(rep(1, npoints(data1)), levels = c(1, 2))
  marks(data2) <- factor(rep(2, npoints(data2)), levels = c(1, 2))
  
  data_all <- superimpose(data1, data2, W = win)
  marks(data_all) <- factor(marks(data_all), levels = c(1, 2))
  
  return(data_all)
}

generate_Thomas <- function(win = square(1), kappa = c(10, 10), mu = c(20, 20), scale = c(0.05, 0.05), dependant = FALSE){
  data1 <- rThomas(kappa[1], scale[1], mu[1], win,
                   nsim=1, drop=TRUE, 
                   saveLambda=FALSE, expand = 4*scale, 
                   poisthresh=1e-6, saveparents= TRUE)
  if(dependant == FALSE){
    data2 <- rThomas(kappa[2], scale[2], mu[2], win,
                     nsim=1, drop=TRUE, 
                     saveLambda=FALSE, expand = 4*scale, 
                     poisthresh=1e-6, saveparents= TRUE)
  }else{
    data2 <- data1
    runif_sh <- runifdisc(data1$n, radius = 0.2)
    data2$x <- data1$x + runif_sh$x
    data2$y <- data1$y + runif_sh$y
    
    keep <- inside.owin(data2$x, data2$y, win)
    data2 <- data2[keep]
  }
  marks(data1) <- factor(rep(1, npoints(data1)), levels = c(1, 2))
  marks(data2) <- factor(rep(2, npoints(data2)), levels = c(1, 2))
  
  data_all <- superimpose(data1, data2, W = win)
  marks(data_all) <- factor(marks(data_all), levels = c(1, 2))
  
  return(data_all)
}