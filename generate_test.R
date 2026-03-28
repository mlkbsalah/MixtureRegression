library(mvtnorm)
library(DirichletReg)
source("rnormix.R")
library(ggplot2)
generate_test <- function(n_obs, k_var, add_intercept=TRUE, seed=1999, plot_data=TRUE){
  set.seed(seed)
  X <- matrix(data = NA, nrow = n_obs, ncol = k_var)
  X <- rmvnorm(n_obs, mean = c(1,4,6,2))
  if (add_intercept) {
    X <- cbind(rep(1,n_obs), X)
    beta <- matrix(c(3,19,4,9,12),ncol = 1)
  } else {
    beta <- matrix(c(19,4,9,12),ncol = 1)
  }
  mu <- X %*% beta
  etaphi <- c(0.149, 0.406, sqrt(1-0.149^2-0.406^2))
  eta1 <- etaphi[1]
  eta2 <- etaphi[2]
  phi <- etaphi[3]
  
  p <- 0.2
  sigma <- matrix(data=rep(100, n_obs),ncol=1)
  mu1 <- mu - sigma * phi * (sqrt(1-p)/sqrt(p))
  mu2 <- mu + sigma * phi * (sqrt(p)/sqrt(1-p))
  sigma1 <- sigma * (eta1/ sqrt(p))
  sigma2 <- sigma * (eta2/sqrt(1-p))
  
  Y <- rnormix(n_obs, p, mu1, sigma1, mu2, sigma2)$samples

  return(list(X=X, Y=Y, beta=beta, etaphi=etaphi))
}