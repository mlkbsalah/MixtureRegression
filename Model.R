Log_lik <- function(p, etaPhi, sigma, mu_XB, Y){
  # Precompute sqrt terms
  sqrt_p <- sqrt(p)
  sqrt_1_p <- sqrt(1 - p)
  
  # sd and mean
  mu1 <- mu_XB - sigma * etaPhi[3] * (sqrt_1_p / sqrt_p)
  mu2 <- mu_XB + sigma * etaPhi[3] * (sqrt_p / sqrt_1_p)
  sigma1 <- sigma * (etaPhi[1] / sqrt_p)
  sigma2 <- sigma * (etaPhi[2] / sqrt_1_p)
  
  # Vectorized calculation of log likelihood
  res <- sum(log(p * dnorm(Y, mu1, sigma1) + (1 - p) * dnorm(Y, mu2, sigma2)))
  
  log_res1 <- log(p) + dnorm(Y, mu1, sigma1, log = TRUE)
  log_res2 <- log1p(-p) + dnorm(Y, mu2, sigma2, log = TRUE)
  
  # Apply the log-sum-exp trick
  max_log_res <- pmax(log_res1, log_res2)
  res <- sum(max_log_res + log(exp(log_res1 - max_log_res) + exp(log_res2 - max_log_res)))
  
  return(res)
}

Log_post_p <- function(p_val, etaPhi_val, sigma_val, mu_XB, Tau, Y_val) {
  return(Log_lik(p_val, etaPhi_val, sigma_val, mu_XB, Y_val) + dbeta(p_val, Tau$tau1, Tau$tau2, log = TRUE))
}

Log_post_etaPhi <- function(p_val, etaPhi_val, sigma_val, mu_XB, Alpha, Y_val){
  return(Log_lik(p_val, etaPhi_val, sigma_val, mu_XB, Y_val) + ddirichlet(matrix(etaPhi_val^2, nrow = 1), Alpha, log = TRUE))
}

Log_post_sigma <- function(p_val, etaPhi_val, sigma_val, mu_XB, ksi, Y_val){
  return(Log_lik(p_val, etaPhi_val, sigma_val, mu_XB, Y_val) + dinvgamma(sigma_val^2, ksi$shape, ksi$rate, log = TRUE))
}

Log_post_beta <- function(p_val, etaPhi_val, sigma_val, beta_val, beta_hat, X_val, XXt, Y_val){
  mu_XB <- X_val %*% beta_val 
  return(Log_lik(p_val, etaPhi_val, sigma_val, mu_XB, Y_val) - (dim(beta_val)[1]/2 + 3/2 +1) * log(t(beta_val-beta_hat) %*% XXt %*% (beta_val-beta_hat)))
}

Log_post_chain <- function(p_val, etaPhi_val, sigma_val, beta_val, Tau, Alpha, ksi, beta_hat, X_val, Y_val){
  mu_XB <- X_val %*% beta_val 
  return(Log_lik(p_val, etaPhi_val, sigma_val, mu_XB, Y_val) + dbeta(p_val, Tau$tau1, Tau$tau2, log = TRUE) + ddirichlet(matrix(etaPhi_val^2, nrow = 1), Alpha, log = TRUE) + dinvgamma(sigma_val^2, ksi$shape, ksi$rate, log = TRUE) + (dim(beta_val)[2]/2 + 3/2 +1) * log(t(beta_val-beta_hat) %*% t(X_val) %*% X_val %*% (beta_val-beta_hat)))
}

