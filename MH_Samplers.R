source("Model.R")

# Function to sample the mixture weight
sample_mixture_weight <- function(p_t, eps_p, temp, etaPhi_t, sigma_sqrd_t, mu_XB, Tau, Y) {
  p_prop <- rbeta(1, (p_t * eps_p) + 1, ((1 - p_t) * eps_p) + 1)
  
  log_acceptance <- (temp * Log_post_p(p_prop, etaPhi_t, sqrt(sigma_sqrd_t), mu_XB, Tau, Y) +
                       dbeta(p_t, p_prop * eps_p + 1, (1 - p_prop) * eps_p + 1, log = TRUE)) -
    (temp * Log_post_p(p_t, etaPhi_t, sqrt(sigma_sqrd_t), mu_XB, Tau, Y) +
       dbeta(p_prop, p_t * eps_p + 1, (1 - p_t) * eps_p + 1, log = TRUE))
  
  s <-(runif(1) < min(1, exp(log_acceptance)))
  return(list(p_next = ifelse(s, p_prop, p_t), accept=s))
}

# Function to sample etaPhi
sample_etaPhi <- function(etaPhi_t, eps_etaPhi_sqrd, temp, p_next, sigma_sqrd_t, mu_XB, Alpha, Y) {
  etaPhi_sqrd_t <- etaPhi_t^2
  etaPhi_sqrd_prop <-
    rdirichlet(1,
               c(
                 etaPhi_sqrd_t[1] * eps_etaPhi_sqrd+1,
                 etaPhi_sqrd_t[2] * eps_etaPhi_sqrd+1,
                 etaPhi_sqrd_t[3] * eps_etaPhi_sqrd+1
               ))
  coeff <- c(1,1,sample(c(-1,1),size = 1))
  etaPhi_prop <- coeff * sqrt(etaPhi_sqrd_prop)
  tmp1 <-
    c(
      etaPhi_sqrd_prop[1] * eps_etaPhi_sqrd + 1,
      etaPhi_sqrd_prop[2] * eps_etaPhi_sqrd + 1,
      etaPhi_sqrd_prop[3] * eps_etaPhi_sqrd + 1
    )
  tmp2 <-
    c(etaPhi_sqrd_t[1] * eps_etaPhi_sqrd + 1,
      etaPhi_sqrd_t[2] * eps_etaPhi_sqrd + 1,
      etaPhi_sqrd_t[3] * eps_etaPhi_sqrd + 1)
  
  log_acceptance <-
    (
      temp*Log_post_etaPhi(p_next, etaPhi_prop, sqrt(sigma_sqrd_t), mu_XB, Alpha, Y) + ddirichlet(matrix(etaPhi_sqrd_t, nrow = 1), tmp1 , log = TRUE)
    ) - (temp*Log_post_etaPhi(p_next, etaPhi_t, sqrt(sigma_sqrd_t), mu_XB, Alpha, Y) + ddirichlet(matrix(etaPhi_sqrd_prop, nrow = 1), tmp2, log = TRUE))
  s <-(runif(1) < min(1, exp(log_acceptance)))
  if (s) etaPhi_next <- etaPhi_prop
  else etaPhi_next <- etaPhi_t
  
  return (list(etaPhi_next = etaPhi_next, accept=s))
}

# Function to sample standard deviation
sample_sigma <- function(sigma_sqrd_t, eps_sigma_sqrd, temp, p_next, etaPhi_next, mu_XB, ksi, Y) {
  sigma_sqrd_prop <- rinvgamma(1, 2 + ((sigma_sqrd_t^2) / eps_sigma_sqrd), sigma_sqrd_t * (1 + (sigma_sqrd_t^2) / eps_sigma_sqrd))
  
  log_acceptance <- (temp * Log_post_sigma(p_next, etaPhi_next, sqrt(sigma_sqrd_prop), mu_XB, ksi, Y) +
                       dinvgamma(sigma_sqrd_t, 2 + ((sigma_sqrd_prop^2) / eps_sigma_sqrd), sigma_sqrd_prop * (1 + (sigma_sqrd_prop^2) / eps_sigma_sqrd), log = TRUE)) -
    (temp * Log_post_sigma(p_next, etaPhi_next, sqrt(sigma_sqrd_t), mu_XB, ksi, Y) +
       dinvgamma(sigma_sqrd_prop, 2 + ((sigma_sqrd_t^2) / eps_sigma_sqrd), sigma_sqrd_t * (1 + (sigma_sqrd_t^2) / eps_sigma_sqrd), log = TRUE))
  
  s <- (runif(1) < min(1, exp(log_acceptance)))

  return(list(sigma_sqrd_next=(ifelse(s, sigma_sqrd_prop, sigma_sqrd_t)), accept=s))
}

# Function to sample beta
sample_beta <- function(beta_t, eps_beta, temp, p_next, etaPhi_next, sigma_sqrd_t, beta_hat, beta_corr, X, XXt, Y) {
  beta_prop <-matrix(rmvnorm(1,
                             mean = beta_t,
                             sigma = eps_beta * beta_corr),
                     ncol = 1)
  
  log_acceptance <-
    (temp*Log_post_beta(p_next, etaPhi_next, sqrt(sigma_sqrd_t), beta_prop, beta_hat, X, XXt, Y) + dmvnorm(t(beta_t), mean = beta_prop, sigma = eps_beta * beta_corr, log = TRUE)) -
    (temp*Log_post_beta(p_next, etaPhi_next, sqrt(sigma_sqrd_t), beta_t   , beta_hat, X, XXt, Y) + dmvnorm(t(beta_prop), mean = beta_t, sigma = eps_beta * beta_corr, log = TRUE))
  s <-(runif(1) < min(1, exp(log_acceptance)))
  
  if (s) beta_next <- beta_prop
  else beta_next <- beta_t
  return (list(beta_next=beta_next, accept=s))
}

# Function to swap states between chains
swap_states <- function(p_t, etaPhi_t, sigma_sqrd_t, Beta_chains_t, Tau, Alpha, ksi, beta_hat, X, Y, N_chains, Temperature){
  swap <- matrix(sample(1:(2*N_chains%/%2), replace = FALSE), ncol = 2)
  accept <- numeric(dim(swap)[1])
  for (c in 1:dim(swap)[1]){
    r <- swap[c,1]
    s <- swap[c,2]
    log_swap_acceptance <-
      Temperature[s] * Log_post_chain(p_t[r], etaPhi_t[r,], sqrt(sigma_sqrd_t[r]), matrix(Beta_chains_t[r,],ncol = 1), Tau, Alpha, ksi, beta_hat, X, Y) +
      Temperature[r] * Log_post_chain(p_t[s], etaPhi_t[s,], sqrt(sigma_sqrd_t[s]), matrix(Beta_chains_t[s,],ncol = 1), Tau, Alpha, ksi, beta_hat, X, Y) -
      Temperature[r] * Log_post_chain(p_t[r], etaPhi_t[r,], sqrt(sigma_sqrd_t[r]), matrix(Beta_chains_t[r,],ncol = 1), Tau, Alpha, ksi, beta_hat, X, Y) -
      Temperature[s] * Log_post_chain(p_t[s], etaPhi_t[s,], sqrt(sigma_sqrd_t[s]), matrix(Beta_chains_t[s,],ncol = 1), Tau, Alpha, ksi, beta_hat, X, Y)
    
    accept[c] <- (min(1, exp(log_swap_acceptance)) > runif(1))
    if (accept[c]){
      # Swap p
      p_tmp <- p_t[r]
      p_t[r] <- p_t[s]
      p_t[s] <- p_tmp
      # cat("p swapped", p_t,"\n")
      
      
      # Swap etaPhi
      etaPhi_tmp <- etaPhi_t[r,]
      etaPhi_t[r,] <- etaPhi_t[s,]
      etaPhi_t[s,] <- etaPhi_tmp
      # cat("etaPhi swapped", etaPhi_t,"\n")
      
      
      # Swap sigma_sqrd
      sigma_sqrd_tmp <- sigma_sqrd_t[r]
      sigma_sqrd_t[r] <- sigma_sqrd_t[s]
      sigma_sqrd_t[s] <- sigma_sqrd_tmp
      # cat("sigma_sqrd swapped", sigma_sqrd_t,"\n")
      
      # Swap Beta_chains
      Beta_chains_tmp <- Beta_chains_t[r,]
      Beta_chains_t[r,] <- Beta_chains_t[s,]
      Beta_chains_t[s,] <- Beta_chains_tmp
    }
    
    return(list(p_next = p_t, etaPhi_next = etaPhi_t, sigma_sqrd_next = sigma_sqrd_t, Beta_chains_next = Beta_chains_t, accept = accept))
  }
}
