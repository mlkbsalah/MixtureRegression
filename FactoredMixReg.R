# Load required libraries
library(MASS)
library(coda)
library(progress)
source("MH_Samplers.R")
source("Model.R")
library(invgamma)

check_adaptation <- function(i, j, acceptance_rate , adaptation_config){
  ac_rate_lower <- adaptation_config$ac_rate_lower
  ac_rate_upper <- adaptation_config$ac_rate_upper
  cycle_period <- adaptation_config$cycle_period
  if ((i%%cycle_period==0) && ((mean(acceptance_rate[j,(i-cycle_period/2):i])>ac_rate_upper) || (mean(acceptance_rate[j,(i-cycle_period/2):i])<ac_rate_lower))) return(TRUE)
  else return(FALSE)
}


MixReg_factor <- function(T, X, Y, Tau, Alpha, ksi, beta_hat, core_temperature) {
  pb <- progress_bar$new(format = ":percent |:bar|  :current/:total  [:elapsed <:eta,  :tick_rateit/s]",
                         total = T,
                         complete = "█",   # Completion bar character
                         incomplete = "░", # Incomplete bar character
                         current = "▒",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar
  pb$tick()
  
  T_adapt <- 100 # Length of the adaptation cycle
  
  Temperature  <- rep(core_temperature, 5) # Temperatures at which to run the chains
  N_chains <- length(Temperature) # Set the number of chains to the number of temperatures
  XXt <- t(X) %*% X # Precompute X'X
  
  # Mixture weights p
  eps_p <- rep(1, N_chains) # Proposition distribution variance
  P <- matrix(NA, nrow = N_chains, ncol = T)
  P[,1] <- rep(0.5, N_chains)
  P_adapt <- numeric(length = T_adapt)
  
  # Initialization etaPhi/etaPhi_sqrd
  eps_etaPhi_sqrd <- rep(1, N_chains) # The proposition is made on \eta_1, \eta_2, \Phi all squared
  etaPhi_sqrd <- array(data = NA, dim = c(N_chains, T ,3))
  etaPhi_sqrd[,1,] <- rdirichlet(N_chains, c(0.5,0.5,0.5))
  
  etaPhi <- array(data = NA, dim = c(N_chains, T ,3))
  coeff <-  cbind( matrix(1, nrow = N_chains, ncol = 2), sample(c(-1, 1), N_chains, replace = TRUE))
  etaPhi[,1,] <- coeff * sqrt(etaPhi_sqrd[,1,])
  etaPhi_adapt <- matrix(data = NA, nrow = T_adapt , ncol = 3)                  
  
  # Mixture std sigma
  eps_sigma_sqrd <- rep(1, N_chains)
  sigma_sqrd <- matrix(nrow = N_chains, ncol = T)
  sigma_sqrd[,1] <- rep(var(Y), N_chains)
  sigma_sqrd_adapt <- matrix(nrow = N_chains, ncol = T_adapt)
  
  # Initialization beta
  eps_beta <- rep(1, N_chains)
  beta_corr <- diag(nrow = ncol(X), ncol = ncol(X))
  Beta_chains <- array(data = NA, dim = c(N_chains, T, k))
  Beta_chains[,1,] <- rmvnorm(N_chains, mean = rep(0,ncol(X)), sigma = beta_corr)
  beta_adapt <- matrix(data = NA, nrow = T_adapt, ncol = k)
  
  # Algorithm vitals
  ac_rate_p <- matrix(NA, nrow = N_chains, ncol = T-1)
  adapt_eps_p <- FALSE
  ac_rate_etaPhi_sqrd <- matrix(NA, nrow = N_chains, ncol = T-1)
  adapt_eps_etaPhi_sqrd <- FALSE
  ac_rate_sigma_sqrd <- matrix(NA, nrow = N_chains, ncol = T-1)
  adapt_eps_sigma_sqrd <- FALSE
  ac_rate_beta <- matrix(NA, nrow = N_chains, ncol = T-1)
  adapt_eps_beta <- FALSE

  
  adaptation_config <- list(ac_rate_lower = .3, 
                            ac_rate_upper = .4, 
                            cycle_period = 50,
                            ac_rate_span = 10)
  
  swap_config <- list(cycle_period = 1)
  
  
  for (i in 1:(T-1)) {
    for (j in 1:N_chains) {
      temp <- Temperature[j]
      mu_XB <- X %*% Beta_chains[j, i, ]
      
      
      # Sample mixture weight
      result <- sample_mixture_weight(P[j, i], eps_p[j], temp, etaPhi[j, i, ], sigma_sqrd[j, i], mu_XB, Tau, Y)
      ac_rate_p[j,i] <- (if (i == 1) result$accept else (ac_rate_p[j,i-1]*(i-1)) + result$accept)/i
      P[j, i + 1] <- result$p_next
      if (check_adaptation(i, j, ac_rate_p, adaptation_config)){
        adapt_eps_p <- TRUE
        ac_rate_window <- numeric(adaptation_config$ac_rate_span)
        P_adapt[1] <- result$p_next
        for (t in 1:(T_adapt-1)){
          result_adapt <- sample_mixture_weight(P_adapt[t], eps_p[j], temp, etaPhi[j, i, ], sigma_sqrd[j, i], mu_XB, Tau, Y)
          P_adapt[t+1] <- result_adapt$p_next
          ac_rate_window[t %% adaptation_config$ac_rate_span + 1] <- result_adapt$accept
          if (t >= adaptation_config$ac_rate_span){ # Check if enough samples were accumulated
            current_ac_rate <- sum(ac_rate_window) / adaptation_config$ac_rate_span
            if (current_ac_rate < adaptation_config$ac_rate_lower) eps_p[j] <- eps_p[j]*1.1
            else if (current_ac_rate > adaptation_config$ac_rate_upper) eps_p[j] <- eps_p[j]*0.8
            }
        }
        adapt_eps_p <- FALSE
      }
      
      # Sample etaPhi
      result <- sample_etaPhi(etaPhi[j, i, ], eps_etaPhi_sqrd[j], temp, P[j, i + 1], sigma_sqrd[j, i], mu_XB, Alpha, Y)
      ac_rate_etaPhi_sqrd[j,i] <- (if (i == 1) result$accept else (ac_rate_etaPhi_sqrd[j,i-1]*(i-1)) + result$accept) / i
      etaPhi[j, i + 1, ] <- result$etaPhi_next
      if (check_adaptation(i, j, ac_rate_etaPhi_sqrd, adaptation_config)){
        adapt_eps_etaPhi_sqrd <- TRUE
        ac_rate_window <- numeric(adaptation_config$ac_rate_span)
        etaPhi_adapt[1, ] <- result$etaPhi_next
        for (t in 1:(T_adapt-1)){
          result_adapt <- sample_etaPhi(etaPhi_adapt[t,], eps_etaPhi_sqrd[j], temp, P[j, i + 1], sigma_sqrd[j, i], mu_XB, Alpha, Y)
          etaPhi_adapt[t+1,] <- result_adapt$etaPhi_next
          ac_rate_window[t %% adaptation_config$ac_rate_span + 1] <- result_adapt$accept
          if (t >= adaptation_config$ac_rate_span){
            current_ac_rate <- sum(ac_rate_window) / adaptation_config$ac_rate_span
            if (current_ac_rate < adaptation_config$ac_rate_lower) eps_etaPhi_sqrd[j] <- eps_etaPhi_sqrd[j] + 100
            else if (current_ac_rate > adaptation_config$ac_rate_upper) eps_etaPhi_sqrd[j] <- max(1, eps_etaPhi_sqrd[j]-100)
          }
        }
        adapt_eps_etaPhi_sqrd <- FALSE
      }
      
      # Sample standard deviation
      result <- sample_sigma(sigma_sqrd[j, i], eps_sigma_sqrd[j], temp, P[j, i + 1], etaPhi[j, i + 1, ], mu_XB, ksi, Y)
      sigma_sqrd[j, i + 1] <- result$sigma_sqrd_next
      ac_rate_sigma_sqrd[j,i] <- (if (i == 1) result$accept else (ac_rate_sigma_sqrd[j,i-1]*(i-1)) + result$accept) / i
      if (check_adaptation(i, j, ac_rate_sigma_sqrd, adaptation_config)){
        adapt_eps_sigma_sqrd <- TRUE
        ac_rate_window <- numeric(adaptation_config$ac_rate_span)
        sigma_sqrd_adapt[j,1] <- result$sigma_sqrd_next
        for (t in 1:(T_adapt-1)){
          result_adapt <- sample_sigma(sigma_sqrd_adapt[j,t], eps_sigma_sqrd[j], temp, P[j, i + 1], etaPhi[j, i + 1, ], mu_XB, ksi, Y)
          sigma_sqrd_adapt[j,t+1] <- result_adapt$sigma_sqrd_next
          ac_rate_window[t %% adaptation_config$ac_rate_span + 1] <- result_adapt$accept
          if (t >= adaptation_config$ac_rate_span){
            current_ac_rate <- sum(ac_rate_window) / adaptation_config$ac_rate_span
            if (current_ac_rate < adaptation_config$ac_rate_lower) eps_sigma_sqrd[j] <- eps_sigma_sqrd[j]*0.95
            else if (current_ac_rate > adaptation_config$ac_rate_upper) eps_sigma_sqrd[j] <- max(1, eps_sigma_sqrd[j]*1.05)
          }
        }
        adapt_eps_sigma_sqrd <- FALSE
      }
      
      # Sample beta
      result <- sample_beta(matrix(Beta_chains[j, i, ]), eps_beta[j], temp, P[j, i + 1], etaPhi[j, i + 1, ], sigma_sqrd[j, i + 1], beta_hat, beta_corr, X, XXt, Y)
      Beta_chains[j, i + 1, ] <- result$beta_next
      ac_rate_beta[j,i] <- (if (i == 1) result$accept else (ac_rate_beta[j,i-1]*(i-1)) + result$accept) / i
      if (check_adaptation(i, j, ac_rate_beta, adaptation_config)){
        adapt_eps_beta <- TRUE
        ac_rate_window <- numeric(adaptation_config$ac_rate_span)
        beta_adapt[1,] <- result$beta_next
        for (t in 1:(T_adapt-1)){
          result_adapt <- sample_beta(matrix(beta_adapt[t,]), eps_beta[j], temp, P[j, i + 1], etaPhi[j, i + 1, ], sigma_sqrd[j, i + 1], beta_hat, beta_corr, X, XXt, Y)
          beta_adapt[t+1,] <- result_adapt$beta_next
          ac_rate_window[t %% adaptation_config$ac_rate_span + 1] <- result_adapt$accept
          if (t >= adaptation_config$ac_rate_span){
            current_ac_rate <- sum(ac_rate_window) / adaptation_config$ac_rate_span
            if (current_ac_rate < adaptation_config$ac_rate_lower) eps_beta[j] <- eps_beta[j]*1
            else if (current_ac_rate > adaptation_config$ac_rate_upper) eps_beta[j] <- eps_beta[j]*1
          }
        }
        adapt_eps_beta <- FALSE
      }
    }
    
    # Swap states between chains
    if (i%%swap_config$cycle_period==0){
      swap_result <- swap_states(P[,i+1], etaPhi[,i+1,], sigma_sqrd[,i+1], Beta_chains[,i+1,], Tau, Alpha, ksi, beta_hat, X, Y, N_chains, Temperature)
      P[,i+1] <- swap_result$p_next
      etaPhi[,i+1,] <- swap_result$etaPhi_next
      sigma_sqrd[,i+1] <- swap_result$sigma_sqrd_next
      Beta_chains[,i+1,] <- swap_result$Beta_chains_next
    }
    
    pb$tick()  # Update progress bar
  }
  
  # Return results
  return(list(P = P, etaPhi = etaPhi, sigma_sqrd = sigma_sqrd, Beta_chains = Beta_chains, 
              ac_rate_p = ac_rate_p, ac_rate_etaPhi_sqrd = ac_rate_etaPhi_sqrd, ac_rate_sigma_sqrd = ac_rate_sigma_sqrd, ac_rate_beta = ac_rate_beta, 
              N_chains = N_chains, Temperature = Temperature))
}
