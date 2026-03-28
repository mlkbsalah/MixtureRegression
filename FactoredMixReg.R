# Load required libraries
library(MASS)
library(coda)
library(progress)
source("MH_Samplers.R")
source("Model.R")
library(invgamma)
library(mvtnorm)
library(DirichletReg)
source("utilities.R")

options(warn=2)

MixReg_factor <- function(T, X, Y, Tau, Alpha, ksi, beta_hat, core_temperature, N_chains) {
  pb <- initialize_progress_bar(T)
  update_progress_bar(pb)
  T_adapt <- 100
  adaptation_config <- list(ac_rate_lower = .3, 
                            ac_rate_upper = .4, 
                            cycle_period = 300,
                            ac_rate_span = 50,
                            T_adapt = 100) # Length of the adaptation cycle
  
  Temperature  <- rep(core_temperature, N_chains) # Temperatures at which to run the chains
  XXt <- t(X) %*% X # Precompute X'X
  k <- ncol(X) # Number of predictors 
  n_obs <- nrow(X) # Number of observations
  
  # Mixture weights p
  eps_p <- rep(1, N_chains) # Proposition distribution variance
  P <- matrix(runif(N_chains), nrow = N_chains, ncol = T)
  
  # Initialization etaPhi/etaPhi_sqrd
  eps_etaPhi_sqrd <- rep(1, N_chains) # The proposition is made on \eta_1, \eta_2, \Phi all squared
  etaPhi_sqrd <- array(data = NA, dim = c(N_chains, T ,3))
  etaPhi_sqrd[,1,] <- rdirichlet(N_chains, c(0.5,0.5,0.5))
  
  etaPhi <- array(data = NA, dim = c(N_chains, T ,3))
  coeff <-  cbind( matrix(1, nrow = N_chains, ncol = 2), sample(c(-1, 1), N_chains, replace = TRUE))
  etaPhi[,1,] <- coeff * sqrt(etaPhi_sqrd[,1,])
  
  # Mixture std sigma
  eps_sigma_sqrd <- rep(1, N_chains)
  sigma_sqrd <- matrix(nrow = N_chains, ncol = T)
  sigma_sqrd[,1] <- rep(var(Y), N_chains)
  
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
  
  
  
  swap_config <- list(cycle_period = 100)
  
  conv_config <- list(
    check_cycle = 250,  # Check every 250 iterations
    window_size = 1000, # Use last 1000 iterations for calculation
    start_iter = 1000    # Start checking after 500 iterations
  )
  gelman <- list()
  
  for (i in 1:(T-1)) {
    for (j in 1:N_chains) {
      temp <- Temperature[j]
      mu_XB <- X %*% Beta_chains[j, i, ]
      
      # Sample mixture weight
      result <- sample_mixture_weight(P[j, i], etaPhi[j, i, ], sigma_sqrd[j, i], mu_XB, Y, temp, eps_p[j], Tau)
      ac_rate_p[j,i] <- (if (i == 1) result$accept else (ac_rate_p[j,i-1]*(i-1)) + result$accept) / i
      P[j, i + 1] <- result$next_val
      # cat("ac_rate_p[",j,",",i,"]: ", ac_rate_p[j,i], "\n")
      if (check_adaptation(i, j, ac_rate_p, adaptation_config)){
        eps_p[j] = adapt_proposition(i=0,P[j, i+1], sample_mixture_weight, list(p_t=P[j,i+1],etaPhi_t=etaPhi[j, i, ], sigma_sqrd_t=sigma_sqrd[j, i], mu_XB=mu_XB, Y=Y, temp=temp, eps=eps_p[j]), "p_t", adaptation_config,Tau = Tau)
      }
      
      
      # Sample etaPhi
      result <- sample_etaPhi(P[j, i + 1], etaPhi[j, i, ], sigma_sqrd[j, i], mu_XB, Y, temp, eps_etaPhi_sqrd[j], Alpha)
      ac_rate_etaPhi_sqrd[j,i] <- (if (i == 1) result$accept else (ac_rate_etaPhi_sqrd[j,i-1]*(i-1)) + result$accept) / i
      etaPhi[j, i + 1, ] <- result$next_val
      if (check_adaptation(i, j, ac_rate_etaPhi_sqrd, adaptation_config)){
        eps_etaPhi_sqrd[j] = adapt_proposition(i=0,etaPhi[j, i+1,], sample_etaPhi, list(p_t=P[j,i+1],etaPhi_t=etaPhi[j, i, ], sigma_sqrd_t=sigma_sqrd[j, i], mu_XB=mu_XB, Y=Y, temp=temp, eps=eps_etaPhi_sqrd[j]), "etaPhi_t", adaptation_config,Alpha = Alpha)
      }
      
      # Check that etaphi squqred sum adds to 1
      if (abs(sum(etaPhi[j, i + 1, ]^2)-1) > 0.0001) {
        cat(sum(etaPhi[j, i + 1, ]^2), "\n")
        cat("etaPhi[j, i + 1, ]", etaPhi[j, i + 1, ], "\n")
        stop("Sum of etaPhi squared is greater than 1")
      }
      
      
      # Sample standard deviation
      result <- sample_sigma(P[j, i + 1], etaPhi[j, i + 1, ], sigma_sqrd[j, i], mu_XB, Y, temp, eps_sigma_sqrd[j], ksi)
      sigma_sqrd[j, i + 1] <- result$next_val
      ac_rate_sigma_sqrd[j,i] <- (if (i == 1) result$accept else (ac_rate_sigma_sqrd[j,i-1]*(i-1)) + result$accept) / i
      if (check_adaptation(i, j, ac_rate_sigma_sqrd, adaptation_config)){
        eps_sigma_sqrd[j] = adapt_proposition(i=i,sigma_sqrd[j, i+1], sample_sigma, list(p_t=P[j,i+1],etaPhi_t=etaPhi[j, i + 1, ], sigma_sqrd_t=sigma_sqrd[j, i], mu_XB=mu_XB, Y=Y, temp=temp, eps=eps_sigma_sqrd[j]), "sigma_sqrd_t", adaptation_config, ksi = ksi)
      }
      
      
      # Sample beta
      result <- sample_beta(matrix(Beta_chains[j, i, ]), eps_beta[j], temp, P[j, i + 1], etaPhi[j, i + 1, ], sigma_sqrd[j, i + 1], beta_hat, beta_corr, X, XXt, Y)
      Beta_chains[j, i + 1, ] <- result$next_val
      ac_rate_beta[j,i] <- (if (i == 1) result$accept else (ac_rate_beta[j,i-1]*(i-1)) + result$accept) / i
      if (check_adaptation(i, j, ac_rate_beta, adaptation_config)){
        adapt_eps_beta <- TRUE
        ac_rate_window <- numeric(adaptation_config$ac_rate_span)
        beta_adapt[1,] <- result$next_val
        for (t in 1:(T_adapt-1)){
          result_adapt <- sample_beta(matrix(beta_adapt[t,]), eps_beta[j], temp, P[j, i + 1], etaPhi[j, i + 1, ], sigma_sqrd[j, i + 1], beta_hat, beta_corr, X, XXt, Y)
          beta_adapt[t+1,] <- result_adapt$next_val
          ac_rate_window[t %% adaptation_config$ac_rate_span + 1] <- result_adapt$accept
          if (t >= adaptation_config$ac_rate_span){
            current_ac_rate <- sum(ac_rate_window) / adaptation_config$ac_rate_span
            if (current_ac_rate < adaptation_config$ac_rate_lower) eps_beta[j] <- eps_beta[j]*1
            else if (current_ac_rate > adaptation_config$ac_rate_upper) eps_beta[j] <- eps_beta[j]*1
          }
        }
        adapt_eps_beta <- FALSE
      }
      
      
      
      # # Sample the latent variable using the multinomial
      # zeta <- sample(1:2, n_obs, prob = c(P[j, i + 1], 1 - P[j, i + 1]), replace = TRUE)
      # X1 <- X[zeta == 1, ]
      # X2 <- X[zeta == 2, ]
      # Y1 <- as.matrix(Y[zeta == 1])
      # Y2 <- as.matrix(Y[zeta == 2])
      # X1tX1 <- (t(X1) %*% X1)
      # f1 <- P[j, i + 1]/etaPhi[j, i + 1, 1]^2
      # X2tX2 <-  (t(X2) %*% X2)
      # f2 <- (1 - P[j, i + 1])/etaPhi[j, i + 1, 2]^2
      # 
      # cov_beta_1 <- (1/(2*sigma_sqrd[j, i])) * f1*X1tX1 + f2*X2tX2
      # cov_beta <- solve(cov_beta_1)
      # sq_p <- sqrt(P[j, i + 1])
      # sq_1_p <- sqrt(1 - P[j, i + 1])
      # y11 <- t(sqrt(sigma_sqrd[j, i]) * (-etaPhi[j, i + 1, 3]) * sq_1_p/sq_p * matrix(1,nrow=nrow(Y1),ncol=ncol(Y1)) - Y1)  %*% X1
      # y22 <- t(sqrt(sigma_sqrd[j, i]) * etaPhi[j, i + 1, 3] * sq_p/sq_1_p * matrix(1,nrow=nrow(Y2),ncol=ncol(Y2)) - Y2)  %*% X2
      # mu_bet <- t(-(1/sigma_sqrd[j, i]) * (y11 + y22) %*% cov_beta)
      # Beta_chains[j, i + 1, ] <- rmvnorm(1, mean = mu_bet, sigma = cov_beta)
      # 
      # # Sample standard deviation
      # beta <- as.matrix(Beta_chains[j, i+1, ])
      # sigma_sqrd[j, i + 1] <- rinvgamma(1, shape = n_obs/2, rate = t(beta - mu_bet) %*% cov_beta_1 %*% (beta - mu_bet))
    }
    
    
    
    
    # Swap states between chains
    if (i%%swap_config$cycle_period==0){
      swap_result <- swap_states(P[,i+1], etaPhi[,i+1,], sigma_sqrd[,i+1], Beta_chains[,i+1,], Tau, Alpha, ksi, beta_hat, X, Y, N_chains, Temperature)
      P[,i+1] <- swap_result$p_next
      etaPhi[,i+1,] <- swap_result$etaPhi_next
      sigma_sqrd[,i+1] <- swap_result$sigma_sqrd_next
      Beta_chains[,i+1,] <- swap_result$Beta_chains_next
    }
    
    
    
    chains_list <- list(
      P = P,
      etaPhi = etaPhi^2,
      sigma_sqrd = sigma_sqrd,
      Beta_chains = Beta_chains,
      N_chains = N_chains
    )
    
    gr_result <- compute_gelman_rubin(
      current_iter = i,
      chains = chains_list,
      check_cycle = conv_config$check_cycle,
      window_size = conv_config$window_size,
      start_iter = conv_config$start_iter
    )
    
    if (!is.null(gr_result)) {
      gelman[['iteration']] <- c(gelman[['iteration']], i)
      gelman[['max_psrf']] <- c(gelman[['max_psrf']], gr_result$max_psrf)
      gelman[['mpsrf']] <- c(gelman[['mpsrf']], gr_result$mpsrf)
    }
    update_progress_bar(pb)  # Update progress bar
  }
  
  #last convergence check
  chains_list <- list(
    P = P,
    etaPhi = etaPhi^2,
    sigma_sqrd = sigma_sqrd,
    Beta_chains = Beta_chains,
    N_chains = N_chains
  )
  
  gr_result <- compute_gelman_rubin(
    current_iter = T,
    chains = chains_list,
    check_cycle = conv_config$check_cycle,
    window_size = conv_config$window_size,
    start_iter = conv_config$start_iter
  )
  
  if (!is.null(gr_result)) {
    gelman[['iteration']] <- c(gelman[['iteration']], T)
    gelman[['max_psrf']] <- c(gelman[['max_psrf']], gr_result$max_psrf)
    gelman[['mpsrf']] <- c(gelman[['mpsrf']], gr_result$mpsrf)
  }
  
  
  # Return results
  return(list(P = P, etaPhi = etaPhi, sigma_sqrd = sigma_sqrd, Beta_chains = Beta_chains, 
              ac_rate_p = ac_rate_p, ac_rate_etaPhi_sqrd = ac_rate_etaPhi_sqrd, ac_rate_sigma_sqrd = ac_rate_sigma_sqrd, ac_rate_beta = ac_rate_beta, 
              N_chains = N_chains, Temperature = Temperature, gelman = gelman))
}