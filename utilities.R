library(ggplot2)
library(tidyr)
library(dplyr)

check_adaptation <- function(i, j, acceptance_rate , adaptation_config){
  ac_rate_lower <- adaptation_config$ac_rate_lower
  ac_rate_upper <- adaptation_config$ac_rate_upper
  cycle_period <- adaptation_config$cycle_period
  
  if ((i%%cycle_period==0) && ((mean(acceptance_rate[j,(i-cycle_period/2):i])>ac_rate_upper) || (mean(acceptance_rate[j,(i-cycle_period/2):i])<ac_rate_lower))) return(TRUE)
  else return(FALSE)
}

# Function to initialize the progress bar
initialize_progress_bar <- function(total_steps) {
  progress_bar$new(
    format = ":percent |:bar|  :current/:total  [:elapsed <:eta,  :tick_rateit/s]",
    total = total_steps,
    complete = "█",   # Completion bar character
    incomplete = "░", # Incomplete bar character
    current = "▒",    # Current bar character
    clear = FALSE,    # If TRUE, clears the bar when finished
    width = 100       # Width of the progress bar
  )
}

# Function to update the progress bar
update_progress_bar <- function(pb) {
  pb$tick()
}

adapt_proposition <- function(i=NA, last_result, sampling_fn, params, param_name, adaptation_config, ...) {
  adapt_eps <- TRUE
  extra_args <- list(...)

  
  ac_rate_window <- numeric(adaptation_config$ac_rate_span)
  eps <- params$eps
  # Make the param_adapt an array of size len(last_result)*T_adapt
  param_adapt <- matrix(0, nrow = adaptation_config$T_adapt, ncol = length(last_result))
  # cat(dim(param_adapt), "\n")
  param_adapt[1, ] <- last_result
  
  # cat("Ready to adapt\n")
  for (t in 1:(adaptation_config$T_adapt - 1)) {
    # params[[param_name]] <- param_adapt[t,]
    # # if ((param_name == "sigma_sqrd_t") && (i > 5100)){
    # #   cat("Adapting at iteration : ", i, " and t: ", t, "\n")
    # #   print(params)
    # # }
    sampling_params <- c(params, extra_args)
    result_adapt <- do.call(sampling_fn, sampling_params)
    param_adapt[t + 1,] <- result_adapt$next_val  # Update with the next value of the adapted eps
    ac_rate_window[t %% adaptation_config$ac_rate_span + 1] <- result_adapt$accept
    # if ((param_name == "sigma_sqrd_t") && (i > 5100)){
    #   cat("Adapting at iteration : ", i, " and t: ", t, "\n")
    #   cat("--> Current accepted: ", result_adapt$accept, "\n")
    #   cat("--> Next param value: ", result_adapt$next_val, "\n")
    #   cat("--> ac_rate_window : ", ac_rate_window, "\n")
    # }
    
    if (t >= adaptation_config$ac_rate_span) {  # Check if enough samples were accumulated
      current_ac_rate <- sum(ac_rate_window) / adaptation_config$ac_rate_span
      if (current_ac_rate < adaptation_config$ac_rate_lower) {
        eps <- eps * 1.1
      } else if (current_ac_rate > adaptation_config$ac_rate_upper) {
        eps <- eps * 0.9
      }
    }
    params$eps <- eps
  }
  
  adapt_eps <- FALSE
  # cat("Adapted <<<<<<<<<<<<<<<<<<<\n")
  return(params$eps)  # Ensure to return the correct value
}

mcmc_from_array <- function(P, etaPhi, sigma_sqrd, Beta_chains, N_chains) {
  chain_list <- list()
  for (i in 1:N_chains) {
    chain_matrix <- cbind(
      P[i,],
      etaPhi[i,,]^2,
      sigma_sqrd[i,],
      Beta_chains[i,,]
    )

    chain_list[[i]] <- coda::mcmc(chain_matrix)
  }
  mcmc_list <- coda::mcmc.list(chain_list)
  return(mcmc_list)
}

compute_gelman_rubin <- function(current_iter, chains, check_cycle, window_size, start_iter) {
  
  if (current_iter < start_iter || current_iter %% check_cycle != 0) {
    return(NULL)
  }
  
  start_index <- max(1, current_iter - window_size + 1)
  chain_list <- list()
  
  for (chain_i in 1:chains$N_chains) {
    # Extract parameter values for this chain
    chain_matrix <- cbind(
      chains$P[chain_i, start_index:current_iter],
      chains$etaPhi[chain_i, start_index:current_iter, ],
      sqrt(chains$sigma_sqrd[chain_i, start_index:current_iter]),
      chains$Beta_chains[chain_i, start_index:current_iter, ]
    )
    
    # Convert to mcmc object
    chain_list[[chain_i]] <- coda::mcmc(chain_matrix)
  }
  
  # Calculate Gelman-Rubin diagnostic
  mcmc_list <- coda::mcmc.list(chain_list)
  gelman_result <- coda::gelman.diag(mcmc_list, multivariate = TRUE, autoburnin = TRUE, transform = TRUE)
  
  return(list(
    max_psrf = max(gelman_result$psrf[, "Point est."]),
    mpsrf = gelman_result$mpsrf, 
    iteration = current_iter
  ))
}

ggplot_overlapping_histograms <- function(matrix,
                                          color_palette = "Spectral",
                                          n_breaks = 30,
                                          title = "Overlapping Histograms",
                                          alpha = 0.5) {
  
  # Convert matrix to long format
  df <- as.data.frame(t(matrix)) %>%
    gather(key = "Row", value = "Value") %>%
    mutate(Row = factor(Row, levels = unique(Row)))
  
  ggplot(df, aes(x = Value, fill = Row)) +
    geom_histogram(position = "identity",
                   alpha = alpha,
                   bins = n_breaks) +
    scale_fill_brewer(palette = color_palette) +
    labs(title = title,
         x = "Value",
         y = "Count") +
    theme_minimal() +
    theme(legend.position = "none")  # Remove legend
}

ggplot_overlapping_lines <- function(matrix,
                                     color_palette = "Spectral",
                                     title = "Overlapping Lines") {
  
  # Convert matrix to long format with index
  df <- as.data.frame(t(matrix)) %>%
    gather(key = "Row", value = "Value") %>%
    group_by(Row) %>%
    mutate(Index = row_number()) %>%
    ungroup() %>%
    mutate(Row = factor(Row, levels = unique(Row)))
  
  ggplot(df, aes(x = Index, y = Value, color = Row)) +
    geom_line() +
    scale_color_brewer(palette = color_palette) +
    labs(title = title,
         x = "Index",
         y = "Value") +
    theme_minimal() +
    theme(legend.position = "none")  # Remove legend
}


