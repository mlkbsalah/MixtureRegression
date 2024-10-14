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

adapt_proposition <- function(last_result, sampling_fn, params, param_name, adaptation_config, ...) {
  # cat("Adapting >>>>>>>>>>>>>>>>>>>>\n")
  adapt_eps <- TRUE
  extra_args <- list(...)
  # for (name in names(params)) {
  #   if (name !="mu_XB" && name!="Y"){
  #     cat(name, " : ", params[[name]], "\n")
  #   }
  # }
  
  
  ac_rate_window <- numeric(adaptation_config$ac_rate_span)
  eps <- params$eps
  # Make the param_adapt an array of size len(last_result)*T_adapt
  param_adapt <- matrix(0, nrow = adaptation_config$T_adapt, ncol = length(last_result))
  # cat(dim(param_adapt), "\n")
  param_adapt[1, ] <- last_result
  
  # cat("Ready to adapt\n")
  for (t in 1:(adaptation_config$T_adapt - 1)) {
    params[[param_name]] <- param_adapt[t,]
    
    sampling_params <- c(params, extra_args)
    result_adapt <- do.call(sampling_fn, sampling_params)
    param_adapt[t + 1,] <- result_adapt$next_val  # Update with the next value of the adapted eps
    ac_rate_window[t %% adaptation_config$ac_rate_span + 1] <- result_adapt$accept
    
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
  return(params$eps)  # Ensure to return the correct value
}