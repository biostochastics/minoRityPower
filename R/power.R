# R/power.R

#' Generate Trial Data for Simulation
#' @keywords internal
generate_trial_data <- function(N, tau2_facility, effect_size, mean_rate) {
  # Generate facility random effects
  n_facilities <- N * 2  # Allow for multiple facilities per trial
  b0_facility <- rnorm(n_facilities, 0, sqrt(tau2_facility))

  # Create trial data frame
  trials <- data.frame(
    trial_id = rep(1:(2*N), each = 2),
    facility_id = 1:n_facilities,
    intervention = rep(rep(c(0,1), each = N), each = 2)
  )

  # Add enrollment rates
  trials$enrollment_rate <- rpois(
    n = nrow(trials),
    lambda = mean_rate * exp(b0_facility + log(effect_size) * trials$intervention)
  )

  # Ensure no zero rates
  trials$enrollment_rate[trials$enrollment_rate == 0] <- 1

  return(trials)
}

#' Power Analysis for Cancer Trial Enrollment Intervention
#'
#' @description
#' Simulates power for detecting effects of interventions on enrollment rates
#' in cancer clinical trials. Uses empirically-derived parameters from
#' ClinicalTrials.gov cancer trial data.
#'
#' @param effect_size Expected intervention effect size (rate ratio)
#' @param N Number of trials per arm
#' @param tau2_facility Facility-level variance (default from empirical data)
#' @param mean_rate Baseline enrollment rate (default from empirical data)
#' @param R_boot Number of bootstrap iterations
#' @return List containing power estimates and simulation details
#' @export
simple_power <- function(effect_size = 1.25,
                  N = 50,
                  tau2_facility = 0.64,
                  mean_rate = 5,
                  R_boot = 1000) {

  # Input validation
  if(N < 10) stop("N must be at least 10")
  if(mean_rate <= 0) stop("mean_rate must be positive")
  if(tau2_facility < 0) stop("tau2_facility must be non-negative")
  if(effect_size <= 0) stop("effect_size must be positive")

  # Initialize results storage
  results <- data.frame(
    iteration = 1:R_boot,
    converged = FALSE,
    p_value = NA_real_,
    estimate = NA_real_
  )

  for(i in 1:R_boot) {
    # Generate trial data
    sim_data <- generate_trial_data(
      N = N,
      tau2_facility = tau2_facility,
      effect_size = effect_size,
      mean_rate = mean_rate
    )

    # Try fitting model
    fit <- try({
      glmer.nb(
        enrollment_rate ~ intervention + (1 | facility_id),
        data = sim_data,
        control = glmerControl(optimizer = "bobyqa")
      )
    }, silent = TRUE)

    # Extract results if model converged
    if(!inherits(fit, "try-error")) {
      coef_summary <- summary(fit)$coefficients
      results$converged[i] <- TRUE
      results$p_value[i] <- coef_summary["intervention", "Pr(>|z|)"]
      results$estimate[i] <- coef_summary["intervention", "Estimate"]
    }
  }

  # Compute summaries
  converged <- results$converged
  power_estimate <- mean(results$p_value[converged] < 0.05, na.rm = TRUE)
  convergence_rate <- mean(converged)

  if(convergence_rate == 0) {
    warning("No simulations converged. Try adjusting parameters.")
  }

  return(list(
    power = power_estimate,
    convergence_rate = convergence_rate,
    results = results[converged, ],
    parameters = list(
      effect_size = effect_size,
      N = N,
      tau2_facility = tau2_facility,
      mean_rate = mean_rate
    )
  ))
}

#' Grid Power Analysis for Cancer Trial Enrollment Intervention
#'
#' @description
#' Performs a grid search of power analyses across different parameter combinations
#' and saves results to disk for later analysis and visualization.
#'
#' @param effect_sizes1 Vector of primary effect sizes to test
#' @param effect_sizes2 Vector of secondary effect sizes to test (for interaction)
#' @param Ns Vector of sample sizes (trials per arm) to test
#' @param tau2_facility Facility-level variance (default from empirical data)
#' @param mean_rate Baseline enrollment rate (default from empirical data)
#' @param R_boot Number of bootstrap iterations
#' @param output_dir Directory to save results (defaults to current directory)
#' @return Invisibly returns a list of all results
#' @export
grid_power <- function(effect_sizes1 = c(1.01, 1.05, 1.10, 1.25, 1.50, 2.00),
                      effect_sizes2 = effect_sizes1,
                      Ns = c(15, 30, 45, 60),
                      tau2_facility = 0.64,
                      mean_rate = 5,
                      R_boot = 1000,
                      output_dir = getwd()) {
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Initialize results storage
  all_results <- list()
  
  # Run grid search
  for(es1 in effect_sizes1) {
    for(es2 in effect_sizes2) {
      for(n in Ns) {
        # Run power analysis
        res <- power(
          effect_size = es1,
          N = n,
          tau2_facility = tau2_facility,
          mean_rate = mean_rate,
          R_boot = R_boot
        )
        
        # Save detailed results
        filename <- file.path(output_dir, 
                            sprintf("p_values_%g_%g_%d.txt", es1, es2, n))
        
        # Save p-values and sample sizes
        write.table(
          data.frame(
            p_value = res$results$p_value,
            N = rep(n, nrow(res$results))
          ),
          file = filename,
          row.names = FALSE,
          col.names = FALSE
        )
        
        all_results[[length(all_results) + 1]] <- list(
          effect_size1 = es1,
          effect_size2 = es2,
          N = n,
          results = res
        )
      }
    }
  }
  
  # Save summary results
  summary_df <- do.call(rbind, lapply(all_results, function(x) {
    data.frame(
      effect_size1 = x$effect_size1,
      effect_size2 = x$effect_size2,
      N = x$N,
      power = x$results$power,
      convergence_rate = x$results$convergence_rate
    )
  }))
  
  write.csv(
    summary_df,
    file = file.path(output_dir, "power_summary.csv"),
    row.names = FALSE
  )
  
  invisible(all_results)
}
