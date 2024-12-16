#' Grid Bootstrap Power Analysis
#'
#' @description
#' Performs a grid search of bootstrap-based power analyses for patient-level data
#' across different parameter combinations and saves results to disk.
#'
#' @param effect_sizes1 Vector of main effect sizes to test
#' @param effect_sizes_int Vector of interaction effect sizes to test
#' @param Ns Vector of sample sizes (trials per arm) to test
#' @param tau2_trial Trial-level variance (default from empirical data)
#' @param K_trial Total number of original trials
#' @param R_boot Number of bootstrap iterations
#' @param n_cores Number of cores for parallel processing
#' @param output_dir Directory to save results (defaults to current directory)
#' @return Invisibly returns a list of all results
#' @export
grid_bootstrap_power <- function(effect_sizes1 = c(1.01, 1.05, 1.10, 1.25, 1.50, 2.00),
                               effect_sizes_int = effect_sizes1,
                               Ns = c(15, 30, 45, 60),
                               tau2_trial = 1.19,
                               K_trial = 360,
                               R_boot = 1000,
                               n_cores = parallel::detectCores() - 1,
                               output_dir = getwd()) {
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Initialize results storage
  all_results <- list()
  
  # Generate seeds for reproducibility
  set.seed(12345)
  seed_vec <- sample.int(1e6, R_boot)
  
  # Setup parallel backend
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterExport(cl, c("generate_interventions_data", 
                               "simulate_participants",
                               "assign_outcomes",
                               "sample_data_for_analysis",
                               "fit_models_and_get_pvalue"))
  
  # Run grid search
  for(es1 in effect_sizes1) {
    for(es_int in effect_sizes_int) {
      for(n in Ns) {
        # Run bootstrap simulations in parallel
        boot_results <- parallel::parLapply(cl, 1:R_boot, function(i) {
          bootstrap_power_sim(
            boot_index = i,
            effect_size1 = es1,
            effect_size_int = es_int,
            N = n,
            seed_vec = seed_vec,
            K_trial = K_trial,
            tau2_trial = tau2_trial
          )
        })
        
        # Extract p-values
        p_values <- sapply(boot_results, function(x) x$p_value)
        
        # Save detailed results
        filename <- file.path(output_dir, 
                            sprintf("bootstrap_p_values_%g_%g_%d.txt", es1, es_int, n))
        
        write.table(
          data.frame(
            p_value = p_values,
            N = rep(n, length(p_values))
          ),
          file = filename,
          row.names = FALSE,
          col.names = FALSE
        )
        
        # Calculate power
        power_05 <- mean(p_values < 0.05, na.rm = TRUE)
        power_01 <- mean(p_values < 0.01, na.rm = TRUE)
        power_005 <- mean(p_values < 0.005, na.rm = TRUE)
        
        all_results[[length(all_results) + 1]] <- list(
          effect_size1 = es1,
          effect_size_int = es_int,
          N = n,
          power_05 = power_05,
          power_01 = power_01,
          power_005 = power_005,
          convergence_rate = mean(!is.na(p_values))
        )
      }
    }
  }
  
  parallel::stopCluster(cl)
  
  # Save summary results
  summary_df <- do.call(rbind, lapply(all_results, as.data.frame))
  write.csv(
    summary_df,
    file = file.path(output_dir, "bootstrap_power_summary.csv"),
    row.names = FALSE
  )
  
  invisible(all_results)
}

#' Single Iteration Bootstrap Simulation
#' @keywords internal
#' @param boot_index Integer, bootstrap iteration index
#' @param effect_size1 Numeric, main effect size for Arm1
#' @param effect_size_int Numeric, interaction effect size
#' @param N Total number of trials selected 
#' @param seed_vec Vector of seeds for reproducibility
#' @param K_trial Total number of original trials
#' @param tau2_trial Variance of random intercept
#' @return A list with p_value for this iteration
bootstrap_power_sim <- function(boot_index,
                                effect_size1,
                                effect_size_int,
                                N,
                                seed_vec,
                                K_trial = 360,
                                tau2_trial = 1.19) {

  interventions <- generate_interventions_data(K_trial)

  # Generate participants
  all_records <- simulate_participants(interventions, tau2_trial, seed = seed_vec[boot_index])
  b0_trial <- attr(all_records, "b0_trial")

  # Assign outcomes
  all_records <- assign_outcomes(all_records, b0_trial, effect_size1, effect_size_int)

  # Sample data
  dat_rr <- sample_data_for_analysis(all_records, N, interventions, seed_vec, boot_index)

  # Fit models
  p_val <- fit_models_and_get_pvalue(dat_rr)

  return(list(p_value = p_val))
}
