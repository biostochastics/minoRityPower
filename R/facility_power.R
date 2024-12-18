#' @import dplyr
#' @import ggplot2
#' @import ggsci
#' @import reshape2
#' @import lme4
NULL

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
#' ClinicalTrials.gov cancer trial data. Uses a literature-based estimate for tau between facilities
#'
#' @param effect_size Expected intervention effect size (rate ratio)
#' @param N Number of trials per arm
#' @param tau2_facility Facility-level variance
#' @param mean_rate Baseline enrollment rate (default from empirical data)
#' @param R_boot Number of bootstrap iterations
#' @return List containing power estimates and simulation details
#' @export
simple_power <- function(effect_size = 1.25,
                         N = 50,
                         tau2_facility = 0.64,
                         mean_rate = 11,
                         R_boot = 100) {

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
#' @param tau2_facility Facility-level variance
#' @param mean_rate Baseline enrollment rate (default from empirical data)
#' @param R_boot Number of bootstrap iterations
#' @param output_dir Directory to save results (defaults to current directory)
#' @return Invisibly returns a list of all results
#' @export
grid_power <- function(effect_sizes1 = c(1.01, 1.05, 1.10, 1.25, 1.50, 2.00),
                       effect_sizes2 = effect_sizes1,
                       Ns = c(15, 30, 45, 60),
                       tau2_facility = 0.64,
                       mean_rate = 11,
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
        res <- simple_power(
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

#' Plot Power Analysis Results
#'
#' @description
#' Creates visualizations of power analysis results showing power to detect
#' intervention effects across different effect sizes and sample sizes.
#'
#' @param results_dir Directory containing power analysis results
#' @param analysis_type Character, either "trial" or "participant"
#' @param alpha_levels Numeric vector of significance levels to analyze
#'
#' @return A ggplot object visualizing power analysis results
#'
#' @examples
#' \dontrun{
#' plot <- plot_power_results(
#'   "power_analysis_results",
#'   analysis_type = "participant",
#'   alpha_levels = c(0.05, 0.01, 0.005)
#' )
#' print(plot)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point theme_bw
#' @importFrom dplyr filter mutate select arrange
#' @importFrom reshape2 melt
#' @export
plot_power_results <- function(results_dir,
                               analysis_type = c("trial", "participant"),
                               alpha_levels = c(0.05, 0.01, 0.005)) {

  analysis_type <- match.arg(analysis_type)

  # Find all result files based on analysis type
  pattern <- if(analysis_type == "trial") {
    "p_values_.*\\.txt$"
  } else {
    "bootstrap_p_values_.*\\.txt$"
  }

  filenames <- list.files(results_dir, pattern = pattern, full.names = TRUE)

  if(length(filenames) == 0) {
    warning("No power analysis results found in directory: ", results_dir)
    return(NULL)
  }

  calc_avg_pval <- function(filename) {
    tryCatch({
      datme <- read.table(filename)
      if(nrow(datme) == 0) {
        warning("Empty data file: ", filename)
        return(NULL)
      }
      datme <- datme %>% filter(!is.na(V1))

      # Extract parameters from filename
      parameters <- basename(filename) %>%
        gsub(ifelse(analysis_type == "trial", "p_values_", "bootstrap_p_values_"), "", .) %>%
        gsub(".txt", "", .) %>%
        strsplit("_", fixed = TRUE) %>%
        unlist()

      if(length(parameters) < 3) {
        warning("Invalid filename format: ", filename)
        return(NULL)
      }

      effect_size1 <- as.numeric(parameters[1])
      effect_size2 <- as.numeric(parameters[2])
      N <- as.numeric(parameters[3])

      # Calculate p-value proportions with error checking
      p_values_05 <- mean(datme$V1 < 0.05, na.rm = TRUE)
      p_values_01 <- mean(datme$V1 < 0.01, na.rm = TRUE)
      p_values_005 <- mean(datme$V1 < 0.005, na.rm = TRUE)

      data.frame(
        effect_size1 = effect_size1,
        effect_size2 = effect_size2,
        N = N,
        avg_p_val_05 = p_values_05,
        avg_p_val_01 = p_values_01,
        avg_p_val_005 = p_values_005
      )
    }, error = function(e) {
      warning("Error processing file ", filename, ": ", e$message)
      return(NULL)
    })
  }

  # Apply the function to each file and bind the results into a data frame
  results_list <- lapply(filenames, calc_avg_pval)
  results_list <- results_list[!sapply(results_list, is.null)]

  if(length(results_list) == 0) {
    warning("No valid results found after processing files")
    return(NULL)
  }

  results <- do.call(rbind, results_list)

  # Create plot data
  df <- results
  df$effect_size1_percent <- paste0((df$effect_size1 - 1) * 100, "%")
  df$effect_size2_percent <- paste0((df$effect_size2 - 1) * 100, "%")
  df <- df %>% filter(effect_size2 < 4)

  if(nrow(df) == 0) {
    warning("No data left after filtering")
    return(NULL)
  }

  df$effect_size1_percent <- factor(df$effect_size1_percent,
                                    levels = unique(df$effect_size1_percent)[order(as.numeric(sub("%", "", unique(df$effect_size1_percent))))])
  df$effect_size2_percent <- factor(df$effect_size2_percent,
                                    levels = unique(df$effect_size2_percent)[order(as.numeric(sub("%", "", unique(df$effect_size2_percent))))])

  # Melt the dataframe for plotting
  df_melt <- reshape2::melt(df,
                            id.vars = c("effect_size1", "effect_size2", "N",
                                        "effect_size1_percent", "effect_size2_percent"),
                            variable.name = "significance_level",
                            value.name = "avg_p_val")

  # Convert significance_level to factor with nice labels
  df_melt$significance_level <- factor(df_melt$significance_level,
                                       levels = c("avg_p_val_05", "avg_p_val_01", "avg_p_val_005"),
                                       labels = c("p<0.05", "p<0.01", "p<0.005"))

  # Create main plot
  title_prefix <- if(analysis_type == "trial") {
    "Trial-level Analysis"
  } else {
    "Participant-level Analysis"
  }

  x_label <- if(analysis_type == "trial") {
    "N trials per arm"
  } else {
    "N participants per trial"
  }

  main_plot <- ggplot(df_melt, aes(x = N, y = avg_p_val, color = significance_level)) +
    geom_line() +
    geom_point(alpha = 0.9, size = 0.85) +
    geom_hline(yintercept = 0.8, linetype = 2, linewidth = 0.5, color = 'black', alpha = 0.7) +
    coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    labs(title = sprintf('%s: Power to Detect CARE Intervention Effects', title_prefix),
         x = x_label,
         y = "Power") +
    theme_bw() +
    theme(
      title = element_text(size = 8),
      text = element_text(size = 8),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 6),
      strip.text = element_text(size = 8)
    ) +
    facet_grid(rows = vars(effect_size2_percent),
               cols = vars(effect_size1_percent),
               scales = "free",
               labeller = label_both) +
    scale_color_jco() +
    guides(color = guide_legend(
      title = "Type I error rate (P)",
      override.aes = list(size = 5, shape = 15),
      keywidth = 1,
      keyheight = 1
    ))

  return(main_plot)
}

#' Run Facility-Level Power Analysis
#'
#' @description
#' Performs power analysis for detecting intervention effects on facility-level
#' enrollment rates in clinical trials. Developed for ARPA-H program evaluation
#' of healthcare system interventions.
#'
#' @details
#' Uses negative binomial mixed-effects models to analyze:
#' \itemize{
#'   \item Intervention effects on enrollment rates
#'   \item Facility-level random effects
#'   \item Overdispersion in enrollment counts
#' }
#'
#' @param effect_sizes Numeric vector of intervention effect sizes (rate ratios)
#' @param Ns Vector of sample sizes (number of facilities per arm)
#' @param tau2_facility Facility-level variance component (default: 0.64)
#' @param mean_rate Mean monthly enrollment rate (default: 11)
#' @param R_boot Number of bootstrap iterations (default: 1000)
#' @param n_cores Number of cores for parallel processing (default: all cores - 1)
#' @param results_dir Directory to save results (default: "power_analysis_results")
#' @param seed Random seed for reproducibility (default: NULL)
#'
#' @return A list containing:
#' \describe{
#'   \item{power_table}{Data frame of power estimates by effect size and sample size}
#'   \item{power_plot}{ggplot object visualizing power analysis results}
#'   \item{full_results}{Complete simulation results and convergence information}
#' }
#'
#' @examples
#' \dontrun{
#' # Basic analysis
#' results <- run_power_facility(
#'   effect_sizes = c(1.25, 1.50),
#'   Ns = c(30, 45),
#'   tau2_facility = 0.64,
#'   R_boot = 100
#' )
#'
#' # View results
#' print(results$power_table)
#' print(results$power_plot)
#' }
#'
#' @seealso \code{\link{run_power_participant}} for participant-level analysis
#' @export
run_power_facility <- function(
    effect_sizes = c(1.01, 1.05, 1.10, 1.25, 1.50, 2.00),
    Ns = c(15, 30, 45, 60),
    tau2_facility = 0.64,
    mean_rate = 5,
    R_boot = 1000,
    results_dir = file.path(getwd(), "power_analysis_results"),
    seed = NULL) {

  # Input validation
  if (!is.numeric(effect_sizes) || any(effect_sizes <= 0)) {
    stop("effect_sizes must be positive numeric values")
  }
  if (!is.numeric(Ns) || any(Ns < 10)) {
    stop("Ns must be numeric values >= 10")
  }
  if (tau2_facility < 0) {
    stop("tau2_facility must be non-negative")
  }
  if (mean_rate <= 0) {
    stop("mean_rate must be positive")
  }

  # Set random seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Create results directory
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

  # Run grid power analysis
  full_results <- grid_power(
    effect_sizes1 = effect_sizes,
    effect_sizes2 = effect_sizes,
    Ns = Ns,
    tau2_facility = tau2_facility,
    mean_rate = mean_rate,
    R_boot = R_boot,
    output_dir = results_dir
  )

  # Extract power summary directly from full_results
  power_summary <- do.call(rbind, lapply(full_results, function(x) {
    data.frame(
      effect_size1 = x$effect_size1,
      N = x$N,
      power = x$results$power
    )
  }))

  # Create formatted power table
  power_table <- tryCatch({
    power_summary %>%
      tidyr::pivot_wider(
        names_from = N,
        values_from = power,
        names_prefix = "N="
      ) %>%
      mutate(
        effect_size_pct = paste0((effect_size1 - 1) * 100, "%")
      ) %>%
      select(effect_size_pct, starts_with("N=")) %>%
      arrange(effect_size_pct)
  }, error = function(e) {
    warning("Error creating power table: ", e$message)
    return(NULL)
  })

  # Write power table to CSV if successfully created
  if (!is.null(power_table)) {
    tryCatch({
      write.csv(
        power_table,
        file = file.path(results_dir, "power_table_summary.csv"),
        row.names = FALSE
      )
    }, error = function(e) {
      warning("Error writing power table to CSV: ", e$message)
    })
  }

  # Generate power plot
  power_plot <- tryCatch({
    plot_power_results(
      results_dir = results_dir,
      analysis_type = "trial",
      alpha_levels = c(0.05, 0.01, 0.005)
    )
  }, error = function(e) {
    warning("Error creating power plot: ", e$message)
    return(NULL)
  })

  # Save plot if successfully created
  if (!is.null(power_plot)) {
    tryCatch({
      ggsave(
        filename = file.path(results_dir, "power_analysis_plot.pdf"),
        plot = power_plot,
        width = 10,
        height = 8,
        units = "in",
        device = "pdf"
      )
    }, error = function(e) {
      warning("Error saving power plot: ", e$message)
    })
  }

  # Return results
  return(list(
    power_table = power_table,
    power_plot = power_plot,
    full_results = full_results
  ))
}
