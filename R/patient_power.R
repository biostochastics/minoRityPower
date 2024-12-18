# Load required packages
library(dplyr)
library(ggplot2)
library(lme4)
library(reshape2)
library(ggsci)
library(parallel)

#' Generate Interventions Data
#' @keywords internal
generate_interventions_data <- function(K_trial) {
  interventions <- data.frame(
    trial_id = 1:K_trial,
    x1 = factor(rep(c("Arm0", "Arm1"), times = c(floor(K_trial / 2),
                                                 K_trial - floor(K_trial / 2))))
  )
  return(interventions)
}

#' Simulate Participant-Level Data
#' @keywords internal
simulate_participants <- function(interventions, tau2_trial, minority_rate = 0.2, seed) {
  set.seed(seed)
  K_trial <- nrow(interventions)
  b0_trial <- rnorm(n = K_trial, mean = 0, sd = sqrt(tau2_trial))
  list_records <- vector("list", K_trial)
  for (i in seq_len(K_trial)) {
    eligible <- max(rnorm(1, mean = 1000, sd = 300), 100)
    minority <- rbinom(eligible, size = 1, prob = minority_rate)
    records <- data.frame(
      trial_id = rep(i, eligible),
      x1 = rep(interventions$x1[i], eligible),
      eligible_index = seq_len(eligible),
      minority = minority
    )
    list_records[[i]] <- records
  }
  all_records <- do.call(rbind, list_records)
  attr(all_records, "b0_trial") <- b0_trial
  return(all_records)
}

#' Assign Outcomes
#' @keywords internal
assign_outcomes <- function(all_records, b0_trial, effect_size1, effect_size_int) {
  baseline_rate <- ifelse(all_records$minority == 1, 0.01, 0.08)
  log_odds <- log(baseline_rate) +
    b0_trial[all_records$trial_id] +
    ifelse(all_records$x1 == "Arm1", log(effect_size1), 0) +
    ifelse(all_records$x1 == "Arm1" & all_records$minority == 1, log(effect_size_int), 0)
  prob <- plogis(log_odds)
  y <- rbinom(n = nrow(all_records), size = 1, prob = prob)
  all_records$y <- y
  return(all_records)
}

#' Sample Participants and Trials
#' @keywords internal
sample_data_for_analysis <- function(all_records, N, interventions, seed_vec, boot_index) {
  dat_sub <- all_records %>%
    group_by(trial_id) %>%
    sample_frac(0.75) %>%
    ungroup()

  N2 <- N / 2
  set.seed(seed_vec[boot_index])

  sampled_trials <- c(
    sample(interventions %>% filter(x1 == "Arm0") %>% pull(trial_id),
           size = N2, replace = TRUE),
    sample(interventions %>% filter(x1 == "Arm1") %>% pull(trial_id),
           size = N2, replace = TRUE)
  )

  dat_rr <- dat_sub[dat_sub$trial_id %in% sampled_trials, ]
  dat_rr$trial_id <- factor(dat_rr$trial_id)
  return(dat_rr)
}

#' Fit Models and Extract P-values
#' @keywords internal
fit_models_and_get_pvalue <- function(dat_rr) {
  tryCatch({
    fit_full <- glmer(y ~ x1 * minority + (1 | trial_id),
                      family = binomial, data = dat_rr)
    fit_null <- glmer(y ~ x1 + minority + (1 | trial_id),
                      family = binomial, data = dat_rr)

    if (isSingular(fit_full) || AIC(fit_full) == Inf) {
      return(NA_real_)
    }

    anova_res <- anova(fit_null, fit_full)
    return(anova_res$`Pr(>Chisq)`[2])
  }, error = function(e) {
    return(NA_real_)
  })
}

#' Single Bootstrap Simulation
#' @keywords internal
bootstrap_power_sim <- function(boot_index,
                                effect_size1,
                                effect_size2,
                                N,
                                seed_vec,
                                K_trial,
                                tau2_trial,
                                minority_rate) {

  interventions <- generate_interventions_data(K_trial)
  all_records <- simulate_participants(interventions, tau2_trial,
                                       minority_rate = minority_rate,
                                       seed = seed_vec[boot_index])
  b0_trial <- attr(all_records, "b0_trial")
  all_records <- assign_outcomes(all_records, b0_trial, effect_size1, effect_size2)
  dat_rr <- sample_data_for_analysis(all_records, N, interventions, seed_vec, boot_index)
  p_val <- fit_models_and_get_pvalue(dat_rr)

  return(list(p_value = p_val))
}

#' Plot Power Analysis Results
#' @keywords internal
plot_power_results <- function(results_dir,
                               analysis_type = c("trial", "participant"),
                               alpha_levels = c(0.05, 0.01, 0.005)) {

  analysis_type <- match.arg(analysis_type)

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

  calc_power <- function(filename) {
    datme <- read.table(filename) %>% filter(!is.na(V1))

    if(nrow(datme) == 0) {
      warning("No valid data in file: ", filename)
      return(NULL)
    }

    parameters <- basename(filename) %>%
      gsub(ifelse(analysis_type == "trial", "p_values_", "bootstrap_p_values_"), "", .) %>%
      gsub(".txt", "", .) %>%
      strsplit("_", fixed = TRUE) %>%
      unlist()

    effect_size1 <- as.numeric(parameters[1])
    effect_size2 <- as.numeric(parameters[2])
    N <- as.numeric(parameters[3])

    power_05 <- mean(datme$V1 < 0.05, na.rm = TRUE)
    power_01 <- mean(datme$V1 < 0.01, na.rm = TRUE)
    power_005 <- mean(datme$V1 < 0.005, na.rm = TRUE)

    return(data.frame(
      effect_size1 = effect_size1,
      effect_size2 = effect_size2,
      N = N,
      power_05 = power_05,
      power_01 = power_01,
      power_005 = power_005
    ))
  }

  results_list <- lapply(filenames, calc_power)
  results_list <- results_list[!sapply(results_list, is.null)]

  if(length(results_list) == 0) {
    warning("No valid results after processing")
    return(NULL)
  }

  df <- do.call(rbind, results_list)
  df <- df %>%
    mutate(
      effect_size1_label = sprintf("Main Effect: %+.0f%%", (effect_size1 - 1) * 100),
      effect_size2_label = sprintf("Minority Effect: %+.0f%%", (effect_size2 - 1) * 100)
    ) %>%
    arrange(effect_size1, effect_size2)

  df$effect_size1_label <- factor(df$effect_size1_label,
                                  levels = unique(df$effect_size1_label))
  df$effect_size2_label <- factor(df$effect_size2_label,
                                  levels = unique(df$effect_size2_label))

  df_melt <- reshape2::melt(df,
                            id.vars = c("effect_size1", "effect_size2", "N",
                                        "effect_size1_label", "effect_size2_label"),
                            variable.name = "significance_level",
                            value.name = "power")

  df_melt$significance_level <- factor(df_melt$significance_level,
                                       levels = c("power_05", "power_01", "power_005"),
                                       labels = c("α = 0.05", "α = 0.01", "α = 0.005"))

  title_prefix <- if(analysis_type == "trial") {
    "Trial-level Analysis"
  } else {
    "Participant-level Analysis"
  }

  x_label <- if(analysis_type == "trial") {
    "Number of Trials per Arm"
  } else {
    "Number of Participants per Trial"
  }

  main_plot <- ggplot(df_melt, aes(x = N, y = power, color = significance_level)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = 0.8, linetype = 2, linewidth = 0.5,
               color = 'black', alpha = 0.7) +
    coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    scale_y_continuous(breaks = seq(0, 1, 0.2),
                       labels = scales::percent,
                       minor_breaks = seq(0.1, 0.9, 0.1)) +
    labs(title = sprintf('%s: Power Analysis Results', title_prefix),
         subtitle = "Interaction Effect Power by Effect Size Combinations",
         x = x_label,
         y = "Statistical Power") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 10),
      text = element_text(size = 9),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8),
      strip.text = element_text(size = 9),
      strip.background = element_rect(fill = "grey95"),
      panel.grid.minor = element_line(color = "grey90"),
      panel.grid.major = element_line(color = "grey85"),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(t = 0)
    ) +
    facet_grid(rows = vars(effect_size2_label),
               cols = vars(effect_size1_label)) +
    scale_color_jco() +
    guides(color = guide_legend(
      title = "Significance Level",
      nrow = 1,
      override.aes = list(size = 3, linewidth = 1.5)
    ))

  return(main_plot)
}

#' Run Participant-Level Power Analysis
#'
#' @description
#' Performs power analysis for detecting intervention effects on minority enrollment
#' in clinical trials using participant-level simulations. Developed for ARPA-H
#' program evaluation of healthcare system interventions.
#'
#' @details
#' Uses mixed-effects logistic regression models to analyze:
#' \itemize{
#'   \item Main intervention effects on enrollment probability
#'   \item Interaction effects between intervention and minority status
#'   \item Clustering of participants within facilities
#' }
#'
#' @param effect_sizes1 Numeric vector of main intervention effect sizes (rate ratios)
#' @param effect_sizes2 Numeric vector of minority interaction effect sizes (rate ratios)
#' @param Ns Vector of sample sizes (number of facilities per arm)
#' @param tau2_trial Trial-level variance component (default: 1.19)
#' @param K_trial Total number of trials to simulate (default: 360)
#' @param minority_rate Expected proportion of minority participants (default: 0.2)
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
#' results <- run_power_participant(
#'   effect_sizes1 = c(1.25, 1.50),
#'   effect_sizes2 = c(1.75, 2.00),
#'   Ns = c(30, 45),
#'   minority_rate = 0.3,
#'   R_boot = 100
#' )
#'
#' # View results
#' print(results$power_table)
#' print(results$power_plot)
#' }
#'
#' @seealso \code{\link{run_power_facility}} for facility-level analysis
#' @export
#'
run_power_participant <- function(
    effect_sizes1 = c(1.01, 1.05, 1.10, 1.25, 1.50, 2.00),
    effect_sizes2 = effect_sizes1,
    Ns = c(15, 30, 45, 60),
    tau2_trial = 1.19,
    K_trial = 360,
    minority_rate = 0.2,
    R_boot = 1000,
    n_cores = parallel::detectCores() - 1,
    results_dir = file.path(getwd(), "power_analysis_results"),
    seed = NULL) {

  # Input validation
  if (!is.numeric(effect_sizes1) || any(effect_sizes1 <= 0)) {
    stop("effect_sizes1 must be positive numeric values")
  }
  if (!is.numeric(effect_sizes2) || any(effect_sizes2 <= 0)) {
    stop("effect_sizes2 must be positive numeric values")
  }
  if (!is.numeric(Ns) || any(Ns < 10)) {
    stop("Ns must be numeric values >= 10")
  }
  if (tau2_trial < 0) {
    stop("tau2_trial must be non-negative")
  }
  if (K_trial < max(Ns)) {
    stop("K_trial must be at least as large as maximum N")
  }
  if (minority_rate <= 0 || minority_rate >= 1) {
    stop("minority_rate must be between 0 and 1")
  }

  # Set random seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Create results directory
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

  # Setup cluster
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl))

  # Generate seeds for reproducibility
  set.seed(12345)
  seed_vec <- sample.int(1e6, R_boot)

  # Export necessary functions and objects to cluster
  parallel::clusterExport(cl, c(
    "generate_interventions_data",
    "simulate_participants",
    "run_power_participant",
    "assign_outcomes",
    "sample_data_for_analysis",
    "fit_models_and_get_pvalue",
    "bootstrap_power_sim",
    "seed_vec",
    "minority_rate",
    "tau2_trial",
    "K_trial"
  ), envir = environment())

  # Load required packages on cluster nodes
  parallel::clusterEvalQ(cl, {
    library(dplyr)
    library(lme4)
  })

  # Initialize results storage
  message("Starting grid bootstrap power analysis...")
  all_results <- list()

  # Run grid search
  for(es1 in effect_sizes1) {
    for(es2 in effect_sizes2) {
      message(sprintf("Processing effect sizes: main = %g, minority = %g", es1, es2))
      for(n in Ns) {
        message(sprintf("  Sample size: N = %d", n))

        # Create simulation environment for parallel execution
        sim_env <- environment()
        sim_env$es1 <- es1
        sim_env$es2 <- es2
        sim_env$n <- n

        parallel::clusterExport(cl, c("es1", "es2", "n"), envir = sim_env)

        # Run bootstrap simulations in parallel
        boot_results <- parallel::parLapply(cl, 1:R_boot, function(i) {
          bootstrap_power_sim(
            boot_index = i,
            effect_size1 = es1,
            effect_size2 = es2,
            N = n,
            seed_vec = seed_vec,
            K_trial = K_trial,
            tau2_trial = tau2_trial,
            minority_rate = minority_rate
          )
        })

        # Extract p-values
        p_values <- sapply(boot_results, function(x) x$p_value)

        # Save detailed results
        filename <- file.path(results_dir,
                              sprintf("bootstrap_p_values_%g_%g_%d.txt", es1, es2, n))

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
        valid_results <- !is.na(p_values)
        power_05 <- mean(p_values[valid_results] < 0.05, na.rm = TRUE)
        power_01 <- mean(p_values[valid_results] < 0.01, na.rm = TRUE)
        power_005 <- mean(p_values[valid_results] < 0.005, na.rm = TRUE)

        all_results[[length(all_results) + 1]] <- list(
          effect_size1 = es1,
          effect_size2 = es2,
          N = n,
          power_05 = power_05,
          power_01 = power_01,
          power_005 = power_005,
          convergence_rate = mean(valid_results)
        )
      }
    }
  }

  # Create summary dataframe
  summary_df <- do.call(rbind, lapply(all_results, as.data.frame))
  write.csv(
    summary_df,
    file = file.path(results_dir, "bootstrap_power_summary.csv"),
    row.names = FALSE
  )

  # Create power table
  power_table <- summary_df %>%
    select(effect_size1, effect_size2, N, power_05) %>%
    tidyr::pivot_wider(
      names_from = N,
      values_from = power_05,
      names_prefix = "N="
    ) %>%
    mutate(
      effect_size1_pct = paste0((effect_size1 - 1) * 100, "%"),
      effect_size2_pct = paste0((effect_size2 - 1) * 100, "%")
    ) %>%
    select(effect_size1_pct, effect_size2_pct, starts_with("N=")) %>%
    arrange(effect_size1_pct, effect_size2_pct)

  # Generate power plot
  message("Generating power analysis plot...")
  power_plot <- plot_power_results(
    results_dir = results_dir,
    analysis_type = "participant",
    alpha_levels = c(0.05, 0.01, 0.005)
  )

  # Save plot if successfully created
  if (!is.null(power_plot)) {
    ggsave(
      filename = file.path(results_dir, "power_analysis_plot.pdf"),
      plot = power_plot,
      width = 10,
      height = 8,
      units = "in",
      device = "pdf"
    )
  }

  message("Analysis complete. Results saved in: ", results_dir)
  return(list(
    power_table = power_table,
    power_plot = power_plot,
    full_results = all_results
  ))
}
