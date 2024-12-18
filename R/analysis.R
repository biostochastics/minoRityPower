#' Analyze Clinical Trials Data
#'
#' @description
#' Main function for processing and analyzing ClinicalTrials.gov data.
#' This function focuses on extracting and summarizing real-world trial data,
#' which can then be used to inform power analysis parameters.
#'
#' @importFrom dplyr %>% mutate filter group_by summarize n_distinct left_join rename
#' @param studies Data frame of ClinicalTrials.gov studies
#' @param conditions Data frame of trial conditions
#' @param facilities Data frame of facility information
#' @param min_duration Minimum trial duration in months
#' @param min_enrollment Minimum enrollment count
#' @return List containing processed trial data and comprehensive statistics
#' @seealso \code{\link{simple_power}} for using these results in power analysis
#' @export
analyze_trials <- function(studies,
                         conditions,
                         facilities,
                         min_duration = 12,
                         min_enrollment = 10) {

  # Initial validation
  if (!all(c("nct_id", "phase", "overall_status", "start_date", "completion_date") %in% names(studies))) {
    stop("Required columns missing from studies data frame")
  }
  if (!all(c("nct_id", "id") %in% names(facilities))) {
    stop("Required columns missing from facilities data frame")
  }

  # Process dates and durations
  studies <- studies %>%
    mutate(
      start_date = as.Date(start_date),
      completion_date = as.Date(completion_date)
    ) %>%
    filter(!is.na(start_date), !is.na(completion_date)) %>%
    mutate(
      duration = as.numeric(difftime(completion_date, start_date, units = "days")) / 30
    )

  # Filter trials
  completed_trials <- studies %>%
    filter(
      phase %in% c("Phase 1/Phase 2", "Phase 2",
                  "Phase 2/Phase 3", "Phase 3", "Phase 4"),
      overall_status == "Completed",
      enrollment_type == "Actual",
      study_type == "Interventional",
      duration >= min_duration,
      enrollment >= min_enrollment
    )

  # Join with facilities
  trials_with_facilities <- completed_trials %>%
    left_join(facilities, by = "nct_id") %>%
    rename(facility_id = id)

  # Calculate trial-level metrics
  trial_metrics <- trials_with_facilities %>%
    group_by(nct_id) %>%
    summarize(
      enrollment_rate = first(enrollment) / first(duration),
      n_facilities = n_distinct(facility_id, na.rm = TRUE),
      .groups = 'drop'
    )

  # Calculate comprehensive statistics
  stats <- list(
    trial_counts = list(
      total_trials = nrow(completed_trials),
      total_facilities = n_distinct(facilities$id),
      trials_per_phase = table(completed_trials$phase)
    ),

    enrollment = list(
      mean_enrollment = mean(completed_trials$enrollment, na.rm = TRUE),
      median_enrollment = median(completed_trials$enrollment, na.rm = TRUE),
      enrollment_quartiles = quantile(completed_trials$enrollment,
                                    probs = c(0.25, 0.5, 0.75),
                                    na.rm = TRUE),
      enrollment_by_phase = tapply(completed_trials$enrollment,
                                 completed_trials$phase,
                                 mean, na.rm = TRUE)
    ),

    duration = list(
      mean_duration = mean(completed_trials$duration, na.rm = TRUE),
      median_duration = median(completed_trials$duration, na.rm = TRUE),
      duration_quartiles = quantile(completed_trials$duration,
                                  probs = c(0.25, 0.5, 0.75),
                                  na.rm = TRUE)
    ),

    facilities = list(
      facilities_per_trial = table(trial_metrics$n_facilities),
      mean_facilities = mean(trial_metrics$n_facilities, na.rm = TRUE),
      median_facilities = median(trial_metrics$n_facilities, na.rm = TRUE)
    ),

    enrollment_rates = list(
      mean_rate = mean(trial_metrics$enrollment_rate, na.rm = TRUE),
      median_rate = median(trial_metrics$enrollment_rate, na.rm = TRUE),
      rate_sd = sd(trial_metrics$enrollment_rate, na.rm = TRUE),
      rate_quartiles = quantile(trial_metrics$enrollment_rate,
                              probs = c(0.25, 0.5, 0.75),
                              na.rm = TRUE)
    )
  )

  # Print key findings
  cat("\nTrial Analysis Summary:\n")
  cat("----------------------\n")
  cat(sprintf("Total Trials Analyzed: %d\n", stats$trial_counts$total_trials))
  cat(sprintf("Mean Enrollment Rate: %.1f patients/month\n",
              stats$enrollment_rates$mean_rate))
  cat(sprintf("Mean Duration: %.1f months\n", stats$duration$mean_duration))
  cat(sprintf("Average Facilities per Trial: %.1f\n",
              stats$facilities$mean_facilities))
  cat("\nPhase Distribution:\n")
  print(stats$trial_counts$trials_per_phase)

  return(list(
    trials = trials_with_facilities,
    metrics = trial_metrics,
    facilities = facilities,
    stats = stats
  ))
}

#' Estimate Parameters for Power Analysis
#'
#' @description
#' Analyzes enrollment patterns from real trial data to estimate parameters
#' needed for power analysis simulations. These parameters can be directly
#' used in \code{simple_power()} or \code{grid_bootstrap_power()}.
#'
#' @details
#' Estimates key parameters including:
#' \itemize{
#'   \item Mean enrollment rate
#'   \item Facility-level variance
#'   \item Trial-level variance
#'   \item Facility count distributions
#' }
#'
#' @param trial_data Output from analyze_trials function
#' @param min_rate Minimum acceptable enrollment rate
#' @param max_rate Maximum acceptable enrollment rate
#' @return List containing parameter estimates for power analysis
#' @seealso
#'   \code{\link{simple_power}} for trial-level power analysis
#'   \code{\link{grid_bootstrap_power}} for participant-level analysis
#' @export
#' @examples
#' \dontrun{
#' # Workflow:
#' # 1. Analyze real trial data
#' trial_results <- analyze_trials(studies, conditions, facilities)
#'
#' # 2. Estimate parameters
#' params <- estimate_power_parameters(trial_results)
#'
#' # 3. Use in power analysis
#' # For trial-level analysis:
#' power_results <- simple_power(
#'   effect_size = 1.25,
#'   N = 50,
#'   tau2_facility = params$parameters$log_rate_variance,
#'   mean_rate = params$parameters$mean_rate
#' )
#'
#' # For participant-level analysis:
#' bootstrap_results <- grid_bootstrap_power(
#'   effect_sizes1 = c(1.25, 1.5),
#'   tau2_trial = params$parameters$log_rate_variance
#' )
#' }
analyze_enrollment <- function(trial_data,
                                   min_rate = 1,
                                   max_rate = 200) {

  # Calculate trial-level enrollment rates
  enrollment_data <- trial_data$trials %>%
    group_by(nct_id) %>%
    summarize(
      enrollment_rate = first(enrollment) / first(duration),
      enrollment = first(enrollment),
      duration = first(duration),
      n_facilities = n_distinct(facility_id)
    ) %>%
    filter(
      enrollment_rate >= min_rate,
      enrollment_rate <= max_rate
    )

  cat(sprintf("\nAnalyzing %d trials with enrollment rates between %.1f and %.1f\n",
              nrow(enrollment_data), min_rate, max_rate))

  # Calculate overall variance in log enrollment rates
  log_rate_variance <- var(log(enrollment_data$enrollment_rate))

  params <- list(
    mean_rate = mean(enrollment_data$enrollment_rate, na.rm = TRUE),
    median_rate = median(enrollment_data$enrollment_rate, na.rm = TRUE),
    rate_variance = var(enrollment_data$enrollment_rate, na.rm = TRUE),
    log_rate_variance = log_rate_variance,
    facility_counts = table(enrollment_data$n_facilities)
  )

  cat(sprintf("\nEstimated parameters:\n"))
  cat(sprintf("Mean enrollment rate: %.2f\n", params$mean_rate))
  cat(sprintf("Overall log rate variance: %.2f\n", params$log_rate_variance))
  cat(sprintf("Median facilities per trial: %.1f\n",
              median(enrollment_data$n_facilities)))

  return(list(
    parameters = params,
    distribution = list(
      rate_quantiles = quantile(enrollment_data$enrollment_rate,
                              probs = c(0.25, 0.5, 0.75)),
      facility_distribution = table(enrollment_data$n_facilities)
    ),
    data = enrollment_data
  ))
}


