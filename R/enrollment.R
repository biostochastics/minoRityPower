#' Analyze Enrollment Patterns and Estimate Parameters
#'
#' @description
#' Analyzes enrollment patterns and estimates parameters for power analysis
#' using trial-level data from ClinicalTrials.gov. Note that facility-level
#' variance cannot be estimated directly as enrollment is only reported at
#' trial level.
#'
#' @param trial_data Output from analyze_trials function
#' @param min_rate Minimum acceptable enrollment rate
#' @param max_rate Maximum acceptable enrollment rate
#' @return List containing enrollment analysis
#' @export
analyze_enrollment_parameters <- function(trial_data,
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
