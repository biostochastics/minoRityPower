#' Generate Interventions Data
#' @keywords internal
#' @param K_trial Number of trials
#' @return A data frame with trial_id and assigned arm (x1).
generate_interventions_data <- function(K_trial) {
  # Assign half to Arm0 and half to Arm1 (or as even as possible)
  interventions <- data.frame(
    trial_id = 1:K_trial,
    x1 = factor(rep(c("Arm0", "Arm1"), times = c(floor(K_trial / 2),
                                                 K_trial - floor(K_trial / 2))))
  )
  return(interventions)
}


#' Simulate Participant-Level Data
#' @keywords internal
#' @param interventions A data frame from generate_interventions_data
#' @param tau2_trial Variance of random trial intercepts
#' @param seed Single integer seed (this will be offset by boot index externally)
#' @return A data frame of all participant-level records before sampling
simulate_participants <- function(interventions, tau2_trial, seed) {
  set.seed(seed)

  K_trial <- nrow(interventions)
  b0_trial <- rnorm(n = K_trial, mean = 0, sd = sqrt(tau2_trial))

  list_records <- vector("list", K_trial)
  for (i in seq_len(K_trial)) {
    eligible <- max(rnorm(1, mean = 1000, sd = 300), 100) # Ensure >=100
    minority <- rbinom(eligible, size = 1, prob = 0.2)
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
#' @param all_records Data frame of participant-level data
#' @param b0_trial Vector of random intercepts for each trial
#' @param effect_size1 Numeric, main effect (Arm1)
#' @param effect_size_int Numeric, interaction effect (Arm1 & minority)
#' @return Data frame with outcome column `y` added
assign_outcomes <- function(all_records, b0_trial, effect_size1, effect_size_int) {
  # Baseline rates differ by minority status
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
