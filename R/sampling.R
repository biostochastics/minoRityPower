#' Sample Participants and Trials
#' @keywords internal
#' @param all_records Full participant-level data
#' @param N Total number of trials to select (half control, half intervention)
#' @param interventions Data frame of trial assignments
#' @param seed_vec Vector of seeds, indexed by bootstrap iteration
#' @param boot_index The current bootstrap iteration index
#' @return A data frame `dat_rr` with the sampled trials and participants
sample_data_for_analysis <- function(all_records, N, interventions, seed_vec, boot_index) {
  # Sample 75% of participants from each trial
  # Using dplyr; ensure dplyr is imported or fully qualified
  dat_sub <- all_records %>%
    dplyr::group_by(trial_id) %>%
    dplyr::sample_frac(0.75) %>%
    dplyr::ungroup()

  N2 <- N / 2
  set.seed(seed_vec[boot_index])

  # Sample trials: N/2 from Arm0 and N/2 from Arm1
  sampled_trials <- c(
    sample(interventions %>% dplyr::filter(x1 == "Arm0") %>% dplyr::pull(trial_id),
           size = N2, replace = TRUE),
    sample(interventions %>% dplyr::filter(x1 == "Arm1") %>% dplyr::pull(trial_id),
           size = N2, replace = TRUE)
  )

  dat_rr <- dat_sub[dat_sub$trial_id %in% sampled_trials, ]
  dat_rr$trial_id <- factor(dat_rr$trial_id)
  return(dat_rr)
}
