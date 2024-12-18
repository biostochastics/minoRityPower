# tests/testthat/test-patient-power.R

test_that("generate_interventions_data creates valid structure", {
  K_trial <- 10
  interventions <- generate_interventions_data(K_trial)
  
  expect_equal(nrow(interventions), K_trial)
  expect_true(all(c("trial_id", "x1") %in% names(interventions)))
  expect_true(all(interventions$x1 %in% c("Arm0", "Arm1")))
  expect_equal(sum(interventions$x1 == "Arm0"), floor(K_trial / 2))
})

test_that("simulate_participants generates valid data", {
  K_trial <- 10
  interventions <- generate_interventions_data(K_trial)
  minority_rate <- 0.2
  tau2_trial <- 0.5
  
  set.seed(123)
  participants <- simulate_participants(
    interventions = interventions,
    tau2_trial = tau2_trial,
    minority_rate = minority_rate,
    seed = 42
  )
  
  expect_true(all(c("trial_id", "x1", "eligible_index", "minority") %in% names(participants)))
  expect_true(all(participants$minority %in% c(0, 1)))
  expect_true(abs(mean(participants$minority) - minority_rate) < 0.1)  # Approximate check
})

test_that("assign_outcomes produces valid probabilities", {
  K_trial <- 10
  interventions <- generate_interventions_data(K_trial)
  set.seed(123)
  participants <- simulate_participants(
    interventions = interventions,
    tau2_trial = 0.5,
    minority_rate = 0.2,
    seed = 42
  )
  b0_trial <- attr(participants, "b0_trial")
  
  result <- assign_outcomes(
    all_records = participants,
    b0_trial = b0_trial,
    effect_size1 = 1.25,
    effect_size_int = 1.5
  )
  
  expect_true("y" %in% names(result))
  expect_true(all(result$y %in% c(0, 1)))
})

test_that("sample_data_for_analysis maintains data structure", {
  K_trial <- 20
  N <- 10
  interventions <- generate_interventions_data(K_trial)
  set.seed(123)
  all_records <- simulate_participants(
    interventions = interventions,
    tau2_trial = 0.5,
    minority_rate = 0.2,
    seed = 42
  )
  
  seed_vec <- 1:100
  boot_index <- 1
  
  sampled_data <- sample_data_for_analysis(
    all_records = all_records,
    N = N,
    interventions = interventions,
    seed_vec = seed_vec,
    boot_index = boot_index
  )
  
  expect_true(is.factor(sampled_data$trial_id))
  expect_equal(length(unique(sampled_data$trial_id)), N)
  expect_true(all(table(sampled_data$x1) > 0))  # Both arms present
})

test_that("fit_models_and_get_pvalue handles convergence", {
  # Create small dataset that should converge
  K_trial <- 10
  interventions <- generate_interventions_data(K_trial)
  set.seed(123)
  participants <- simulate_participants(
    interventions = interventions,
    tau2_trial = 0.5,
    minority_rate = 0.2,
    seed = 42
  )
  b0_trial <- attr(participants, "b0_trial")
  
  data <- assign_outcomes(
    all_records = participants,
    b0_trial = b0_trial,
    effect_size1 = 2.0,  # Large effect for clear signal
    effect_size_int = 2.0
  )
  
  p_value <- fit_models_and_get_pvalue(data)
  
  expect_true(!is.na(p_value))
  expect_true(p_value >= 0 && p_value <= 1)
})

test_that("bootstrap_power_sim produces valid output", {
  result <- bootstrap_power_sim(
    boot_index = 1,
    effect_size1 = 1.25,
    effect_size2 = 1.5,
    N = 30,
    seed_vec = 1:100,
    K_trial = 100,
    tau2_trial = 0.5,
    minority_rate = 0.2
  )
  
  expect_type(result, "list")
  expect_true("p_value" %in% names(result))
  expect_true(is.numeric(result$p_value))
  expect_true(is.na(result$p_value) || (result$p_value >= 0 && result$p_value <= 1))
})

test_that("run_power_participant handles basic case", {
  temp_dir <- tempdir()
  
  result <- run_power_participant(
    effect_sizes1 = c(1.25),
    effect_sizes2 = c(1.5),
    Ns = c(30),
    tau2_trial = 0.5,
    K_trial = 100,
    minority_rate = 0.2,
    R_boot = 3,  # Small for testing
    results_dir = file.path(temp_dir, "participant_power")
  )
  
  expect_true(!is.null(result))
  expect_true(file.exists(file.path(temp_dir, "participant_power", "power_summary.csv")))
})
