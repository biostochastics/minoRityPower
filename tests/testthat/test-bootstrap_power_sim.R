# tests/testthat/test-bootstrap_power_sim.R

test_that("grid_bootstrap_power handles basic case correctly", {
  # Create temp directory for output
  temp_dir <- tempdir()
  
  result <- grid_bootstrap_power(
    effect_sizes1 = c(1.1, 1.2),
    effect_sizes_int = c(1.1, 1.2),
    Ns = c(30, 40),
    tau2_trial = 0.64,
    K_trial = 100,
    output_dir = file.path(temp_dir, "test_bootstrap"),
    R_boot = 3  # Small for testing
  )
  
  # Check output files exist
  expect_true(file.exists(file.path(temp_dir, "test_bootstrap", "power_results.rds")))
  expect_true(file.exists(file.path(temp_dir, "test_bootstrap", "parameters.rds")))
  
  # Check results structure
  expect_type(result, "list")
  expect_equal(nrow(result$power_grid), 4)  # 2x2 parameter combinations
  expect_true(all(result$power_grid$power >= 0 & result$power_grid$power <= 1))
})

test_that("grid_bootstrap_power respects parameter bounds", {
  temp_dir <- tempdir()
  
  expect_error(
    grid_bootstrap_power(
      effect_sizes1 = c(0.5, 1.0),  # Invalid effect size
      effect_sizes_int = c(1.1, 1.2),
      Ns = c(30, 40),
      tau2_trial = 0.64,
      K_trial = 100,
      output_dir = file.path(temp_dir, "test_bootstrap")
    ),
    "Effect sizes must be greater than 1"
  )
  
  expect_error(
    grid_bootstrap_power(
      effect_sizes1 = c(1.1, 1.2),
      effect_sizes_int = c(1.1, 1.2),
      Ns = c(5, 10),  # Invalid N
      tau2_trial = 0.64,
      K_trial = 100,
      output_dir = file.path(temp_dir, "test_bootstrap")
    ),
    "N must be at least"
  )
})

test_that("bootstrap simulation generates valid data", {
  sim_data <- simulate_participant_data(
    N = 30,
    K = 100,
    tau2_trial = 0.64,
    effect_size = 1.25,
    effect_size_int = 1.5,
    prob_minority = 0.3
  )
  
  # Check structure
  expect_true(all(c("trial_id", "participant_id", "intervention",
                    "minority_status", "enrolled") %in% names(sim_data)))
  
  # Check data properties
  expect_equal(length(unique(sim_data$trial_id)), 60)  # 2*N trials
  expect_true(all(sim_data$enrolled %in% c(0,1)))
  expect_true(all(sim_data$minority_status %in% c(0,1)))
  expect_true(mean(sim_data$minority_status) > 0.2 && 
              mean(sim_data$minority_status) < 0.4)  # Roughly 30% minority
})

test_that("bootstrap power analysis handles parallel processing", {
  skip_on_cran()  # Skip on CRAN to avoid heavy computation
  
  temp_dir <- tempdir()
  
  # Test with different numbers of cores
  result_serial <- grid_bootstrap_power(
    effect_sizes1 = c(1.1),
    effect_sizes_int = c(1.1),
    Ns = c(30),
    tau2_trial = 0.64,
    K_trial = 50,
    output_dir = file.path(temp_dir, "test_serial"),
    R_boot = 3,
    n_cores = 1
  )
  
  result_parallel <- grid_bootstrap_power(
    effect_sizes1 = c(1.1),
    effect_sizes_int = c(1.1),
    Ns = c(30),
    tau2_trial = 0.64,
    K_trial = 50,
    output_dir = file.path(temp_dir, "test_parallel"),
    R_boot = 3,
    n_cores = 2
  )
  
  # Results should be similar (not exactly equal due to random generation)
  expect_equal(
    result_serial$power_grid$power,
    result_parallel$power_grid$power,
    tolerance = 0.1
  )
})
