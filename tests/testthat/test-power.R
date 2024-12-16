# tests/testthat/test-power.R

test_that("power produces valid output structure", {
  result <- power(
    effect_size = 1.5,
    N = 30,
    R_boot = 5
  )

  # Check structure only
  expect_type(result, "list")
  expect_named(result, c("power", "convergence_rate", "results", "parameters"))
  expect_type(result$power, "double")
  expect_true(!is.null(result$parameters$effect_size))
})

test_that("power preserves input parameters", {
  params <- list(
    effect_size = 1.5,
    N = 30,
    tau2_facility = 0.3,
    mean_rate = 10
  )

  result <- do.call(power, c(params, list(R_boot = 5)))

  # Check parameter storage
  expect_equal(result$parameters$effect_size, params$effect_size)
  expect_equal(result$parameters$N, params$N)
})

test_that("generate_trial_data creates valid structure", {
  data <- generate_trial_data(
    N = 10,
    tau2_facility = 0.3,
    effect_size = 1.25,
    mean_rate = 10
  )

  # Check structure
  expect_true(all(c("trial_id", "facility_id", "intervention",
                    "enrollment_rate") %in% names(data)))
  expect_true(all(data$intervention %in% c(0,1)))
  expect_equal(length(unique(data$trial_id)), 20)
})

test_that("power handles invalid inputs appropriately", {
  expect_error(power(N = 5))
  expect_error(power(mean_rate = 0))
  expect_error(power(tau2_facility = -1))
  expect_error(power(effect_size = 0))
})

test_that("simple_power produces valid output structure", {
  result <- simple_power(
    effect_size = 1.5,
    N = 30,
    R_boot = 5
  )

  # Check structure
  expect_type(result, "list")
  expect_named(result, c("power", "convergence_rate", "results", 
                        "parameters", "diagnostics"))
  expect_type(result$power, "double")
  expect_true(!is.null(result$parameters$effect_size))
  expect_true(!is.null(result$diagnostics$test_data_summary))
})

test_that("simple_power preserves input parameters", {
  params <- list(
    effect_size = 1.5,
    N = 30,
    tau2_facility = 0.3,
    mean_rate = 10
  )

  result <- do.call(simple_power, c(params, list(R_boot = 5)))

  # Check parameter storage
  expect_equal(result$parameters$effect_size, params$effect_size)
  expect_equal(result$parameters$N, params$N)
  expect_equal(result$parameters$tau2_facility, params$tau2_facility)
  expect_equal(result$parameters$mean_rate, params$mean_rate)
})

test_that("generate_trial_data creates valid structure", {
  data <- generate_trial_data(
    N = 10,
    tau2_facility = 0.3,
    effect_size = 1.25,
    mean_rate = 10
  )

  # Check structure
  expect_true(all(c("trial_id", "facility_id", "intervention",
                    "enrollment_rate") %in% names(data)))
  expect_true(all(data$intervention %in% c(0,1)))
  expect_equal(length(unique(data$trial_id)), 20)  # 2*N trials
  expect_true(all(data$enrollment_rate > 0))  # No zero rates
})

test_that("simple_power handles invalid inputs appropriately", {
  expect_error(simple_power(N = 5), "must be at least 10")
  expect_error(simple_power(mean_rate = 0), "must be positive")
  expect_error(simple_power(tau2_facility = -1), "must be non-negative")
  expect_error(simple_power(effect_size = 0), "must be positive")
})

test_that("grid_power handles basic case correctly", {
  # Create temp directory for output
  temp_dir <- tempdir()
  
  result <- grid_power(
    effect_sizes1 = c(1.1, 1.2),
    effect_sizes2 = c(1.1, 1.2),
    Ns = c(30, 40),
    output_dir = file.path(temp_dir, "test_output"),
    R_boot = 3  # Small for testing
  )
  
  # Check output files exist
  expect_true(file.exists(file.path(temp_dir, "test_output", "power_results.rds")))
  expect_true(file.exists(file.path(temp_dir, "test_output", "parameters.rds")))
  
  # Check results structure
  expect_type(result, "list")
  expect_equal(nrow(result$power_grid), 4)  # 2x2 parameter combinations
  expect_true(all(result$power_grid$power >= 0 & result$power_grid$power <= 1))
})

test_that("grid_power respects parameter bounds", {
  temp_dir <- tempdir()
  
  expect_error(
    grid_power(
      effect_sizes1 = c(0.5, 1.0),  # Invalid effect size
      Ns = c(30, 40),
      output_dir = file.path(temp_dir, "test_output")
    ),
    "Effect sizes must be greater than 1"
  )
  
  expect_error(
    grid_power(
      effect_sizes1 = c(1.1, 1.2),
      Ns = c(5, 10),  # Invalid N
      output_dir = file.path(temp_dir, "test_output")
    ),
    "N must be at least"
  )
})
