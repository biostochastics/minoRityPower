# tests/testthat/test-analysis.R

# Create sample data for tests
create_test_data <- function() {
  list(
    studies = data.frame(
      nct_id = 1:3,
      phase = rep("Phase 2", 3),
      overall_status = "Completed",
      enrollment_type = "Actual",
      study_type = "Interventional",
      start_date = as.Date("2020-01-01"),
      completion_date = as.Date("2021-01-01"),
      enrollment = c(100, 150, 200)
    ),
    facilities = data.frame(
      nct_id = rep(1:3, each = 2),
      facility_id = 1:6,
      name = paste0("Facility", 1:6)
    ),
    conditions = data.frame(
      nct_id = 1:3,
      condition = "Cancer"
    )
  )
}

test_that("analyze_trials handles basic case correctly", {
  test_data <- create_test_data()
  
  result <- analyze_trials(
    studies = test_data$studies,
    conditions = test_data$conditions,
    facilities = test_data$facilities
  )

  # Check structure
  expect_type(result, "list")
  expect_named(result, c("trials", "metrics", "facilities"))
  
  # Check content
  expect_equal(nrow(result$trials), 3)
  expect_true(all(result$metrics$n_facilities == 2))
  expect_true(all(result$trials$enrollment_rate > 0))
})

test_that("analyze_trials respects minimum requirements", {
  test_data <- create_test_data()
  test_data$studies$enrollment <- c(10, 150, 200)  # One below min
  
  result <- analyze_trials(
    studies = test_data$studies,
    conditions = test_data$conditions,
    facilities = test_data$facilities,
    min_enrollment = 100
  )
  
  expect_equal(nrow(result$trials), 2)  # Should exclude the low enrollment trial
})

test_that("estimate_power_parameters produces valid output", {
  test_data <- create_test_data()
  trial_results <- analyze_trials(
    studies = test_data$studies,
    conditions = test_data$conditions,
    facilities = test_data$facilities
  )
  
  params <- estimate_power_parameters(
    trial_results,
    min_rate = 1,
    max_rate = 200
  )
  
  # Check structure
  expect_type(params, "list")
  expect_named(params, c("parameters", "distribution", "data"))
  
  # Check parameter values
  expect_true(params$parameters$log_rate_variance > 0)
  expect_true(params$parameters$mean_rate > 0)
  
  # Check data filtering
  expect_true(all(params$data$enrollment_rate >= 1))
  expect_true(all(params$data$enrollment_rate <= 200))
})

test_that("estimate_power_parameters handles edge cases", {
  test_data <- create_test_data()
  trial_results <- analyze_trials(
    studies = test_data$studies,
    conditions = test_data$conditions,
    facilities = test_data$facilities
  )
  
  # Test with very restrictive rate bounds
  expect_warning(
    estimate_power_parameters(trial_results, min_rate = 100, max_rate = 110),
    "Few trials"
  )
  
  # Test with invalid bounds
  expect_error(
    estimate_power_parameters(trial_results, min_rate = -1),
    "must be positive"
  )
  expect_error(
    estimate_power_parameters(trial_results, min_rate = 200, max_rate = 100),
    "must be greater than"
  )
})

test_that("analyze_trials handles missing columns gracefully", {
  test_data <- create_test_data()
  test_data$studies$nct_id <- NULL
  
  expect_error(
    analyze_trials(
      studies = test_data$studies,
      conditions = test_data$conditions,
      facilities = test_data$facilities
    ),
    "Required columns missing"
  )
})

test_that("analyze_trials calculates duration correctly", {
  test_data <- create_test_data()
  test_data$studies$start_date <- as.Date("2020-01-01")
  test_data$studies$completion_date <- as.Date("2021-01-01")
  
  result <- analyze_trials(
    studies = test_data$studies,
    conditions = test_data$conditions,
    facilities = test_data$facilities
  )
  
  expect_equal(unique(result$trials$duration), 12)  # 12 months duration
})

test_that("analyze_trials handles NA dates appropriately", {
  test_data <- create_test_data()
  test_data$studies$start_date[1] <- NA
  
  result <- analyze_trials(
    studies = test_data$studies,
    conditions = test_data$conditions,
    facilities = test_data$facilities
  )
  
  expect_equal(nrow(result$trials), 2)  # Should exclude the NA date trial
})

test_that("analyze_trials calculates enrollment rate correctly", {
  test_data <- create_test_data()
  test_data$studies$enrollment <- c(120, 240, 360)  # 10, 20, 30 per month
  test_data$studies$start_date <- as.Date("2020-01-01")
  test_data$studies$completion_date <- as.Date("2021-01-01")
  
  result <- analyze_trials(
    studies = test_data$studies,
    conditions = test_data$conditions,
    facilities = test_data$facilities
  )
  
  expect_equal(result$metrics$enrollment_rate, c(10, 20, 30))
})
