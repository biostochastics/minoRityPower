# tests/testthat/test-enrollment.R

test_that("analyze_cancer_enrollment handles typical cancer trial data", {
  # Create realistic test data
  trials <- data.frame(
    nct_id = rep(1:3, each = 2),
    facility_id = 1:6,
    enrollment = rep(c(100, 150, 200), each = 2),
    duration = rep(12, 6),
    phase = rep(c("Phase 2", "Phase 3", "Phase 2"), each = 2)
  )

  result <- analyze_cancer_enrollment(trials)

  # Test structure
  expect_type(result, "list")
  expect_named(result, c("enrollment_data", "facility_patterns",
                         "summary", "distribution"))

  # Test rate calculations
  expect_equal(
    result$enrollment_data$enrollment_rate,
    c(8.3, 12.5, 16.7),
    tolerance = 0.1
  )

  # Test facility patterns
  expect_equal(nrow(result$facility_patterns), 6)
  expect_true(!is.null(result$facility_patterns$mean_rate))
})

test_that("estimate_enrollment_parameters calculates variance components", {
  trials <- data.frame(
    nct_id = rep(1:10, each = 3),
    facility_id = rep(1:5, 6),
    enrollment_rate = rpois(30, lambda = 10)
  )

  params <- estimate_enrollment_parameters(trials, facilities = NULL)

  expect_type(params, "list")
  expect_named(params, c("mean_rate", "tau2_facility",
                         "dispersion", "empirical_variance"))
  expect_true(params$tau2_facility > 0)
  expect_true(params$mean_rate > 0)
})
