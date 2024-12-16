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
