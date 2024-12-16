# tests/testthat/test-plots.R

test_that("plot_power_results handles trial-level analysis", {
  # Create sample results
  temp_dir <- tempdir()
  dir.create(file.path(temp_dir, "test_trial"))
  
  power_grid <- data.frame(
    effect_size1 = rep(c(1.1, 1.2), each = 2),
    effect_size2 = rep(c(1.1, 1.2), 2),
    N = 30,
    power = runif(4),
    convergence_rate = rep(0.95, 4)
  )
  
  saveRDS(power_grid, file.path(temp_dir, "test_trial", "power_results.rds"))
  saveRDS(list(tau2_facility = 0.64), file.path(temp_dir, "test_trial", "parameters.rds"))
  
  plots <- plot_power_results(
    file.path(temp_dir, "test_trial"),
    analysis_type = "trial"
  )
  
  # Check plot objects
  expect_type(plots, "list")
  expect_true("power_heatmap" %in% names(plots))
  expect_s3_class(plots$power_heatmap, "ggplot")
})

test_that("plot_power_results handles participant-level analysis", {
  # Create sample results
  temp_dir <- tempdir()
  dir.create(file.path(temp_dir, "test_participant"))
  
  power_grid <- data.frame(
    effect_size1 = rep(c(1.1, 1.2), each = 2),
    effect_size_int = rep(c(1.1, 1.2), 2),
    N = 30,
    power = runif(4),
    convergence_rate = rep(0.95, 4)
  )
  
  saveRDS(power_grid, file.path(temp_dir, "test_participant", "power_results.rds"))
  saveRDS(list(tau2_trial = 0.64), file.path(temp_dir, "test_participant", "parameters.rds"))
  
  plots <- plot_power_results(
    file.path(temp_dir, "test_participant"),
    analysis_type = "participant"
  )
  
  # Check plot objects
  expect_type(plots, "list")
  expect_true("power_heatmap" %in% names(plots))
  expect_s3_class(plots$power_heatmap, "ggplot")
})

test_that("plot_power_results handles missing files gracefully", {
  temp_dir <- tempdir()
  dir.create(file.path(temp_dir, "empty_dir"))
  
  expect_error(
    plot_power_results(
      file.path(temp_dir, "empty_dir"),
      analysis_type = "trial"
    ),
    "Could not find results"
  )
})

test_that("plot_power_results validates analysis_type", {
  temp_dir <- tempdir()
  
  expect_error(
    plot_power_results(
      temp_dir,
      analysis_type = "invalid"
    ),
    "analysis_type must be"
  )
})
