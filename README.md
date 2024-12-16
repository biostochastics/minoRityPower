# minoRityPower

**minoRityPower** is an R package designed to estimate statistical power for detecting effects in multi-center clinical trials, with a particular focus on minority participant enrollment and intervention effects. The package provides both trial-level and participant-level analyses through mixed-effects models and bootstrap-based simulations.

## Overview

Many clinical trials, particularly those aiming to improve minority inclusion or evaluate specialized interventions, involve hierarchical data structures (e.g., patients nested within facilities, facilities nested within trials). This package provides tools to:

1. Analyze real ClinicalTrials.gov data to estimate empirical parameters
2. Estimate power for detecting intervention effects on overall enrollment rates
3. Analyze interaction effects between interventions and minority status
4. Support both trial-level and participant-level analyses
5. Visualize power across different effect sizes and sample sizes

## Installation

```r
# Install from GitHub
devtools::install_github("biostochastics/minoRityPower")
```

## Usage

The package offers three main workflows:

### 1. Analyzing Real Trial Data

Start by analyzing real ClinicalTrials.gov data to estimate parameters for power analysis:

```r
library(minoRityPower)

# Analyze trial data
trial_results <- analyze_trials(
  studies = studies_data,
  conditions = conditions_data,
  facilities = facilities_data,
  min_duration = 11,    # Minimum trial duration in months
  min_enrollment = 24   # Minimum enrollment count
)

# Estimate parameters for power analysis
params <- estimate_power_parameters(
  trial_results,
  min_rate = 1,     # Minimum acceptable enrollment rate
  max_rate = 200    # Maximum acceptable enrollment rate
)

# Use these parameters in power analysis
power_results <- simple_power(
  effect_size = 1.25,
  N = 50,
  tau2_facility = params$parameters$log_rate_variance,
  mean_rate = params$parameters$mean_rate
)
```

### 2. Trial-Level Analysis

Focuses on facility-level enrollment counts within trial arms:

```r
# Single power analysis
results <- simple_power(
  effect_size = 1.25,  # 25% increase in enrollment
  N = 50,              # 50 trials per arm
  tau2_facility = 0.64 # Facility-level variance
)

# Grid search across parameters
grid_results <- grid_power(
  effect_sizes1 = c(1.01, 1.05, 1.10, 1.25, 1.50, 2.00),
  effect_sizes2 = c(1.01, 1.05, 1.10, 1.25, 1.50, 2.00),
  Ns = c(15, 30, 45, 60),
  output_dir = "results/trial_level"
)

# Visualize results
trial_plots <- plot_power_results(
  "results/trial_level", 
  analysis_type = "trial"
)
```

### 3. Participant-Level Analysis

Simulates individual-level enrollment decisions and outcomes:

```r
# Grid search with bootstrap simulations
bootstrap_results <- grid_bootstrap_power(
  effect_sizes1 = c(1.01, 1.05, 1.10, 1.25, 1.50, 2.00),
  effect_sizes_int = c(1.01, 1.05, 1.10, 1.25, 1.50, 2.00),
  Ns = c(15, 30, 45, 60),
  tau2_trial = params$parameters$log_rate_variance, # Use estimated parameter
  K_trial = 360,
  output_dir = "results/participant_level"
)

# Visualize results
participant_plots <- plot_power_results(
  "results/participant_level", 
  analysis_type = "participant"
)
```

## Key Features

- **Real Data Analysis**: Analyze ClinicalTrials.gov data to estimate empirical parameters
- **Flexible Parameter Exploration**: Test multiple effect sizes, sample sizes, and variance components
- **Two Analysis Levels**: Choose between trial-level or participant-level analysis
- **Parallel Processing**: Efficient computation for bootstrap simulations
- **Result Storage**: Save detailed results for further analysis
- **Publication-Ready Plots**: Generate high-quality visualizations of power analyses

## Output

The package provides comprehensive output at each stage:

### Data Analysis Output
- Summary statistics of real trial data
- Estimated parameters for power analysis
- Facility and trial-level variance components
- Enrollment rate distributions

### Power Analysis Output
- Power estimates at multiple significance levels (0.05, 0.01, 0.005)
- Convergence rates for model fitting
- Detailed simulation results saved to disk
- Summary statistics in CSV format
- Visualizations through `ggplot2`

## References

If you use this package in your research, please cite:

```bibtex
@software{kornilov2024minoritypower,
  author       = {Kornilov, Sergey},
  title        = {minoRityPower: Power Analysis Tools for Minority Recruitment Interventions},
  year         = {2024},
  publisher    = {GitHub},
  version      = {0.1.0},
  url          = {https://github.com/biostochastics/minoRityPower}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

We welcome contributions! Please submit issues and pull requests through GitHub.
