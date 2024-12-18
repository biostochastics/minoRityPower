# minoRityPower

**minoRityPower** is an R package designed to estimate statistical power for detecting effects of healthcare system-level interventions on clinical trial enrollment, with a particular focus on accelerating minority participant recruitment. Developed for the ARPA-H program application, this package provides both facility-level and participant-level analyses through mixed-effects models and bootstrap-based simulations.

## Overview

Clinical trials often face challenges in recruiting minority participants, leading to underrepresentation in medical research. Healthcare system-level interventions aim to address this issue by implementing systematic changes to improve minority enrollment. This package provides tools to:

1. Estimate power for detecting intervention effects on facility-level enrollment rates
2. Analyze participant-level enrollment probabilities and intervention interactions with minority status
3. Support comprehensive power analyses at both facility and individual levels
4. Account for the hierarchical nature of clinical trial data
5. Visualize power across different effect sizes and sample sizes

## Installation

```r
# Install from GitHub
devtools::install_github("biostochastics/minoRityPower")
```

## Analysis Approaches

The package implements two complementary approaches to power analysis:

### 1. Facility-Level Analysis

This approach models the rate of enrollment at the facility level, focusing on how interventions affect the number of participants enrolled per facility:

- Uses negative binomial mixed-effects models
- Accounts for facility-level random effects
- Models enrollment counts directly
- Suitable for analyzing system-level effects on recruitment capacity

```r
# Facility-level power analysis
facility_results <- run_power_facility(
  effect_sizes1 = c(1.25, 1.50),    # 25% and 50% increase in enrollment rate
  Ns = c(30, 45, 60),               # Number of facilities per arm
  tau2_facility = 0.64,             # Facility-level variance
  mean_rate = 11                    # Average monthly enrollment rate
)
```

### 2. Participant-Level Analysis

This approach models individual enrollment probabilities, focusing on how interventions affect the likelihood of enrollment for minority participants:

- Uses logistic mixed-effects models
- Incorporates interaction between intervention and minority status
- Models individual enrollment decisions
- Suitable for analyzing differential intervention effects

```r
# Participant-level power analysis
participant_results <- run_power_participant(
  effect_sizes1 = c(1.25, 1.50),    # Main intervention effects
  effect_sizes2 = c(1.75, 2.00),    # Minority-specific effects
  Ns = c(30, 45, 60),               # Facilities per arm
  minority_rate = 0.3,              # Expected minority enrollment proportion
  tau2_trial = 1.19                 # Trial-level variance
)
```

## Key Features

### Facility-Level Analysis Features
- Models enrollment rates using negative binomial distribution
- Accounts for overdispersion in enrollment counts
- Incorporates facility-random effects
- Estimates power for detecting rate ratios

### Participant-Level Analysis Features
- Models enrollment probabilities using logistic regression
- Tests interaction effects with minority status
- Accounts for clustering within facilities
- Supports varying minority enrollment rates

### Shared Features
- Parallel processing for computational efficiency
- Bootstrap-based simulation approaches
- Multiple significance level testing (0.05, 0.01, 0.005)
- Publication-ready visualizations
- Comprehensive output storage

## Output

The package provides detailed output for both analysis types:

### Numerical Output
- Power estimates at multiple significance levels
- Convergence rates and model diagnostics
- Effect size estimates and confidence intervals
- Summary tables in CSV format

### Graphical Output
- Power curves across sample sizes
- Faceted plots by effect size combinations
- Interaction effect visualizations
- Publication-ready figures in PDF format

## Example Workflow

```r
library(minoRityPower)

# 1. Facility-level analysis
facility_results <- run_power_facility(
  effect_sizes1 = c(1.25, 1.50),
  Ns = c(30, 45, 60),
  tau2_facility = 0.64,
  R_boot = 1000
)

# 2. Participant-level analysis
participant_results <- run_power_participant(
  effect_sizes1 = c(1.25, 1.50),    # Main effects
  effect_sizes2 = c(1.75, 2.00),    # Minority interaction effects
  Ns = c(30, 45, 60),
  minority_rate = 0.3,
  R_boot = 1000
)

# 3. View results
print(facility_results$power_table)
print(participant_results$power_table)

# 4. Display plots
print(facility_results$power_plot)
print(participant_results$power_plot)
```

## References

If you use this package, please cite:

```bibtex
@software{your_citation_2024,
  author       = {Sergey Kornilov},
  title        = {minoRityPower: Power Analysis for Healthcare System Interventions in Clinical Trial Enrollment},
  year         = {2024},
  publisher    = {GitHub},
  version      = {0.1.1},
  url          = {https://github.com/biostochastics/minoRityPower}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please submit issues and pull requests through GitHub.

## Acknowledgments

This work was developed to support the evaluation of healthcare system-level interventions for the ARPA-H program application, focusing on accelerating clinical trial enrollment among minority participants.
