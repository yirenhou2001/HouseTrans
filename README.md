# HouseTrans

<!-- badges: start -->
[![R-CMD-check](https://github.com/yirenhou2001/HouseTrans/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yirenhou2001/HouseTrans/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

A Bayesian framework for simulating household infection dynamics and estimating transmission parameters using Stan. The model incorporates viral load dynamics, seasonal forcing, role-specific susceptibility/infectivity, covariates (e.g., vaccination status), and support for reinfections with waning immunity.

## Core Workflow

1. **`GenSyn()`**: Simulation and validation
   - Simulates household structures with viral load trajectories
   - Estimates parameters using Bayesian Stan
   - Validates results by comparing estimates against known ground truth

2. **`TransmissionChainAnalysis()`**: User data analysis
   - Estimates transmission parameters from your own observational data
   - Accepts either long-format testing data or per-person episode tables
   
## Features

- **Household transmission modeling** with role-specific parameters (adult, infant, toddler, elderly)
- **Viral load dynamics** using double-exponential (log10) or piecewise linear (Ct) trajectories
- **Covariate support** for susceptibility and infectivity modifiers (e.g., vaccination efficacy)
- **Reinfection modeling** with gamma-distributed waning immunity
- **Interval sampling** for infection timing with recovery tail modeling
- **Flexible Bayesian priors** (Normal, Uniform, LogNormal)
- **Pre-compiled Stan model** for faster fitting
- **Covariate-aware transmission chain reconstruction** from posterior estimates
- **Comprehensive visualization** tools

## Installation

```r
# Install from GitHub
devtools::install_github("yirenhou2001/HouseTrans")
```

## Quick Start

### 1. Simulate and Estimate (GenSyn)

```r
library(HouseTrans)

# Basic simulation with estimation
result <- GenSyn(
  n_households = 50,
  start_date = "2024-07-01",
  end_date = "2025-06-30",
  stan_chains = 2,
  stan_iter = 1000,
  stan_warmup = 500,
  seed = 123
)

# View results
print(result)

# Plot posterior distributions
plot(result, which = "posterior")
```

### 2. With Covariates (e.g., Vaccination)

```r
# Define vaccination covariate
vacc_config <- list(
  list(
    name = "vacc_status",
    efficacy = 0.8,
    effect_on = "both",  # affects susceptibility and infectivity
    coverage = list(infant = 0.8, toddler = 0, adult = 0, elderly = 0)
  )
)

result <- GenSyn(
  n_households = 50,
  covariates_config = vacc_config,
  covariates_susceptibility = "vacc_status",
  covariates_infectivity = "vacc_status",
  stan_chains = 2,
  stan_iter = 1000
)

# View covariate effects
plot(result, which = "covariate_effects")
```

### 3. Analyze Your Own Data (TransmissionChainAnalysis)

```r
# Per-person episode format
df_person <- data.frame(
  hh_id = c("HH1","HH1","HH1","HH2","HH2","HH2"),
  person_id = c(1, 2, 3, 1, 2, 3),
  role = c("adult","infant","elderly","adult","infant","elderly"),
  infection_time = c(2, 4, NA, 1, 3, NA),
  infectious_start = c(3, 6, NA, 2, 5, NA),
  infectious_end = c(8, 9, NA, 7, 9, NA),
  infection_resolved = c(9, 10, NA, 8, 10, NA)
)

result <- TransmissionChainAnalysis(
  user_data = df_person,
  max_days = 30,
  stan_chains = 2,
  stan_iter = 1000,
  stan_warmup = 500
)

print(result)

# Plot epidemic curve of your data
plot(result, which = "epidemic_curve")
```

## Input Data Formats

### Per-person Episode Table (recommended)

| Column | Description |
|--------|-------------|
| `hh_id` | Household identifier |
| `person_id` | Individual ID within household |
| `role` | Family role: "adult", "infant", "toddler", "elderly" |
| `infection_time` | Day of infection onset (NA if not infected) |
| `infectious_start` | Day infectiousness begins |
| `infectious_end` | Day infectiousness ends |
| `infection_resolved` | Day infection resolves |

### Long-format Testing Table

| Column | Description |
|--------|-------------|
| `HH` | Household identifier |
| `individual_ID` | Individual ID within household |
| `role` | Family role |
| `test_date` | Integer day or Date object |
| `infection_status` | Binary (0/1) indicating positive test |

## Surveillance Data

**Surveillance data** is real-world disease monitoring data from public health agencies (e.g., CDC, local health departments). It provides daily or weekly case counts representing the actual epidemic in the community.

### Purpose

1. **Seasonal forcing**: Tells the model when the virus was circulating in the community. High case counts = higher community transmission risk.
2. **Validation**: For `GenSyn()`, overlay simulated infections against real surveillance to validate your simulation matches real-world patterns.
3. **Context**: For `TransmissionChainAnalysis()`, optionally overlay your household data against population-level trends.

### Format

```r
surveillance_df <- data.frame(
  date = seq(as.Date("2024-07-01"), as.Date("2025-06-30"), by = "day"),
  cases = c(10, 12, 15, 20, 25, 30, ...)  # daily case counts
)
```

### Usage

```r
# With GenSyn - surveillance affects simulation AND creates comparison plot
result <- GenSyn(
  n_households = 50,
  surveillance_df = surveillance_df,  # Used for seasonal forcing
  start_date = "2024-07-01",
  end_date = "2025-06-30"
)
plot(result, which = "epidemic_curve")  # Compares simulation to surveillance

# With TransmissionChainAnalysis - optional overlay on your data
result <- TransmissionChainAnalysis(
  user_data = df_person,
  surveillance_df = surveillance_df  # Optional: for comparison plot
)
plot(result, which = "epidemic_curve")  # Shows your data, optionally with surveillance
```

## Advanced Configuration

### Prior Configuration

Customize priors using a named list:

```r
my_priors <- list(
  beta1 = list(dist = "normal", params = c(-5, 1)),
  beta2 = list(dist = "normal", params = c(-5, 1)),
  alpha = list(dist = "normal", params = c(-6, 2)),
  covariates = list(dist = "normal", params = c(0, 1)),
  gen_shape = list(dist = "lognormal", params = c(1.5, 0.5)),
  gen_rate = list(dist = "lognormal", params = c(0.0, 0.5)),
  ct50 = list(dist = "normal", params = c(35.0, 3.0)),
  slope = list(dist = "lognormal", params = c(0.4, 0.5))
)

result <- GenSyn(n_households = 50, priors = my_priors)
```

### Recovery Parameters

Customize the immunity waning tail by role:

```r
recovery_params <- list(
  adult   = list(shape = 2, scale = 3),
  infant  = list(shape = 2, scale = 3),
  toddler = list(shape = 2, scale = 3),
  elderly = list(shape = 2, scale = 3)
)

result <- TransmissionChainAnalysis(
  user_data = df_person,
  recovery_params = recovery_params
)
```

### Household Profile Configuration

Customize household composition probabilities:

```r
household_profile <- list(
  prob_adults   = c(0, 0, 1),      # Probability of 0/1/2 adults (default: always 2)
  prob_infant   = 1.0,             # Probability of having an infant (default: always 1)
  prob_siblings = c(0, 0.8, 0.2),  # Probability of 0/1/2 toddlers/siblings
  prob_elderly  = c(0.7, 0.1, 0.2) # Probability of 0/1/2 elderly
)

result <- GenSyn(
  n_households = 100,
  household_profile_list = household_profile
)
```

**Role naming:**
| Profile Parameter | Concept | Role Generated |
|------------------|---------|----------------|
| `prob_adults` | Parents | `"adult"` |
| `prob_infant` | Baby/index child | `"infant"` |
| `prob_siblings` | Older siblings | `"toddler"` |
| `prob_elderly` | Grandparents | `"elderly"` |

### Viral Curve Imputation Parameters

Customize mechanistic viral curve imputation:

```r
# For Ct data
imputation_params <- list(
  adult   = list(Cpeak = 33, r = 1.5, d = 1.2, t_peak = 5),
  infant  = list(Cpeak = 33.3, r = 2.11, d = 1.38, t_peak = 5.06),
  toddler = list(Cpeak = 34, r = 1.26, d = 1.27, t_peak = 4.75),
  elderly = list(Cpeak = 33, r = 1.49, d = 1.22, t_peak = 5.14)
)

# For log10 viral load data
imputation_params <- list(
  adult   = list(v_p = 4.14, t_p = 5.09, lambda_g = 2.31, lambda_d = 2.71),
  infant  = list(v_p = 5.84, t_p = 4.09, lambda_g = 2.82, lambda_d = 1.01),
  toddler = list(v_p = 5.84, t_p = 4.09, lambda_g = 2.82, lambda_d = 1.01),
  elderly = list(v_p = 2.95, t_p = 5.1, lambda_g = 3.15, lambda_d = 0.87)
)
```

## Available Plots

Use `plot(result, which = "...")` with:

### For GenSynResult

| Plot | Description |
|------|-------------|
| `"posterior"` | Posterior distributions of phi/kappa by role |
| `"covariate_effects"` | Forest plot of covariate coefficients |
| `"epidemic_curve"` | Simulated infections vs surveillance data |
| `"transmission_chains"` | Transmission link probabilities for a household |
| `"all"` | All available plots |

### For TransmissionChainResult

| Plot | Description |
|------|-------------|
| `"posterior"` | Posterior distributions of phi/kappa by role |
| `"covariate_effects"` | Forest plot of covariate coefficients |
| `"epidemic_curve"` | User data infections (with optional surveillance overlay) |
| `"transmission_chains"` | Transmission link probabilities for a household |
| `"all"` | All available plots |

### Plot Options

```r
# Transmission chain plot for specific household
plot(result, which = "transmission_chains", hh_id = 3, prob_cutoff = 0.1)

# Epidemic curve with custom bin width (default: 7 days)
plot(result, which = "epidemic_curve", bin_width = 14)
```

## Model Details

The Stan model estimates:

- **`phi[role]`**: Role-specific susceptibility multipliers
- **`kappa[role]`**: Role-specific infectivity multipliers
- **`beta1`**: Baseline transmission coefficient
- **`beta2`**: Viral load contribution to transmission
- **`alpha_comm`**: Community transmission rate
- **`gen_shape`, `gen_rate`**: Generation interval parameters
- **`Ct50`, `slope_ct`**: Viral load infectivity parameters
- **`beta_susc`, `beta_inf`**: Covariate coefficients (when covariates specified)

## Output Structure

### GenSynResult

- `$call`: The matched function call
- `$n_households`: Number of households simulated
- `$simulation`: Raw simulation output (`hh_df`, `diagnostic_df`)
- `$surveillance_df`: Surveillance data (if provided)
- `$start_date`, `$end_date`: Study period
- `$stan_data`: Data prepared for Stan
- `$fit`: The stanfit object
- `$postprocessing`: Tidy posterior summary
- `$attack_rates`: Primary attack rate and reinfection summaries
- `$transmission_chains`: Reconstructed transmission links

### TransmissionChainResult

- `$call`: The matched function call
- `$user_data`: The processed input data
- `$surveillance_df`: Surveillance data (if provided)
- `$start_date`, `$end_date`: Study period
- `$stan_data`: Data prepared for Stan
- `$fit`: The stanfit object
- `$postprocessing`: Tidy posterior summary
- `$transmission_chains`: Reconstructed transmission links
