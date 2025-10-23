
### {Household.Transmission.Chain.Data.Analysis}: Household Transmission Chain Simulation and Estimation

{Household.Transmission.Chain.Data.Analysis} helps users do two things
in a streamlined process:

1.  **Generate synthetic household data** and run the **full estimation
    pipeline** end-to-end with `GenSyn()`.
2.  **Estimate transmission parameters from user data** with
    `TransmissionChainAnalysis()`.

Within this pipeline, the package:

- summarizes individuals and imputes infection timelines (latent
  periods, reporting delays, infectious periods);
- builds a person-day table suitable for likelihood-based estimation;
- fits penalized models for community and household risks with optional
  covariates (shared or role-specific);
- provides post-processing that reports mean estimates, uncertainty
  (standard error), bias, and relative bias;
- provides the option to summarize data by infection episodes per
  individual.

## Package purpose

This package lets user simulate household transmission data, estimate
parameters of interest, and compare results either against known “true”
values (specified for synthetic data simulation and analysis) or compare
parameter estimates of user’s own data (which may have missing values or
limited covariates) with parameter estimates of simulated synthetic data
for population of interest.

- **One call, full workflow**: simulate data/user data -\> summarize
  individuals -\> impute infection time -\> estimate parameters -\>
  post-process of estimates.
- **Emulates** user study’s structure (dates, testing frequency,
  covariates, missingness) to see how coefficient estimates shift
  relative to the “ideal” simulated situation.
- **Compares** to ground truth for simulations, post-processing reports
  bias and relative bias versus the known parameters.

## Quick start

``` r
library(Household.Transmission.Chain.Data.Analysis)
# 1) Simulate and estimate
out1 <- GenSyn(
  n_households = 10,
  n_runs       = 10,
)
```

    ## Initialized start_par of length 8 (based on available covariates).

``` r
# 2) Estimate from your own long-format data
HH = c(rep(1L, 6), rep(2L, 6))
individual_ID = c(1,1,2,2,3,3, 1,1,2,2,3,3)
role <- c("infant","infant","adult","adult","sibling","sibling",
"infant","infant","adult","adult","elder","elder")
test_date = c(1,8,1,8,1,8, 1,8,1,8,1,8) 
infection_status = c(0,1,0,0,0,0, 0,0,0,1,0,0) 
community_risk = rep(0.001, length(HH))

df = data.frame(HH, individual_ID, role, test_date, infection_status, community_risk)

out2 <- TransmissionChainAnalysis(
  user_data = df,                 # see required columns in ?TransmissionChainAnalysis
  n_runs    = 20
)
```

    ## Initialized start_par of length 8 (based on available covariates).

## Inputs

- `GenSyn()` — for synthetic data only (errors if synthetic_data =
  FALSE).
- `TransmissionChainAnalysis()` — for user data with required columns
  (HH, individual_ID, role, test_date, infection_status,
  community_risk). Additional columns can be mapped as covariates.

## Outputs

Both functions return a list:

- Results: raw simulations (if synthetic), summaries, person–day table,
  and estimates;
- Postprocessing: a comparison table of mean estimates versus ground
  truth values (for synthetic data) or summary statistics (for user
  data)
