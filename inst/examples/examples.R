library(devtools)

remove.packages("HouseTrans")
.rs.restartR()
devtools::document()
meg = devtools::check()
meg$notes
meg$warnings
meg$errors

devtools::install()
devtools::build_manual()


####################################### Examples #######################################
library(HouseTrans)

# Basic simulation with estimation
result <- GenSyn(
  n_households = 50,
  start_date = "2024-07-01",
  end_date = "2025-06-30",
  stan_chains = 2,
  stan_iter = 800,
  stan_warmup = 10,
  seed = 123
)

# View results
print(result)

# Plot posterior distributions
plot(result, which = "posterior")



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
  stan_chains = 1,
  stan_iter = 200,
  stan_warmup = 10
)

# View covariate effects
plot(result, which = "covariate_effects")

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
  stan_iter = 200,
  stan_warmup = 10
)

print(result)
plot(result, which = "all")
