library(devtools)
library(tidyverse)
#library(pbapply)

remove.packages("Household.Transmission.Chain.Data.Analysis")
.rs.restartR()
devtools::document()
meg = devtools::check()
meg$notes
meg$warnings
meg$errors

devtools::install()
devtools::build_manual()
library(Household.Transmission.Chain.Data.Analysis)

result_example1 <- GenSyn(
  n_households = 10,
  n_runs       = 10,
  engine = "legacy",
  estimation_method = "mle"
)
result_example1$postprocessing


result_example2 <- GenSyn(
  n_households = 10,
  n_runs       = 10,
  data_summary = TRUE,

  engine = "legacy",
  estimation_method = "mle"
)
result_example2$data_summary
result_example2$postprocessing
print(result_example2)


seasonal_forcing_list <- readRDS("seasonal_forcing_list.rds")
result_example3 <- GenSyn(
  n_households = 5,
  plots = c("daily","weekly","timeline","sar"),

  engine = "rsv_vl",
  estimation_method = "stan",
  seasonal_forcing_list = seasonal_forcing_list,

  stan_chains = 1,
  stan_iter   = 800,
  stan_warmup = 400,
  stan_control = list(adapt_delta = 0.98, max_treedepth = 12),
  stan_refresh = 25,
  stan_cores   = 1
)
print(result_example3)
result_example3$plots
result_example3$plot_list[1]
result_example3$plot_list[2]
result_example3$plot_list[3]
result_example3$plot_list[4]





HH = c(rep(1L, 6), rep(2L, 6))
individual_ID = c(1,1,2,2,3,3, 1,1,2,2,3,3)
role <- c("infant","infant","adult","adult","sibling","sibling",
          "infant","infant","adult","adult","elder","elder")
test_date = c(1,8,1,8,1,8, 1,8,1,8,1,8)
infection_status = c(0,1,0,0,0,0, 0,0,0,1,0,0)
community_risk = rep(0.001, length(HH))

df = data.frame(HH, individual_ID, role, test_date, infection_status, community_risk)

result_example4 <- TransmissionChainAnalysis(
  user_data = df,
  n_runs    = 20,
  estimation_method = "mle"
)
print(result_example4)




T_max <- 12

df_person <- rbind(
  data.frame(
    hh_id             = "HH1",
    person_id         = 1:3,
    role              = c("adult","child","elderly"),
    infection_time    = c( 2,  4, NA),
    infectious_start  = c( 3,  6, NA),
    infectious_end    = c( 8,  9, NA),
    infection_resolved= c( 9, 10, NA),
    stringsAsFactors  = FALSE
  ),
  data.frame(
    hh_id             = "HH2",
    person_id         = 1:3,
    role              = c("adult","child","elderly"),
    infection_time    = c( 1,  3, NA),
    infectious_start  = c( 2,  5, NA),
    infectious_end    = c( 7,  9, NA),
    infection_resolved= c( 8, 10, NA),
    stringsAsFactors  = FALSE
  )
)

# Optional list-cols (safe to omit)
df_person$vl_full_trajectory    <- vector("list", nrow(df_person))
df_person$viral_loads_test_days <- vector("list", nrow(df_person))

# Flat seasonal forcing for all roles (length must match max_days)
seasonal_forcing_list <- list(
  adult   = rep(1, T_max),
  child   = rep(1, T_max),
  elderly = rep(1, T_max),
  toddler = rep(1, T_max)
)


result_example5 <- TransmissionChainAnalysis(
  user_data         = df_person,    # <-- per-person episodes format
  estimation_method = "stan",

  # RSV/VL + Stan knobs
  seasonal_forcing_list = seasonal_forcing_list,
  max_days              = T_max,

  # keep Stan light so it finishes quickly
  stan_chains  = 1,
  stan_iter    = 300,
  stan_warmup  = 150,
  stan_control = list(adapt_delta = 0.98, max_treedepth = 12),
  stan_init    = "random",
  stan_refresh = 50,
  stan_cores   = 1
)
print(result_example5)
result_example5$plots
