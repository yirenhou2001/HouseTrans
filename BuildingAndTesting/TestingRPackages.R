library(devtools)
library(tidyverse)
#library(pbapply)

remove.packages("Household.Transmission.Chain.Data.Analysis")
.rs.restartR()
#devtools::document()
meg = devtools::check()
meg$notes
meg$warnings
meg$errors

devtools::install()
#devtools::build_manual()

library(Household.Transmission.Chain.Data.Analysis)
# 1) Simulate and estimate
out1 <- GenSyn(
  n_households = 10,
  n_runs       = 10,
)
out1

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
out2
