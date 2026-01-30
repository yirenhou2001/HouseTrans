library(devtools)
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

