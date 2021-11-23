rm(list = ls())

# source function
source("0_sources.R")

PARAMS = matrix(nrow = length(FX), ncol = 8)

# Differential Evolution
for (i in 1:length(FX)) {
  cur = i
  source("2_1_Bates_DE_algorithm.R")
  PARAMS[cur,] = DEparam
}

colnames(PARAMS) = c("kappa", "theta", "sigmav", "v0", "rho", "lambdaJ", "muJ", "sigmaJ")
rownames(PARAMS) = FX

rm(list = setdiff(ls(), "PARAMS"))

# store back into option_data.RData
load("option_data.RData")
save.image("option_data.RData")
