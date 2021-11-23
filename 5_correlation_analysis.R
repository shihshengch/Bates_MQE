rm(list = ls())

load("MQEdata.RData")

source("0_sources.R")

no_cores = detectCores()-1
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = ls(), envir = environment())
clusterEvalQ(cl = cl, source("0_sources.R"))

corr.sim = rowMeans(pbreplicate(1000, 
                                BatesMQE(model_corr, PARAMS, NT, NS, 
                                         T[nT], Nasset)$R.sim, cl = cl), 
                    dims = 2)

stopCluster(cl)

corr
corr.sim
model_corr

mean(((corr - corr.sim)^2)[upper.tri(corr)])

model.corr.full = FullCorrelation(PARAMS[,"rho"], model_corr)
model.corr.full
