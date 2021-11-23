rm(list = ls())

load("MQEdata.RData")

source("0_sources.R")

NS = 5000
NB = 30
testK = 1

BOprice = WOprice = rep(0,4)
names(BOprice) = paste0("T = ", T*365, " days")
WOprice = BOprice

no_cores = detectCores()-1
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = ls(), envir = environment())
clusterEvalQ(cl = cl, source("0_sources.R"))

BOprice[1] = mean(pbreplicate(NB, best.of.Call(PARAMS, testK, r, T = T[1]), cl = cl))
BOprice[2] = mean(pbreplicate(NB, best.of.Call(PARAMS, testK, r, T = T[2]), cl = cl))
BOprice[3] = mean(pbreplicate(NB, best.of.Call(PARAMS, testK, r, T = T[3]), cl = cl))
BOprice[4] = mean(pbreplicate(NB, best.of.Call(PARAMS, testK, r, T = T[4]), cl = cl))

WOprice[1] = mean(pbreplicate(NB, worst.of.Call(PARAMS, testK, r, T = T[1]), cl = cl))
WOprice[2] = mean(pbreplicate(NB, worst.of.Call(PARAMS, testK, r, T = T[2]), cl = cl))
WOprice[3] = mean(pbreplicate(NB, worst.of.Call(PARAMS, testK, r, T = T[3]), cl = cl))
WOprice[4] = mean(pbreplicate(NB, worst.of.Call(PARAMS, testK, r, T = T[4]), cl = cl))

stopCluster(cl)

BOprice
WOprice