rm(list = ls())

# source function
source("0_sources.R")

Nasset = length(FX)

nT = 2
NT = T[nT] * 365
NS = 1000

# lower and upper bound for rho
lb = rep(-0.999, Nasset*(Nasset-1)/2)
ub = -lb

### Find the Differential Evolution parameters
no_cores = detectCores()-1
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = ls(), envir = environment())
clusterEvalQ(cl = cl, source("0_sources.R"))

rho.param = DEoptim(fn = corr_adjustment, lower = lb, upper = ub,
                    control = DEoptim.control(cluster = cl, NP = 120, itermax = 200),
                    his_corr = corr, PARAMS = PARAMS, 
                    NT = NT, NS = NS, T = T[nT], Nasset = Nasset, ZJ = TRUE)

stopCluster(cl)

rho.cal = rho.param$optim$bestmem

model_corr = diag(nrow = Nasset, ncol = Nasset)
model_corr[upper.tri(model_corr)] = rho.cal
model_corr[lower.tri(model_corr)] = t(model_corr)[lower.tri(model_corr)]
rownames(model_corr) = colnames(model_corr) = FX
corr.sim = rowMeans(replicate(2, BatesMQE(model_corr, PARAMS, NT, NS, T[nT], Nasset)$R.sim), dims = 2)

corr.full = FullCorrelation(PARAMS[,"rho"], corr)
model.corr.full = FullCorrelation(PARAMS[,"rho"], model_corr)

save.image("MQEdata.RData")
