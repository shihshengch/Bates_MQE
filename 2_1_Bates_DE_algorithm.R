# lower and upper bound for Bates parameter
l = 1e-7
u = 5
#      kappa theta sigmav  v0    rho lambdaJ  muJ sigmaJ
lb = c(    l,    l,     l,  l, -.999,      l,  -u,     l)
ub = c(   22,    u,     u,  u,  .999,      u,   u,     u)

### Find the Differential Evolution parameters
no_cores = detectCores()-1
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = ls(), envir = environment())
clusterEvalQ(cl = cl, source("0_sources.R"))

DEparam = DEoptim(fn = BatesIVLoss, lower = lb, upper = ub,
                  control = DEoptim.control(cluster = cl, NP = 120, itermax = 500),
                  S0 = S0[cur], r = r, q = q[cur], T = T, K = K[[cur]], 
                  MktPrice = MktPrice[[cur]], MktIV = MktIV[[cur]])

stopCluster(cl)

DEparam = DEparam$optim$bestmem