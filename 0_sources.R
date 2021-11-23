# library
if (!require("statmod"))  install.packages("statmod")
if (!require("Matrix"))   install.packages("Matrix")
if (!require("DEoptim"))  install.packages("DEoptim")
if (!require("pbapply"))  install.packages("pbapply")
if (!require("parallel")) install.packages("parallel")
if (!require("purrr"))    install.packages("purrr")
if (!require("abind"))    install.packages("abind")
if (!require("abind"))    install.packages("abind")
if (!require("devtools")) install.packages("devtools")
if (!require("implvol"))  devtools::install_github("cmhein/implvol")

library(statmod)
library(implvol)
library(Matrix)
library(DEoptim)
library(pbapply)
library(parallel)
library(purrr)
library(abind)

# load data
load("option_data.RData")

# Merton Characteristic Function
MertonCF = function(u, params, S0, r, q, T){
  lambdaJ  = params[1]
  muJ      = params[2]
  sigmaJ   = params[3]
  sigma    = params[4]
  
  x = log(S0)
  muM = (r - q - lambdaJ*muJ) - sigma^2/2
  
  MertonCF = exp(1i*u*x + muM*1i*u*T - (sigma^2 * u^2)/2 * T
                 + T*lambdaJ*((1+muJ)^(1i*u)
                              *exp(0.5*sigmaJ^2*1i*u*(1i*u-1))-1))
  
  return(MertonCF)
}

# Jump Characteristic Function
JumpCF = function(u, params, T){
  lambdaJ  = params[6]
  muJ      = params[7]
  sigmaJ   = params[8]
  
  JumpCF = exp(- lambdaJ*muJ*1i*u*T 
               + T*lambdaJ*((1+muJ)^(1i*u)
                            *exp(0.5*sigmaJ^2*1i*u*(1i*u-1))-1))
  
  return(JumpCF)
}

# Heston Characteristic Function
HestonCF = function(u, params, S0, r, q, T){
  kappa  = params[1]
  theta  = params[2]
  sigmav = params[3]
  v0     = params[4]
  rho    = params[5]
  lambda = 0
  
  x = log(S0)
  
  a = kappa*theta
  b = kappa + lambda
  
  d = sqrt((rho*sigmav*1i*u - b)^2 + sigmav^2*(1i*u + u^2))
  g = (b - rho*sigmav*1i*u + d) / (b - rho*sigmav*1i*u - d)
  c = 1/g
  
  C = (r-q)*1i*u*T + 
    a/sigmav^2*((b - rho*sigmav*1i*u - d)*T
                - 2*log((1 - c*exp(-d*T))/(1-c)))
  D = (b- rho*sigmav*1i*u - d)/sigmav^2*((1-exp(-d*T))/(1-c*exp(-d*T)))
  
  HestonCF = exp(C + D*v0 + 1i*u*x)
  
  return(HestonCF)
}

# Bates characteristic function
BatesCF = function(u, params, S0, r, q, T){
  BatesCF = HestonCF(u, params, S0, r, q, T) * JumpCF(u, params, T) 
  return(BatesCF)
}

# Bates integrand
BatesIntegrand = function(u, params, S0, K, r, q, T, Pnum){
  if (Pnum == 2) {
    Chi = BatesCF(u, params, S0, r, q, T)
    BatesIntegrand = Re( exp(-1i*u*log(K))*Chi / (1i*u) )
  }
  else if (Pnum==1) {
    Chinum = BatesCF(u-1i, params, S0, r, q, T)
    Chiden = BatesCF(   -1i, params, S0, r, q, T)
    BatesIntegrand = Re( exp(-1i*u*log(K))*Chinum / (1i*u*Chiden) )
  } 
  return(BatesIntegrand)
}

# Bates SVJD option pricing model
BatesPrice = function(params, S0, r, q, T, K, ...){
  
  # Gauss-Laguerre nodes and weights
  n = 30
  x = gauss.quad(n,kind="laguerre")$nodes
  w = gauss.quad(n,kind="laguerre")$weights * exp(x)
  
  int1 = c()
  int2 = c()
  for (j in 1:length(x)) {
    int1[j] = w[j] * BatesIntegrand(x[j], params, S0, K, r, q, T, 1)
    int2[j] = w[j] * BatesIntegrand(x[j], params, S0, K, r, q, T, 2)
  }
  
  P1 = 1/2 + 1/pi * sum(int1)
  P2 = 1/2 + 1/pi * sum(int2)
  
  CallPrice = S0*exp(-q*T)*P1 - K*exp(-r*T)*P2
  
  return(CallPrice)
}

# Loss function for Bates' IV using MSE
BatesIVLoss = function(params, S0, r, q, T, K, MktPrice, MktIV, ...){
  
  # Gauss-Laguerre nodes and weights
  n = 30
  x = gauss.quad(n,kind="laguerre")$nodes
  w = gauss.quad(n,kind="laguerre")$weights * exp(x)
  
  NK = dim(MktPrice)[1]
  NT = dim(MktPrice)[2]
  Chi = Chinum = Chiden = int1 = int2 = c()
  ModelPrice = ModelIV = error = matrix(nrow = NK, ncol = NT)
  for (t in 1:NT) {
    for (j in 1:length(x)) {
      u = x[j]
      Chi[j]    = BatesCF(u   ,params,S0,r,q,T[t])
      Chinum[j] = BatesCF(u-1i,params,S0,r,q,T[t])
      Chiden[j] = BatesCF(   -1i,params,S0,r,q,T[t])
    }
    for (k in 1:NK) {
      for (j in 1:length(x)) {
        u = x[j]
        int2[j] = w[j] * Re( exp(-1i*u*log(K[k]))*Chi[j] 
                             / (1i*u) )
        int1[j] = w[j] * Re( exp(-1i*u*log(K[k]))*Chinum[j] 
                             / (1i*u*Chiden[j]) )
      }
      
      P1 = 1/2 + 1/pi*sum(int1)
      P2 = 1/2 + 1/pi*sum(int2)
      
      CallPrice = S0*exp(-q*T[t])*P1 - K[k]*exp(-r*T[t])*P2
      
      ModelPrice[k,t] = CallPrice
      
      ModelIV[k,t] = implvol(price = ModelPrice[k,t], forward = S0*exp((r-q)*T[t]), 
                             strike = K[k], time = T[t], call = TRUE)
      
      error[k,t] = (ModelIV[k,t] - MktIV[k,t])^2
    }
  }
  loss = sum(error) / (NT*NK)
  if (is.na(loss)) loss = Inf
  return(loss)
}

# Single-asset Bates QE scheme
BatesQE = function(params, S0, r, q, T, NT, NS, dZ = NULL, dZJ = NULL) {
  
  kappa   = params[1]
  theta   = params[2]
  sigmav  = params[3]
  v0      = params[4]
  rho     = params[5]
  lambdaJ = params[6]
  muJ     = params[7]
  sigmaJ  = params[8]
  
  psiC   = 1.5
  gamma1 = gamma2 = 0.5
  
  # path of stock price, variance, and log price
  S.path = V.path = matrix(rep(0, NS*(NT+1)), nrow = NS, ncol = NT+1)
  lnS    = V.temp = matrix(rep(0, NS*(NT+1)), nrow = NS, ncol = NT+1)
  
  # time step
  dt = T/NT
  
  lnS[,1]    = log(S0*exp(-q*T))
  V.temp[,1] = v0
  
  alpha = log(1+muJ) - sigmaJ^2/2
  
  k1 = exp(-kappa*dt)
  k2 = sigmav^2*k1/kappa * (1-k1)
  k3 = theta/2 * k2/k1 * (1-k1)
  
  K1 = gamma1*dt*(kappa*rho/sigmav - 0.5) - rho/sigmav
  K2 = gamma2*dt*(kappa*rho/sigmav - 0.5) + rho/sigmav
  K3 = gamma1*dt*(1-rho^2)
  K4 = gamma2*dt*(1-rho^2)
  
  A = K2 + 0.5*K4
  
  if (is.null(dZ)) {
    Zs = matrix(rnorm(NS*NT), nrow = NS, ncol = NT)
    Zv = matrix(rnorm(NS*NT), nrow = NS, ncol = NT)
  } else {
    Zs = t(dZ[1,,])
    Zv = t(dZ[2,,])
  }
  Uv = pnorm(Zv)
  
  if (!is.null(dZJ)) Zj = t(dZJ)
  
  lnM = rep(0,NS)
  
  for (i in 2:(NT+1)) {
    
    m   = theta + (V.temp[,i-1] - theta) * k1
    s2  = V.temp[,i-1]*k2 + k3
    psi = s2/m^2
    
    b2 = 2/psi - 1 + sqrt(2/psi * (2/psi - 1))
    a  = m/(1+b2)
    
    p = (psi-1) / (psi+1)
    beta = (1-p)/m
    
    I1 = which(psi <= psiC)
    I2 = !I1
    
    if (length(I1) > 0) {
      V.temp[I1, i] = a[I1]*(sqrt(b2[I1]) + Zv[I1, i-1])^2
    }
    lnM[I1] = log( exp(A*b2[I1]*a[I1]/(1-2*A*a[I1])) / sqrt(1-2*A*a[I1]) )
    
    
    V.temp[(Uv[,i-1] <= p) & (psi > psiC), i] = 0
    I1b = which((Uv[,i-1] > p) & (psi > psiC))
    
    if (length(I1b) > 0) {
      V.temp[I1b, i] = log( (1-p[I1b]) / (1-Uv[I1b,i-1]) ) / beta[I1b]
    }
    lnM[I2] = log( p[I2] + beta[I2]*(1-p[I2])/(beta[I2]-A) )
    
    P = rpois(NS, lambdaJ*dt)
    
    if (is.null(dZJ)) {
      lnJ = alpha*P + sigmaJ*sqrt(P)*rnorm(NS)
    } else {
      lnJ = alpha*P + sigmaJ*sqrt(P)*Zj[,i-1]
    }
    
    lnS[,i] = lnS[,i-1] + (r-q)*dt - lnM - (K1+0.5*K3)*V.temp[,i-1] + K1*V.temp[,i-1] +
      K2*V.temp[,i] + sqrt(K3*V.temp[,i-1] + K4*V.temp[,i])*Zs[,i-1] +
      lnJ - dt*lambdaJ*muJ
  }
  
  S.path = exp(lnS)
  V.path = V.temp
  
  return(list(S.path = S.path, V.path = V.path))
}

# Function for constructing Full correlation matrix C
FullCorrelation = function(rho, corr) {
  
  Nasset = length(rho)
  
  # correlation for the Brownian motions driving 
  # the 2 * Nasset vector (S1,V1,S2,V2,...,Sn,Vn)
  C = matrix(0, ncol = 2 * Nasset, nrow = 2 * Nasset)
  
  ## correlation matrix C(i,k) representation function
  cor_block = function(rho, corr, i, j) {
    C = matrix(ncol = 2, nrow = 2)
    if (i == j) {
      C = matrix(c(1,rho[i],rho[i],1), ncol = 2, byrow = TRUE)
    }
    else {
      C = corr[i,j]*matrix(c(1,rho[j],rho[i],rho[i]*rho[j]), ncol = 2, byrow = TRUE)
    }
    return(C)
  }
  
  ## full correlation matrix of (S1,V1,S2,V2,...,Sn,Vn)
  for(i in 1:Nasset){
    for(j in 1:Nasset){
      C[2*(i-1)+ c(1,2), 2*(j-1)+ c(1,2)] = cor_block(rho, corr, i, j)
    }
  }
  
  return(C)
}

# Function for constructing decorrelator Q
decorr = function(rho) {
  
  Nasset = length(rho)
  q = q.inv = list()
  
  for (i in 1:Nasset) {
    q[[i]] = matrix(c(sqrt(1-rho[i]^2), rho[i], 0, 1), ncol = 2, byrow = TRUE)
    q.inv[[i]] = solve(q[[i]])
  }
  
  Q.inv = bdiag(q.inv)
  return(Q.inv)
}

# Function for constructing asset-asset correlation matrix R
construct.R = function(model_corr_ij) {
  
  # Nasset*(Nasset-1)/2 = length(model_corr_ij)
  Nasset = sqrt(2*length(model_corr_ij)+0.25) + 0.5
  
  I_N = diag(nrow = Nasset, ncol = Nasset)
  varrho = matrix(0, nrow = Nasset, ncol = Nasset)
  varrho[upper.tri(varrho)] = model_corr_ij
  
  R = I_N + varrho + t(varrho)
  return(R)
}

# Function for checking if positive-semidefinite
is.psd = function(x) {
  eigenvalues = eigen(x, only.values = TRUE)$values
  n = nrow(x)
  for (i in 1:n) {
    if (abs(eigenvalues[i]) < 1e-08) {
      eigenvalues[i] = 0
    }
  }
  if (any(eigenvalues < 0)) {
    return(FALSE)
  }
  return(TRUE)
}

# objective function of multi-asset Bates model calibration using correlation adjustment 
corr_adjustment = function(model_corr_ij, his_corr, PARAMS, NT, NS, T, Nasset, ZJ = TRUE) {
  
  R = construct.R(model_corr_ij)
  
  # Bates MQE scheme
  C = FullCorrelation(rho = PARAMS[, "rho"], corr = R)
  if (!is.psd(C)) return(Inf)
  L = t(chol(C))
  
  Q.inv = as.matrix(decorr(rho = PARAMS[, "rho"]))
  
  dZ = replicate(NS, Q.inv %*% L %*% matrix(rnorm(2*Nasset*NT), nrow = 2*Nasset))
  
  LJ = t(chol(R))
  if (ZJ) {
    dZJ = replicate(NS, LJ %*% matrix(rnorm(Nasset*NT), nrow = Nasset))
    
    SV = list()
    
    for (cur in 1:Nasset) {
      SV[[cur]] = BatesQE(PARAMS[cur,], S0[cur], r, q[cur], T = T, 
                          NT = NT, NS = NS, dZ = dZ[(2*cur-1):(2*cur),,], dZJ = dZJ[cur,,])
    }
  } else {
    SV = list()
    
    for (cur in 1:Nasset) {
      SV[[cur]] = BatesQE(PARAMS[cur,], S0[cur], r, q[cur], T = T, 
                          NT = NT, NS = NS, dZ = dZ[(2*cur-1):(2*cur),,])
    }
  }
  
  S.sim = map(SV, 1)
  
  S.sim = abind(S.sim, along = 3)
  
  Return.sim = apply(S.sim, c(1,3), function(y) diff(log(y)))
  
  R.sim = rowMeans(apply(Return.sim, 2, cor))
  
  rho.sim = R.sim[upper.tri(R)]
  rho.his = his_corr[upper.tri(his_corr)]
  
  loss = mean((rho.sim - rho.his)^2)
  
  return(loss)
}

# Bates MQE simulation scheme
BatesMQE = function(model_corr, PARAMS, NT, NS, T, Nasset, ZJ = TRUE, LJcor = NULL) {
  
  C = FullCorrelation(rho = PARAMS[, "rho"], corr = model_corr)
  L = t(chol(C))
  
  Q.inv = as.matrix(decorr(rho = PARAMS[, "rho"]))
  
  dZ = replicate(NS, Q.inv %*% L %*% matrix(rnorm(2*Nasset*NT), nrow = 2*Nasset))
  
  if (!is.null(LJcor)) {
    LJ = t(chol(diag(nrow = Nasset, ncol = Nasset) + LJcor - 
                  diag(LJcor, nrow = Nasset, ncol = Nasset)))
  } else {
    LJ = t(chol(model_corr))
  }
  
  if (ZJ) {
    dZJ = replicate(NS, LJ %*% matrix(rnorm(Nasset*NT), nrow = Nasset))
    
    SV = list()
    
    for (cur in 1:Nasset) {
      SV[[cur]] = BatesQE(PARAMS[cur,], S0[cur], r, q[cur], T = T, 
                          NT = NT, NS = NS, dZ = dZ[(2*cur-1):(2*cur),,], dZJ = dZJ[cur,,])
    }
  } else {
    SV = list()
    
    for (cur in 1:Nasset) {
      SV[[cur]] = BatesQE(PARAMS[cur,], S0[cur], r, q[cur], T = T, 
                          NT = NT, NS = NS, dZ = dZ[(2*cur-1):(2*cur),,])
    }
  }
  
  S.sim = map(SV, 1)
  V.sim = map(SV, 2)
  
  S.sim = abind(S.sim, along = 3)
  V.sim = abind(V.sim, along = 3)
  
  Return.sim = apply(S.sim, c(1,3), function(y) diff(log(y)))
  
  R.sim = matrix(rowMeans(apply(Return.sim, 2, cor)), nrow = Nasset)
  
  rownames(R.sim) = colnames(R.sim) = rownames(model_corr)
  
  return(list(S.sim = S.sim, V.sim = V.sim, R.sim = R.sim))
}

# Function for pricing the best-of-four call option
best.of.Call = function(P, K, r, T, ZJ = TRUE, LJcor = NULL) {
  MQE_S.sim = BatesMQE(model_corr, P, NT, NS, T, Nasset, ZJ = ZJ, LJcor = LJcor)$S.sim
  ST = MQE_S.sim[,ncol(MQE_S.sim),]
  S0 = MQE_S.sim[,1,]
  
  SimPrice = c()
  for (i in 1:length(K)) {
    Payoff = pmax(apply(ST/S0, 1 , max) - K[i], 0)
    SimPrice[i] = exp(-r*T)*mean(Payoff)
  }

  return(SimPrice)
}

# Function for pricing the worst-of-four call option
worst.of.Call = function(P, K, r, T, ZJ = TRUE, LJcor = NULL) {
  MQE_S.sim = BatesMQE(model_corr, P, NT, NS, T, Nasset, ZJ = ZJ, LJcor = LJcor)$S.sim
  ST = MQE_S.sim[,ncol(MQE_S.sim),]
  S0 = MQE_S.sim[,1,]
  
  SimPrice = c()
  for (i in 1:length(K)) {
    Payoff = pmax(apply(ST/S0, 1 , min) - K[i], 0)
    SimPrice[i] = exp(-r*T)*mean(Payoff)
  }
  
  return(SimPrice)
}