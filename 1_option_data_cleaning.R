# Data and Parameter Requirement
## 1. Call option data
##    source: https://www.cmegroup.com/tools-information/quikstrike/option-settlement.html
## 2. Expiry date
##    source: https://www.cmegroup.com/tools-information/quikstrike/option-settlement.html
## 3. Domestic risk-free rate (3 Month ICE LIBOR)
##    source: https://www.theice.com/marketdata/reports/170
## 4. Foreign risk-free rate (3 Month ICE LIBOR)
##    source: https://www.theice.com/marketdata/reports/170
## 5. Historical price
##    source: library(priceR) from European Central Bank Statistical Data Warehouse

rm(list = ls())

# read library
if (!require("readr"))       install.packages("reader")
if (!require("readxl"))      install.packages("readxl")
if (!require("xts"))         install.packages("xts")
if (!require("dplyr"))       install.packages("dplyr")
if (!require("priceR"))      install.packages("priceR")
if (!require("matrixStats")) install.packages("matrixStats")

library(readr)
library(readxl)
library(xts)
library(dplyr)
library(priceR)
library(matrixStats)

# data path
option.path = list.files(pattern = "_rawdata")

# currency name
FX = substring(option.path, 1, 3)

# read option price
call = Call = iv = IV = K = MktPrice = MktIV = list()
for (j in option.path) {
  
  sheetname = excel_sheets(j)
  
  call[[j]] = iv[[j]] = list()
  for (i in sheetname) {
    call[[j]][[i]] = data.frame(read_xlsx(j, col_names = TRUE, sheet = i))[,c(2,1)]
    iv[[j]][[i]] = data.frame(read_xlsx(j, col_names = TRUE, sheet = i))[,c(2,3)]
    names(call[[j]][[i]]) = names(iv[[j]][[i]]) = c("K", i)
  }
  Call[[j]] = Reduce(full_join, call[[j]])
  Call[[j]] = na.omit(Call[[j]][order(Call[[j]]$K),])
  
  IV[[j]] = Reduce(full_join, iv[[j]])
  IV[[j]] = na.omit(IV[[j]][order(IV[[j]]$K),])
  
  K[[j]]  = Call[[j]]$K
  
  MktPrice[[j]] = Call[[j]][,-1]
  MktIV[[j]] =   IV[[j]][,-1]/100
}

names(K) = names(MktPrice) = names(MktIV) = FX

# set maturity
Expiry = as.Date(sheetname)
today  = as.Date("2020-11-20")
Maturity = as.numeric(Expiry - today)
T  = Maturity/365

# domestic risk-free rate (USD)
r = 0.20488/100

# foreign risk-free rate (CHF, EUR, GBP, JPY)
q = c(-0.76780, -0.54386, 0.05100, -0.10200)/100

# historical price and return
price = matrix(nrow = 366, ncol = length(FX))
n = 0
for (i in FX) {
  n = n + 1
  temp = historical_exchange_rates(from = i, to = "USD",
                                   start_date = today - 365, end_date = today)
  price[,n] = temp[,2]
}
S = xts(price, order.by = as.Date(temp[,1]))
colnames(S) = FX
S0 = price[nrow(price),]

Returns = diff(log(S))
R = na.omit(Returns)

# return standard deviation
sigma = colSds(R)

# correlation matrix
corr = cor(xts::last(R, paste(Maturity[ceiling(length(T)/2)], "days")), method = "pearson")

# keep necessary data
rm(list = setdiff(ls(), c("S", "S0", "R", "r", "q", "MktPrice", "MktIV", "K", "T", 
                          "sigma", "corr", "FX", "Maturity")))

# save data
save.image("option_data.RData")
