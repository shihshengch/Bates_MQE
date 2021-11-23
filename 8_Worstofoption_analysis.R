rm(list = ls())

load("MQEdata.RData")

source("0_sources.R")

nT = 2
NT = T[nT]*365
NS = 5000
NB = 30
testK = seq(0.9,1.3,0.01)

PARAMSJ1 = PARAMSJ2 = PARAMSJ3 = PARAMS

PARAMSJ1[,6] = 0.05
PARAMSJ2[,6] = 0.1
PARAMSJ3[,6] = 0.2

WOprice0 = WOprice1 = WOprice2 = WOprice3 = matrix(nrow = length(testK), ncol = 6)
colnames(WOprice0) = c("original", "Jump Corr = 0", "Jump Corr = 0.25", 
                       "Jump Corr = 0.5", "Jump Corr = 0.75", "Jump Corr = 0.999")

WOprice1 = WOprice2 = WOprice3 = WOprice0
no_cores = detectCores()-1
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = ls(), envir = environment())
clusterEvalQ(cl = cl, source("0_sources.R"))

WOprice0[,1] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMS, testK, r, T = T[nT]), cl = cl))
WOprice0[,2] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMS, testK, r, T = T[nT], ZJ = FALSE), cl = cl))
WOprice0[,3] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMS, testK, r, T = T[nT], LJcor = 0.25), cl = cl))
WOprice0[,4] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMS, testK, r, T = T[nT], LJcor = 0.5 ), cl = cl))
WOprice0[,5] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMS, testK, r, T = T[nT], LJcor = 0.75), cl = cl))
WOprice0[,6] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMS, testK, r, T = T[nT], LJcor = 0.999), cl = cl))

WOprice1[,1] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMSJ1, testK, r, T = T[nT]), cl = cl))
WOprice1[,2] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMSJ1, testK, r, T = T[nT], ZJ = FALSE), cl = cl))
WOprice1[,3] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMSJ1, testK, r, T = T[nT], LJcor = 0.25), cl = cl))
WOprice1[,4] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMSJ1, testK, r, T = T[nT], LJcor = 0.5 ), cl = cl))
WOprice1[,5] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMSJ1, testK, r, T = T[nT], LJcor = 0.75), cl = cl))
WOprice1[,6] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMSJ1, testK, r, T = T[nT], LJcor = 0.999), cl = cl))

WOprice2[,1] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMSJ2, testK, r, T = T[nT]), cl = cl))
WOprice2[,2] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMSJ2, testK, r, T = T[nT], ZJ = FALSE), cl = cl))
WOprice2[,3] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMSJ2, testK, r, T = T[nT], LJcor = 0.25), cl = cl))
WOprice2[,4] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMSJ2, testK, r, T = T[nT], LJcor = 0.5 ), cl = cl))
WOprice2[,5] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMSJ2, testK, r, T = T[nT], LJcor = 0.75), cl = cl))
WOprice2[,6] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMSJ2, testK, r, T = T[nT], LJcor = 0.999), cl = cl))

WOprice3[,1] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMSJ3, testK, r, T = T[nT]), cl = cl))
WOprice3[,2] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMSJ3, testK, r, T = T[nT], ZJ = FALSE), cl = cl))
WOprice3[,3] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMSJ3, testK, r, T = T[nT], LJcor = 0.25), cl = cl))
WOprice3[,4] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMSJ3, testK, r, T = T[nT], LJcor = 0.5 ), cl = cl))
WOprice3[,5] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMSJ3, testK, r, T = T[nT], LJcor = 0.75), cl = cl))
WOprice3[,6] = rowMeans(pbreplicate(NB, worst.of.Call(PARAMSJ3, testK, r, T = T[nT], LJcor = 0.999), cl = cl))

stopCluster(cl)

WOprice = list(cbind(testK,WOprice0), cbind(testK,WOprice1), cbind(testK,WOprice2), cbind(testK,WOprice3))

# Plot
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("reshape")) install.packages("reshape")
if (!require("ggpubr"))  install.packages("ggpubr")

library(ggplot2)
library(reshape)
library(ggpubr)

dpi = 600
tiff(paste0("Figure_7.tif"), width=7*dpi, height=6*dpi, res = dpi)

plots = list()
LAMBDAJ = c("original", "0.05", "0.1", "0.2")
for (t in 1:4) {
  df = as.data.frame(WOprice[[t]])
  WOdf = melt(df, id.vars = "testK")
  plots[[t]] =
    ggplot(WOdf, aes(x = testK, y = value
                     , color = variable, linetype = variable)) +
    geom_line(aes(lty = variable)) +
    scale_linetype_manual(values=c(1:6)) +
    labs(title = paste0("Jump intensity = ", LAMBDAJ[t]),
         x = "Strike Price", y = "Option Price")+
    theme(legend.title = element_blank()) +
    theme(text = element_text(family = "serif"))
}

ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]],ncol=2, nrow =2,
          common.legend = TRUE, legend = "bottom")

model.plot = recordPlot()
dev.off()

dpi = 600
tiff(paste0("Figure_8.tif"), width=10*dpi, height=6*dpi, res = dpi)

plots = list()
JumpCorr = c("Asset Corr", "0", "0.25", "0.5", "0.75", "0.999")
for (p in 1:6) {
  df = cbind(testK, WOprice0[,p],WOprice1[,p],WOprice2[,p],WOprice3[,p])
  colnames(df) = c("testK", paste0("Jump intensity = ", LAMBDAJ))
  df = as.data.frame(df)
  WOdf = melt(df, id.vars = "testK")
  plots[[p]] =
    ggplot(WOdf, aes(x = testK, y = value
                     , color = variable, linetype = variable)) +
    geom_line(aes(lty = variable)) +
    scale_linetype_manual(values=c(1:4)) +
    labs(title = paste0("Jump Corr = ", JumpCorr[p]),
         x = "Strike Price", y = "Option Price")+
    theme(legend.title = element_blank()) +
    theme(text = element_text(family = "serif"))
}

ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], 
          ncol=3, nrow =2,
          common.legend = TRUE, legend = "bottom")

model.plot = recordPlot()
dev.off()