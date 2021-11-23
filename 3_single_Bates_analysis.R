# Parameter Estimation
rm(list = ls())

# source function
source("0_sources.R")

if (!require("ggplot2")) install.packages("ggplot2")
if (!require("reshape")) install.packages("reshape")
if (!require("ggpubr"))  install.packages("ggpubr")

library(ggplot2)
library(reshape)
library(ggpubr)

# change cur to plot a different currency
cur = 4

# recording the plot and save the plot in working directory
dpi = 600
tiff(paste0("Figure_", cur, ".tif"), width=12*dpi, height=3*dpi, res = dpi)
plots = list()

DEparam = PARAMS[cur,]
# Compute MSE Loss
BatesIVLoss(PARAMS[cur,], S0[cur], r, q[cur], T, K[[cur]], MktPrice[[cur]], MktIV[[cur]])

# compute implied volatilities
NK = dim(MktPrice[[cur]])[1]
NT = dim(MktPrice[[cur]])[2]
BPrice = BatesIV = matrix(nrow = NK, ncol = NT)
for (k in 1:NK) {
  for (t in 1:NT) {
    BPrice[k,t]   = BatesPrice(DEparam, S0 = S0[cur], r = r, q = q[cur], T = T[t], K = K[[cur]][k])
    BatesIV[k,t]  = implvol(price = BPrice[k,t], forward = S0[cur]*exp((r-q[cur])*T[t]), strike = K[[cur]][k],
                            time = T[t], call = TRUE)
  }
}

# Plot the implied volatilities

NT = length(T)

for (t in 1:NT) {
  df = data.frame("K" = K[[cur]], "Market" = MktIV[[cur]][,t], "Model" = BatesIV[,t])
  IVdf = melt(df, id.vars = "K")
  plots[[t]] =
    ggplot(IVdf, aes(x = K, y = value*100
                     , color = variable, linetype = variable)) +
    geom_smooth(se = FALSE, method = "loess") +
    scale_linetype_manual(values=c(2:1)) +
    labs(title = paste(as.character(T[t]*365), "days to maturity"),
         x = "Strike Price K")+
    theme(legend.title = element_blank(), axis.title.y = element_blank()) +
    theme(text = element_text(family = "serif"))
}

plotmix = ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]],ncol=4, nrow =1,
                    common.legend = TRUE, legend = "bottom")

annotate_figure(plotmix,
                left = text_grob("Implied Volatility (%)", rot = 90, size = 12,
                                 family = "serif", hjust = 0.3))

model.plot = recordPlot()
dev.off()