library(zinbwave)
library(SummarizedExperiment)
library(ggplot2)
library(cowplot)
library(gridExtra)

# test
if(TRUE){
  set.seed(64839)
  m <- zinbModel(n=100, J=500, K=2)
  x <- zinbSim(m)
  counts = x$counts[rowSums(x$counts) != 0 , ]
  counts = counts[, colSums(counts) != 0]

  fit <- lapply(0:1, function(k){
    print(k)
    print(system.time(ff<-zinbFit(counts, K = k)))
    ff
  })

  plot(sapply(fit, zinb.AIC, t(counts)))
  plot(sapply(fit, zinb.BIC, t(counts)))
  plot(sapply(fit, loglik, t(counts)))
}


