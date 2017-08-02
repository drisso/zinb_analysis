library(zinbwave)
library(SummarizedExperiment)
library(ggplot2)
library(cowplot)
library(gridExtra)


# Functions
nParams <- function(model){
  X_pi <- getX_pi(model)
  X_mu <- getX_mu(model)
  V_pi <- getV_pi(model)
  V_mu <- getV_mu(model)
  disp <- get
  
  n <- nSamples(model)
  J <- nFeatures(model)
  K <- nFactors(model)
  M_pi <- NCOL(X_pi)
  M_mu <- NCOL(X_mu)
  L_pi <- NCOL(V_pi)
  L_mu <- NCOL(V_mu)
  ndisp <- length(unique(getZeta(model)))
  
  total <- J * (M_mu + M_pi) + n * (L_mu + L_pi) + 2 * K * J + n * K + ndisp
  return(total)
}

zinb.AIC <- function(model, x){
  if ((nSamples(model) != nrow(x))|(nFeatures(model) != ncol(x))) {
    stop("x and model should have the same dimensions!")
  }
  stopifnot()
  k <- nParams(model)
  ll <- loglik(model, x)
  return(2*k - 2*ll)
}

zinb.BIC <- function(model, x){
  n <- nSamples(model)
  if ((n != nrow(x))|(nFeatures(model) != ncol(x))) {
    stop("x and model should have the same dimensions!")
  }
  k <- nParams(model)
  ll <- loglik(model, x)
  return(log(n)*k - 2*ll)
}
  
# test
if(FALSE){
  set.seed(64839)
  m <- zinbModel(n=100, J=500, K=2)
  x <- zinbSim(m)
  counts = x$counts[rowSums(x$counts) != 0 , ]
  counts = counts[, colSums(counts) != 0]
  
  fit <- lapply(0:5, function(k){
    print(k)
    print(system.time(ff<-zinbFit(counts, K = k)))
    ff
  })
  
  plot(sapply(fit, zinb.AIC, t(counts)))
  plot(sapply(fit, zinb.BIC, t(counts)))
  plot(sapply(fit, loglik, t(counts)))
}







