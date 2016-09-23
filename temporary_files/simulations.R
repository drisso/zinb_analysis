## simulate data

simulateSimple <- function(J=2000, n0=50, n1=50, nde=0, phi=1, seed=123, 
                           fixedPhi=TRUE) {
  set.seed(seed)
  
  n <- n0 + n1
  group <- as.factor(c(rep(0, n0), rep(1, n1)))
  
  loglambda <- rnorm(J, mean=2.5, sd=1.3)
  log2foldchange <- rep(0, J)
  logOR <- rep(0, J)
  
  if(nde>0) {
    log2foldchange[1:nde] <- sample(1:4, nde, replace=TRUE) * sample(c(1, -1), nde, replace=TRUE)
  }

  logfoldchange <- log(2) * log2foldchange
  
  lambda0 <- exp(loglambda - 0.5 * logfoldchange)
  lambda1 <- exp(loglambda + 0.5 * logfoldchange)
  
  if(!fixedPhi) {
    phi <- exp(rnorm(J, mean=log(phi), sd=0.5))
  }
  
  alpha0 <- rnorm(n, mean=2, sd=2)
  alpha1 <- rnorm(n, mean=-1, sd=0.3)
  
  pi <- sapply(1:n , function(i) {
    eta <- alpha0[i] + alpha1[i] * loglambda
    return(logistic(eta))
  })
  
  Z <- matrix(rbinom(n*J, 1, pi), ncol=n)
  mu <- matrix(data=0, ncol=n, nrow=J)
  mu[,group==0] <- lambda0
  mu[,group==1] <- lambda1
  
  Y <- matrix(data=0, ncol=n, nrow=J)
  Y[Z==0] <- rnbinom(sum(Z==0), mu = mu[Z==0], size = 1)
  head(Y)
  
  rownames(Y) <- rownames(Z) <- rownames(pi) <- paste0("gene", 1:J)
  colnames(Y) <- colnames(Z) <- colnames(pi) <- names(group) <- paste0("sample", 1:n)
  
  if(nde>0) {
    pos <- rownames(Y)[1:nde]
  } else {
    pos <- NULL
  }
  neg <- rownames(Y)[(nde+1):J]
  
  return(list(Y=Y, Z=Z, pi=pi, pos=pos, neg=neg, x=group, logfoldchange=logfoldchange,
              phi=phi, lambda0=lambda0, lambda1=lambda1, alpha0=alpha0, alpha1=alpha1))
  
}

sim <- simulateSimple()

## First, let's use the (unknown) true log abundance as covariate
parameters <- c(log(sim$lambda0), sim$alpha0, sim$alpha1, log(sim$phi))
n <- 100
J <- 20000

X <- matrix(rep(1, n), ncol=1)
dim(X)
W <- model.matrix(~log(sim$lambda0))
dim(W)

nbeta <- J
nalpha <- n*2

loglik_small(parameters, sim$Y, sim$Y>0, X, W, nbeta, nalpha, 0, 0, binomial())
grad_small(parameters, sim$Y, sim$Y>0, X, W, nbeta, nalpha, 0, 0, binomial())

## optimize the likelihood with optim, starting from the true values
system.time(fit <- optim(fn = loglik_small, gr = grad_small, par = parameters, Y=sim$Y, Y1=sim$Y>0, X=X, W=W, 
                         kx=nbeta, kw=nalpha, offsetx=0, offsetw=0, linkobj=binomial(),
                         hessian = FALSE, method = "BFGS", control=list(fnscale=-1)))

betahat <- matrix(fit$par[1:nbeta], ncol=J, nrow=ncol(X))
muhat <- t(exp(X %*% betahat))[,1]

plot(log(sim$lambda0), log(muhat) - log(sim$lambda0))
lines(lowess(log(sim$lambda0), log(muhat) - log(sim$lambda0)), col=2, lwd=2)
abline(h=0, lty=2)

alphahat <- matrix(fit$par[(nbeta + 1):(nbeta + nalpha)], nrow=ncol(W), ncol=n, byrow=TRUE)
etahat <- W %*% alphahat
pihat <- logistic(etahat)

plot(sim$alpha0, alphahat[1,] - sim$alpha0)
lines(lowess(sim$alpha0, alphahat[1,] - sim$alpha0), col=2, lwd=2)
abline(h=0, lty=2)

plot(sim$alpha1, alphahat[2,] - sim$alpha1)
lines(lowess(sim$alpha1, alphahat[2,] - sim$alpha1), col=2, lwd=2)
abline(h=0, lty=2)

plot(sort(sim$pi[,1]), pihat[order(sim$pi[,1]),1], type='l', ylim=c(0, 1), xlim=c(0, 1))
for(i in 1:n) {
  lines(sort(sim$pi[,i]), pihat[order(sim$pi[,i]),i])
}
lines(sort(rowMeans(sim$pi)), rowMeans(pihat)[order(rowMeans(sim$pi))], col=2, lwd=2)

plot(sim$pi[,1], pihat[,1])
abline(0, 1)

plot(colMeans(sim$Y), colMeans(sim$pi-pihat))
plot(rowMeans(sim$pi), rowMeans(sim$pi-pihat))

thetahat <- exp(fit$par[nbeta + nalpha + 1])
thetahat

## optimize the likelihood with optim, starting from "random" values
parameters <- c(rowMeans(log(sim$Y+1)), rep(0, 2*n), log(1))

system.time(fit <- optim(fn = loglik_small, gr = grad_small, par = parameters, Y=sim$Y, Y1=sim$Y>0, X=X, W=W, 
                         kx=nbeta, kw=nalpha, offsetx=0, offsetw=0, linkobj=binomial(),
                         hessian = FALSE, method = "BFGS", control=list(fnscale=-1)))

betahat <- matrix(fit$par[1:nbeta], ncol=J, nrow=ncol(X))
muhat <- t(exp(X %*% betahat))[,1]

plot(log(sim$lambda0), log(muhat) - log(sim$lambda0))
lines(lowess(log(sim$lambda0), log(muhat) - log(sim$lambda0)), col=2, lwd=2)
abline(h=0, lty=2)

alphahat <- matrix(fit$par[(nbeta + 1):(nbeta + nalpha)], nrow=ncol(W), ncol=n, byrow=TRUE)
etahat <- W %*% alphahat
pihat <- logistic(etahat)

plot(sim$alpha0, alphahat[1,] - sim$alpha0)
lines(lowess(sim$alpha0, alphahat[1,] - sim$alpha0), col=2, lwd=2)
abline(h=0, lty=2)

plot(sim$alpha1, alphahat[2,] - sim$alpha1)
lines(lowess(sim$alpha1, alphahat[2,] - sim$alpha1), col=2, lwd=2)
abline(h=0, lty=2)

plot(sort(sim$pi[,1]), pihat[order(sim$pi[,1]),1], type='l', ylim=c(0, 1), xlim=c(0, 1))
for(i in 1:n) {
  lines(sort(sim$pi[,i]), pihat[order(sim$pi[,i]),i])
}
lines(sort(rowMeans(sim$pi)), rowMeans(pihat)[order(rowMeans(sim$pi))], col=2, lwd=2)

plot(sim$pi[,1], pihat[,1])
abline(0, 1)

plot(colMeans(sim$Y), colMeans(sim$pi-pihat))
plot(rowMeans(sim$pi), rowMeans(sim$pi-pihat))

thetahat <- exp(fit$par[nbeta + nalpha + 1])
thetahat

####
parameters <- c(rowMeans(log(sim$Y+1)), rep(1, n), rep(-1, n), log(1))

system.time(fit <- optim(fn = loglik_small, gr = grad_small, par = parameters, Y=sim$Y, Y1=sim$Y>0, X=X, W=W, 
                         kx=nbeta, kw=nalpha, offsetx=0, offsetw=0, linkobj=binomial(),
                         hessian = FALSE, method = "BFGS", control=list(fnscale=-1)))

betahat <- matrix(fit$par[1:nbeta], ncol=J, nrow=ncol(X))
muhat <- t(exp(X %*% betahat))[,1]

plot(log(sim$lambda0), log(muhat) - log(sim$lambda0))
lines(lowess(log(sim$lambda0), log(muhat) - log(sim$lambda0)), col=2, lwd=2)
abline(h=0, lty=2)

alphahat <- matrix(fit$par[(nbeta + 1):(nbeta + nalpha)], nrow=ncol(W), ncol=n, byrow=TRUE)
etahat <- W %*% alphahat
pihat <- logistic(etahat)

plot(sim$alpha0, alphahat[1,] - sim$alpha0)
lines(lowess(sim$alpha0, alphahat[1,] - sim$alpha0), col=2, lwd=2)
abline(h=0, lty=2)

plot(sim$alpha1, alphahat[2,] - sim$alpha1)
lines(lowess(sim$alpha1, alphahat[2,] - sim$alpha1), col=2, lwd=2)
abline(h=0, lty=2)

plot(sort(sim$pi[,1]), pihat[order(sim$pi[,1]),1], type='l', ylim=c(0, 1), xlim=c(0, 1))
for(i in 1:n) {
  lines(sort(sim$pi[,i]), pihat[order(sim$pi[,i]),i])
}
lines(sort(rowMeans(sim$pi)), rowMeans(pihat)[order(rowMeans(sim$pi))], col=2, lwd=2)

plot(sim$pi[,1], pihat[,1])
abline(0, 1)

plot(colMeans(sim$Y), colMeans(sim$pi-pihat))
plot(rowMeans(sim$pi), rowMeans(sim$pi-pihat))

thetahat <- exp(fit$par[nbeta + nalpha + 1])
thetahat

## Plug-in average of Y as log abundance covariate
X <- matrix(rep(1, n), ncol=1)
dim(X)

Y <- sim$Y[rowMeans(sim$Y)>0,]
Y0 <- Y
Y0[Y==0] <- NA
J <- nrow(Y)
W <- model.matrix(~log(rowMeans(Y0, na.rm=TRUE)+1))
dim(W)

nbeta <- J
nalpha <- n*2

parameters <- c(rowMeans(log(Y+1)), rep(1, n), rep(-1, n), log(1))

system.time(fit <- optim(fn = loglik_small, gr = grad_small, par = parameters, Y=Y, Y1=Y>0, X=X, W=W, 
                         kx=nbeta, kw=nalpha, offsetx=0, offsetw=0, linkobj=binomial(),
                         hessian = FALSE, method = "BFGS", control=list(fnscale=-1)))

betahat <- matrix(fit$par[1:nbeta], ncol=J, nrow=ncol(X))
muhat <- t(exp(X %*% betahat))[,1]

plot(log(sim$lambda0), log(muhat) - log(sim$lambda0))
lines(lowess(log(sim$lambda0), log(muhat) - log(sim$lambda0)), col=2, lwd=2)
abline(h=0, lty=2)

alphahat <- matrix(fit$par[(nbeta + 1):(nbeta + nalpha)], nrow=ncol(W), ncol=n, byrow=TRUE)
etahat <- W %*% alphahat
pihat <- logistic(etahat)

plot(sim$alpha0, alphahat[1,] - sim$alpha0)
lines(lowess(sim$alpha0, alphahat[1,] - sim$alpha0), col=2, lwd=2)
abline(h=0, lty=2)

plot(sim$alpha1, alphahat[2,] - sim$alpha1)
lines(lowess(sim$alpha1, alphahat[2,] - sim$alpha1), col=2, lwd=2)
abline(h=0, lty=2)

pi <- sim$pi
pi <- pi[rownames(Y),]
plot(sort(pi[,1]), pihat[order(pi[,1]),1], type='l', ylim=c(0, 1), xlim=c(0, 1))
for(i in 1:n) {
  lines(sort(pi[,i]), pihat[order(pi[,i]),i])
}
lines(sort(rowMeans(pi)), rowMeans(pihat)[order(rowMeans(pi))], col=2, lwd=2)

plot(pi[,1], pihat[,1])
abline(0, 1)

plot(colMeans(Y), colMeans(pi-pihat))
plot(rowMeans(pi), rowMeans(pi-pihat))

thetahat <- exp(fit$par[nbeta + nalpha + 1])
thetahat
