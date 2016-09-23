source("../R/standard.R")
source("../R/reg_smallmatrix.R")

## true expression: low, high, medium genes.
set.seed(12324)
n <- 5
J <- 9
mu_pop <- rep(c(10, 50, 100), each=3)

## homogeneous population == same mu across all samples
mu <- matrix(data=rep(mu_pop, n), ncol=n, nrow=J)

## sample from real GCC and length
info <- read.table("~/git/zinbAnalysis/data/gene.info", row.names=1)
l <- sample(info[,1], J)
g <- sample(info[,2], J)

## values based on BRAIN data
alpha0 <- rnorm(n, mean=2, sd=2)
alpha1 <- rnorm(n, mean=-1, sd=0.3)
alpha2 <- rnorm(n, mean=0, sd=0.15)
alpha3 <- rnorm(n, mean=1.15, sd=2)

logistic <- binomial()$linkinv

## generate pi = prob. of "detection"
pi <- sapply(1:n , function(i) {
  eta <- alpha0[i] + alpha1[i] * log(mu_pop) + alpha2[i] * log(l) + alpha3[i] * g
  return(logistic(eta))
})

## generate Z, indicator of expression
Z <- matrix(rbinom(n*J, 1, pi), ncol=n)
Z

## generate Y, read counts (same dispersion for all genes)
Y <- matrix(data=0, ncol=n, nrow=J)
Y[Z==0] <- rnbinom(sum(Z==0), mu = mu[Z==0], size = 1)
Y

## compute likelihood at the true parameter value
parameters <- c(log(mu_pop), alpha0, alpha1, alpha2, alpha3, log(1))

y <- as.vector(Y)
y0 <- y <= 0
y1 <- y > 0

a0 <- as.factor(as.vector(col(Y)))
b0 <- as.factor(as.vector(row(Y)))
a1 <- model.matrix(~a0-1) * log(mu_pop)
colnames(a1) <- paste("a1_", 1:n, sep="")
a2 <- model.matrix(~a0-1) * log(l)
colnames(a2) <- paste("a2_", 1:n, sep="")
a3 <- model.matrix(~a0-1) * g
colnames(a3) <- paste("a3_", 1:n, sep="")

X <- model.matrix(~b0 - 1)
dim(X)
W <- model.matrix(~a0 + a1 + a2 + a3 - 1)
dim(W)

loglik(parameters, y, y0, y1, X, W, ncol(X), ncol(W), 0, 0, binomial())
grad(parameters, y, y0, y1, X, W, ncol(X), ncol(W), 0, 0, binomial())

## optimize the likelihood with optim, starting from the true values
system.time(fit1 <- optim(fn = loglik, gr = grad, par = parameters, Y=y, Y0=y0, Y1=y1, X=X, W=W, 
                          kx=ncol(X), kw=ncol(W), offsetx=0, offsetz=0, linkobj=binomial(),
                          hessian = FALSE, method = "BFGS", control=list(fnscale=-1)))
fit1$value
fit1$par

## optimize the likelihood with optim, starting from random values
system.time(fit2 <- optim(fn = loglik, gr = grad, par = c(rep(5, J), rnorm(n*4), 0), Y=y, Y0=y0, Y1=y1, X=X, W=W, 
                         kx=ncol(X), kw=ncol(W), offsetx=0, offsetz=0, linkobj=binomial(),
                         hessian = FALSE, method = "BFGS", control=list(fnscale=-1)))
fit2$value
fit2$par

## compute likelihood at the true parameter value (with the memory efficient function)
parameters <- c(log(mu_pop), alpha0, alpha1, alpha2, alpha3, log(1))

X <- matrix(rep(1, n), ncol=1)
dim(X)
W <- model.matrix(~log(mu_pop) + log(l) + g)
dim(W)

nbeta <- J
nalpha <- n*4

loglik_small(parameters, Y, Y>0, X, W, nbeta, nalpha, 0, 0, binomial())
grad_small(parameters, Y, Y>0, X, W, nbeta, nalpha, 0, 0, binomial())

## optimize the likelihood with optim, starting from the true values
system.time(fit3 <- optim(fn = loglik_small, gr = grad_small, par = parameters, Y=Y, Y1=Y>0, X=X, W=W, 
                          kx=nbeta, kw=nalpha, offsetx=0, offsetw=0, linkobj=binomial(),
                          hessian = FALSE, method = "BFGS", control=list(fnscale=-1)))
fit3$value
fit3$par

