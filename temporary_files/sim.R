source("../R/standard.R")
source("../R/reg_smallmatrix.R")
source("../R/em.R")

## true expression: from BRAIN data
set.seed(12324)
library(matrixStats)
library(brainData)
data(brain01_counts)

filteredCells <- read.table("~/git/brainAnalysis/brain/data/brain01_filtered_cells.txt", as.is=TRUE)

exclude <- !(brain01_counts$Barcode %in% filteredCells[,1])
brain01 <- brain01_counts[,!exclude]

counts <- exprs(brain01)
type <- as.factor(brain01$"Cell.type")

genes <- rownames(brain01)[fData(brain01)[,3]==1]
genes <- intersect(genes, rownames(na.omit(counts)))

filtered <- filterCount(counts[genes,], nRead=10, nCell=10)

geneMap <- read.table("~/git/brainAnalysis/brain/data/gene.map", row.names=1)

ens <- rownames(geneMap)
names(ens) <- geneMap[,1]
allgenes <- ens[rownames(filtered)]

allgenes <- na.omit(allgenes)
info <- read.table("~/git/brainAnalysis/brain/data/gene.info", row.names=1)
len <- info[allgenes,1]
gcc <- info[allgenes,2]
names(len) <- names(gcc) <- names(allgenes)

allgenes <- allgenes[names(na.omit(len))]

len <- len[names(allgenes)]
gcc <- gcc[names(allgenes)]
filtered <- filtered[names(allgenes),]

filtered0 <- filtered
filtered0[filtered==0] <- NA

means <- round(rowMedians(filtered0, na.rm=TRUE))
names(means) <- rownames(filtered0)

J <- 100
n <- 10

mu_pop <- sample(means, J)

## homogeneous population == same mu across all samples
mu <- matrix(data=rep(mu_pop, n), ncol=n, nrow=J)

## sample from real GCC and length
l <- len[names(mu_pop)]
g <- gcc[names(mu_pop)]

## values based on BRAIN data
alpha0 <- rnorm(n, mean=2, sd=2)
alpha1 <- rnorm(n, mean=-1, sd=0.3)
alpha2 <- rnorm(n, mean=0, sd=0.15)
alpha3 <- rnorm(n, mean=1.15, sd=2)

logistic <- binomial()$linkinv

## generate pi = prob. of "detection"
pi <- sapply(1:n , function(i) {
  eta <- alpha0[i] + alpha1[i] * log(mu_pop) #+ alpha2[i] * log(l) + alpha3[i] * g
  return(logistic(eta))
})

## generate Z, indicator of expression
Z <- matrix(rbinom(n*J, 1, pi), ncol=n)
head(Z)

## generate Y, read counts (same dispersion for all genes)
Y <- matrix(data=0, ncol=n, nrow=J)
Y[Z==0] <- rnbinom(sum(Z==0), mu = mu[Z==0], size = 1)
head(Y)

## compute likelihood at the true parameter value (with the memory efficient function)
parameters <- c(log(mu_pop), alpha0, alpha1, log(1))

X <- matrix(rep(1, n), ncol=1)
dim(X)
W <- model.matrix(~log(mu_pop))
dim(W)

nbeta <- J
nalpha <- n*2

loglik_small(parameters, Y, Y>0, X, W, nbeta, nalpha, 0, 0, binomial())
grad_small(parameters, Y, Y>0, X, W, nbeta, nalpha, 0, 0, binomial())

## optimize the likelihood with optim, starting from the true values
system.time(fit <- optim(fn = loglik_small, gr = grad_small, par = parameters, Y=Y, Y1=Y>0, X=X, W=W, 
                          kx=nbeta, kw=nalpha, offsetx=0, offsetw=0, linkobj=binomial(),
                          hessian = FALSE, method = "BFGS", control=list(fnscale=-1)))
fit$value

## compute likelihood at the true parameter value (with the standard function)
y <- as.vector(Y)
y0 <- y <= 0
y1 <- y > 0

a0 <- as.factor(as.vector(col(Y)))
b0 <- as.factor(as.vector(row(Y)))
a1 <- model.matrix(~a0-1) * log(mu_pop)
colnames(a1) <- paste("a1_", 1:n, sep="")

X <- model.matrix(~b0 - 1)
dim(X)
W <- model.matrix(~a0 + a1 - 1)
dim(W)

loglik(parameters, y, y0, y1, X, W, ncol(X), ncol(W), 0, 0, binomial())
grad(parameters, y, y0, y1, X, W, ncol(X), ncol(W), 0, 0, binomial())

## optimize the likelihood with optim, starting from the true values
system.time(fit1 <- optim(fn = loglik, gr = grad, par = parameters, Y=y, Y0=y0, Y1=y1, X=X, W=W, 
                          kx=ncol(X), kw=ncol(W), offsetx=0, offsetz=0, linkobj=binomial(),
                          hessian = FALSE, method = "BFGS", control=list(fnscale=-1)))
fit1$value
## the memory efficient implementation is 18 time fasters!

### Estimate evaluation (note that I'm cheating because I'm using the true parameter values as starting point for optim)
J <- 10000
n <- 1000

mu_pop <- sample(means, J, replace=TRUE)

## homogeneous population == same mu across all samples
mu <- matrix(data=rep(mu_pop, n), ncol=n, nrow=J)

## sample from real GCC and length
l <- len[names(mu_pop)]
g <- gcc[names(mu_pop)]

## values based on BRAIN data
alpha0 <- rnorm(n, mean=2, sd=2)
alpha1 <- rnorm(n, mean=-1, sd=0.3)
alpha2 <- rnorm(n, mean=0, sd=0.15)
alpha3 <- rnorm(n, mean=1.15, sd=2)

logistic <- binomial()$linkinv

## generate pi = prob. of "detection"
pi <- sapply(1:n , function(i) {
  eta <- alpha0[i] + alpha1[i] * log(mu_pop) #+ alpha2[i] * log(l) + alpha3[i] * g
  return(logistic(eta))
})

## generate Z, indicator of expression
Z <- matrix(rbinom(n*J, 1, pi), ncol=n)
head(Z)

## generate Y, read counts (same dispersion for all genes)
Y <- matrix(data=0, ncol=n, nrow=J)
Y[Z==0] <- rnbinom(sum(Z==0), mu = mu[Z==0], size = 1)
head(Y)

## compute likelihood at the true parameter value (with the memory efficient function)
parameters <- c(log(mu_pop), alpha0, alpha1, log(1))

X <- matrix(rep(1, n), ncol=1)
dim(X)
W <- model.matrix(~log(mu_pop))
dim(W)

nbeta <- J
nalpha <- n*2

loglik_small(parameters, Y, Y>0, X, W, nbeta, nalpha, 0, 0, binomial())
grad_small(parameters, Y, Y>0, X, W, nbeta, nalpha, 0, 0, binomial())

## optimize the likelihood with optim, starting from the true values
system.time(fit <- optim(fn = loglik_small, gr = grad_small, par = parameters, Y=Y, Y1=Y>0, X=X, W=W, 
                         kx=nbeta, kw=nalpha, offsetx=0, offsetw=0, linkobj=binomial(),
                         hessian = FALSE, method = "BFGS", control=list(fnscale=-1)))
fit$value
## for 1,000 genes and 100 cells, optim takes 26 seconds
## for 10,000 genes and 1,000 cells, optim takes 34 minutes

betahat <- matrix(fit$par[1:nbeta], ncol=J, nrow=ncol(X))
muhat <- t(exp(X %*% betahat))[,1]

plot(log(mu_pop), log(muhat) - log(mu_pop))
lines(lowess(log(mu_pop), log(muhat) - log(mu_pop)), col=2, lwd=2)
abline(h=0, lty=2)

alphahat <- matrix(fit$par[(nbeta + 1):(nbeta + nalpha)], nrow=ncol(W), ncol=n, byrow=TRUE)
etahat <- W %*% alphahat
pihat <- logistic(etahat)

plot(alpha0, alphahat[1,] - alpha0)
lines(lowess(alpha0, alphahat[1,] - alpha0), col=2, lwd=2)
abline(h=0, lty=2)

plot(alpha1, alphahat[2,] - alpha1)
lines(lowess(alpha1, alphahat[2,] - alpha1), col=2, lwd=2)
abline(h=0, lty=2)


plot(sort(pi[,1]), pihat[order(pi[,1]),1], type='l', ylim=c(0, 1), xlim=c(0, 1))
for(i in 1:n) {
  lines(sort(pi[,i]), pihat[order(pi[,i]),i])
}
lines(sort(rowMeans(pi)), rowMeans(pihat)[order(rowMeans(pi))], col=2, lwd=2)

plot(pi[,1], pihat[,1])
abline(0, 1)

plot(pi[,4], pihat[,4])
abline(0, 1)

plot(colMeans(Y), colMeans(pi-pihat))
plot(rowMeans(pi), rowMeans(pi-pihat))

thetahat <- exp(fit$par[nbeta + nalpha + 1])
thetahat

