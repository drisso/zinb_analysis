logit <- binomial()$linkfun
logistic <- binomial()$linkinv

#' Log-likelihood function of a zero-inflated negative binomial model
#' 
#' This function computes the log-likelihood of a zero-inflated negative binomial model
#' 
#' This simplified model assumes that that pi depends only on log(mu) and not on other technical covariates.
#' It should be easy to extend this model to other technical covariates, but I suspect that it will make the
#' optimization slower.
#'  
#' @param parms a vector of parameters: should contain the values of beta, followed by those of alpha, followed by the 1/phi
#' @param Y the data matrix (genes in rows, cells in columns)
#' @param X the design matrix for the regression on mu (n x k_X)
#' @param kx the number of beta parameters
#' @param kw the number of alpha parameters
loglik <- function(parms, Y, X, kx, kw) {

  beta <- parms[1:kx]
  mu <- exp(X %*% beta)
  W <- model.matrix(~beta)
  
  alpha <- matrix(parms[(kx + 1):(kw+kx)], ncol=ncol(W), nrow=n)
  eta <- alpha %*% t(W)
  pi <- logistic(eta)

  theta <- parms[(kw + kx + 1):(kw+2*kx)]
  
  
  loglik0 <- log(pi + exp(log(1-pi) + dnbinom(0, size = theta, mu = mu, log = TRUE)))
  loglik1 <- log(1-pi) + dnbinom(t(Y), size = theta, mu = mu, log = TRUE)
  return(sum(loglik0[which(Y==0)]) + sum(loglik1[which(Y>0)]))
}

#' Log-likelihood function of a zero-inflated negative binomial model
#' 
#' This function computes the log-likelihood of a zero-inflated negative binomial model
#' 
#' This model let the user specify both a regression on mu and on pi. Note that pi will always depend on beta0
#' in addition to the covariates specified by the user. For this to be meaningful, we should consider a parametrization
#' of X for which beta0 is the overall mean (e.g., contrasts = "contr.sum").
#'  
#' @param parms a vector of parameters: should contain the values of beta, followed by those of alpha, followed by those of phi
#' @param Y the data matrix (genes in rows, cells in columns)
#' @param X the design matrix for the regression on mu (n x k_X)
#' @param W the design matrix for the regression on pi (J x k_W)
#' @param kx the number of beta parameters
#' @param kw the number of alpha parameters
loglikfull <- function(parms, Y, X, W, kx, kw) {
  
  beta <- matrix(parms[1:kx], ncol=ncol(X), nrow=J)
  mu <- exp(X %*% t(beta))
  
  W <- cbind(W, beta[,1])
  
  alpha <- matrix(parms[(kx + 1):(kw+kx)], ncol=ncol(W), nrow=n)
  eta <- alpha %*% t(W)
  pi <- logistic(eta)
  
  phi <- parms[(kw + kx + 1):(kw+2*kx)]
  theta <- matrix(1/theta, ncol=J, nrow=4, byrow=TRUE)
  
  loglik0 <- log(pi + exp(log(1-pi) + dnbinom(0, size = theta, mu = mu, log = TRUE)))
  loglik1 <- log(1-pi) + dnbinom(t(Y), size = theta, mu = mu, log = TRUE)
  return(sum(loglik0[which(Y==0)]) + sum(loglik1[which(Y>0)]))
}

## true expression: low, high, medium genes.
set.seed(12324)
n <- 4
J <- 9
mu_pop <- rep(c(10, 50, 100), each=3)
beta <- log(mu_pop)

## homogeneous population == same mu across all samples
mu <- matrix(data=rep(mu_pop, n), ncol=n, nrow=J)

## values based on BRAIN data
alpha0 <- rnorm(n, mean=2, sd=2)
alpha1 <- rnorm(n, mean=-1, sd=0.3)

## generate pi = prob. of "detection"
pi <- sapply(1:n , function(i) {
  eta <- alpha0[i] + alpha1[i] * beta
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
parameters <- c(beta, alpha0, alpha1, rep(1, J))

X <- matrix(rep(1, n), ncol=1)
dim(X)
W <- model.matrix(~log(mu_pop))
dim(W)
