#' Log-likelihood function of the zero-inflated negative binomial model
#' 
#' This function computes the log-likelihood of a standard regression model
#' 
#' This is a (hopefully) memory-efficient implementation of the log-likelihood of a 
#' zero-inflated negative binomial regression model.
#' In this attempt, the design matrices don't have n*J rows, but n and J, respectively.
#' The computation is a bit slower, but the memory usage should be much smaller for
#' large J and n.
#' 
#' @param parms a vector of parameters: should contain the values of beta, followed by those of alpha, followed by the log(1/phi)
#' @param Y the data matrix (genes in rows, cells in columns)
#' @param Y1 a logical indicator of Y>0
#' @param X the design matrix for the regression on mu (n x k_X)
#' @param W the design matrix for the regression on pi (J x k_W)
#' @param kx the number of beta parameters
#' @param kw the number of alpha parameters
#' @param offsetx the offset for the regression on X
#' @param offsetw the offset for the regression on W
#' @param linkobj the link function object for the regression on pi (typically the result of binomial())
loglik_small <- function(parms, Y, Y1, X, W, kx, kw, offsetx, offsetw, linkobj) {
  
  J <- nrow(Y)
  n <- ncol(Y)
  
  beta <- matrix(parms[1:kx], ncol=J, nrow=ncol(X))
  mu <- t(exp(X %*% beta))
  
  alpha <- matrix(parms[(kx + 1):(kw+kx)], nrow=ncol(W), ncol=n, byrow=TRUE)
  eta <- W %*% alpha
  pi <- logistic(eta)
    
  # for now just one phi
  theta <- exp(parms[(kw + kx) + 1])
  #theta <- exp(parms[(kw + kx + 1):(kw + kx + J)])
  
  loglik0 <- log(pi + exp(log(1 - pi) + suppressWarnings(dnbinom(0, size = theta, mu = mu, log = TRUE))))
  loglik1 <- log(1 - pi) + suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))
  
  return(sum(loglik0[which(Y==0)]) + sum(loglik1[which(Y>0)]))
}

#' Gradient function of the zero-inflated negative binomial model
#' 
#' This function computes the gradient of a standard regression model
#' 
#' This is a standard implementation of the gradient of a 
#' zero-inflated negative binomial regression model.
#' The problem with this implementation is that, when applied to our model, it requires
#' the design matrices to be (n*J) x p, where p = n * (k_W + 1) + J * (k_X + 1) + J,
#' n is the number of cells, J the number of genes, k_W is the number of covariates
#' of the logistic regression, k_X is the number of covariates in the log-linear regression.
#' 
#' @param parms a vector of parameters: should contain the values of beta, followed by those of alpha, followed by the log(1/phi)
#' @param Y the data matrix (genes in rows, cells in columns)
#' @param Y1 a logical indicator of Y>0
#' @param X the design matrix for the regression on mu (n x k_X)
#' @param W the design matrix for the regression on pi (J x k_W)
#' @param kx the number of beta parameters
#' @param kw the number of alpha parameters
#' @param offsetx the offset for the regression on X
#' @param offsetw the offset for the regression on W
#' @param linkobj the link function object for the regression on pi (typically the result of binomial())
grad_small <- function(parms, Y, Y1, X, W, kx, kw, offsetx, offsetw, linkobj) {
  
  J <- nrow(Y)
  n <- ncol(Y)

  beta <- matrix(parms[1:kx], ncol=J, nrow=ncol(X))
  mu <- t(exp(X %*% beta))
  
  alpha <- matrix(parms[(kx + 1):(kw+kx)], nrow=ncol(W), ncol=n, byrow=TRUE)
  eta <- W %*% alpha
  pi <- logistic(eta)
  
  # for now just one phi
  theta <- exp(parms[(kw + kx) + 1])
  #theta <- exp(parms[(kw + kx + 1):(kw + kx + J)])

  clogdens0 <- dnbinom(0, size = theta, mu = mu, log = TRUE)
  dens0 <- pi * (1 - as.numeric(Y1)) + exp(log(1 - pi) + clogdens0)
  
  wres_count <- matrix(ifelse(Y1, Y - mu * (Y + theta)/(mu + theta), 
                       -exp(-log(dens0) + log(1 - pi) + clogdens0 + log(theta) - log(mu + theta) + log(mu))),
                       ncol=n)
  
  wres_zero <- t(matrix(ifelse(Y1, -1/(1 - pi) * linkobj$mu.eta(eta), 
                      (linkobj$mu.eta(eta) - exp(clogdens0) * linkobj$mu.eta(eta))/dens0),
                      ncol=n))
  
  wres_theta <- theta * ifelse(Y1, digamma(Y + theta) - digamma(theta) + log(theta) - log(mu + theta)
                               + 1 - (Y + theta)/(mu + theta), exp(-log(dens0) + log(1 - pi) + clogdens0)
                               * (log(theta) - log(mu + theta) + 1 - theta/(mu + theta)))
  
  return(c(as.vector(wres_count %*% X), as.vector(wres_zero %*% W), sum(wres_theta)))
}



#' Log-likelihood function of a zero-inflated negative binomial model
#' 
#' This function computes the log-likelihood of a zero-inflated negative binomial model
#' 
#' This model let the user specify both a regression on mu and on pi. Note that pi will always depend on beta0
#' in addition to the covariates specified by the user. For this to be meaningful, we should consider a parametrization
#' of X for which beta0 is the overall mean (e.g., contrasts = "contr.sum").
#'  
#' @param parms a vector of parameters: should contain the values of beta, followed by those of alpha, followed by those of log(1/phi)
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
  
  theta0 <- exp(parms[(kw + kx + 1):(kw+2*kx)])
  theta <- matrix(theta0, ncol=J, nrow=4, byrow=TRUE)
  
  loglik0 <- log(pi + exp(log(1 - pi) + suppressWarnings(dnbinom(0, size = theta, mu = mu, log = TRUE))))
  loglik1 <- log(1 - pi) + suppressWarnings(dnbinom(t(Y), size = theta, mu = mu, log = TRUE))
  
  return(sum(t(loglik0)[which(Y==0)]) + sum(t(loglik1)[which(Y>0)]))
}
