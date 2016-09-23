logit <- binomial()$linkfun
logistic <- binomial()$linkinv

#' Log-likelihood function of the zero-inflated negative binomial model
#' 
#' This function computes the log-likelihood of a standard regression model
#' 
#' This is a standard implementation of the log-likelihood of a 
#' zero-inflated negative binomial regression model.
#' The problem with this implementation is that, when applied to our model, it requires
#' the design matrices to be (n*J) x k_W and (n*J) x k_X, where p = n * k_W + J * k_X + J,
#' n is the number of cells, J the number of genes, k_W is the number of covariates
#' of the logistic regression, k_X is the number of covariates in the log-linear regression.
#' Moreover, this model assumes only one global dispersion parameter, but this should be easy to change.
#' 
#' @param parms a vector of parameters: should contain the values of beta, followed by those of alpha, followed by the log(1/phi)
#' @param Y the data matrix (genes in rows, cells in columns)
#' @param Y0 an indicator of Y==0
#' @param Y1 an indicator of Y==1
#' @param X the design matrix for the regression on mu
#' @param W the design matrix for the regression on pi
#' @param kx the number of beta parameters
#' @param kw the number of alpha parameters
#' @param offsetx the offset for the regression on X
#' @param offsetw the offset for the regression on W
#' @param linkobj the link function object for the regression on pi (typically the result of binomial())
loglik <- function(parms, Y, Y0, Y1, X, W, kx, kw, offsetx, offsetz, linkobj) {
  mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
  phi <- as.vector(linkobj$linkinv(W %*% parms[(kx + 1):(kx + kw)] + offsetz))
  theta <- exp(parms[(kx + kw) + 1])
  loglik0 <- log(phi + exp(log(1 - phi) + suppressWarnings(dnbinom(0, size = theta, mu = mu, log = TRUE))))
  loglik1 <- log(1 - phi) + suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))
  loglik <- sum(loglik0[Y0]) + sum(loglik1[Y1])
  return(loglik)
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
#' @param Y0 an indicator of Y==0
#' @param Y1 an indicator of Y==1
#' @param X the design matrix for the regression on mu
#' @param W the design matrix for the regression on pi
#' @param kx the number of beta parameters
#' @param kw the number of alpha parameters
#' @param offsetx the offset for the regression on X
#' @param offsetw the offset for the regression on W
#' @param linkobj the link function object for the regression on pi (typically the result of binomial())
grad <- function(parms, Y, Y0, Y1, X, W, kx, kw, offsetx, offsetz, linkobj) {
  eta <- as.vector(X %*% parms[1:kx] + offsetx)
  mu <- exp(eta)
  etaz <- as.vector(W %*% parms[(kx + 1):(kx + kw)] + offsetz)
  muz <- linkobj$linkinv(etaz)
  theta <- exp(parms[(kx + kw) + 1])
  clogdens0 <- dnbinom(0, size = theta, mu = mu, log = TRUE)
  dens0 <- muz * (1 - as.numeric(Y1)) + exp(log(1 - muz) + clogdens0)
  wres_count <- ifelse(Y1, Y - mu * (Y + theta)/(mu + theta), 
                       -exp(-log(dens0) + log(1 - muz) + clogdens0 + log(theta) - log(mu + theta) + log(mu)))
  wres_zero <- ifelse(Y1, -1/(1 - muz) * linkobj$mu.eta(etaz), 
                      (linkobj$mu.eta(etaz) - exp(clogdens0) * linkobj$mu.eta(etaz))/dens0)
  wres_theta <- theta * ifelse(Y1, digamma(Y + theta) - digamma(theta) + log(theta) - log(mu + theta)
                               + 1 - (Y + theta)/(mu + theta), exp(-log(dens0) + log(1 - muz) + clogdens0)
                               * (log(theta) - log(mu + theta) + 1 - theta/(mu + theta)))
  return(colSums(cbind(wres_count * X, wres_zero * W, wres_theta)))
}
