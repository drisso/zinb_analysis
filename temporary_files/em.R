#' EM algorithm for zero-inflated negative binomial model
#' 
#' This function implements an EM algorithm for a zero-inflated negative binomial model
#' 
#' This function implements an expectation-maximization algorithm for a zero-inflated
#' model, with a constant gene-specific mean parameter, gene-specific dispersion and
#' a mixing probability that depends solely on the mean.
#' 
#' @param Y the data matrix (genes in rows, cells in columns)
#' 
#' @return a list with the estimates of mu, pi, theta and with the coefficients of the logistic regression
#' 
em <- function(Y, maxiter=100, verbose=FALSE) {

  n <- ncol(Y)
  J <- nrow(Y)
  
  # 0. Initial estimate of mu, pi and z
  mu0 <- rowMeans(Y)

  pi0 <- apply(Y, 2, function(x) {
    y <- as.numeric(x<=0)
    fit <- glm(y~log(mu0), family = binomial)
    return(fitted.values(fit))
  })

  zhat <- pi0/(pi0 + (1 - pi0) * dnbinom(0, size = 1, mu = matrix(mu0, nrow=J, ncol=n)))
  zhat[Y>0] <- 0
  thetahat <- rep(1, J)

  coefs_mu <- log(mu0)
  coefs_pi <- sapply(seq_len(n), function(i) {
    fit <- suppressWarnings(glm(zhat[,i]~log(mu0), family = binomial(link = logit)))
    return(coefficients(fit))
  })

  X <- matrix(rep(1, n), ncol=1)
  W <- model.matrix(~log(mu0))
  linkobj <- binomial()

  ll_new <- loglik_small(c(coefs_mu, coefs_pi[1,], coefs_pi[2,], log(thetahat)), Y, Y>0, X, W, J, n*2, 0, 0, linkobj)
  ll_old <- 2 * ll_new

  ## EM iteration
  iter <- 0
  while (abs((ll_old - ll_new)/ll_old) > 1e-4 & iter<maxiter) {
      ll_old <- ll_new
      
      fit_mu <- lapply(seq_len(J), function(i) {
        fit <- suppressWarnings(glm.nb(Y[i,] ~ 1, weights = (1 - zhat[i,]), init.theta = thetahat[i], start=coefs_mu[i]))
        return(fit)
      })
      coefs_mu <- sapply(fit_mu, coefficients)
      muhat <- exp(coefs_mu)
      thetahat <- sapply(fit_mu, function(x) x$theta)
    
      fit_pi <- lapply(seq_len(n), function(i) {
        fit <- suppressWarnings(glm(zhat[,i]~log(muhat), family = binomial(link = logit), start=coefs_pi[,i]))
        return(fit)
      })
      coefs_pi <- sapply(fit_pi, coefficients)
      pihat <- sapply(fit_pi, fitted.values)
      
      zhat <- pihat/(pihat + (1 - pihat) * dnbinom(0, size = matrix(thetahat, nrow=J, ncol=n), mu = matrix(muhat, nrow=J, ncol=n)))
      zhat[Y>0] <- 0
    
      W <- model.matrix(~log(muhat))
      ll_new <- loglik_small(c(coefs_mu, coefs_pi[1,], coefs_pi[2,], log(thetahat)), Y, Y>0, X, W, J, n*2, 0, 0, linkobj)
      
      if(verbose) {
        print(ll_new)
      }
      
      iter <- iter + 1
  }
  
  convergence <- 0
  if(iter == maxiter) {
    convergence <- 1
  }
  
  return(list(muhat=muhat, pihat=pihat, thetahat=thetahat, coefs=coefs_pi, convergence=convergence))
}

#' EM algorithm for zero-inflated negative binomial model
#' 
#' This function implements an EM algorithm for a zero-inflated negative binomial model
#' 
#' This function implements an expectation-maximization algorithm for a zero-inflated
#' model, with a constant gene-specific mean parameter, gene-specific dispersion and
#' a mixing probability that depends solely on the mean.
#' 
#' @param Y the data matrix (genes in rows, cells in columns)
#' @param offset the offset (useful for normalization)
#' 
#' @return a list with the estimates of mu, pi, theta and with the coefficients of the logistic regression
#' 
em_offset <- function(Y, offset, maxiter=100, verbose=FALSE) {
  
  n <- ncol(Y)
  J <- nrow(Y)
  
  # 0. Initial estimate of mu, pi and z
  mu0 <- rowMeans(Y)
  
  pi0 <- apply(Y, 2, function(x) {
    y <- as.numeric(x<=0)
    fit <- glm(y~log(mu0), family = binomial)
    return(fitted.values(fit))
  })
  
  zhat <- pi0/(pi0 + (1 - pi0) * dnbinom(0, size = 1, mu = matrix(mu0, nrow=J, ncol=n)))
  zhat[Y>0] <- 0
  thetahat <- rep(1, J)
  
  coefs_mu <- log(mu0)
  coefs_pi <- sapply(seq_len(n), function(i) {
    fit <- suppressWarnings(glm(zhat[,i]~log(mu0), family = binomial(link = logit)))
    return(coefficients(fit))
  })
  
  X <- matrix(rep(1, n), ncol=1)
  W <- model.matrix(~log(mu0))
  linkobj <- binomial()
  
  ll_new <- loglik_small(c(coefs_mu, coefs_pi[1,], coefs_pi[2,], log(thetahat)), Y, Y>0, X, W, J, n*2, offset, 0, linkobj)
  ll_old <- 2 * ll_new
  
  ## EM iteration
  iter <- 0
  while (abs((ll_old - ll_new)/ll_old) > 1e-4 & iter<maxiter) {
    ll_old <- ll_new
    
    fit_mu <- lapply(seq_len(J), function(i) {
      tryCatch(fit <- suppressWarnings(glm.nb(Y[i,] ~ 1 + offset(offset), weights = (1 - zhat[i,]), init.theta = thetahat[i], start=coefs_mu[i])),
               error=function(e){cat(paste("Gene", i, "failed\n"))})
      if(any(class(fit)=="glm")) {
        return(fit)
      } else {
        return(NULL)
      }
    })
    
    coefs_mu <- sapply(fit_mu, coefficients)
    muhat <- exp(coefs_mu)
    thetahat <- sapply(fit_mu, function(x) x$theta)
    
    fit_pi <- lapply(seq_len(n), function(i) {
      fit <- suppressWarnings(glm(zhat[,i]~log(muhat), family = binomial(link = logit), start=coefs_pi[,i]))
      return(fit)
    })
    coefs_pi <- sapply(fit_pi, coefficients)
    pihat <- sapply(fit_pi, fitted.values)
    
    zhat <- pihat/(pihat + (1 - pihat) * dnbinom(0, size = matrix(thetahat, nrow=J, ncol=n), mu = matrix(muhat, nrow=J, ncol=n)))
    zhat[Y>0] <- 0
    
    W <- model.matrix(~log(muhat))
    ll_new <- loglik_small(c(coefs_mu, coefs_pi[1,], coefs_pi[2,], log(thetahat)), Y, Y>0, X, W, J, n*2, offset, 0, linkobj)
    
    if(verbose) {
      print(ll_new)
    }
    
    iter <- iter + 1
  }
  
  convergence <- 0
  if(iter == maxiter) {
    convergence <- 1
  }
  
  return(list(muhat=muhat, pihat=pihat, thetahat=thetahat, coefs=coefs_pi, convergence=convergence))
}
