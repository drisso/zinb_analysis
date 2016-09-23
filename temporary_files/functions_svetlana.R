
# TODO: dangerous to define global variables since they may be overriden by other parts. Better to put them as function arguments when you need them

linkstr <- "logit"
linkobj <- make.link(linkstr)
linkinv <- linkobj$linkinv


# TODO: to document a function, do as in the 'standard.R' file of Davide, so that the documentation will be automatically created by roxygen 


#' Log-likelihood function of the zero-inflated negative binomial model for the estimation of "right parts" with 
#' log(M)=X.M*alpha.M+U*V; logit(Pi)=X.pi*alpha.pi+U*W; M=E[Y] and Y data matrix
#' @param parms a vector of parameters: j-th column of alpha.M, followed by j-th column of V, followed by j-th column of alpha.pi, followed by j-th column of W, followed by the log(1/phi)
#' @param j is the number of column (gene) to be optimized over
#' @param Y the data matrix (cells in rows, genes in columns)
#' @param X.M is the design matrix for M 
#' @param X.pi is the design matrix for Pi.
#' @param U is the current value for the latent factor matrix
#' @param offsetx the offset for the regression on M
#' @param offsetw the offset for the regression on Pi
#' @param linkobj the link function object for the regression on pi (typically the result of binomial())
#' @export


ziNegBin <- function(parms,j,X.M=NULL,X.pi=NULL,U,Y,epsilon=0,offsetx=0,offsetz=0) {
    X=cbind(X.M,U)
    Z=cbind(X.pi,U)
    kx=ncol(X)
    kz=ncol(Z)

    mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
    phi <- as.vector(linkinv(Z %*% parms[(kx + 1):(kx + kz)] + 
                                 offsetz))
    theta <- exp(parms[(kx + kz) + 1])
    loglik0 <- log(phi + exp(log(1 - phi) + suppressWarnings(dnbinom(0, 
                                                                     size = theta, mu = mu, log = TRUE))))
    loglik1 <- log(1 - phi) + suppressWarnings(dnbinom(Y[,j], 
                                                       size = theta, mu = mu, log = TRUE))
    Y0 <- Y <= 0 #==1 if counts is 0
    Y1 <- Y > 0 #==1 if count is not 0
    loglik <- sum(loglik0[Y0[,j]]) + sum(loglik1[Y1[,j]])
    loglik<-loglik-sum(parms[(kx + 1):(kx + kz)]^2)*epsilon
}

#' Gradient with respect to "right parts" of log-likelihood function of the zero-inflated negative binomial model with 
#' log(M)=X.M*alpha.M+U*V; logit(Pi)=X.pi*alpha.pi+U*W; M=E[Y] and Y data matrix
#' @param parms a vector of parameters: j-th column of alpha.M, followed by j-th column of V, followed by j-th column of alpha.pi, followed by j-th column of W, followed by the log(1/phi)
#' @param j is the number of column (gene) to be optimized over
#' @param Y the data matrix (cells in rows, genes in columns)
#' @param X.M is the design matrix for M 
#' @param X.pi is the design matrix for Pi.
#' @param U is the current value for the latent factor matrix
#' @param offsetx the offset for the regression on M
#' @param offsetw the offset for the regression on Pi
#' @param linkobj the link function object for the regression on pi (typically the result of binomial())
#' @export
gradNegBin <- function(parms,j,X.M=NULL,X.pi=NULL,U,Y,epsilon=0,offsetx=0,offsetz=0) {
    X=cbind(X.M,U)
    Z=cbind(X.pi,U)
    kx=ncol(X)
    kz=ncol(Z)
    eta <- as.vector(X %*% parms[1:kx] + offsetx)
    mu <- exp(eta)
    etaz <- as.vector(Z %*% parms[(kx + 1):(kx + kz)] + offsetz)
    muz <- linkinv(etaz)
    theta <- exp(parms[(kx + kz) + 1])
    Y0 <- Y <= 0 #==1 if counts is 0
    Y1 <- Y > 0 #==1 if count is not 0
    clogdens0 <- dnbinom(0, size = theta, mu = mu, log = TRUE)
    dens0 <- muz * (1 - as.numeric(Y1[,j])) + exp(log(1 - muz) + 
                                                      clogdens0)
    wres_count <- ifelse(Y1[,j], Y[,j] - mu * (Y[,j] + theta)/(mu + theta), 
                         -exp(-log(dens0) + log(1 - muz) + clogdens0 + log(theta) - 
                                  log(mu + theta) + log(mu)))
    wres_zero <- ifelse(Y1[,j], -1/(1 - muz) * linkobj$mu.eta(etaz), 
                        (linkobj$mu.eta(etaz) - exp(clogdens0) * linkobj$mu.eta(etaz))/dens0)
    wres_theta <- theta * ifelse(Y1[,j], digamma(Y[,j] + theta) - 
                                     digamma(theta) + log(theta) - log(mu + theta) + 1 - 
                                     (Y[,j] + theta)/(mu + theta), exp(-log(dens0) + log(1 - 
                                                                                             muz) + clogdens0) * (log(theta) - log(mu + theta) + 
                                                                                                                      1 - theta/(mu + theta)))
    c(colSums(wres_count * X), colSums(wres_zero * Z)-2*epsilon*parms[(kx + 1):(kx + kz)],sum(wres_theta))
}



######################################################################################
#Optimization with respect to U
######################################################################################

#' Log-likelihood function of the zero-inflated negative binomial model for the estimation of U in the model 
#' log(M)=X.M*alpha.M+U*V; logit(Pi)=X.pi*alpha.pi+U*W; M=E[Y] with Y data matrix
#' @param parms a vector of parameters: i-th line of U
#' @param i is the number of ligne (cell) to be optimized over
#' @param Y the data matrix (cells in rows, genes in columns)
#' @param X.M is the design matrix for M 
#' @param X.pi is the design matrix for Pi.
#' @param V,W,alpha.M,alpha.pi are the "right parts" to be given as arguments
#' @param theta vector of length J with gene specific dispersion parameters
#' @param linkobj the link function object for the regression on pi (typically the result of binomial())
#' @export
ziNegBin.U <- function(parms,i,V,W,X.M=F,X.pi=F,alpha.M=F,alpha.pi=F,theta,Y){
    if((X.M!=F)&(alpha.M!=F)){
        offset.M=(t(alpha.M)%*%t(X.M))[,i]
    }else{offset.M=0}
    if((X.pi!=F)&(alpha.pi!=F)){
        offset.pi=(t(alpha.pi)%*%t(X.pi))[,i]
    }else{offset.pi=0}
    theta=exp(theta)
    Y0 <- Y == 0
    Y1 <- Y != 0
    mu <- as.vector(exp(t(V) %*% parms + offset.M))
    phi <- as.vector(linkinv(t(W) %*% parms + 
                                 offset.pi))
    loglik0 <- log(phi + exp(log(1 - phi) + suppressWarnings(dnbinom(0, 
                                                                     size = theta, mu = mu, log = TRUE))))
    # counts[i,] is i-th column of t(counts)
    loglik1 <- log(1 - phi) + suppressWarnings(dnbinom(Y[i,], 
                                                       size = theta, mu = mu, log = TRUE))
    loglik <- sum(loglik0[Y0[i,]]) + sum(loglik1[Y1[i,]])
    loglik
}
#sum(sapply(seq(1000),function(i){ziNegBin.U(U[i,],i,V,W,t(alpha.M)%*%t(X.M),t(alpha.pi)%*%t(X.pi),counts,thetas)}))
# gradient function for optimization in U
#sum(sapply(seq(20),function(j){ziNegBin(trueparam[,j],j)}))
# gradient function for optimization in U

#' Gradient of the log-likelihood function of the zero-inflated negative binomial model for the estimation of U in the model 
#' log(M)=X.M*alpha.M+U*V; logit(Pi)=X.pi*alpha.pi+U*W; M=E[Y] with Y data matrix
#' @param parms a vector of parameters: i-th line of U
#' @param i is the number of ligne (cell) to be optimized over
#' @param Y the data matrix (cells in rows, genes in columns)
#' @param X.M is the design matrix for M 
#' @param X.pi is the design matrix for Pi.
#' @param V,W,alpha.M,alpha.pi are the "right parts" to be given as arguments
#' @param theta vector of length J with gene specific dispersion parameters
#' @param linkobj the link function object for the regression on pi (typically the result of binomial())
#' @export
gradNegBin.U <- function(parms,i,V,W,X.M=F,X.pi=F,alpha.M=F,alpha.pi=F,Y,theta) {
    Y0 <- Y == 0
    Y1 <- Y != 0
    if((X.M!=F)&(alpha.M!=F)){
        offset.M=(t(alpha.M)%*%t(X.M))[,i]
    }else{offset.M=0}
    if((X.pi!=F)&(alpha.pi!=F)){
        offset.pi=(t(alpha.pi)%*%t(X.pi))[,i]
    }else{offset.pi=0}
    eta <- as.vector(t(V) %*% parms + offset.M)
    mu <- exp(eta)
    etaz <- as.vector(t(W) %*% parms + offset.pi)
    muz <- linkinv(etaz)
    theta <- exp(theta)
    clogdens0 <- dnbinom(0, size = theta, mu = mu, log = TRUE)
    dens0 <- muz * (1 - as.numeric(Y1[i,])) + exp(log(1 - muz) + 
                                                            clogdens0)
    wres_count <- ifelse(Y1[i,],  Y[i,] - mu * (Y[i,] + theta)/(mu + theta), 
                         -exp(-log(dens0) + log(1 - muz) + clogdens0 + log(theta) - 
                                  log(mu + theta) + log(mu)))
    wres_zero <- ifelse(Y1[i,], -1/(1 - muz) * linkobj$mu.eta(etaz), 
                        (linkobj$mu.eta(etaz) - exp(clogdens0) * linkobj$mu.eta(etaz))/dens0)
    colSums(wres_count * t(V) + wres_zero * t(W))
}

#' Estimation of latent factors from a count matrix
#' @param datamatrix the count data matrix (cells in rows, genes in columns)
#' @param alt.number maximum number of iterations (default 25)
#' @export
zinb = function(datamatrix, alt.number=25){
    n=nrow(datamatrix)
    J=ncol(datamatrix)
    PCA.init=prcomp(log(1+datamatrix),center=TRUE,scale.=TRUE)
    #U.0=matrix(1,nrow=n,ncol=2)
    U.0=PCA.init$x[,1:2]
    #V.0=matrix(1,nrow=2,ncol=J)
    V.0=t(PCA.init$rotation[,1:2])
    colnames(U.0)=NULL
    colnames(V.0)=NULL
    W.0=matrix(1,nrow=2,J)
    theta0=rep(1,ncol(datamatrix))
    
   # for(l in 1:J){
    #theta0[l]=fitdistr(datamatrix[,l],densfun="negative binomial")$estimate[1]
    #}
    
    
    total.lik=rep(0,alt.number)
 
    #UNE VERSION ALTERNATIVE AVEC PENALIZATION
    for (alt in 1:alt.number){
        print(alt)
        #evaluate total likelihood before alternation num alt
        total.lik[alt]=sum(sapply(1:n,function(t) ziNegBin.U(U.0[t,],t,V.0,W.0,X.M=F,X.pi=F,alpha.M=F,alpha.pi=F,log(theta0),Y=datamatrix)))
        #if the increase in likelihood is smaller than 0.5%, stop maximization
        #if(alt>1){if(abs((total.lik[alt]-total.lik[alt-1])/total.lik[alt-1])<0.005)break}
        print(total.lik[alt])
        #if the increase in likelihood is smaller than 0.5%, stop maximization
        # if(alt>1){if(abs((total.lik[alt]-total.lik[alt-1])/total.lik[alt-1])<0.005)break}
        if(alt>1){if(abs(total.lik[alt]-total.lik[alt-1])<1)break}
        #optimization for V and W, by gene, theta is optimized also
        for (gene in 1:J){     
            #test=zeroinfl(X1~0+X2+X3,data.frame(cbind(datamatrix[,gene],U.0)),start=list(count=V.0[,gene],zero=W.0[,gene],theta=theta0[gene]),dist="negbin")
            #test=zeroinfl(X1~0+X2+X3,data.frame(cbind(datamatrix[,gene],U.0)),dist="negbin")
            #print(ziNegBin(parms=c(V.0[,gene],W.0[,gene],log(theta0)[gene]),gene,X.M=NULL,X.pi=NULL,U=U.0,Y=datamatrix,epsilon=0,offsetx=0,offsetz=0))
            #print(ziNegBin(parms=estimate,gene,X.M=NULL,X.pi=NULL,U=U.0,Y=datamatrix,epsilon=0,offsetx=0,offsetz=0))
            estimate=optim(fn=ziNegBin,gr=gradNegBin,j=gene,U=U.0,par=c(V.0[,gene],W.0[,gene],log(theta0)[gene]),Y=datamatrix,
                           control=list(fnscale=-1),method="BFGS",epsilon=0)$par 
            V.0[,gene]=estimate[1:length(V.0[,gene])]
            W.0[,gene]=estimate[(length(V.0[,gene])+1):(length(V.0[,gene])+length(W.0[,gene]))]
            theta0[gene]=exp(estimate[length(V.0[,gene])+length(W.0[,gene])+1])    
        }
        print(sum(sapply(1:n,function(t) ziNegBin.U(U.0[t,],t,V.0,W.0,X.M=F,X.pi=F,alpha.M=F,alpha.pi=F,log(theta0),Y=datamatrix))))
        #optimization for U, by cell, V,W and theta are fixed 
        for(cell in 1:n){
            U.0[cell,]=optim(fn=ziNegBin.U,gr=gradNegBin.U,i=cell,par=U.0[cell,],Y=datamatrix,
                             theta=log(theta0),V=V.0,W=W.0,X.M=F,X.pi=F,
                             alpha.M=F,alpha.pi=F,
                             control=list(fnscale=-1),method="BFGS")$par        
        }
    }    
    zinb.result <- list(U=U.0,V=V.0,W=W.0,theta=theta0)
}


zinbloglikelihood <- function(datamatrix, U, V, W, theta) {
    # Compute the loglikelihood for the simple model with only U, V and W. 
    # Choose the fastest direction to sum (rows or columns) depending on the size of the matrix
    # WE DONT USE IT BECAUSE IT IS FASTER TO JUST COMPUTE IT DIRECTLY WITHOUT CALLING THIS FUNCTION
    if (nrow(datamatrix) <= ncol(datamatrix)) {
        return( sum(sapply(1:nrow(datamatrix),function(t) ziNegBin.U(U[t,],t,V,W,X.M=F,X.pi=F,alpha.M=F,alpha.pi=F,log(theta),Y=datamatrix))) )
    } else {
        return( sum(sapply(1:ncol(datamatrix),function(t) ziNegBin(c(V[,t],W[,t],log(theta)[t]), j=t, U=U,Y=datamatrix))) )
    }
}
#' Estimation of latent factors from a count matrix
#' @param datamatrix the count data matrix (cells in rows, genes in columns)
#' @param k number of latent factors (default 2)
#' @param alt.number maximum number of iterations (default 25)
#' @param epsilon regularization parameter (default 0.1)
#' @param verbose print information (default TRUE)
#' @export
zinb2 = function(datamatrix, k=2, alt.number=25, epsilon=0.1, verbose=TRUE, use.old=FALSE){
    n=nrow(datamatrix)
    J=ncol(datamatrix)
    # Initialize U and V by PCA on log(count+1) matrix
    PCA.init=prcomp(log(1+datamatrix),center=TRUE,scale.=TRUE)
    U.0=PCA.init$x[,1:k]
    V.0=t(PCA.init$rotation[,1:k])
    colnames(U.0)=NULL
    colnames(V.0)=NULL
    # Initialize W and theta to 1
    W.0=matrix(0,nrow=k,ncol=J)
    theta0=rep(1,ncol(datamatrix))
    
    # for(l in 1:J){
    #theta0[l]=fitdistr(datamatrix[,l],densfun="negative binomial")$estimate[1]
    #}
    
    
    total.lik=rep(NA,alt.number)
    
    #UNE VERSION ALTERNATIVE AVEC PENALIZATION
    for (alt in 1:alt.number){
        if (verbose) {cat("Iteration ",alt,"\n",sep="")}

        #evaluate total likelihood before alternation num alt
        total.lik[alt]=sum(sapply(1:n,function(t) ziNegBin.U(U.0[t,],t,V.0,W.0,X.M=F,X.pi=F,alpha.M=F,alpha.pi=F,log(theta0),Y=datamatrix)))
        
        
        #if the increase in likelihood is smaller than 0.5%, stop maximization
        #if(alt>1){if(abs((total.lik[alt]-total.lik[alt-1])/total.lik[alt-1])<0.005)break}
        if (verbose) {cat("log-likelihood = ",total.lik[alt],"\n",sep="")}
        #if the increase in likelihood is smaller than 0.5%, stop maximization
        # if(alt>1){if(abs((total.lik[alt]-total.lik[alt-1])/total.lik[alt-1])<0.005)break}
        if(alt>1){if(abs(total.lik[alt]-total.lik[alt-1])<1)break}
        #optimization for V and W, by gene, theta is optimized also
        
        if (use.old) {
            # OLD VERSION (to remove once the new version is validated)
        ptm <- proc.time()
        for (gene in 1:J){     
            #test=zeroinfl(X1~0+X2+X3,data.frame(cbind(datamatrix[,gene],U.0)),start=list(count=V.0[,gene],zero=W.0[,gene],theta=theta0[gene]),dist="negbin")
            #test=zeroinfl(X1~0+X2+X3,data.frame(cbind(datamatrix[,gene],U.0)),dist="negbin")
            #print(ziNegBin(parms=c(V.0[,gene],W.0[,gene],log(theta0)[gene]),gene,X.M=NULL,X.pi=NULL,U=U.0,Y=datamatrix,epsilon=0,offsetx=0,offsetz=0))
            #print(ziNegBin(parms=estimate,gene,X.M=NULL,X.pi=NULL,U=U.0,Y=datamatrix,epsilon=0,offsetx=0,offsetz=0))
            estimate=optim(fn=ziNegBin,gr=gradNegBin,j=gene,U=U.0,par=c(V.0[,gene],W.0[,gene],log(theta0)[gene]),Y=datamatrix,control=list(fnscale=-1),method="BFGS",epsilon=0)$par 
            V.0[,gene]=estimate[1:length(V.0[,gene])]
            W.0[,gene]=estimate[(length(V.0[,gene])+1):(length(V.0[,gene])+length(W.0[,gene]))]
            theta0[gene]=exp(estimate[length(V.0[,gene])+length(W.0[,gene])+1])    
        }
        ptm <- proc.time()-ptm
        
        } else {
            # NEW version, faster and more precise
            ptm <- proc.time()
            estimate <- sapply(seq(J), function(gene) {optim(fn=ziNegBinNEW,gr=gradNegBinNEW,par=c(V.0[,gene],W.0[,gene],log(theta0)[gene]),U=U.0,Y=datamatrix[,gene],control=list(fnscale=-1),method="BFGS",epsilon=epsilon)$par })
            if (verbose) {print(proc.time()-ptm)}
            V.0 <- estimate[1:k,]
            W.0 <- estimate[(k+1):(2*k),]
            theta0 <- exp(estimate[(2*k+1),])
        }
       
        
        if (verbose) {cat("log-likelihood = ",sum(sapply(1:n,function(t) ziNegBin.U(U.0[t,],t,V.0,W.0,X.M=F,X.pi=F,alpha.M=F,alpha.pi=F,log(theta0),Y=datamatrix))),"\n",sep="")}
        
        #optimization for U, by cell, V,W and theta are fixed 
        for(cell in 1:n){
            U.0[cell,]=optim(fn=ziNegBin.U,gr=gradNegBin.U,i=cell,par=U.0[cell,],Y=datamatrix,
                             theta=log(theta0),V=V.0,W=W.0,X.M=F,X.pi=F,
                             alpha.M=F,alpha.pi=F,
                             control=list(fnscale=-1),method="BFGS")$par        
        }
    }    
    zinb.result <- list(U=U.0,V=V.0,W=W.0,theta=theta0)
}




#' Log-likelihood function of the zero-inflated negative binomial model for the estimation of "right parts" with 
#' log(M)=X.M*alpha.M+U*V; logit(Pi)=X.pi*alpha.pi+U*W; M=E[Y] and Y data matrix
#' @param parms a vector of parameters: j-th column of alpha.M, followed by j-th column of V, followed by j-th column of alpha.pi, followed by j-th column of W, followed by the log(1/phi)
#' @param j is the number of column (gene) to be optimized over
#' @param Y the data matrix (cells in rows, genes in columns)
#' @param X.M is the design matrix for M 
#' @param X.pi is the design matrix for Pi.
#' @param U is the current value for the latent factor matrix
#' @param offsetx the offset for the regression on M
#' @param offsetw the offset for the regression on Pi
#' @export
ziNegBinNEW <- function(parms,X.M=NULL,X.pi=NULL,U,Y,epsilon=0,offsetx=0,offsetz=0) {
        X=cbind(X.M,U)
        Z=cbind(X.pi,U)
        #    X=cbind(X.M,U)
#    Z=cbind(X.pi,U)
    kx=ncol(X)
    kz=ncol(Z)
    
    mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
    #    phi <- as.vector(linkinv(Z %*% parms[(kx + 1):(kx + kz)] + 
    #                                 offsetz))
    eta <- as.vector(Z %*% parms[(kx + 1):(kx + kz)] + offsetz)
    
    lognorm <- -eta - copula::log1pexp(-eta)
    
    theta <- exp(parms[(kx + kz) + 1])
    
    logPnb <- suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))
    
    sum(logPnb[Y>0]) + sum(logPnb[Y==0] + copula::log1pexp(eta[Y==0]-logPnb[Y==0])) + sum(lognorm) -sum(parms[(kx + 1):(kx + kz)]^2)*epsilon
    
}



#' Gradient with respect to "right parts" of log-likelihood function of the zero-inflated negative binomial model with 
#' log(M)=X.M*alpha.M+U*V; logit(Pi)=X.pi*alpha.pi+U*W; M=E[Y] and Y data matrix
#' @param parms a vector of parameters: j-th column of alpha.M, followed by j-th column of V, followed by j-th column of alpha.pi, followed by j-th column of W, followed by the log(1/phi)
#' @param j is the number of column (gene) to be optimized over
#' @param Y the data matrix (cells in rows, genes in columns)
#' @param X.M is the design matrix for M 
#' @param X.pi is the design matrix for Pi.
#' @param U is the current value for the latent factor matrix
#' @param offsetx the offset for the regression on M
#' @param offsetw the offset for the regression on Pi
#' @export
gradNegBinNEW <- function(parms,X.M=NULL,X.pi=NULL,U,Y,epsilon=0,offsetx=0,offsetz=0) {
        X=cbind(X.M,U)
        Z=cbind(X.pi,U)
        #    X=cbind(X.M,U)
#    Z=cbind(X.pi,U)
    kx=ncol(X)
    kz=ncol(Z)
    eta <- as.vector(X %*% parms[1:kx] + offsetx)
    mu <- exp(eta)
    etaz <- as.vector(Z %*% parms[(kx + 1):(kx + kz)] + offsetz)
    muz <- linkinv(etaz)
    logtheta <- parms[(kx + kz) + 1]
    theta <- exp(logtheta)
    Y0 <- Y <= 0 #==1 if counts is 0
    Y1 <- Y > 0 #==1 if count is not 0
    clogdens0 <- dnbinom(0, size = theta, mu = mu, log = TRUE)
    dens0 <- muz * (1 - as.numeric(Y1)) + exp(log(1 - muz) + 
                                                      clogdens0)
    wres_count <- ifelse(Y1, Y - mu * (Y + theta)/(mu + theta), 
                         -exp(-log(dens0) + log(1 - muz) + clogdens0 + logtheta - 
                                  log(mu + theta) + log(mu)))
    wres_zero <- ifelse(Y1, -1/(1 - muz) * linkobj$mu.eta(etaz), 
                        (linkobj$mu.eta(etaz) - exp(clogdens0) * linkobj$mu.eta(etaz))/dens0)
    wres_theta <- theta * ifelse(Y1, digamma(Y + theta) - 
                                     digamma(theta) + logtheta - log(mu + theta) + 1 - 
                                     (Y + theta)/(mu + theta), exp(-log(dens0) + log(1 - 
                                    muz) + clogdens0) * (logtheta - log(mu + theta) +
                                    1 - theta/(mu + theta)))
    c(colSums(wres_count * X), colSums(wres_zero * Z)-2*epsilon*parms[(kx + 1):(kx + kz)],sum(wres_theta))
}

