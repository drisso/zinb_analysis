# Simulate data
###############

# TODO: need to set the seed of the random number generator if we want to be able to reproduce results:
set.seed(12324)

# Number of cells
n=100
# Number of genes
J=60
# Number of latent actors
p=2

# Generate the n*p matrix U of p latent factors for the cells
U=cbind(runif(n),runif(n))
# TODO: U has only two columns, however the number of latent factors is p. You should define U with p columns so that the code will work when we change p to some other values. For example:
# U <- matrix(runif(n*p),nrow=n)

#Simulate true U-coefficients for M
# TODO: same comment as for U, V should have p rows
V=rbind(runif(J),runif(J))

#design matrix X.M
X.M=matrix(0,nrow=n,ncol=p)
# TODO: I'm lost, here M has p columns, I thought p was for the latent factors in U and V?
X.M[,1]=c(rep(1,n/2),rep(0,n/2))
X.M[,2]=c(rep(0,n/2),rep(1,n/2))

#overall design matrix for M
X=cbind(X.M,U)

#coefficients for X.M
alpha.M=matrix(0,nrow=p,ncol=J)
alpha.M[1,]=5*runif(J)
alpha.M[2,]=6*runif(J)

#matrix M
M=exp(X.M%*%alpha.M+U%*%V)

#simulate true U-coefficients for Pi
W=-rbind(c(rep(1,J/2),rep(2,J/2)),c(rep(4,J/2),rep(3,J/2)))

#design matrix X.pi, take the same for the moment
X.pi=matrix(0,nrow=n,ncol=p)
X.pi[,1]=c(rep(1,n/2),rep(0,n/2))
X.pi[,2]=c(rep(0,n/2),rep(1,n/2))

#true coefficients for X.pi
alpha.pi=rbind(c(rep(0.01,J/2),rep(0.005,J/2)),c(rep(0.007,J/2),rep(0.02,J/2)))

#matrix of probabilities of dropout
Pi=binomial()$linkinv(X.pi%*%alpha.pi+U%*%W)
min(Pi)

#simulate negative binomial from M (matrix of expressions)
exprs=NULL
for(i in 1:(n*J)){
    exprs[i]=rnbinom(1,mu=as.vector(M)[i],size=log(3))
}
exprs=matrix(exprs,nrow=n)


indic=NULL
for(i in 1:(n*J)){
    indic[i]=1-rbinom(1,size=1,prob=as.vector(Pi)[i])
}
indic=matrix(indic,nrow=n)

Y=matrix(0,nrow=n,ncol=J)
Y[indic==T]=exprs[indic==T]
sum(Y==0)
#which(apply(Y==0,2,sum)==max(apply(Y==0,2,sum)))
#sum(Y[,17]==0)
sum(apply(Y!=0,2,sum)==0)
######################### end of data simulation  ##########################


#########################  beginning of the model estimation  ##############


#X_M design matrix for M, X_pi same for Pi
#U latent factor matrix

#matrix of covariates for M regression 
#(kx=nb of known M-factors + nb of latent factors)
# TODO: X was already defined, no need to redefine it
X=cbind(X.M,U)

#matrix of covariates for Pi regression
#(kz=nb of known Pi-factors + nb of latent factors)
# TODO: to be consistant, Z should be defined in the same time as X
Z=cbind(X.pi,U)
#Z=X.pi

#for each gene, parms is a vector of size kx+kz+1 
#(total nb of known + nb of latent factors + theta)

# TODO: the 2 lines below do not seem to be needed (not used later)
offsetx=0 #no offset for the moment
offsetz=0 #no offset for the moment

#Y is matrix of counts n times J
## TODO: the 5 lines below do not seem to be needed (either already defined, or not used later)
n <- nrow(Y) #number of cells
kx <- NCOL(X) #number of M-factors
kz <- NCOL(Z) #number of Pi-factors
Y0 <- Y <= 0 #==1 if counts is 0
Y1 <- Y > 0 #==1 if count is not 0

#test of the first part : optimization only wrt the right part
trueparam=rbind(alpha.M,V,alpha.pi,W,log(6))
optim(fn=ziNegBin,gr=gradNegBin,j=10,epsilon=0,X=X,Z=Z,Y=Y,par=c(1,1,0.5,0.1,1,1,1,1,1),control=list(fnscale=-1),method="BFGS")$par
trueparam[,10]


#optim(fn=ziNegBin.U,gr=gradNegBin.U,i=1,par=c(0.7,0.5),counts=Y,
#     thetas=rep(0.3266343,20),V=V,W=W,
#   offset.M=t(alpha.M)%*%t(X.M),offset.pi=t(alpha.pi)%*%t(X.pi),
#   control=list(fnscale=-1,trace=1),method="BFGS")

#Put everything together
#alt.number - number of alternations
alt.number=50
#initialization
p=2 # TODO: no need to redefine p?
U0=matrix(3,nrow=n,ncol=p)
V0=matrix(5,nrow=p,ncol=J)
alpha.M0=matrix(1,nrow=p,ncol=J)
alpha.pi0=matrix(1,nrow=p,ncol=J)
W0=matrix(5,nrow=p,ncol=J)
theta0=0.7
kxM=ncol(X.M)
kxV=p
kxPi=ncol(X.pi)
kxW=p
par.est=matrix(0,nrow=kxM+kxPi+2*p+1,ncol=J)
for (alt in 1:alt.number){
    for (gene in 1:J){
        par.est[,gene]=optim(fn=ziNegBin,Y=Y,j=gene,gr=gradNegBin,
                             par=c(alpha.M0[,gene],V0[,gene],alpha.pi0[,gene],W0[,gene],theta0),
                             X=cbind(X.M,U0),Z=cbind(X.pi,U0),
                             control=list(fnscale=-1),method="BFGS")$par
    }
    alpha.M0=par.est[1:kxM,]
    V0=par.est[(kxM+1):(kxM+p),]
    alpha.pi0=par.est[(kxM+p+1):(kxM+p+kxPi),]
    W0=par.est[(kxM+p+kxPi+1):(kxM+2*p+kxPi),]
    thetas=par.est[kxM+2*p+kxPi+1,]
    for(cell in 1:n){
        U0[cell,]=optim(fn=ziNegBin.U,gr=gradNegBin.U,i=cell,par=U0[cell,],counts=Y,
                        thetas=thetas,V=V0,W=W0,
                        offset.M=t(alpha.M0)%*%t(X.M),offset.pi=t(alpha.pi0)%*%t(X.pi),
                        control=list(fnscale=-1),method="BFGS")$par
    }
}

# TODO 