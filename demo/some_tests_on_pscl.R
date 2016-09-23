library(pscl)
library(zinb)

# goal: simulate data K times with the same value of parameters (we take one gene only) 
# and check the variability of estimators
# 
# # sample size
# n=100  #cell number
# K=500  #number of repeats
# J=100  #number of genes
#     
# # matrix to store the results of simulations
# results=matrix(0,nrow=K,ncol=5)
# X1=2*runif(n)
# U1=3*(1+runif(n))
# for (k in 1:K){
# 
#     X=cbind(X1,U1)
#     
#     alphaM=c(1.5,0.8)
#     alphaPi=c(-0.3,-0.05)
#     
#     M=exp(X%*%alphaM)
#     Pi=binomial()$linkinv(X%*%alphaPi)
#     gene.exp=NULL
#     for(i in 1:n){
#         gene.exp[i]=rnbinom(1,mu=M[i],size=1)
#     }
#     gene.is1=NULL
#     for(i in 1:n){
#         gene.is1[i]=1-rbinom(1,size=1,prob=as.vector(Pi)[i])
#     }
#     gene.exp.zeroinf=gene.exp
#     gene.exp.zeroinf[gene.is1==0]=0
#     
#     tab=cbind(gene.exp.zeroinf,X)
#     colnames(tab)=c("expr","X1","X2")
#     testdata=data.frame(tab)
#     
# #     tab2=cbind(gene.exp,X)
# #     colnames(tab2)=c("expr","X1","X2")
# #     testdata2=data.frame(tab2)
#     
#     estimate=zeroinfl(expr~X1+X2-1|X1+X2-1,data=testdata,dist="negbin",link="logit")
#     results[k,1:4]=c(estimate$coefficients$count,estimate$coefficients$zero)
#     results[k,5]=estimate$theta
#     #glm.nb(expr~X1+X2-1,data=testdata2)
# }
# 
# #plot histograms for estimators of each parameter
# hist(results[,1],breaks=50,xlab="alpha_M (true = 1.5)",main="sample size n=150, 500 simulations")
# hist(results[,2],breaks=50)
# hist(results[,3],breaks=50,xlab="alpha_Pi (true = 0.3)",main="sample size n=150, 500 simulations")
# hist(results[,4],breaks=50,xlab="alpha_Pi (true = 0.05)",main="sample size n=150, 500 simulations")
# hist(results[,5],breaks=50,xlab="theta (true = 1.38)",main="sample size n=150, 500 simulations")
# 
# #########################################################################################################
# #   ESTIMATE U*V and U*W without design matrices
# #   data simulated from one factor model where U is a column, V and W are two lines, all of them will be supposed unknown
# #########################################################################################################
# 
# U=matrix(runif(n)*3+1,ncol=1)
# V=matrix(runif(J)*2,nrow=1)
# W=-matrix(runif(J)*(2.5+2*runif(J)),nrow=1)
# theta=1+runif(J)  #a known theta same for all genes for the moment 
# theta.matr=matrix(0,nrow=n,ncol=J)
# for(j in 1:length(theta)){
# theta.matr[,j]=theta[j]
# }
# #mean of NB
# M=exp(U%*%V)
# #zero inflation probabilities
# Pi=binomial()$linkinv(U%*%W)
# 
# gene.exp=NULL
# 
# for(i in 1:(n*J)){
#     gene.exp[i]=rnbinom(1,mu=as.vector(M)[i],size=as.vector(theta.matr)[i])
# }
# gene.exp=matrix(gene.exp,ncol=J)
# 
# gene.is1=NULL
# for(i in 1:(n*J)){
#     gene.is1[i]=1-rbinom(1,size=1,prob=as.vector(Pi)[i])
# }
# gene.is1=matrix(gene.is1,ncol=J)
# 
# gene.exp.zeroinf=gene.exp
# gene.exp.zeroinf[gene.is1==0]=0


# #initialization
# U.0=matrix(1,nrow=n,ncol=1)
# V.0=matrix(1,nrow=1,ncol=J)
# W.0=matrix(1,nrow=1,ncol=J)
# theta0=rep(1,J)
# #for this simulation : X.M=0, X.pi=0 alpha.M=0 alpha.pi=0
# 
# alt.number=30
# for (alt in 1:alt.number){
#     for (gene in 1:J){
#         if(min(gene.exp.zeroinf[,gene])==0){
#             estimate=zeroinfl(X1~X2-1|X2-1,data=data.frame(cbind(gene.exp.zeroinf[,gene],U.0)),dist="negbin",link="logit")
#             V.0[,gene]=estimate$coefficients$count
#             W.0[,gene]=estimate$coefficients$zero
#             theta0[gene]=estimate$theta
#             #if no zeros, fit negative binomial (would be better to include a test of zero inflation)
#         }else{
#             estimate=glm.nb(X1~X2-1,data=data.frame(cbind(gene.exp.zeroinf[,gene],U.0)),link="log")
#             V.0[,gene]=estimate$coefficients
#             W.0[,gene]=0
#             theta0[gene]=estimate$theta
#         }
#     }
# 
#     for(cell in 1:n){
#         U.0[cell,]=optim(fn=ziNegBin.U,gr=gradNegBin.U,i=cell,par=U.0[cell,],Y=gene.exp.zeroinf,
#                         theta=log(theta0),V=V.0,W=W.0,X.M=F,X.pi=F,
#                         alpha.M=F,alpha.pi=F,
#                         control=list(fnscale=-1),method="BFGS")$par        
#     }
# }


########################################################################
#      31 JANUARY PCA with rank =2
########################################################################
Nsim=50
zero.fraction=rep(0,Nsim)
L2.zinb=rep(0,Nsim)
L2.pca=rep(0,Nsim)

for(k in 1:Nsim){
n=100
J=100

U=matrix(0,nrow=n,ncol=2)
U[,1]=c(rnorm(n/2)*0.6+2,rnorm(n/2)*0.6+4)
#c(rep(1,n/2)+runif(n/2)*1.5,2.5+runif(n/2)*1.5)
U[,2]=c(rnorm(n/2)*0.6+3,rnorm(n/2)*0.6+5)
#c(1+runif(n/2)*1.5,rep(3,n/2)+runif(n/2))
#U=U/mean(U)
#plot(U)
V=matrix(runif(2*J)*2.2,nrow=2)
#W=-0.8*matrix(runif(2*J)*(2+2*runif(2*J)),nrow=2) #22% of zeros
#W=-45*matrix(runif(2*J)*(2+2*runif(2*J)),nrow=2) #10% of zeros
W=-0.4*matrix(runif(2*J)*(2+2*runif(2*J)),nrow=2)
theta=1+runif(J)  #a known theta same for all genes for the moment 
theta.matr=matrix(0,nrow=n,ncol=J)
for(j in 1:length(theta)){
    theta.matr[,j]=theta[j]
}
#mean of NB

M=exp(U%*%V)
#zero inflation probabilities
perturb=matrix(rnorm(n*2),ncol=2)
Pi=binomial()$linkinv((U+perturb)%*%W)

gene.exp=NULL

for(i in 1:(n*J)){
    gene.exp[i]=rnbinom(1,mu=as.vector(M)[i],size=as.vector(theta.matr)[i])
}
gene.exp=matrix(gene.exp,ncol=J)

gene.is1=NULL
for(i in 1:(n*J)){
    gene.is1[i]=1-rbinom(1,size=1,prob=as.vector(Pi)[i])
}
gene.is1=matrix(gene.is1,ncol=J)

gene.exp.zeroinf=gene.exp
gene.exp.zeroinf[gene.is1==0]=0

zero.fraction[k]=sum(gene.exp.zeroinf==0)/(n*J)*100
#sizes=1+runif(n)/2
#sizemat=matrix(rep(sizes,J),ncol=J)
#gene.exp.zeroinf=round(gene.exp.zeroinf*sizemat)

#initialization at random
#U.0=matrix(1+runif(2*n),nrow=n,ncol=2)
#V.0=matrix(1+runif(2*J),nrow=2,ncol=J)
#W.0=matrix(1,nrow=2,ncol=J)


#initialization with PCA
PCA.init=prcomp(log(1+gene.exp.zeroinf),center=TRUE,scale.=TRUE)
U.0=PCA.init$x[,1:2]
V.0=t(PCA.init$rotation[,1:2])
colnames(U.0)=NULL
colnames(V.0)=NULL
W.0=matrix(1,nrow=2,ncol=J)
theta0=rep(1,J)

#for this simulation : X.M=0, X.pi=0 alpha.M=0 alpha.pi=0

alt.number=15 #max number of alternations
total.lik=rep(0,alt.number)


# for (alt in 1:alt.number){
#     for (gene in 1:J){
#  #       if(sum(gene.exp.zeroinf[,gene]==0)>6){
#             estimate=zeroinfl(X1~X2+X3-1|X2+X3-1,data=data.frame(cbind(gene.exp.zeroinf[,gene],U.0)),dist="negbin",link="logit")
#             V.0[,gene]=estimate$coefficients$count
#             W.0[,gene]=estimate$coefficients$zero
#             theta0[gene]=estimate$theta
#             #if no zeros, fit negative binomial (would be better to include a test of zero inflation)
#  #       }else{
#   #          estimate=glm.nb(X1~X2+X3-1,data=data.frame(cbind(gene.exp.zeroinf[,gene],U.0)),link="log")
#    #         V.0[,gene]=estimate$coefficients
#     #        W.0[,gene]=0
#      #       theta0[gene]=estimate$theta
# #        }
#     if(gene %in% c(1:9)*1000){paste(gene)}
#     }
# 
#     for(cell in 1:n){
#         U.0[cell,]=optim(fn=ziNegBin.U,gr=gradNegBin.U,i=cell,par=U.0[cell,],Y=gene.exp.zeroinf,
#                          theta=log(theta0),V=V.0,W=W.0,X.M=F,X.pi=F,
#                          alpha.M=F,alpha.pi=F,
#                          control=list(fnscale=-1),method="BFGS")$par        
#     }
#     total.lik[alt]=sum(sapply(1:n,function(t) ziNegBin.U(U.0[t,],t,V.0,W.0,X.M=F,X.pi=F,alpha.M=F,alpha.pi=F,theta0,Y=gene.exp.zeroinf)))
#     alt=alt+1
# }

#UNE VERSION ALTERNATIVE AVEC PENALIZATION
for (alt in 1:alt.number){
    print(alt)
    #evaluate total likelihood before alternation num alt
    total.lik[alt]=sum(sapply(1:n,function(t) ziNegBin.U(U.0[t,],t,V.0,W.0,X.M=F,X.pi=F,alpha.M=F,alpha.pi=F,theta0,Y=gene.exp.zeroinf)))
    #if the increase in likelihood is smaller than 0.5%, stop maximization
    #if(alt>1){if(abs((total.lik[alt]-total.lik[alt-1])/total.lik[alt-1])<0.005)break}
    print(total.lik[alt])
    #if the increase in likelihood is smaller than 0.5%, stop maximization
    # if(alt>1){if(abs((total.lik[alt]-total.lik[alt-1])/total.lik[alt-1])<0.005)break}
    if(alt>1){if(abs(total.lik[alt]-total.lik[alt-1])<1)break}
    #optimization for V and W, by gene, theta is optimized also
    for (gene in 1:J){            
        estimate=optim(fn=ziNegBin,gr=gradNegBin,j=gene,U=U.0,par=c(V.0[,gene],W.0[,gene],log(theta0)[gene]),Y=gene.exp.zeroinf,
                       control=list(fnscale=-1),method="BFGS",epsilon=0.001)$par 
        V.0[,gene]=estimate[1:length(V.0[,gene])]
        W.0[,gene]=estimate[(length(V.0[,gene])+1):(length(V.0[,gene])+length(W.0[,gene]))]
        theta0[gene]=min(exp(estimate[length(V.0[,gene])+length(W.0[,gene])+1]),10)       
    }
    #optimization for U, by cell, V,W and theta are fixed 
    for(cell in 1:n){
        U.0[cell,]=optim(fn=ziNegBin.U,gr=gradNegBin.U,i=cell,par=U.0[cell,],Y=gene.exp.zeroinf,
                         theta=log(theta0),V=V.0,W=W.0,X.M=F,X.pi=F,
                         alpha.M=F,alpha.pi=F,
                         control=list(fnscale=-1),method="BFGS")$par        
    }
    alt=alt+1
}

#true matrix of log expressions
ref=U%*%V
#matrix of log expressions reconstructed by zinb
logM.zinb=U.0%*%V.0
#L2error of zinb
L2.zinb[k]=sqrt(sum((ref-logM.zinb)^2))

#do PCA
pca.reconstruct=prcomp(log(1+gene.exp.zeroinf),center=TRUE,scale.=TRUE)
#scaled matrix of log expressions reconstructed by first two PCs
logM.pca=pca.reconstruct$x[,1:2]%*%t(pca.reconstruct$rotation[,1:2])
#means of initial log-counts
means=apply(log(gene.exp.zeroinf+1),2,mean)
#standard deviations of initial log-counts
sds=apply(log(gene.exp.zeroinf+1),2,sd)
#go back to unscaled matrix of log counts: multiply by sd and add mean to each column
for(j in 1:ncol(gene.exp.zeroinf)){
    logM.pca[,j]=logM.pca[,j]*sds[j]+means[j]
}
#L2error of pca
L2.pca[k]=sqrt(sum((ref-logM.pca)^2))
paste(k)
}
#################################################################################
#######         ACCESS PERFORMANCE IN TERMS OF CLUSTERING QUALITY
#################################################################################
#function which performs kmeans and returns rand and silhouette indices
partition.quality=function(truelabels,true.U,nclusters,pca.result,zinb.result){
    pca.partition=kmeans(pca.result,centers=nclusters)
    zinb.partition=kmeans(zinb.result,centers=nclusters)
    #compute matrices of pairwise distances for each method and for true U
    pca.pairwisedist=dist(pca.result,method="euclidean")
    zinb.pairwisedist=dist(zinb.result,method="euclidean")
    true.pairwisedist=dist(true.U,method="euclidean")
    #compute rand index of similarity between the true and the obtained partition 
    require("mclust")
    rand.pca=adjustedRandIndex(truelabels, pca.partition$cluster)
    rand.zinb=adjustedRandIndex(truelabels, zinb.partition$cluster)
    #compute silhouette statistics for each cluster and overall
    #silhouette s(i) of each observation i measures how well it is clustered
    #silhouette of cluster is mean s(i) over all points in cluster
    require("cluster")
    #compute silhouette objects for zinb, pca and true U
    silhouette.pca=silhouette(pca.partition$cluster,pca.pairwisedist)
    silhouette.zinb=silhouette(zinb.partition$cluster,zinb.pairwisedist)
    silhouette.bench=silhouette(truelabels,true.pairwisedist)
    quality.measures=list(rand.pca=rand.pca,rand.zinb=rand.zinb,silhouette.pca=summary(silhouette.pca),
                          silhouette.zinb=summary(silhouette.zinb),silhouette.bench=summary(silhouette.bench))
}

#################################################################################
#######         ACCESS PERFORMANCE IN TERMS OF PAIRWISE MATRICES
#################################################################################
#function which returns correlations of pairwise distances for each method and of "true distances"
pairwise.quality=function(true.U,pca.result,zinb.result){
    pca.pairwisedist=dist(pca.result,method="euclidean")
    zinb.pairwisedist=dist(zinb.result,method="euclidean")
    true.pairwisedist=dist(true.U,method="euclidean")
    pca.cor=cor(as.vector(pca.pairwisedist),as.vector(true.pairwisedist))
    zinb.cor=cor(as.vector(zinb.pairwisedist),as.vector(true.pairwisedist))
}

    
zero.fraction.22=zero.fraction
L2.pca.22=L2.pca
L2.zinb.22=L2.zinb

zero.fraction.10=zero.fraction
L2.pca.10=L2.pca
L2.zinb.10=L2.zinb

zero.fraction.30=zero.fraction
L2.pca.30=L2.pca
L2.zinb.30=L2.zinb

hist(L2.pca.30/L2.zinb.30, xlab="Ratio of PCA L2 error over ZINB L2 error")
hist(zero.fraction.30)

frame10=data.frame(val=c(L2.zinb.10,L2.pca.10),lev=c(rep("ZINB",50),rep("PCA",50)))
bp=boxplot(val ~ lev, frame10, range = 0.3, varwidth = TRUE,
           col = "skyblue",main="L2 errors, 10% of zeros")

frame20=data.frame(val=c(L2.zinb.22,L2.pca.22),lev=c(rep("ZINB",50),rep("PCA",50)))
bp=boxplot(val ~ lev, frame20, range = 0.3, varwidth = TRUE,
           col = "skyblue",main="L2 errors, 20% of zeros")

frame30=data.frame(val=c(L2.zinb.30,L2.pca.30),lev=c(rep("ZINB",50),rep("PCA",50)))
bp=boxplot(val ~ lev, frame30, range = 0.3, varwidth = TRUE,
           col = "skyblue",main="L2 errors, 30% of zeros")

framebis=data.frame(val=c(L2.pca.10/L2.zinb.10,L2.pca.22/L2.zinb.22,L2.pca.30/L2.zinb.30),lev=c(rep("10%",50),rep("20%",50),rep("30%",50)))
bp=boxplot(val ~ lev, framebis, range = 0.3, varwidth = TRUE,
           col = "skyblue",main="Ratio of L2 errors: PCA over ZINB",xlab="fraction of zeros")

plot(c(zero.fraction.10,zero.fraction.22,zero.fraction.30),c(L2.pca.10/L2.zinb.10,L2.pca.22/L2.zinb.22,L2.pca.30/L2.zinb.30))
c(zero.fraction.10,zero.fraction.22,zero.fraction.30),c(L2.pca.10/L2.zinb.10,L2.pca.22/L2.zinb.22,L2.pca.30/L2.zinb.30)
fr <- data.frame(val = log10(1+Mtimedata[gene,]),
                 lev = factor(daybis))
bp <- boxplot(val ~ lev, fr, range = 0.3, varwidth = TRUE,
              col = "skyblue", ylab = "boxplots of expression values",xlab=("Days"), ylim = c(min(fr$val)-0.1, max(fr$val) + 0.1), 
              main = gene)
plot(main="dd",c(zero.fraction.10,zero.fraction.22,zero.fraction.30),c(L2.pca.10/L2.zinb.10,L2.pca.22/L2.zinb.22,L2.pca.30/L2.zinb.30))


#########################################################################################################
################    FLUIDIGM data set
#########################################################################################################

install.packages("/Users/svetlanagribkova/GitHub/scRNAseq",repos=NULL,type="source")
library(scRNAseq)
data("fluidigm")
fluidigm
#object consisting of several data sets
assayNames(fluidigm) #"counts" and "fpkm"
dim(fluidigm)
#?assays gives information on mode of accession of data
#assays(...) to pick assays and then $... to pick the data I need

names(colData(fluidigm))
length(colData(fluidigm)$Biological_Condition)
sum(colData(fluidigm)$Coverage_Type=="High")#65 cells
test=assays(fluidigm)$counts[,which(colData(fluidigm)$Coverage_Type=="High")]
colData(fluidigm)$Coverage_Type=="High"

ncol(test)
nrow(test)#26262
#find genes to be filtered out
all0=apply(test>30,1,sum)<40
sum(!all0) #6995
#exclude genes with too low expressions
fluidigm2=test[!all0,]

######## SIZE FACTORS ESTIMATION  ############
#to estimate size factors, find genes never 0
all1=(apply(fluidigm2==0,1,sum)==0) 
sum(all1) #81
#estimate size factors
sizeFgenes=fluidigm2[all1==1,]
#calculate geometric means of lines
geom.means.lines=apply(sizeFgenes^(1/ncol(sizeFgenes)),1,prod)
#divide each line by its geom mean
for(i in 1:nrow(sizeFgenes)){
    sizeFgenes[i,]=sizeFgenes[i,]/geom.means.lines[i]
}
#compute size factors
size.factors=apply(sizeFgenes,2,median)
###############################################

#initialization with PCA
gene.exp.zeroinf=t(fluidigm2)#[,500:1000]
n=nrow(gene.exp.zeroinf)
J=ncol(gene.exp.zeroinf)
PCA.init=prcomp(log(1+gene.exp.zeroinf),center=TRUE,scale.=TRUE)
plot(PCA.init$x[,1:2])
U.0=PCA.init$x[,1:2]
V.0=t(PCA.init$rotation[,1:2])
colnames(U.0)=NULL
colnames(V.0)=NULL
W.0=matrix(1,nrow=2,J)
theta0=rep(1,ncol(gene.exp.zeroinf))

alt.number=25 #max number of alternations
total.lik=rep(0,alt.number)#vector to stock likelihoods

#OPTIMIZATION PROCEDURE
for (alt in 1:alt.number){
    print(alt)
    #evaluate total likelihood before alternation num alt
    total.lik[alt]=sum(sapply(1:n,function(t) ziNegBin.U(U.0[t,],t,V.0,W.0,X.M=F,X.pi=F,alpha.M=F,alpha.pi=F,theta0,Y=gene.exp.zeroinf)))
    print(total.lik[alt])
    #if the increase in likelihood is smaller than 0.5%, stop maximization
   # if(alt>1){if(abs((total.lik[alt]-total.lik[alt-1])/total.lik[alt-1])<0.005)break}
   if(alt>1){if(abs(total.lik[alt]-total.lik[alt-1])<1)break}
    #optimization for V and W, by gene, theta is optimized also
    for (gene in 1:J){        
        estimate=optim(fn=ziNegBin,gr=gradNegBin,j=gene,U=U.0,par=c(V.0[,gene],W.0[,gene],log(theta0)[j]),Y=gene.exp.zeroinf,
                       control=list(fnscale=-1),method="BFGS",epsilon=0.001)$par 
        V.0[,gene]=estimate[1:length(V.0[,gene])]
        W.0[,gene]=estimate[(length(V.0[,gene])+1):(length(V.0[,gene])+length(W.0[,gene]))]
        #in some cases theta is extremely high, put the limitation
        theta0[gene]=exp(estimate[length(V.0[,gene])+length(W.0[,gene])+1])
        #theta0[gene]=min(exp(estimate[length(V.0[,gene])+length(W.0[,gene])+1]),10^2)
    }
    #optimization for U, by cell, V,W and theta are fixed 
    print("Step2")
    for(cell in 1:n){
        U.0[cell,]=optim(fn=ziNegBin.U,gr=gradNegBin.U,i=cell,par=U.0[cell,],Y=gene.exp.zeroinf,
                         theta=log(theta0),V=V.0,W=W.0,X.M=F,X.pi=F,
                         alpha.M=F,alpha.pi=F,
                         control=list(fnscale=-1),method="BFGS")$par        
    }
    alt=alt+1
}


PCA.final=prcomp(U.0%*%V.0,center=TRUE,scale.=TRUE)
plot(PCA.final$x[,1:2])
plot(U.0)
which(U.0[,1]>10)














#center and standartize the true matrix of log expressions
ref=scale(U%*%V,center=TRUE,scale=TRUE)
#reconstruction by zinb of matrix of log expressions
logM.zinb=U.0%*%V.0
#center and standartize zinb reconstruction of matrix of log expressions
logM.zinb=scale(logM.zinb,center=TRUE,scale=TRUE)
#calculate PCA reconstruction based on first two PCs
pca.reconstruct=prcomp(log(1+gene.exp.zeroinf),center=TRUE,scale.=TRUE)
logM.pca=pca.reconstruct$x[,1:2]%*%t(pca.reconstruct$rotation[,1:2])
#calculate L2 error
L2.zinb=sqrt(sum((ref-logM.zinb)^2))
#calculate L2 error
L2.pca=sqrt(sum((ref-logM.pca)^2))

#another way: find the initial non scaled matrix 
means=apply(log(gene.exp.zeroinf+1),2,mean)
sds=apply(log(gene.exp.zeroinf+1),2,sd)
for(j in 1:ncol(gene.exp.zeroinf)){
    logM.pca[,j]=logM.pca[,j]*sds[j]+means[j]
}

L2.zinb.2=sqrt(sum((U%*%V-U.0%*%V.0)^2))
L2.pca.2=sqrt(sum((U%*%V-logM.pca)^2))





#compare distances between points
#put all on the same scale (of log)
dist.ref.zinb=dist(U%*%V)
dist.ref.zinb2=dist(scale(log(1+U%*%V)))
dist.ref.pca=dist(scale(U%*%V))

dist.zinb=dist(U.0)
dist.zinb2=dist(scale(log(1+U.0%*%V.0)))
dist.pca=dist(U.0.0)

cor(dist.ref.pca,dist.pca,method="spearman")
cor(dist.ref.zinb,dist.zinb,method="spearman")

plot(dist.ref,dist.zinb)
plot(dist.ref.pca,dist.pca)





plot(as.vector(U%*%V),as.vector(U.0%*%V.0))
sum(sapply(1:n,function(t) ziNegBin.U(U[t,],t,V,W,X.M=F,X.pi=F,alpha.M=F,alpha.pi=F,theta0,Y=gene.exp.zeroinf)))
sum(sapply(1:n,function(t) ziNegBin.U(U.0[t,],t,V.0,W.0,X.M=F,X.pi=F,alpha.M=F,alpha.pi=F,theta0,Y=gene.exp.zeroinf)))

plot(V,V.0)
plot(U,U.0)
plot(W,W.0)
mean(V[1,]/V.0[1,])
mean(U.0[1,]/U[1,])
mean(W.0[1,]/W[1,])
plot(as.vector(U%*%V),as.vector(U.0%*%V.0))
plot(as.vector(binomial()$linkinv(U%*%W)),as.vector(binomial()$linkinv(U.0%*%W.0)))

# COMPARE WITH PCA
#lines of U = loadings of cells on lines of V
plot(U.0[,1],U.0[,2])
plot(U[,1],U[,2])
pca.sim=prcomp(log(1+gene.exp.zeroinf),scale.=TRUE,center=TRUE)
plot(pca.sim$x[,c(1,2)])

points(pca.sim$x[U[,2]<2.5,c(1,2)],col="red")
#kmeans criterion
kPCA=kmeans(pca.sim$x[,c(1,2)],2)
kZINB=kmeans(U.0,2)
plot(pca.sim$x[,c(1,2)])
points(pca.sim$x[kPCA$cluster==1,c(1,2)],col="red")




sqrt(sum((U%*%V-U.0%*%V.0)^2))#error of zinb comapring to the true matrix
sqrt(sum((U%*%V)^2))
hist(U%*%V-U.0%*%V.0,main="U*V-its zinb estimation")
pca.reconstruct=pca.sim$x[,1:2]%*%t(pca.sim$rotation[,1:2])
#pca.reconstruct reconstructs the scaled centered matrix...

standard=function(matr){
for(j in 1:ncol(M)){
    matr[,j]=(matr[,j]-mean(matr[,j]))/sd(matr[,j])
}
return(matr)
}

hist(pca.reconstruct-standard(U%*%V),main="PCA recon. of stand. U*V - st. U*V")
mean(pca.reconstruct-standard(log(M)))
median(pca.reconstruct-standard(log(M)))

#test if PCA does correct thing
hist(pca.reconstruct-standard(log(1+gene.exp.zeroinf))) #OK
sqrt(sum((pca.reconstruct-standard(U%*%V))^2)) #error of PCA with 2 PCs = 74


#real data
test=read.table("rawcounts.txt",header=F)
test=test[2:nrow(test),3:98]
exp.mat=matrix(as.numeric(t(test)),nrow=nrow(t(test)))
exp.mat=exp.mat[,apply(exp.mat!=0,2,sum)>20]
dim(exp.mat)
exp.mat.pca=prcomp(exp.mat,scale.=TRUE,center=TRUE)
n=nrow(exp.mat)
J=ncol(exp.mat)
#normaliser
#m1=median(U.0[,1])
#m2=median(U.0[,2])
#sfcts=apply(cbind(U.0[,1]/m1,U.0[,2]/m2),1,mean)
#plot(U.0[,1]/sfcts,U.0[,2]/sfcts)

# old version (to keep) 23 January 2016
# alt.number=30
# for (alt in 1:alt.number){
#     for (gene in 1:J){
#         #       V.0[,gene]=optim(fn=ziNegBin,gr=gradNegBin,j=gene,
#         #       par=c(V[,gene]+0.3,W[,gene]+0.4,theta+0.1),
#         #       X.M=NULL,U=U,X.pi=NULL,Y=gene.exp.zeroinf,
#         #       control=list(fnscale=-1),method="BFGS")$par
#         if(min(gene.exp.zeroinf[,gene])==0){
#             estimate=zeroinfl(X1~X2-1|X2-1,data=data.frame(cbind(gene.exp.zeroinf[,gene],U.0)),dist="negbin",link="logit")
#             V.0[,gene]=estimate$coefficients$count
#             W.0[,gene]=estimate$coefficients$zero
#             theta0[gene]=theta[gene]#estimate$theta
#             #if no zeros, fit negative binomial (would be better to include a test of zero inflation)
#         }else{
#             estimate=glm.nb(X1~X2-1,data=data.frame(cbind(gene.exp.zeroinf[,gene],U.0)),link="log")
#             V.0[,gene]=estimate$coefficients
#             W.0[,gene]=0
#             theta0[gene]=estimate$theta
#         }
#     }
#     #    alpha.M0=par.est[1:kxM,]
#     #    V0=par.est[(kxM+1):(kxM+p),]
#     #    alpha.pi0=par.est[(kxM+p+1):(kxM+p+kxPi),]
#     #    W0=par.est[(kxM+p+kxPi+1):(kxM+2*p+kxPi),]
#     #    thetas=par.est[kxM+2*p+kxPi+1,]
#     for(cell in 1:n){
#         #    U0[cell,]=optim(fn=ziNegBin.U,gr=gradNegBin.U,i=cell,par=U0[cell,],counts=Y,
#         #    thetas=thetas,V=V0,W=W0,
#         #    offset.M=t(alpha.M0)%*%t(X.M),offset.pi=t(alpha.pi0)%*%t(X.pi),
#         #    control=list(fnscale=-1),method="BFGS")$par
#         U.0[cell,]=optim(fn=ziNegBin.U,gr=gradNegBin.U,i=cell,par=U.0[cell,],Y=gene.exp.zeroinf,
#                          theta=theta0,V=V.0,W=W.0,X.M=F,X.pi=F,
#                          alpha.M=F,alpha.pi=F,
#                          control=list(fnscale=-1),method="BFGS")$par
#         
#     }
# }





#test that my functions give the same result on these data
# matrix with one gene only
#Y=matrix(gene.exp.zeroinf,ncol=1)
#j=1
#same predictors for M and for Pi
#Z=X

#optim(fn=ziNegBin,gr=gradNegBin,j=j,X.M=matrix(X1,ncol=1),X.pi=matrix(X1,ncol=1),U=matrix(U1,ncol=1),Y=Y,par=c(1,1,0.5,0.1,1),control=list(fnscale=-1),method="BFGS")
#optim(fn=ziNegBin,gr=gradNegBin,j=j,X.M=NULL,X.pi=NULL,U=cbind(X1,U1),Y=Y,par=c(1,1,0.5,0.1,1),control=list(fnscale=-1),method="BFGS")


