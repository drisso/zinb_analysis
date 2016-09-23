source("../R/standard.R")
source("../R/reg_smallmatrix.R")

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

## real data (we cannot use the fq normalization because otherwise the zeros are the same across all cells)
library(EDASeq)
uq <- betweenLaneNormalization(filtered, which="upper")
x <- as.factor(pData(brain01_counts[,colnames(filtered)])$Cell.type)
X <- model.matrix(~x)
dim(X)
uq0 <- uq
uq0[uq==0] <- NA

W <- model.matrix(~log(rowMeans(uq0, na.rm=TRUE)))
dim(W)

parameters <- c(log(rowMeans(uq0, na.rm=TRUE)), rep(0, nrow(filtered) + ncol(filtered)*2), log(1))

nbeta <- nrow(filtered)*2
nalpha <- ncol(filtered)*2

loglik_small(parameters, uq, uq>0, X, W, nbeta, nalpha, 0, 0, binomial())
grad_small(parameters, uq, uq>0, X, W, nbeta, nalpha, 0, 0, binomial())

## optimize the likelihood with optim (~20 mins)
system.time(fit_real <- optim(fn = loglik_small, gr = grad_small, par = parameters, Y=uq, Y1=uq>0, X=X, W=W, 
                              kx=nbeta, kw=nalpha, offsetx=0, offsetw=0, linkobj=binomial(),
                              hessian = FALSE, method = "BFGS", control=list(fnscale=-1)))

betahat_real <- matrix(fit_real$par[1:nbeta], ncol=nrow(uq), nrow=ncol(X))
head(t(betahat_real))

plot(density(betahat_real[1,]))
plot(density(betahat_real[2,]))

muhat_real <- t(exp(X %*% betahat_real))

alphahat_real <- matrix(fit_real$par[(nbeta + 1):(nbeta + nalpha)], nrow=ncol(W), ncol=ncol(uq), byrow=TRUE)
etahat_real <- W %*% alphahat_real
pihat_real <- logistic(etahat_real)

plot(log(colMeans(uq)), colMeans(pihat_real))
plot(log(colMeans(uq)), log(colMeans(muhat_real)))

thetahat <- exp(fit_real$par[nbeta + nalpha + 1])
1/thetahat
