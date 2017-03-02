library(zinb)
library(EDASeq)

COUNT_GEN <- readRDS("sims/figures/lunSims/function.rds")
ngenes <- 1000
B <- 10

for (ncells in c(100, 1000, 10000)){
  for (x in c(0, .33, .67)){
    n1 = n2 = round(ncells/3)
    n3 = ncells - n1 * 2
    labels = c(rep(1, n1), rep(2, n2), rep(3, n3))
    simData = lapply(seq_len(B), function(i){
      seed = i
      set.seed(seed)
      COUNT_GEN(labels, ngenes, zinb=TRUE, zi.add = x)
    })
    fileName = sprintf('simLun_%s_ziadd%s.rda', ncells, x)
    save(simData, labels, file = fileName)
  }
}

cc = simData[[1]]$counts
sum(cc == 0)/(ncol(cc)*nrow(cc))

# lib size
barplot(colSums(cc), col = labels)

# pca versus zinb
ccNo = cc[rowSums(cc) != 0, ]
fq = betweenLaneNormalization(ccNo, which="full")
pca = prcomp(t(log1p(fq)))
zinb = zinbFit(ccNo, K=2, ncores = 2, commondispersion = F, epsilon = nrow(cc))
par(mfrow = c(1, 2))
plot(pca$x[,1:2], col = labels)
plot(zinb@W, col = labels)

# zero fraction
lapply(c(0,.33,.67), function(x){
  load(sprintf('simLun_10000_ziadd%s.rda', x))
  sapply(1:10, function(i){
    cc = simData[[i]]$counts
    sum(cc == 0)/(ncol(cc)*nrow(cc))
  })
})





