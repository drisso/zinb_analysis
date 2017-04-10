library(zinbwave)
library(EDASeq)

COUNT_GEN <- readRDS("fig6e-g/function.rds")
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
    fileName = sprintf('fig6e-g/simLun_%s_ziadd%s.rda', ncells, x)
    save(simData, labels, file = fileName)
  }
}

# zero fraction
cc = simData[[1]]$counts
sum(cc == 0)/(ncol(cc)*nrow(cc))

# zero fraction
lapply(c(0,.33,.67), function(x){
  load(sprintf('fig6e-g/simLun_10000_ziadd%s.rda', x))
  sapply(1:10, function(i){
    cc = simData[[i]]$counts
    sum(cc == 0)/(ncol(cc)*nrow(cc))
  })
})





