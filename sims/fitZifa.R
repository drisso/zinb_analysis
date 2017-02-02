library(cluster)
library(parallel)
library(zinb)
library(EDASeq)
library(digest)
library(edgeR)
library(DESeq2)

wrapRzifa <- function(Y, block = T){
  # wrapper R function for ZIFA.
  # md5 hashing and temporary files are used not to re-run zifa 
  # if it has already be run on this computer.
  d = digest(Y, "md5")
  tmp = paste0(tempdir(), '/', d)
  write.csv(Y, paste0(tmp, '.csv'))
  
  if (!file.exists(paste0(tmp, '_zifa.csv'))){
    print('run ZIFA')
    bb = ifelse(block, '-b ', '')
    cmd = sprintf('python real_data/run_zifa.py %s%s.csv %s_zifa.csv', bb, tmp, tmp)
    system(cmd)
  }
  read.csv(sprintf("%s_zifa.csv", tmp), header=FALSE)
}

for (nc in c(1000)){
  for (offs in c(3.5)){
    mclapply(c(1, .85), function(aa){
      pref = sprintf('sims/datasets/simAllen_%s_a%s_offs%s_seed9128', nc, aa, offs)
      print(nc)
      print(aa)
      print(offs)
      load(paste0(pref, '.rda'))
      
      print('zifa raw')
      zifa = lapply(1:10, function(i){
        print(i)
        Y = t(simData[[i]]$counts)
        Y = Y[rowSums(Y) != 0, ]
        Y = log1p(Y)
        wrapRzifa(Y)
      })
      save(zifa, file = paste0(pref, '_zifa.rda'))
      
      print('zifa tc')
      zifaTC = lapply(1:10, function(i){
        print(i)
        Y = t(simData[[i]]$counts)
        Y = Y[rowSums(Y) != 0, ]
        mult = sum(Y) / (ncol(Y) * nrow(Y))
        fact = colSums(Y)
        Y = mult * (t(Y) / fact)
        Y = log1p(t(Y))
        wrapRzifa(Y)
      })
      save(zifaTC, file = paste0(pref, '_zifaTC.rda'))
      
      print('zifa fq')
      zifaFQ = lapply(1:10, function(i){
        print(i)
        Y = t(simData[[i]]$counts)
        Y = Y[rowSums(Y) != 0, ]
        fq <- betweenLaneNormalization(Y, which="full")
        Y = log1p(fq)
        wrapRzifa(Y)
      })
      save(zifaFQ, file = paste0(pref, '_zifaFQ.rda'))
      
      print('zifa tmm')
      zifaTMM = lapply(1:10, function(i){
        print(i)
        counts = t(simData[[i]]$counts)
        counts = counts[rowSums(counts) != 0, ]
        y = DGEList(counts)
        y = calcNormFactors(y, method="TMM")
        tmm = t(t(counts) / (y$samples$lib.size * y$samples$norm.factors))
        Y = log1p(tmm)
        wrapRzifa(Y)
      })
      save(zifaTMM, file = paste0(pref, '_zifaTMM.rda'))
      
      print('zifa deseq2')
      zifaDeseq2 = lapply(1:10, function(i){
        print(i)
        counts = t(simData[[i]]$counts)
        counts = counts[rowSums(counts) != 0, ]
        condition = factor(rep(1, ncol(counts)))
        ## remove very high counts, otherwise DESeq2 creates an error
        ## only one or two very high count per dataset in Allen simulations
        counts[which(counts > 1000000000)] = max(counts[which(counts < 1000000000)])
        dds = DESeqDataSetFromMatrix(counts, DataFrame(condition), ~ 1)
        dds = estimateSizeFactors(dds)
        Y = log1p(counts(dds, normalized=TRUE))
        wrapRzifa(Y)
      })
      save(zifaDeseq2, file = paste0(pref, '_zifaDeseq2.rda'))
      
    }, mc.cores = 2)
  }
}



