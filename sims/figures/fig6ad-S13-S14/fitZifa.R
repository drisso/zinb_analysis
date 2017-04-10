library(cluster)
library(parallel)
library(zinbwave)
library(EDASeq)
library(digest)
library(edgeR)

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
    cmd = sprintf('python ../../real_data/run_zifa.py %s%s.csv %s_zifa.csv', bb, tmp, tmp)
    system(cmd)
  }
  read.csv(sprintf("%s_zifa.csv", tmp), header=FALSE)
}


lun = T
if (lun){
  for (nc in c(100, 1000, 10000)){
    for (ziadd in c('_ziadd0', '_ziadd0.33', '_ziadd0.67')){
      pref = sprintf('fig6e-g/simLun_%s%s', nc, ziadd)
      load(paste0(pref, '.rda'))
      
      print('zifa raw')
      zifa = lapply(1:length(simData), function(i){
        print(i)
        if (ziadd == ''){
          Y = simData[[i]] 
        }else{
          Y = simData[[i]]$counts
        }
        Y = Y[rowSums(Y) != 0, ]
        Y = log1p(Y)
        zz <- tryCatch({z = wrapRzifa(Y)},
                       error = function(err) return(NULL), 
                       finally = function() return(z))
        zz
      })
      save(zifa, file = paste0(pref, '_zifa.rda'))
      
      print('zifa tc')
      zifaTC = lapply(1:length(simData), function(i){
        print(i)
        if (ziadd == ''){
          Y = simData[[i]] 
        }else{
          Y = simData[[i]]$counts
        }
        Y = Y[rowSums(Y) != 0, ]
        mult = sum(Y) / (ncol(Y) * nrow(Y))
        fact = colSums(Y)
        Y = mult * (t(Y) / fact)
        Y = log1p(t(Y))
        zz <- tryCatch({z = wrapRzifa(Y)},
                       error = function(err) return(NULL), 
                       finally = function() return(z))
        zz
      })
      save(zifaTC, file = paste0(pref, '_zifaTC.rda'))
      
      print('zifa fq')
      zifaFQ = lapply(1:length(simData), function(i){
        print(i)
        if (ziadd == ''){
          Y = simData[[i]] 
        }else{
          Y = simData[[i]]$counts
        }
        Y = Y[rowSums(Y) != 0, ]
        fq <- betweenLaneNormalization(Y, which="full")
        Y = log1p(fq)
        zz <- tryCatch({z = wrapRzifa(Y)},
                       error = function(err) return(NULL), 
                       finally = function() return(z))
        zz
      })
      save(zifaFQ, file = paste0(pref, '_zifaFQ.rda'))
      
      print('zifa tmm')
      zifaTMM = lapply(1:length(simData), function(i){
        print(i)
        if (ziadd == ''){
          Y = simData[[i]] 
        }else{
          Y = simData[[i]]$counts
        }
        counts = Y[rowSums(Y) != 0, ]
        y = DGEList(counts)
        y = calcNormFactors(y, method="TMM")
        tmm = t(t(counts) / (y$samples$lib.size * y$samples$norm.factors))
        Y = log1p(tmm)
        zz <- tryCatch({z = wrapRzifa(Y)},
                       error = function(err) return(NULL), 
                       finally = function() return(z))
        zz
      })
      save(zifaTMM, file = paste0(pref, '_zifaTMM.rda'))
    }
  }
}


allen = F
if (allen){
  for (nc in c(100, 1000)){
    for (b2 in c(1, 5, 50)){
      for (offs in c(0, 2, 5)){
        pref = sprintf('fig6ad-S13-S14/simAllen_nc%s_ratio%s_offs%s', nc, b2, offs)
        load(paste0(pref, '.rda'))
        
        print('zifa raw')
        zifa = lapply(1:length(simData), function(i){
          print(i)
          Y = t(simData[[i]]$counts)
          Y = Y[rowSums(Y) != 0, ]
          Y = log1p(Y)
          zz <- tryCatch({z = wrapRzifa(Y)},
                         error = function(err) return(NULL), 
                         finally = function() return(z))
          zz
        })
        save(zifa, file = paste0(pref, '_zifa.rda'))
        
        print('zifa tc')
        zifaTC = lapply(1:length(simData), function(i){
          print(i)
          Y = t(simData[[i]]$counts)
          Y = Y[rowSums(Y) != 0, ]
          mult = sum(Y) / (ncol(Y) * nrow(Y))
          fact = colSums(Y)
          Y = mult * (t(Y) / fact)
          Y = log1p(t(Y))
          zz <- tryCatch({z = wrapRzifa(Y)},
                         error = function(err) return(NULL), 
                         finally = function() return(z))
          zz
        })
        save(zifaTC, file = paste0(pref, '_zifaTC.rda'))
        
        print('zifa fq')
        zifaFQ = lapply(1:length(simData), function(i){
          print(i)
          Y = t(simData[[i]]$counts)
          Y = Y[rowSums(Y) != 0, ]
          fq <- betweenLaneNormalization(Y, which="full")
          Y = log1p(fq)
          zz <- tryCatch({z = wrapRzifa(Y)},
                         error = function(err) return(NULL), 
                         finally = function() return(z))
          zz
        })
        save(zifaFQ, file = paste0(pref, '_zifaFQ.rda'))
        
        print('zifa tmm')
        zifaTMM = lapply(1:length(simData), function(i){
          print(i)
          counts = t(simData[[i]]$counts)
          counts = counts[rowSums(counts) != 0, ]
          y = DGEList(counts)
          y = calcNormFactors(y, method="TMM")
          tmm = t(t(counts) / (y$samples$lib.size * y$samples$norm.factors))
          Y = log1p(tmm)
          zz <- tryCatch({z = wrapRzifa(Y)},
                         error = function(err) return(NULL), 
                         finally = function() return(z))
          zz
        })
        save(zifaTMM, file = paste0(pref, '_zifaTMM.rda'))
      }
    }
  }
}
