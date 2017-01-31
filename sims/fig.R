require(dplyr)
require(reshape)
require(ggplot2)
setwd('~/Documents/BRAIN/gitrepo/zinb_analysis/sims/datasets/data/data')

corr = lapply(c(100, 1000), function(nc){
  lapply(c(1, .85), function(aa){
    lapply(c(-3.5, 0), function(offs){
      pp = sprintf('simAllen_%s_a%s_offs%s_seed9128', nc, aa, offs)
      load(paste0(pp, '_dist.rda'))
      cc <- lapply(1:13, function(i) rowMeans(sapply(res, function(x) x[[i]])))
      do.call(cbind, cc)
    })
  })
})
corr = do.call(rbind, do.call(rbind, do.call(rbind, corr)))
corr = data.frame(corr, stringsAsFactors = F)
colnames(corr) <- c("true W", paste0("zinb k=", 1:4), "PCA (k=2)",
                    "PCA TC (k=2)", "PCA TMM (k=2)", "PCA FQ (k=2)",
                    'ZIFA (k=2)', 'ZIFA TC (k=2)', 'ZIFA TMM (k=2)',
                    'ZIFA FQ (k=2)')

corr$pzero = c(rep(rep(c(.25,.45), each = 100), 2),
               rep(rep(c(.25,.45), each = 1000), 2))
corr$var = c(rep(paste('Clustering', 1:2), each = 200),
             rep(paste('Clustering', 1:2), each = 2000))
corr$nc = c(rep('100', 400), rep('1000', 4000))
corrMolten = melt(corr, id.vars = c('pzero', 'var', 'nc'))
corrSum = corrMolten %>% group_by(pzero, var, variable, nc) %>%
  summarize(mean = mean(value), sd = sd(value)) %>% ungroup() %>%
  as.data.frame()
corrSum$pzero = factor(corrSum$pzero)

c1 = corrSum[corrSum$variable != 'true W',]
c1$method = sapply(strsplit(as.vector(c1$variable), ' '), '[[', 1)
ggplot(c1, aes(x = pzero, y = mean, col = variable, group = variable)) +
  geom_point() + geom_line() + labs(col='') + 
  theme_bw() + xlab('Zero Fraction') + facet_grid(nc ~ var) +
  ylab('Correlation between true and estimated sample distances in W')