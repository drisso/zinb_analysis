## simulate data with Lun model

library(zinbwave)

COUNT_GEN <- readRDS("../figures/lunSims/function.rds")
ngenes <- 1000
B <- 10
x <- .33
nc <- c(50, 100, 500, 1000, 5000, 10000)

for (ncells in nc){
  n1 = n2 = round(ncells/3)
  n3 = ncells - n1 * 2
  labels = c(rep(1, n1), rep(2, n2), rep(3, n3))
  simData = lapply(seq_len(B), function(i){
    seed = i
    set.seed(seed)
    COUNT_GEN(labels, ngenes, zinb=TRUE, zi.add = x)
  })

  fileName = sprintf('simLun_%s.rda', ncells)
  save(simData, labels, file = fileName)
}

cpuTime = lapply(nc, function(ncells){
  fileName = sprintf('simLun_%s.rda', ncells)
  load(fileName)
  tt = lapply(1:B, function(j){
    counts = simData[[j]]$counts
    counts = counts[rowSums(counts) > 5, colSums(counts) > 5 ]
    system.time(zinbFit(counts, K = 2, commondispersion = TRUE,
                        epsilon = 1000, ncores = 7))
  })
  tt
})
save(cpuTime, file = 'cpuTime.rda')

###########
# plot
###########
library(ggplot2)
library(cowplot)
library(RColorBrewer)
mycol = c(brewer.pal(11,"RdYlGn")[c(8:11, 1:4)], brewer.pal(11,"RdYlBu")[8:11])
load('cpuTime.rda')
elapse = lapply(1:length(cpuTime), function(nc){
  print(nc)
  sapply(1:length(cpuTime[[1]]), function(i){
    print(i)
    cpuTime[[nc]][[i]][[3]]
  })
})
mean = sapply(elapse, median)
sd = sapply(elapse, mad)
nc = c(50, 100, 500, 1000, 5000, 10000)
dfTime = data.frame(mean = mean, sd = sd, nc = nc)
pd <- position_dodge(0.1)
cpu = ggplot(dfTime, aes(x = nc, y = mean, group = 1)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .1, position = pd) +
  geom_line(position=pd) +
  geom_point(position=pd) + background_grid(major = 'xy', minor = 'xy') +
  xlab('Number of cells') + ylab('CPU time (sec.)') +
  coord_trans(y = "log10", x = 'log10') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(breaks = c(25, 50, 100, 250, 500, 1000, 2000)) +
  scale_x_continuous(breaks = nc)
cpu
ggsave(filename="cpu.png", plot = cpu, device = 'png')
