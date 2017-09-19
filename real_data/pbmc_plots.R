run_analysis <- FALSE
if(run_analysis) {
  set.seed(8283)

  library(cellrangerRkit)
  pipestance_path <- "tenx_pbmc4k"
  if(!file.exists(paste0(pipestance_path, "/outs"))) {
    download_sample(sample_name="pbmc4k",sample_dir=pipestance_path,
                    host="http://cf.10xgenomics.com/samples/cell-exp/1.3.0/")
  }
  pbmc <- load_cellranger_matrix(pipestance_path)

  use_genes <- get_nonzero_genes(pbmc)
  raw <- as.matrix(exprs(pbmc[use_genes, ]))

  vars <- rowVars(log1p(raw))
  names(vars) <- rownames(raw)
  vars <- sort(vars, decreasing = TRUE)
  vargenes <- names(vars)[1:1000]

  dense <- raw[vargenes,]

  #ZINB
  library(BiocParallel)
  library(doParallel)
  registerDoParallel(6)
  register(DoparParam())

  library(zinbwave)
  system.time(zinb_res <- zinbFit(dense, K=2, epsilon=1e2))

  # PCA
  totalcount = function (ei)
  {
    sums = colSums(ei)
    eo = t(t(ei)*mean(sums)/sums)
    return(eo)
  }

  tc <- totalcount(raw)
  fq <- FQT_FN(raw)
  tmm <- TMM_FN(raw)
  deseq <- DESEQ_FN(raw)

  library(rARPACK)
  fastpca <- function(expr, scale=FALSE) {
    svd_raw <- svds(scale(t(expr), center=TRUE, scale=scale), k=3, nu=3, nv=0)
    pc_raw <- svd_raw$u %*% diag(svd_raw$d[1:3])
    return(pc_raw)
  }

  pc_raw <- fastpca(log1p(raw[vargenes,]))
  pc_tc <- fastpca(log1p(tc[vargenes,]))
  pc_fq <- fastpca(log1p(fq[vargenes,]))
  pc_tmm <- fastpca(log1p(tmm[vargenes,]))
  pc_deseq <- fastpca(log1p(deseq[vargenes,]))

  # ZIFA
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
      cmd = sprintf('python run_zifa.py %s%s.csv %s_zifa.csv', bb, tmp, tmp)
      system(cmd)
    }
    read.csv(sprintf("%s_zifa.csv", tmp), header=FALSE)
  }

  zifa_raw <- wrapRzifa(log1p(raw[vargenes,]))
  zifa_tc <- wrapRzifa(log1p(tc[vargenes,]))
  zifa_tmm <- wrapRzifa(log1p(tmm[vargenes,]))
  zifa_fq <- wrapRzifa(log1p(fq[vargenes,]))

  ## compute QC metrics
  library(scater)
  sceset <- newSCESet(countData = exprs(pbmc[use_genes,]))

  sceset <- calculateQCMetrics(sceset)
  qc <- pData(sceset)[,c(2, 5, 7, 12:15)]
  head(qc)

  save(zinb_res, pc_raw, pc_tc, pc_fq, pc_tmm, pc_deseq,
       zifa_raw, zifa_tc, zifa_fq, zifa_tmm, qc, file="pbmc_4k_res.rda")
}
load("pbmc_4k_res.rda")
load("tenx_pbmc4k/zinb_eps.rda")

clusters <- read.csv("tenx_pbmc4k/clusters.csv")
cl <- as.factor(clusters[,2])
names(cl) <- clusters[,1]
table(cl)

library(cowplot)
library(magrittr)
library(RColorBrewer)
col1 <- brewer.pal(9, "Set1")
col2 <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Set3"), brewer.pal(8, "Set1"))

data.frame(Dim1=pc_fq[,1], Dim2=pc_fq[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=cl)) + geom_point() +
  scale_color_manual(values=col2) -> panel1_pca
data.frame(Dim1=zifa_fq[,1], Dim2=zifa_fq[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=cl)) + geom_point() +
  scale_color_manual(values=col2) -> panel1_zifa
data.frame(Dim1=zinb_res3@W[,1], Dim2=zinb_res3@W[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=cl)) + geom_point() +
  scale_color_manual(values=col2) -> panel1_zinb

p1 <- plot_grid(panel1_pca + theme(legend.position = "none"),
                panel1_zifa + theme(legend.position = "none"),
                panel1_zinb + theme(legend.position = "none"),
                labels=c("a", "c", "e"), align = "h", ncol=3)

legend <- get_legend(panel1_pca)
upper <- plot_grid(p1, legend, rel_widths = c(3, .6))


cors <- lapply(1:2, function(i) abs(cor(pc_fq[,i], qc)))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors, QC=rep(stringr::str_to_lower(colnames(qc)), 2), Dimension=as.factor(rep(1:2, each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) + geom_bar(stat="identity", position='dodge') + scale_fill_manual(values=col2) + ylim(0, 1) -> panel2_pca

cors <- lapply(1:2, function(i) abs(cor(zifa_fq[,i], qc)))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors, QC=rep(stringr::str_to_lower(colnames(qc)), 2), Dimension=as.factor(rep(1:2, each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) + geom_bar(stat="identity", position='dodge') + scale_fill_manual(values=col2) + ylim(0, 1) -> panel2_zifa

cors <- lapply(1:2, function(i) abs(cor(zinb_res3@W[,i], qc)))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors, QC=rep(stringr::str_to_lower(colnames(qc)), 2), Dimension=as.factor(rep(1:2, each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) + geom_bar(stat="identity", position='dodge') + scale_fill_manual(values=col2) + ylim(0, 1) -> panel2_zinb

p2 <- plot_grid(panel2_pca + theme(legend.position = "none"),
                panel2_zifa + theme(legend.position = "none"),
                panel2_zinb + theme(legend.position = "none"),
                labels=c("b", "d", "f"), align = "h", ncol=3)

legend2 <- get_legend(panel2_pca)
lower <- plot_grid(p2, legend2, rel_widths = c(3, 1))

fig1 <- plot_grid(upper, lower, ncol=1, nrow=2)
fig1

save_plot("pbmc_fig1.pdf", fig1,
          ncol = 3,
          nrow = 3,
          base_aspect_ratio = 1.3
)

library(cluster)
mycol <- c(brewer.pal(11,"RdYlGn")[c(8:11, 1:4)], brewer.pal(11,"RdYlBu")[8:11])
colsil <- mycol[c(5:12, 2, 3)]

methods <- list(pc_raw[,1:2], pc_tc[,1:2], pc_tmm[,1:2], pc_fq[,1:2],
                zifa_raw, zifa_tc, zifa_tmm, zifa_fq,
                zinb_res3@W[,1:2])
names(methods) <- c(paste0("PCA_", c("RAW", "TC", "TMM", "FQ")),
                    paste0("ZIFA_", c("RAW", "TC", "TMM", "FQ")),
                    "ZINB-WaVE")
met_type <- as.factor(c(rep(c("PCA", "ZIFA"), each=4), "ZINB-WaVE"))

sil_lay <- sapply(seq_along(methods), function(i) {
  d <- dist(methods[[i]])
  ss <- silhouette(as.numeric(cl), d)
  mean(ss[,3])
})

bars <- data.frame(AverageSilhouette=sil_lay,
                   Method=factor(names(methods),
                                 levels=c("PCA_RAW", "PCA_TC", "PCA_TMM", "PCA_FQ",
                                          "ZIFA_RAW", "ZIFA_TC", "ZIFA_TMM", "ZIFA_FQ",
                                          "ZINB-WaVE")),
                   Type=met_type)

bars %>%
  ggplot(aes(Method, AverageSilhouette, group=Type, fill=Method)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=colsil) + coord_flip() +
  theme(legend.position = "none") -> sil_pbmc

sil_pbmc
