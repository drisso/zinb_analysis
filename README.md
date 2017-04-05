# Data analysis and simulations for the ZINB-WaVE paper

This repository is designed to allow interested people to reproduce the results and figures of the paper:

Risso D, Perraudeau F, Gribkova S, Dudoit S, Vert JP. ZINB-WaVE: A general and flexible method for signal extraction from single-cell RNA-seq data.

## Dependencies

In order to be able to run the code in this repo, it is required to have `R` (>=3.3) and `python` (>=2.7).

### R packages

- [zinbwave (R package)](https://github.com/drisso/zinbwave)
- cluster
- matrixStats
- magrittr
- RColorBrewer
- ggplot2
- reshape
- dplyr
- knitr
- rmarkdown
- mclust
- cowplot
- rARPACK
- Rtsne
- parallel
- digest


### Bioconductor packages

- EDASeq
- biomaRt
- scRNAseq
- SummarizedExperiment
- edgeR
- scran
- scater
- scone
- DESeq2


### python packages

- [ZIFA](https://github.com/epierson9/ZIFA)

## Getting started

### Real data

For each of the real datasets analyzed in the paper, there are a `.Rmd` file and a `.R` file in the `real_data` folder, e.g.,
for the Patel data, the files are [patel_covariates.Rmd](https://github.com/drisso/zinb_analysis/blob/master/real_data/patel_covariates.Rmd)
and [patel_plots.R](https://github.com/drisso/zinb_analysis/blob/master/real_data/patel_plots.R).

One needs to compile the .Rmd file first. This will have two effects: (i) it will create an HTML report with useful analyses of
the dataset; and (ii) it will create a `.rda` file with the results of `zinbwave`, `pca`, and `zifa`. Once this file is generated,
one can use the `.R` file to generate the dataset-specific plots found in the paper.

To generate the plots related to silhouette width, one needs to source the [silhouette.R](https://github.com/drisso/zinb_analysis/blob/master/real_data/silhouette.R)
file.

To generate the plots related to the goodness-of-fit, run the `.Rmd` files in the `real_data` folder starting with `goodness_of_fit`, e.g., for the Patel data, the file is [goodness_of_fit_patel.Rmd].

### Simulations

To create the simulated datasets from the real datasets used in the paper, first run the code in [simFunction.R](https://github.com/drisso/zinb_analysis/blob/master/sims/figures/simFunction.R). Then, run the `.R` files in the folders in `sims/figures`. Finally, run [figuresPaper.Rmd](https://github.com/drisso/zinb_analysis/blob/master/sims/figures/figuresPaper.Rmd).

To simulate the datasets from the Lun & Marioni model, use the steps described in the methods section of the paper.

To fit the simulated datasets with 10,000 cells, run [Makefile](https://github.com/drisso/zinb_analysis/blob/master/sims/figures/fig6ad-S13-S14/Makefile). It launches a job on a servor [fitZinb10000.R](https://github.com/drisso/zinb_analysis/blob/master/sims/figures/fig6ad-S13-S14/fitZinb10000.R) with the arguments in the Makefile. Alternatively, you can just call [fitZinb10000.R](https://github.com/drisso/zinb_analysis/blob/master/sims/figures/fig6ad-S13-S14/fitZinb10000.R) from your terminal with the arguments you want.


