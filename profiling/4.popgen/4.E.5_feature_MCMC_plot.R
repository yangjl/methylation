### Jinliang Yang
### 10-14-2016
### run MCMC for comet on chr10

source("lib/mplots.R")
source("lib/mcmcbc.R")


ob <- load("largedata/sfs_cg_comet_0.33_gene.RData")
sfsplot(res,burnin=0.2,rates=c(1E8,1E8,1E8), sfsplot="plotmean", Ne=150000, k=0:40)
mplot(res, burnin=0.2, rates=c(1E8,1E8,1E8))
# posterior mu [ 1.09605214547552e-06 ], nu [ 4.26242752079954e-08 ] and s [ 1.63352605892327e-05 ]

ob <- load("largedata/sfs_cg_comet_0.33_exon.RData")
sfsplot(res,burnin=0.2,rates=c(1E8,1E8,1E5), sfsplot="plotmean", Ne=150000, k=0:40)
mplot(res, burnin=0.2, rates=c(1E8,1E8,1E5))
# posterior mu [ 4.64966149155177e-06 ], nu [ 2.30681386298064e-05 ] and s [ 2.42128603993789e-10 ]

ob <- load("largedata/sfs_cg_comet_0.33_intron.RData")

mplot(res, burnin=0.2, rates=c(1E8,1E8, 1E5))
sfsplot(res,burnin=0.2,rates=c(1E8,1E8,1E8), sfsplot="plotmean", Ne=150000, k=0:40)



