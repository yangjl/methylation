### Jinliang Yang
### 10-14-2016
### run MCMC for comet on chr10

source("lib/mplots.R")
source("lib/mcmcbc.R")


ob <- load("largedata/sfs_cg_comet_0.33_gene.RData")
sfsplot(res,burnin=0.2,rates=c(1E8,1E8,1E8), sfsplot="plotmean", Ne=150000, k=0:40)
mplot(res, burnin=0.2, rates=c(1E8,1E8,1E8))
# posterior mu [ 4.60934455109854e-07 ], nu [ 6.77645732778247e-07 ] and s [ 1.79033264284426e-06 ]

ob <- load("largedata/sfs_cg_comet_0.33_exon.RData")
sfsplot(res,burnin=0.2,rates=c(1E8,1E8,1E8), sfsplot="plotmean", Ne=150000, k=0:40)
mplot(res, burnin=0.2, rates=c(1E8,1E8,1E5))
# posterior mu [ 6.62599293755031e-07 ], nu [ 1.03146663162443e-07 ] and s [ 7.32085981794232e-06 ]

ob <- load("largedata/sfs_cg_comet_0.33_intron.RData")

mplot(res, burnin=0.2, rates=c(1E8,1E8, 1E5))
sfsplot(res,burnin=0.2,rates=c(1E8,1E8,1E8), sfsplot="plotmean", Ne=150000, k=0:40)

###########
ob <- load("largedata/ngbm_sfs_cg_comet_0.33_exon.RData")
mplot(res, burnin=0.2, rates=c(1E8,1E8, 1E5))
sfsplot(res,burnin=0.2,rates=c(1E8,1E8,1E5), sfsplot="plotmean", Ne=150000, k=0:40)
sfsplot(res,burnin=0.2,rates=c(1E8,1E8,1E5), sfsplot="plotmode", Ne=150000, k=0:40)
#posterior mu [ 5.19604993078061e-07 ], nu [ 1.23984784487354e-07 ] and s [ 4.84685320372646e-06 ]

ob <- load("largedata/ngbm_sfs_cg_comet_0.33_gene.RData")
mplot(res, burnin=0.2, rates=c(1E8,1E8, 1E5))
sfsplot(res,burnin=0.2,rates=c(1E8,1E8,1E5), sfsplot="plotmean", Ne=150000, k=0:40)
sfsplot(res,burnin=0.2,rates=c(1E8,1E8,1E5), sfsplot="plotmode", Ne=150000, k=0:40)
#posterior mu [ 7.65364773455305e-07 ], nu [ 1.15642243453368e-07 ] and s [ 7.40731878049114e-06 ]

ob <- load("largedata/ngbm_sfs_cg_comet_0.33_intron.RData")
mplot(res, burnin=0.2, rates=c(1E8,1E8, 1E5))
sfsplot(res,burnin=0.2,rates=c(1E8,1E8,1E5), sfsplot="plotmean", Ne=150000, k=0:40)
sfsplot(res,burnin=0.2,rates=c(1E8,1E8,1E5), sfsplot="plotmode", Ne=150000, k=0:40)
#posterior mu [ 7.70753203030945e-07 ], nu [ 1.10005583819213e-07 ] and s [ 8.22688696216015e-06 ]

################ gbM
ob <- load("largedata/gbm_sfs_cg_comet_0.33_exon.RData")
mplot(res, burnin=0.2, rates=c(1E8,1E8, 1E5))
sfsplot(res,burnin=0.2,rates=c(1E8,1E8,1E5), sfsplot="plotmean", Ne=150000, k=0:40)
sfsplot(res,burnin=0.2,rates=c(1E8,1E8,1E5), sfsplot="plotmode", Ne=150000, k=0:40)
#posterior mu [ 7.77170188577186e-07 ], nu [ 5.42423118338264e-08 ] and s [ 1.41743515403035e-05 ]

ob <- load("largedata/gbm_sfs_cg_comet_0.33_gene.RData")
mplot(res, burnin=0.2, rates=c(1E8,1E8, 1E5))
sfsplot(res,burnin=0.2,rates=c(1E8,1E8,1E5), sfsplot="plotmean", Ne=150000, k=0:40)
sfsplot(res,burnin=0.2,rates=c(1E8,1E8,1E5), sfsplot="plotmode", Ne=150000, k=0:40)
#posterior mu [ 9.03339445303956e-07 ], nu [ 5.72949539919047e-08 ] and s [ 1.52062951256493e-05 ]

ob <- load("largedata/gbm_sfs_cg_comet_0.33_intron.RData")
mplot(res, burnin=0.2, rates=c(1E8,1E8, 1E5))
sfsplot(res,burnin=0.2,rates=c(1E8,1E8,1E5), sfsplot="plotmean", Ne=150000, k=0:40)
#sfsplot(res,burnin=0.2,rates=c(1E8,1E8,1E5), sfsplot="plotmode", Ne=150000, k=0:40)
#posterior mu [ 7.53712629946411e-07 ], nu [ 4.37777207495372e-08 ] and s [ 1.45222758073796e-05 ]



