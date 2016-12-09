### Jinliang Yang
### 10-14-2016
### run MCMC for comet on chr10

##get command line args
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

JOBID <- as.numeric(as.character(args[1]))
print(JOBID)

source("lib/mplots.R")
source("lib/mcmcbc.R")

sfs <- read.csv("cache/gbm_sfs_genomic.csv")

b <- read.csv("cache/SFS_comet_blocks_CG.csv")
hist(log10(b$length), xlab="COMET Length log10(bp)", main="CG COMETs", col="#cdc0b0")
abline(v=quantile(log10(b$length))[2:4], lwd=3, lty=2 )
# 0%   25%   50%   75%  100% 
# 1    96   309   784 76180 

b1 <- subset(b, length <= 96)
b1 <- subset(b, length > 700)

out <- table(b1$sfs)
sfs <- data.frame(out)
plot(sfs, type="p")

b2 <- subset(b, length > 96 & length <= 309)
out <- table(b2$sfs)
sfs2 <- data.frame(out)
plot(sfs2, type="p")
# If acceptance too high, increase these values to explore wider space. If acceptance too low, decrease.
res2 <- MCMCBC(my_sfs=sfs2$Freq, rates=c(1E8,1E8,1E8), sd=c(0.05,0.05,0.05), k=0:40,
               conditional=FALSE, Ne=150000, ngen=100000, verbose=TRUE)

b3 <- subset(b, length > 309 & length <= 784)
out <- table(b3$sfs)
sfs3 <- data.frame(out)
plot(sfs3, type="p")
# If acceptance too high, increase these values to explore wider space. If acceptance too low, decrease.
res3 <- MCMCBC(my_sfs=sfs3$Freq, rates=c(1E8,1E8,1E8), sd=c(0.05,0.05,0.05), k=0:40,
               conditional=FALSE, Ne=150000, ngen=100000, verbose=TRUE)


b4 <- subset(b, length > 784)
out <- table(b4$sfs)
sfs4 <- data.frame(out)
plot(sfs4, type="p")


# If acceptance too high, increase these values to explore wider space. If acceptance too low, decrease.
res4 <- MCMCBC(my_sfs=sfs4$Freq, rates=c(1E8,1E8,1E8), sd=c(0.05,0.05,0.05), k=0:40,
              conditional=FALSE, Ne=150000, ngen=100000, verbose=TRUE)

#####
out <- gsub("cache", "largedata", files[JOBID])
out <- gsub("csv", "RData", out)
save(list="res", file=out)
### plot trace and posteriors
sfsplot(res4,burnin=0.2,rates=c(1E8,1E8,1E8), sfsplot="plotmean", Ne=150000, k=0:40)
mplot(res, burnin=0.1, rates=c(1E8,1E8,1E8))


