### Jinliang Yang
### Jan 10th, 2017
### run COMET overlap with features and then MCMC

## COMET length
ln <- 0:4
## context
type <- c("CG", "CHG", "CHH")
## gene means 1kb upstream
fs <- c("exon", "intron", "gene", "up1k")

## body methylation
gbody <- c("gbm", "ngbm")

a <- data.frame(length=rep(ln, each=24), type=rep(type, each=8), fs=rep(fs, each=2), gbody=rep(gbody, each=1))

a <- subset(a, type != "CHH")
a$id <- paste(a$type, a$fs, a$length, a$gbody, sep="_")
a$TE <- "no"
write.table(a, "largedata/gbody_type.csv", sep=",", row.names=FALSE, quote=FALSE)

library("farmeR")
run_Rcodes(inputdf=data.frame(file=1:80, out=1:80), outdir="slurm-script", cmdno=1,
           rcodes = "profiling/5.MCMC_scan/5.B.1_gbm_overlap_MCMC.R",
           arrayshid = "slurm-script/run_mcmc_80.sh",
           email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", 1, "8G"))

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemm --mem 8G --ntasks=1 --exclude=bigmem1,bigmem2,bigmem6 --time 24:00:00 slurm-script/run_mcmc_80.sh

