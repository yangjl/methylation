### Jinliang Yang
### Jan 10th, 2017
### run COMET overlap with features and then MCMC

## COMET length
ln <- 1:4
## context
type <- c("CG", "CHG", "CHH")
## gene means 1kb upstream
fs <- c("exon", "intron", "gene", "up1k")

## body methylation
gbody <- c("gbm", "ngbm")

a <- data.frame(length=rep(ln, each=24), fs=rep(fs, each=8), type=rep(type, each=2), gbody=rep(gbody, each=1))

a <- subset(a, type != "CHH")
a$id <- paste(a$type, a$fs, a$length, sep="_")
a$TE <- "no"
write.table(a, "largedata/gbody_type.csv", sep=",", row.names=FALSE, quote=FALSE)

library("farmeR")
run_Rcodes(inputdf=data.frame(file=1:24, out=1:24), outdir="slurm-script", cmdno=1,
           rcodes = "profiling/4.popgen/4.H.3_overlap_MCMC.R",
           arrayshid = "slurm-script/run_mcmc_64.sh",
           email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", 1, "8G"))

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemm --mem 8G --ntasks=1 --exclude=bigmem1,bigmem6 --time 24:00:00 slurm-script/run_mcmc_24.sh

