### Jinliang Yang
### Jan 10th, 2017
### run COMET overlap with features and then MCMC

## COMET length
ln <- 0
## context
type <- c("CG", "CHG", "CHH")
## gene means 1kb upstream
fs <- c("exon", "intron", "gene", "up1k")

## MUST: gset
gset <- c("above", "below")

a <- data.frame(length=rep(ln, each=24), type=rep(type, each=8), fs=rep(fs, each=2), gset=rep(gset, each=1))

a <- subset(a, type != "CHH")
a$id <- paste(a$type, a$fs, a$length, a$gset, sep="_")
a$TE <- "no"

### RPKM mean
write.table(a, "largedata/rpkm_mean_type.csv", sep=",", row.names=FALSE, quote=FALSE)
### RPKM variance
write.table(a, "largedata/rpkm_var_type.csv", sep=",", row.names=FALSE, quote=FALSE)


library("farmeR")
run_Rcodes(inputdf=data.frame(file=1:16, out=1:16), outdir="slurm-script", cmdno=1,
           rcodes = "profiling/5.MCMC_scan/5.C.1_rpkm_MCMC.R",
           arrayshid = "slurm-script/run_mcmc_16.sh",
           email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", 1, "8G"))

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemm --mem 8G --ntasks=1 --exclude=bigmem1,bigmem2,bigmem6 --time 24:00:00 slurm-script/run_mcmc_16.sh



library("farmeR")
run_Rcodes(inputdf=data.frame(file=1:16, out=1:16), outdir="slurm-script", cmdno=1,
           rcodes = "profiling/5.MCMC_scan/5.C.1_rpkm_MCMC_var.R",
           arrayshid = "slurm-script/run_mcmc_var_16.sh",
           email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", 1, "8G"))

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemm --mem 8G --ntasks=1 --exclude=bigmem1,bigmem2,bigmem6 --time 24:00:00 slurm-script/run_mcmc_var_16.sh
