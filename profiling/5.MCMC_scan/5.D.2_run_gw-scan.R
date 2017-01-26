### Jinliang Yang
### Jan 10th, 2017
### run COMET overlap with features and then MCMC

## COMET length
ln <- 0:4
## context
type <- c("CG", "CHG")
## gene means 1kb upstream
fs <- NULL

## MUST: gset
gset <- NULL

a <- data.frame(length=rep(ln, each=2), type=rep(type, each=1))
a$id <- paste(a$type, a$length, sep="_")
a$TE <- "no"

### RPKM mean
write.table(a, "largedata/gs_type.csv", sep=",", row.names=FALSE, quote=FALSE)
### RPKM variance
write.table(a, "largedata/rpkm_var_type.csv", sep=",", row.names=FALSE, quote=FALSE)


library("farmeR")
run_Rcodes(inputdf=data.frame(file=1:10, out=1:10), outdir="slurm-script", cmdno=1,
           rcodes = "profiling/5.MCMC_scan/5.D.1_gw-scan.R",
           arrayshid = "slurm-script/run_mcmc_gs10.sh",
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
