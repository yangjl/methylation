### Jinliang Yang
### 12-08-2016
### purpose: run MCMC for WG CG



library("farmeR")
run_Rcodes(inputdf=data.frame(file=1:10, out=1:10), outdir="slurm-script", cmdno=1,
           rcodes = "profiling/4.popgen/4.H.1_chrall_SFS_MCMC.R",
           arrayshid = "slurm-script/run_mcmc_gc.sh",
           email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", 1, "8G"))

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemm --mem 8G --ntasks=1 --time 24:00:00 --exclude=bigmem1,bigmem6 slurm-script/run_mcmc_gc.sh

