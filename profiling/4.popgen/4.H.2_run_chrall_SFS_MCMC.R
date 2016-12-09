### Jinliang Yang
### 12-08-2016
### purpose: run MCMC for WG CG



library("farmeR")
run_Rcodes(inputdf=data.frame(file=1, out=1), outdir="slurm-script", cmdno=1,
           rcodes = "profiling/4.popgen/4.H.1_chrall_SFS_MCMC.R",
           arrayshid = "slurm-script/run_mcmc_array.sh",
           email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemh", 2, "16G"))

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemm --mem 16G --ntasks=2 --time 48:00:00 slurm-script/run_rcode_array.sh

