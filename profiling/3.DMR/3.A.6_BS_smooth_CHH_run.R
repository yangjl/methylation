### Jinliang Yang
### 10-11-2016
### purpose: run BSmoothing



library("farmeR")
run_Rcodes(inputdf=data.frame(file=1:20, out=1), outdir="slurm-script", cmdno=1,
           rcodes = "profiling/3.DMR/3.A.5_BS_smooth_CHH.R",
           arrayshid = "slurm-script/run_rcode_array_chh.sh",
           email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", 16, "120G"))

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemm --mem 130G --ntasks=16 --time 100:00:00 slurm-script/run_rcode_array_chh.sh
