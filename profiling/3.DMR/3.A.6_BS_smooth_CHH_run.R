### Jinliang Yang
### 10-11-2016
### purpose: run BSmoothing



library("farmeR")
run_Rcodes(inputdf=data.frame(file=1:20, out=1), outdir="slurm-script", cmdno=1,
           rcodes = "profiling/3.DMR/3.A.5_BS_smooth_CHH.R",
           arrayshid = "slurm-script/run_rcode_array_chh.sh",
           email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", 8, "60G"))

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemm --mem 60G --ntasks=8 --time 1000:00:00 slurm-script/run_rcode_array_chh.sh
