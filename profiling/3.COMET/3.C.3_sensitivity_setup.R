### Jinliang Yang
### 10-11-2016
### purpose: run BSmoothing



library("farmeR")
run_Rcodes(inputdf=data.frame(file=1:9, out=1), outdir="slurm-script", cmdno=1,
           rcodes = "profiling/3.DMR/3.C.2_comet2blocks_sensitivity.R",
           arrayshid = "slurm-script/run_rcode_array.sh",
           email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", 4, "30G"))

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemm --mem 30G --ntasks=4 --time 24:00:00 slurm-script/run_rcode_array.sh
