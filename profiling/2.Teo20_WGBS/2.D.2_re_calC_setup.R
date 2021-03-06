### Jinliang Yang
### 10-04-2016
### purpose: re-calibrate the C/T and G/A heterozygote sites



library("farmeR")
run_Rcodes(inputdf=data.frame(file=1:20, out=1), outdir="slurm-script", cmdno=1,
           rcodes = "profiling/2.Teo20_WGBS/2.D.1_re_cal_C_sites.R",
           arrayshid = "slurm-script/run_rcode_array.sh",
           email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", 8, "60G"))

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p jclarke --mem 60G --ntasks=8 --time 100:00:00 slurm-script/run_rcode_array.sh