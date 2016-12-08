### Jinliang Yang
### 10-12-2016
### purpose: run COMETgazer

dirs <- dir(path="largedata/COMET", pattern="JR", full.names=TRUE, recursive=FALSE)
for(i in 1:length(dirs)){
    file.copy(from="profiling/3.DMR/1.COMETgazer.sh", to=dirs[i], overwrite=TRUE)
}


library("farmeR")
run_Rcodes(inputdf=data.frame(file=1:20, out=1), outdir="slurm-script", cmdno=1,
           rcodes = "profiling/3.DMR/3.B.2_COMMETgazer.R",
           arrayshid = "slurm-script/run_rcode_array.sh",
           email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", 2, "16G"))

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemm --mem 60G --ntasks=8 --time 24:00:00 slurm-script/run_rcode_array.sh
####>>> sbatch -p bigmemm --mem 16G --ntasks=2 --time 12:00:00 slurm-script/run_rcode_array.sh

dirs <- dir(path="largedata/COMET_CHG", pattern="JR", full.names=TRUE, recursive=FALSE)
for(i in 1:length(dirs)){
    file.copy(from="profiling/3.DMR/1.COMETgazer.sh", to=dirs[i], overwrite=TRUE)
}


library("farmeR")
run_Rcodes(inputdf=data.frame(file=1:20, out=1), outdir="slurm-script", cmdno=1,
           rcodes = "profiling/3.DMR/3.B.2_COMMETgazer_CHG.R",
           arrayshid = "slurm-script/run_rgazer_array.sh",
           email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", 2, "16G"))

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
####>>> sbatch -p bigmemm --mem 16G --ntasks=2 --time 12:00:00 slurm-script/run_rgazer_array.sh
