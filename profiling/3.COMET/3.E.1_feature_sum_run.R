### Jinliang Yang
### 12-01-2016
### purpose: run BSmoothing

cx <- c("largedata/COMET", "largedata/COMET_CHG", "largedata/COMET_CHH")

df <- data.frame(pwd=rep(cx, each=1), output=c(1,2,3), context=rep(c("CG", "CHG", "CHH"), each=1))
df$output <- paste0("cache/", df$context, "_chr1_gene_fea.csv")
# col, pwd="largedata/COMET"
# col: type="class=I"
# col: output="cache/CG_chr1_TE_class1.csv"
write.csv(df, "largedata/run_gene_df.csv")

library("farmeR")
run_Rcodes(inputdf=data.frame(file=1:3, out=1), outdir="slurm-script", cmdno=1,
           rcodes = "profiling/3.COMET/3.E.1_feature_summary.R",
           arrayshid = "slurm-script/run_fea_array.sh",
           email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", 10, "80G"))

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemm --time 24:00:00 --mem 80G --ntasks=10 --exclude=bigmem1,bigmem2 slurm-script/run_fea_array.sh
