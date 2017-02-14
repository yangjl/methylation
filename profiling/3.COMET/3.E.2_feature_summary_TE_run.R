### Jinliang Yang
### 12-01-2016
### purpose: run BSmoothing

cx <- c("largedata/COMET", "largedata/COMET_CHG", "largedata/COMET_CHH")
cl <- c("class=I", "class=II")

df <- data.frame(pwd=rep(cx, each=2), class=rep(cl), output=c(1,2), context=rep(c("CG", "CHG", "CHH"), each=2))
df$output <- paste0("cache/", df$context, "_chr1_TE_class", df$output, ".csv")
# col, pwd="largedata/COMET"
# col: type="class=I"
# col: output="cache/CG_chr1_TE_class1.csv"
write.table(df, "largedata/run_df.csv")

library("farmeR")
run_Rcodes(inputdf=data.frame(file=1:6, out=1), outdir="slurm-script", cmdno=1,
           rcodes = "profiling/3.COMET/3.E.2_feature_summary_TE.R",
           arrayshid = "slurm-script/run_TE_array.sh",
           email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", 8, "60G"))

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemm --mem 60G --ntasks=8 --time 24:00:00 slurm-script/run_rcode_array.sh
