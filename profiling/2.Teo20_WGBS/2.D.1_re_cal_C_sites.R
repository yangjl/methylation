### Jinliang Yang
### 09-13-2016
### purpose: re-calibrate the C/T and G/A heterozygote sites



########
library(data.table)
files <- list.files(path="largedata/wgbs_smoothed", pattern="cg$", full.names = TRUE)

res <- re_calC(files)

run_Rcodes(inputdf, outdir, cmdno = 100, rcodes = "lib/C_format.R",
           arrayshid = "slurm-script/run_bcf_query_array.sh", email = NULL,
           runinfo = c(FALSE, "bigmemh", 1))
