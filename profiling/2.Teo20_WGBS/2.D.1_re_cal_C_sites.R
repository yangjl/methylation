### Jinliang Yang
### 09-13-2016
### purpose: re-calibrate the C/T and G/A heterozygote sites


##get command line args
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

JOBID <- as.numeric(as.character(args[1]))
print(JOBID)
###########
library(data.table)
source("lib/re_calc.R")

files <- list.files(path="largedata/wgbs_smoothed", pattern="cg$", full.names = TRUE)
out <- re_calc(files[JOBID])

write.table(out, paste0(files[JOBID], ".out"), sep=",", row.names=FALSE, quote=FALSE)
