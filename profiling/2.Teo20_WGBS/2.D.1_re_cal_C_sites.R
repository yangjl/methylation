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

### CG
#file1 <- list.files(path="largedata/wgbs_smoothed", pattern="cg$", full.names = TRUE)
#out1 <- re_calc(file1[JOBID])
#write.table(out1, paste0(file1[JOBID], ".out"), sep=",", row.names=FALSE, quote=FALSE)

## CHG
#file2 <- list.files(path="largedata/wgbs_smoothed", pattern="chg$", full.names = TRUE)
#out2 <- re_calc(file2[JOBID])
#write.table(out2, paste0(file2[JOBID], ".out"), sep=",", row.names=FALSE, quote=FALSE)

## CHH
file3 <- list.files(path="largedata/wgbs_smoothed", pattern="chh$", full.names = TRUE)
out3 <- re_calc(file3[JOBID])
fwrite(out3, paste0(file3[JOBID], ".out"), sep=",", row.names=FALSE, quote=FALSE)

