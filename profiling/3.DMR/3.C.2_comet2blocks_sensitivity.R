### Jinliang Yang
### 10-12-2016
### Chop COMET into shared blocks

##get command line args
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

JOBID <- as.numeric(as.character(args[1]))
print(JOBID)

library("farmeR")
library("plyr")
library("GenomicRanges")
library("data.table")
source("lib/comet2blocks.R")

files <- list.files(path="largedata/COMET/CG_COMET", pattern="csv", full.names=T)
num <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.33, 0.35, 0.4, 0.45)

res <- comet2blocks(files, chri=10, cutoff=c(num[JOBID], 1-num[JOBID]))

out <- paste0("largedata/COMET/CG_COMET/comet_blocks_", num[JOBID], ".csv")
write.table(res, outfile, sep=",", row.names=FALSE, quote=FALSE)
    