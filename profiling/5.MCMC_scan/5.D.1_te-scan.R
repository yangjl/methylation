### Jinliang Yang
### Jan 10th, 2017
### run COMET overlap with features and then MCMC

##get command line args
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

JOBID <- as.numeric(as.character(args[1]))
print(JOBID)

###########
source("lib/runoverlap.R")
source("lib/run_mcmc_bygeneset.R")
source("lib/mplots.R")
source("lib/mcmcbc.R")

library("data.table")
library("GenomicRanges")



set.seed(12345679)
run_mcmc_bygeneset(JOBID, typefile="largedata/te_type.csv", geneset=NULL, 
                   cutoff=median(geneset$value), outid="te_len")

