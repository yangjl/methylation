### Jinliang Yang
### Jan 10th, 2017
### run COMET overlap with features and then MCMC

##get command line args
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

JOBID <- as.numeric(as.character(args[1]))
print(JOBID)

###########

source("lib/runmcmc_te_geneset.R")

library("data.table")
library("GenomicRanges")
library(gsl) #Gnu scientific Library is a collection of numerical routines for scientific computing.
library(coda) #Output Analysis and Diagnostics for MCMC
library(utils)
library(tidyr)
library(dplyr)
library(cowplot)


set.seed(1234579)

runmcmc_te_geneset(JOBID, inputdf="largedata/mcmc_scan_input.csv")
## JOBID: control the array job number. [num, 1]
## inputdf: [data.frame, cols=comet_file:(comet blocks file, i.e. largedata/lcache/SFS_comet_blocks_CG.csv), 
#                    length:(blocks quantile len, 0,1,2,3,4), TE:(yes, no), 
#                    outRD:(chr, "out.RData"),
#                    geneset_file:(chr, csv, cols: geneid, value, i.e. "geneset.tex"),
#                    
#                    optional:
#                    feature:(chr, "exon", "intron", "up1k", "gene"),
#                    gset:(chr, "below" or "above", "wholeset"),
#                    cutoff:(num, num, 0.5)]
## geneset: gene set to determine which set to cal sfs. [data.frame, cols: geneid, value]
## cutoff: cutoff for the values to get up set and down set of the gene. [num, 0.5]

    