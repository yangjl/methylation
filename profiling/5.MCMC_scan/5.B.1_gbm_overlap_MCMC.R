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

source("lib/mplots.R")
source("lib/mcmcbc.R")

library("data.table")
library("GenomicRanges")


run1mcmc <- function(JOBID, typefile="largedata/type.csv"){
    ###### main codes:    
    a <- read.csv(typefile)
    
    mya <- a[JOBID, ]
    
    ### type
    if(mya$type == "CG"){
        df <- fread("largedata/lcache/SFS_comet_blocks_CG.csv", data.table=FALSE)
    }else if(mya$type == "CHG"){
        df <- fread("largedata/lcache/SFS_comet_blocks_CHG.csv", data.table=FALSE)
    }else if(mya$type == "CHH"){
        df <- fread("largedata/lcache/SFS_comet_blocks_CHH.csv", data.table=FALSE)
    }else{
        stop("### type error!")
    }
    df$chr <- gsub("chr", "", df$chr)
    df$cid <- paste(df$chr, df$bid, sep="_")
    
    ### length quantile
    qt <- quantile(df$length)
    if(mya$length == 1){
        df <- subset(df, length <= qt[2])
    }else if(mya$length == 2){
        df <- subset(df, length > qt[2] & length <= qt[3])
    }else if(mya$length == 3){
        df <- subset(df, length > qt[3] & length <= qt[4])
    }else if(mya$length == 4){
        df <- subset(df, length > qt[4])
    }
    
    
    if(mya$TE == "yes"){
        repeats <- fread("~/dbcenter/AGP/AGPv2/repeats/ZmB73_5a_MIPS_repeats.gff", header=FALSE, data.table=FALSE)
        names(repeats) <- c("seqname", "source", "feature", "start", "end", "score",
                            "strand", "frame", "attribute")
        repeats$class <- gsub(".*type=|;name=.*", "", repeats$attribute)
        gff <- repeats
    }else if(mya$TE == "no"){
        ######## format GFF and repeat files
        gff <- fread("~/dbcenter/AGP/AGPv2/ZmB73_5b_FGS.gff", header=TRUE, data.table=FALSE)
        names(gff) <- c("seqname", "source", "feature", "start", "end", "score",
                        "strand", "frame", "attribute")
        
        if(mya$gbody == "gbm"){
            res <- read.csv("cache/stat_exon_mean_var.csv")
            gff$geneid <- gsub(";.*|_.*", "", gff$attribute)
            gff$geneid <- gsub(".*=", "", gff$geneid)
            ####
            gbM <- subset(res, mm > 0.6)
            gff <- subset(gff, geneid %in% gbM$geneid)
            
        }else if(mya$gbody == "ngbm"){
            res <- read.csv("cache/stat_exon_mean_var.csv")
            gff$geneid <- gsub(";.*|_.*", "", gff$attribute)
            gff$geneid <- gsub(".*=", "", gff$geneid)
            ####
            ngbM <- subset(res, mm <= 0.6)
            gff <- subset(gff, geneid %in% gbM$geneid)
        }
        
    }
    
    ##########
    out <- run_overlap_MCMC(df, gff, fea=mya$fs, runid=mya$id, outdir="largedata/lcache/", outid="gbody")
    
}

set.seed(12345679)
run1mcmc(JOBID, typefile="largedata/gbody_type.csv")



    
    
