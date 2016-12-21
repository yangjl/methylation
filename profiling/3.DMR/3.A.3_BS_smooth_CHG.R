### Jinliang Yang
### 10-10-2016
### purpose: BS smooth

##get command line args
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

JOBID <- as.numeric(as.character(args[1]))
print(JOBID)

### start running
library(data.table)
library(bsseq)

runbsmooth <- function(infile="largedata/wgbs_smoothed/JRA1_pe.cg.rc", outdir="largedata/COMET",
                       cores=8, chrs=1:10, BIN){
    
    meth <- fread(infile, data.table=FALSE)
    myid <- gsub(".*/|_.*", "", infile)
    myp <- paste(outdir, myid, sep="/")
    dir.create(path=myp, showWarnings = FALSE)
    
    for(i in chrs){
        #tmp <- dat[dat$strand == "+",]
        message(sprintf("[runbsmooth] Working on [%s] at [chr%s] ...", myid, i))
        
        meth0 <- subset(meth, V2 %in% i)
        meth0$V2 <- paste0("chr", meth0$V2)
       
        #meth0 <- subset(meth0, V7 > 0)
        meth0 <- meth0[order(meth0$V3), ]
        
        tmp1 <- meth0[meth0$V4 == "+",]
        #tmp1 <- tmp1[1:10000,]
        BS.forward <- BSseq(pos = tmp1$V3, chr = tmp1$V2, M = as.matrix(tmp1$V5, ncol = 1),
                            Cov = as.matrix(tmp1$V7, ncol = 1), sampleNames = "forward")
        tmp2 <- meth0[meth0$V4 == "-",]
        #tmp2 <- tmp2[1:12000,]
        BS.reverse <- BSseq(pos = tmp2$V3 - 2L, chr = tmp2$V2, M = as.matrix(tmp2$V5, ncol = 1),
                            Cov = as.matrix(tmp2$V7, ncol = 1), sampleNames = "reverse")
        BS <- combine(BS.forward, BS.reverse)
        BS <- collapseBSseq(BS, columns = c("a", "a"))
        
        tot <- ceiling(length(BS)/BIN)
        for(b in 1:tot){
            sline <- BIN*(b-1) + 1
            eline <- BIN*b
            if(eline > length(BS)){
                eline <- length(BS)
            }
            myBS <- BS[sline:eline, ]
            
            message(sprintf("[runbsmooth] Smoothing [chr%s], bin [%s/%s: %s - %s] ...",
                            i, b, tot, sline, eline))
            res <- BSmooth(myBS, ns = 70, h = 1000, maxGap = 10^8,
                           parallelBy = c("sample", "chromosome"), mc.preschedule = FALSE,
                           mc.cores = 1, keep.se = FALSE, verbose = TRUE)
            
            df <- data.frame(seqnames=seqnames(res),
                             starts=start(res),
                             ends=end(res),
                             scores=getMeth(res, type="smooth"))
            
            if(b == 1){
                write.table(df, paste0(myp, "/", "chr", i, ".txt"), 
                            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
            }else{
                write.table(df, paste0(myp, "/", "chr", i, ".txt"), 
                            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append = TRUE)
            }
        }
    }
}

### test
#runbsmooth(infile="largedata/wgbs_smoothed/JRA1_pe.cg.rc", outdir="largedata/COMET", 
#           cores=1, chrs=10, BIN=1000000)

## CG
file1 <- list.files(path="largedata/wgbs_smoothed", pattern="chg.rc$", full.names = TRUE)
file1 <- file1[c(9,10,16)]
runbsmooth(infile=file1[JOBID], outdir="largedata/COMET_CHG", cores=1, chrs=1:10, BIN=1000000)

