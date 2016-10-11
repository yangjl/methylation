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
                       cores=8){
    
    meth <- fread(infile, data.table=FALSE)
    myid <- gsub(".*/|_.*", "", infile)
    myp <- paste(outdir, myid, sep="/")
    dir.create(path=myp, showWarnings = TRUE)
    for(i in 1:10){
        #tmp <- dat[dat$strand == "+",]
        meth0 <- subset(meth, V2 == i)
        meth0$V2 <- paste0("chr", meth0$V2)
        #meth0 <- subset(meth0, V7 > 0)
        meth0 <- meth0[order(meth0$V3), ]
        
        BS <- BSseq(chr = meth0$V2, pos = meth0$V3,
                    M = as.matrix(meth0$V5, ncol = 1),
                    Cov = as.matrix(meth0$V7, ncol = 1), 
                    sampleNames = myid)
        
        tmp <- meth0[meth0$V4 == "+",]
        BS.forward <- BSseq(pos = tmp$V3, chr = tmp$V2, M = as.matrix(tmp$V5, ncol = 1),
                            Cov = as.matrix(tmp$V7, ncol = 1), sampleNames = "forward")
        tmp <- meth0[meth0$V4 == "-",]
        BS.reverse <- BSseq(pos = tmp$V3, chr = tmp$V2, M = as.matrix(tmp$V5, ncol = 1),
                            Cov = as.matrix(tmp$V7, ncol = 1), sampleNames = "reverse")
        BS <- combine(BS.forward, BS.reverse)
        BS <- collapseBSseq(BS, columns = c("a", "a"))

        
        res <- BSmooth(BS, ns = 70, h = 1000, maxGap = 10^8,
                       parallelBy = c("sample", "chromosome"), mc.preschedule = FALSE,
                       mc.cores = cores, keep.se = FALSE, verbose = TRUE)
        
        message(sprintf("###>>> output [chr%s] in [%s]", i, myp))
        meth0$sm <- getMeth(res, type="smooth")
        write.table( meth0[, c(2,3,3,8)], paste0(myp, "/", "chr", i, ".txt"), 
                     sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    }
}

### test
#runbsmooth(infile="largedata/wgbs_smoothed/JRA1_pe.cg.rc", outdir="largedata/COMET", cores=8)

## CG
file1 <- list.files(path="largedata/wgbs_smoothed", pattern="cg.rc$", full.names = TRUE)
runbsmooth(infile=file1[JOBID], outdir="largedata/COMET", cores=8)

