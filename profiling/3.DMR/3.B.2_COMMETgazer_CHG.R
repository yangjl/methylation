### Jinliang Yang
### 10-12-2016
### using smooth/segmentation approach for DMR
### http://bioconductor.org/packages/devel/bioc/html/bsseq.html
### manual: http://bioconductor.org/packages/devel/bioc/vignettes/bsseq/inst/doc/bsseq.pdf


##get command line args
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

JOBID <- as.numeric(as.character(args[1]))
print(JOBID)

dirs <- dir(path="largedata/COMET_CHG", pattern="JR", full.names=TRUE, recursive=FALSE)
file1 <- file1[c(9,10,16)] #for the unfinished ones

setwd(dirs[JOBID])
system("sh 1.COMETgazer.sh")
