### Jinliang Yang
### 02-03-2017
### purpose: using smoothed data to summarize methylation levels

##get command line args
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

JOBID <- as.numeric(as.character(args[1]))
print(JOBID)


library("data.table")
library("GenomicRanges")
library("tidyr")
source("lib/find_cval_gff.R")

runit_gene <- function(df, id){
    mypwd <- as.character(df$pwd[id])
    #mytype <- as.character(df$type[id])
    myout <- as.character(df$output[id]) 
    
    ######## format GFF and repeat files
    gff <- fread("~/dbcenter/AGP/AGPv2/ZmB73_5b_FGS.gff", header=TRUE, data.table=FALSE)
    names(gff) <- c("seqname", "source", "feature", "start", "end", "score",
                    "strand", "frame", "attribute")
    
    ####### smoothed files CG
    pwd1 <- list.files(path=mypwd, pattern="^J", full.names = TRUE)
    
    res1 <- data.frame()
    for(i in 1:length(pwd1)){
        out1 <- find_cval_gff(infile=paste0(pwd1[i], "/chr1.txt"), gff,
                              features=c("exon", "intron", "up1k", "down1k", "gene", "CDS"))
        res1 <- rbind(res1, out1)
    }
    write.csv(res1, myout)
}


## control elements
df <- read.csv("largedata/run_gene_df.csv")
# col, pwd="largedata/COMET"
# col: output="cache/CG_chr1_TE_class1.csv"
runit_gene(df, id=JOBID)










