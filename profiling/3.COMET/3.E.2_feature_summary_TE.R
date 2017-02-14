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


#########################################
runit <- function(df, id, gff){
    
    mypwd <- df$pwd[id]
    mytype <- df$type[id]
    myout <- df$out[id]
    
    ####### smoothed files CG
    pwd1 <- list.files(path=mypwd, pattern="^J", full.names = TRUE)
    
    res1 <- data.frame()
    for(i in 1:length(pwd1)){
        #out1 <- find_cval_gff(infile=paste0(pwd1[i], "/chr1.txt"), gff)
        out1 <- find_cval_gff(infile=paste0(pwd1[i], "/chr1.txt"), gff=subset(gff, feature %in% mytype), 
                              features=c(mytype, "up1k", "down1k"))
        res1 <- rbind(res1, out1)
    }
    write.csv(res1, myout)
}


######## format GFF and repeat files
# AGPv2 three annotation files

repeats <- fread("~/dbcenter/AGP/AGPv2/repeats/ZmB73_5a_MTEC_repeats.gff", header=FALSE, data.table=FALSE)
names(repeats) <- c("seqname", "source", "feature", "start", "end", "score",
                    "strand", "frame", "attribute")
repeats$class <- gsub(";.*", "", repeats$attribute)
gff <- repeats[, -3]
names(gff)[ncol(gff)] <- "feature"


## control elements
df <- read.csv("lardedata/run_df.csv")
# col, pwd="largedata/COMET"
# col: type="class=I"
# col: output="cache/CG_chr1_TE_class1.csv"
runit(df, id=JOBID, gff)










