### Jinliang Yang
### 02-03-2017
### purpose: using smoothed data to summarize methylation levels

library("data.table")
library("GenomicRanges")
library("tidyr")
source("lib/find_cval_gff.R")

######## format GFF and repeat files
gff <- fread("~/dbcenter/AGP/AGPv2/ZmB73_5b_FGS.gff", header=TRUE, data.table=FALSE)
names(gff) <- c("seqname", "source", "feature", "start", "end", "score",
                "strand", "frame", "attribute")


####### smoothed files CG
pwd1 <- list.files(path="largedata/COMET", pattern="^J", full.names = TRUE)

res1 <- data.frame()
for(i in 1:length(pwd1)){
    out1 <- find_cval_gff(infile=paste0(pwd1[i], "/chr1.txt"), gff)
    res1 <- rbind(res1, out1)
}
write.csv(res1, "cache/CG_chr1_gff_feas.csv")


####### smoothed files CHG
pwd2 <- list.files(path="largedata/COMET_CHG/", pattern="^J", full.names = TRUE)

res2 <- data.frame()
for(i in 1:length(pwd2)){
    out2 <- find_cval_gff(infile=paste0(pwd2[i], "/chr1.txt"), gff)
    res2 <- rbind(res2, out2)
}
write.csv(res2, "cache/CHG_chr1_gff_feas.csv")


####### smoothed files CHG
pwd3 <- list.files(path="largedata/COMET_CHH/", pattern="^J", full.names = TRUE)

res3 <- data.frame()
for(i in 1:length(pwd3)){
    out3 <- find_cval_gff(infile=paste0(pwd3[i], "/chr1.txt"), gff)
    res3 <- rbind(res3, out3)
}
write.csv(res3, "cache/CHH_chr1_gff_feas.csv")












