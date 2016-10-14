### Jinliang Yang
### 10-16-2016
### purpose: get SFS for features

library(GenomicFeatures)
library("data.table")

res <- read.csv("cache/stat_exon_mean_var.csv")


comet <- fread("largedata/COMET/CG_COMET/comet_blocks.csv", data.table=FALSE)


get_comet_sfs <- function(comet, context="CG", cols=3:22){
    
    comet$miss <- apply(comet[, cols], 1, function(x){
        return(sum(x < 0))
    })
    comet <- subset(comet, miss == 0)
    ###########
    
    f <- apply(comet[, cols], 1, function(x){
        n2 <- sum(x == 2)
        n1 <- sum(x == 1)
        return(2*n2 + n1)
    })
    sfs <- table(f)
    return(sfs)
}


sfs <- getsfs(context="CHG", cols=3:22, BINSIZE=100)
write.table(sfs, "cache/sfs_cg_comet.csv", sep=",", row.names=FALSE, quote=FALSE)



