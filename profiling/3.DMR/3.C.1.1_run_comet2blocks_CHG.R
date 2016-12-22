### Jinliang Yang
### 10-12-2016
### Chop COMET into shared blocks

library("farmeR")
library("plyr")
library("GenomicRanges")
library("data.table")
source("lib/comet2blocks.R")
# https://www.r-bloggers.com/fitting-mixture-distributions-with-the-r-package-mixtools/
library(mixtools)


files <- list.files(path="largedata/COMET/CHG_COMET", pattern="COMET.csv", full.names=T)
res <- data.frame()
for(myi in 1:10){
    out <- comet2blocks(files, chri=myi, collapse=TRUE, verbose=T, cutoff=NULL)
    #out <- comet2blocks(files, chri=myi, collapse=TRUE, verbose=T, cutoff=c(0.2, 0.7))
    res <- rbind(out, res)
}

write.table(out, "largedata/COMET/CG_COMET/chrall_comet_blocks.csv", sep=",", row.names=FALSE, quote=FALSE)
    

################ 
res <- fread("largedata/COMET_CHG/CHG_COMET/chrall_comet_blocks.csv", data.table=FALSE)

df <- res[, 1:2]
df$sfs <- apply(res[,-1:-2], 1, sum)
df$start <- as.numeric(as.character(gsub("_.*", "", df$bid)))
df$end <- as.numeric(as.character(gsub(".*_", "", df$bid)))
df$length <- df$end - df$start + 1

write.table(df, "cache/SFS_comet_blocks_CG.csv", sep=",", row.names=FALSE, quote=FALSE)



####### plot
b <- read.csv("cache/SFS_comet_blocks_CG.csv")
hist(log10(b$length), xlab="COMET Length log10(bp)", main="CG COMETs", col="#cdc0b0")
abline(v=quantile(log10(b$length))[2:4], lwd=3, lty=2 )
# 0%   25%   50%   75%  100% 
# 1    96   309   784 76180 





### determining variable or conserved sites
vout <- res[, 1:2]
vout$sites <- unlist(apply(res[, -1:-2], 1, function(x){
    x <- x[x >= 0]
    tab <- table(x)
    if(length(tab) > 1){
        return(3)
    }else if(length(tab) == 1){
        return(names(tab))
    }else{
        return(9)
    }
}))
 
vout$start <- as.numeric(as.character(gsub("_.*", "", vout$bid)))
vout$end <- as.numeric(as.character(gsub(".*_", "", vout$bid)))
vout$bp <- vout$end - vout$start

sum(subset(vout, sites == 2)$bp)
sum(subset(vout, sites == 1)$bp)
sum(subset(vout, sites == 0)$bp)
sum(subset(vout, sites == 3)$bp)



head(subset(res, bid %in% subset(vout, sites==0)$bid))


vout0 <- subset(vout, sites == 0) # low
vout2 <- subset(vout, sites == 2) # high
vout1 <- subset(vout, sites == 1) # med
vout3 <- subset(vout, sites == 3) #variable

save(list=c("vout0", "vout1", "vout2", "vout3"), file="largedata/chr10_comet.RData")

plot(c(1, max(vout$end)), c(0, 4), type="n", pch=16, cex=0.5, col="grey", 
     xlab="Chr10", ylab="Methylation Ratio", main="Chr10 CG methylation (N=20)")
#points(out$start, out$score, pch=16, cex=0.6, col=makeTransparent("darkblue", alpha=0.5))

library(farmeR)
c <- makeTransparent(c("#65534f", "#36566d"), alpha=0.6)
rect(xleft=vout0$start, ybottom=0, xright=vout0$end, ytop=1, col="red")
rect(xleft=vout1$start, ybottom=1, xright=vout1$end, ytop=2, col="#8e2c21")
rect(xleft=vout2$start, ybottom=2, xright=vout2$end, ytop=3, col=c[1], border=NA)
rect(xleft=vout3$start, ybottom=3, xright=vout3$end, ytop=4, col=c[2], border= NA)



  