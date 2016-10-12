### Jinliang Yang
### 10-12-2016
### understand how COMMET works

library("data.table")
library(bsseq)
chr10 <- fread("largedata/COMET/JRA1/data.frame/chr10.txt", data.table=FALSE)
df <- fread("largedata/COMET/JRA1/data.frame/chr10.df.txt", data.table=FALSE)



### original data
getraw <-function(){
    meth <- fread("largedata/wgbs_smoothed/JRA1_pe.cg.rc", data.table=FALSE)
    
    raw <- subset(meth, V2 == 10 & V3 < 1000000)
    
    tmp1 <- raw[raw$V4 == "+",]
    #tmp1 <- tmp1[1:10000,]
    BS.forward <- BSseq(pos = tmp1$V3, chr = tmp1$V2, M = as.matrix(tmp1$V5, ncol = 1),
                        Cov = as.matrix(tmp1$V7, ncol = 1), sampleNames = "forward")
    tmp2 <- raw[raw$V4 == "-",]
    #tmp2 <- tmp2[1:12000,]
    BS.reverse <- BSseq(pos = tmp2$V3 - 1L, chr = tmp2$V2, M = as.matrix(tmp2$V5, ncol = 1),
                        Cov = as.matrix(tmp2$V7, ncol = 1), sampleNames = "reverse")
    BS <- combine(BS.forward, BS.reverse)
    BS <- collapseBSseq(BS, columns = c("a", "a"))
    
    df <- data.frame(seqnames=seqnames(BS), starts=start(BS), ends=end(BS),
                     scores=getMeth(BS, type="raw"))
    return(df)
}

### get chr=10, <1Mb Region
raw <- getraw()

bl <- fread("largedata/COMET/JRA1/data.frame/chr10.blocks.df.txt", data.table=FALSE)
bl <- subset(bl, chr == 10, start < 1000000)

out <- merge(bl, raw[, c("starts", "a")], by.x="start", by.y="starts")
write.table(out, "largedata/chr10_1Mb.csv", sep=",", row.names=FALSE, quote=FALSE)
################################

out <- read.csv("largedata/chr10_1Mb.csv")
chr10 <- read.table("largedata/COMET/JRA1/COMETs/chr10.blocks.verified.txt", header=T)
chr10 <- subset(chr10, stop < 1000000)
low <- subset(chr10, meth <= 0.33)
med <- subset(chr10, meth > 0.33 & meth < 0.66)
high <- subset(chr10, meth >= 0.66)



plot(out$start, out$a, pch=16, cex=0.5, col="grey", xlab="Chr10, 0-1Mb", ylab="Methylation Ratio")
points(out$start, out$score, pch=16, cex=0.6, col=makeTransparent("darkblue", alpha=0.5))

c <- makeTransparent(c("#65534f", "#36566d"), alpha=0.5)
rect(xleft=high$start, ybottom=0, xright=high$stop, ytop=high$meth, col=c[1], border = NA)
rect(xleft=low$start, ybottom=0, xright=low$stop, ytop=low$meth, col="#8e2c21", border=NA)
rect(xleft=med$start, ybottom=0, xright=med$stop, ytop=med$meth, col="#36566d", border=NA)




     